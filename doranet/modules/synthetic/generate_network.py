"""Synthetic network generation code."""

import collections.abc
import dataclasses
import re
import time
import typing
from datetime import datetime

from rdkit import Chem
from rdkit.Chem import rdqueries
from rdkit.Chem.rdMolDescriptors import CalcMolFormula

import doranet as dn
from doranet import interfaces, metadata
from doranet.modules.synthetic.Reaction_Smarts_Forward import op_smarts
from doranet.modules.synthetic.Reaction_Smarts_Retro import op_retro_smarts


@typing.final
@dataclasses.dataclass(frozen=True, slots=True)
class Chem_Rxn_dH_Calculator(metadata.RxnPropertyCalc[float]):
    dH_key: collections.abc.Hashable
    direction: str
    Mole_Hf: typing.Callable[[str], float]

    @property
    def key(self) -> collections.abc.Hashable:
        return self.dH_key

    @property
    def resolver(self) -> metadata.MetaDataResolverFunc[float]:
        return metadata.TrivialMetaDataResolverFunc

    def __call__(
        self,
        data: interfaces.ReactionExplicit,
        prev_value: typing.Optional[float] = None,
    ) -> typing.Optional[float]:
        if prev_value is not None:
            return prev_value
        if data.reaction_meta is not None and self.dH_key in data.reaction_meta:
            return None
        if self.Mole_Hf is None:
            return "No_Thermo"
        dH = 0.0

        if data.operator.meta is None:
            raise RuntimeError("No operator metadata found!")

        for idx, mol in enumerate(data.products):
            if not isinstance(mol.item, interfaces.MolDatRDKit):
                raise NotImplementedError(
                    f"""Calculator only implemented for molecule type \
                        MolDatRDKit, not {type(mol.item)}"""
                )
            pro_Hf = self.Mole_Hf(mol.item.smiles)
            if pro_Hf is None:
                return float("nan")
            dH = dH + pro_Hf * data.operator.meta["products_stoi"][idx]
        for idx, mol in enumerate(data.reactants):
            if not isinstance(mol.item, interfaces.MolDatRDKit):
                raise NotImplementedError(
                    f"""Calculator only implemented for molecule type \
                        MolDatRDKit, not {type(mol.item)}"""
                )
            rea_Hf = self.Mole_Hf(mol.item.smiles)
            if rea_Hf is None:
                return float("nan")
            dH = dH - rea_Hf * data.operator.meta["reactants_stoi"][idx]
        correction = data.operator.meta["enthalpy_correction"]
        if correction is None:
            correction = 0
        if self.direction == "forward":
            dH = (dH + correction) / data.operator.meta["number_of_steps"]
        if self.direction == "retro":
            dH = (-dH + correction) / data.operator.meta["number_of_steps"]
        return round(dH, 4)

    @property
    def meta_required(self) -> interfaces.MetaKeyPacket:
        return interfaces.MetaKeyPacket(
            operator_keys={
                "reactants_stoi",
                "products_stoi",
                "enthalpy_correction",
                "number_of_steps",
            }
        )


@typing.final
@dataclasses.dataclass(frozen=True)
class Rxn_dH_Filter(metadata.ReactionFilterBase):
    __slots__ = ("max_dH", "dH_key")
    max_dH: float
    dH_key: collections.abc.Hashable

    def __call__(self, recipe: interfaces.ReactionExplicit) -> bool:
        # if metadata not available, reject reaction
        # (clearly something has gone wrong here)
        if recipe.reaction_meta is None:
            return False
        dH = recipe.reaction_meta[self.dH_key]
        if dH == "No_Thermo":
            return True
        if dH == float("nan"):
            return False
        return dH < self.max_dH

    @property
    def meta_required(self) -> interfaces.MetaKeyPacket:
        return interfaces.MetaKeyPacket()


@typing.final
@dataclasses.dataclass(frozen=True)
class Ring_Issues_Filter(metadata.ReactionFilterBase):
    def __call__(self, recipe: interfaces.ReactionExplicit) -> bool:
        if recipe.operator.meta is None:
            raise RuntimeError("No operator metadata found!")

        if (
            recipe.operator.meta["ring_issue"] is True
            and recipe.operator.meta["enthalpy_correction"] is None
        ):
            reactants_dict: dict[str, float | int] = dict()
            products_dict: dict[str, float | int] = dict()
            pattern = r"([A-Z][a-z]*)(\d*)"
            for idx, mol in enumerate(recipe.reactants):
                if not isinstance(mol.item, interfaces.MolDatRDKit):
                    raise NotImplementedError(
                        f"""Calculator only implemented for molecule type \
                            MolDatRDKit, not {type(mol.item)}"""
                    )
                smiles = CalcMolFormula(mol.item.rdkitmol)
                matches = re.findall(pattern, smiles)
                for match in matches:
                    element, count = match
                    count = int(count) if count else 1
                    reactants_dict[element] = (
                        reactants_dict.get(element, 0)
                        + count * recipe.operator.meta["reactants_stoi"][idx]
                    )
            for idx, mol in enumerate(recipe.products):
                if not isinstance(mol.item, interfaces.MolDatRDKit):
                    raise NotImplementedError(
                        f"""Calculator only implemented for molecule type \
                            MolDatRDKit, not {type(mol.item)}"""
                    )
                if "." in mol.item.smiles:
                    # if there're fragments in a mol, indicates invalid rxn
                    return False
                smiles = CalcMolFormula(mol.item.rdkitmol)
                matches = re.findall(pattern, smiles)
                for match in matches:
                    element, count = match
                    count = int(count) if count else 1
                    products_dict[element] = (
                        products_dict.get(element, 0)
                        + count * recipe.operator.meta["products_stoi"][idx]
                    )
            if reactants_dict != products_dict:
                return False
        return True

    @property
    def meta_required(self) -> interfaces.MetaKeyPacket:
        return interfaces.MetaKeyPacket(
            operator_keys={
                "reactants_stoi",
                "products_stoi",
                "ring_issue",
                "enthalpy_correction",
            }
        )


@typing.final
@dataclasses.dataclass(frozen=True)
class Retro_Not_Aromatic_Filter(metadata.ReactionFilterBase):
    def __call__(self, recipe: interfaces.ReactionExplicit) -> bool:
        rea_aro_ring_num = 0
        pro_aro_ring_num = 0

        if recipe.operator.meta is None:
            raise RuntimeError("No operator metadata found!")

        if recipe.operator.meta["Retro_Not_Aromatic"] is True:
            for idx, mol in enumerate(recipe.reactants):
                if not isinstance(mol.item, interfaces.MolDatRDKit):
                    raise NotImplementedError(
                        f"""Filter only implemented for molecule type \
                            MolDatRDKit, not {type(mol.item)}"""
                    )
                rea_aro_ring_num += (
                    Chem.rdMolDescriptors.CalcNumAromaticRings(
                        mol.item.rdkitmol
                    )
                    * recipe.operator.meta["reactants_stoi"][idx]
                )
            for idx, mol in enumerate(recipe.products):
                if not isinstance(mol.item, interfaces.MolDatRDKit):
                    raise NotImplementedError(
                        f"""Filter only implemented for molecule type \
                            MolDatRDKit, not {type(mol.item)}"""
                    )
                pro_aro_ring_num += (
                    Chem.rdMolDescriptors.CalcNumAromaticRings(
                        mol.item.rdkitmol
                    )
                    * recipe.operator.meta["products_stoi"][idx]
                )
            if rea_aro_ring_num < pro_aro_ring_num:
                return False
        return True

    @property
    def meta_required(self) -> interfaces.MetaKeyPacket:
        return interfaces.MetaKeyPacket(
            operator_keys={
                "reactants_stoi",
                "products_stoi",
                "Retro_Not_Aromatic",
            }
        )


@typing.final
@dataclasses.dataclass(frozen=True)
class Check_balance_filter(metadata.ReactionFilterBase):
    def __call__(self, recipe: interfaces.ReactionExplicit) -> bool:
        if (
            recipe.operator.meta is not None
            and recipe.operator.meta["enthalpy_correction"] is None
            and recipe.operator.meta["ring_issue"] is False
        ):
            reactants_dict: dict[str, float | int] = dict()
            products_dict: dict[str, float | int] = dict()
            pattern = r"([A-Z][a-z]*)(\d*)"
            for idx, mol in enumerate(recipe.reactants):
                if not isinstance(mol.item, interfaces.MolDatRDKit):
                    raise NotImplementedError(
                        f"""Filter only implemented for molecule type \
                            MolDatRDKit, not {type(mol.item)}"""
                    )
                smiles = CalcMolFormula(mol.item.rdkitmol)
                matches = re.findall(pattern, smiles)
                for match in matches:
                    element, count = match
                    count = int(count) if count else 1
                    reactants_dict[element] = (
                        reactants_dict.get(element, 0)
                        + count * recipe.operator.meta["reactants_stoi"][idx]
                    )
            for idx, mol in enumerate(recipe.products):
                if not isinstance(mol.item, interfaces.MolDatRDKit):
                    raise NotImplementedError(
                        f"""Filter only implemented for molecule type \
                            MolDatRDKit, not {type(mol.item)}"""
                    )
                smiles = CalcMolFormula(mol.item.rdkitmol)
                matches = re.findall(pattern, smiles)
                for match in matches:
                    element, count = match
                    count = int(count) if count else 1
                    products_dict[element] = (
                        products_dict.get(element, 0)
                        + count * recipe.operator.meta["products_stoi"][idx]
                    )
            if reactants_dict != products_dict:
                return False
        return True

    @property
    def meta_required(self) -> interfaces.MetaKeyPacket:
        return interfaces.MetaKeyPacket(
            operator_keys={
                "reactants_stoi",
                "products_stoi",
                "ring_issue",
                "enthalpy_correction",
                "name",
            }
        )


@typing.final
@dataclasses.dataclass(frozen=True)
class Cross_Reaction_Filter(interfaces.RecipeFilter):
    # no rea1 + rea2 or pure helper rxns. rea1 + rea1 or rea1 + helper ok
    # __slots__ = "coreactants"
    coreactants: collections.abc.Container[interfaces.MolIndex]

    def __call__(self, recipe: interfaces.RecipeExplicit) -> bool:
        reactant_set = set()
        for mol in recipe.reactants:
            if mol.i not in self.coreactants:
                reactant_set.add(mol.i)
        return len(set(reactant_set)) == 1

    @property
    def meta_required(self) -> interfaces.MetaKeyPacket:
        return interfaces.MetaKeyPacket()


@typing.final
@dataclasses.dataclass(frozen=True)
class Enol_filter_forward(metadata.ReactionFilterBase):
    # if an enol is in reactants, only allow tautomerization rxn
    def __call__(self, recipe: interfaces.ReactionExplicit) -> bool:
        for mol in recipe.reactants:
            if not isinstance(mol.item, interfaces.MolDatRDKit):
                raise NotImplementedError(
                    f"""Calculator only implemented for molecule type \
                        MolDatRDKit, not {type(mol.item)}"""
                )
            if (
                mol.item.rdkitmol.HasSubstructMatch(
                    Chem.MolFromSmarts("[C]=[C]-[OH]")
                )
                and recipe.operator.meta is not None
                and (
                    recipe.operator.meta["name"]
                    != "Keto-enol Tautomerization Reverse"
                )
            ):
                return False
        return True

    @property
    def meta_required(self) -> interfaces.MetaKeyPacket:
        return interfaces.MetaKeyPacket(operator_keys={"name"})


@typing.final
@dataclasses.dataclass(frozen=True)
class Enol_filter_retro(metadata.ReactionFilterBase):
    # if an enol is in reactants, only allow tautomerization rxn
    def __call__(self, recipe: interfaces.ReactionExplicit) -> bool:
        for mol in recipe.products:
            if not isinstance(mol.item, interfaces.MolDatRDKit):
                raise NotImplementedError(
                    f"""Filter only implemented for molecule type \
                        MolDatRDKit, not {type(mol.item)}"""
                )
            if (
                mol.item.rdkitmol.HasSubstructMatch(
                    Chem.MolFromSmarts("[C]=[C]-[OH]")
                )
                and recipe.operator.meta is not None
                and (
                    recipe.operator.meta["name"]
                    != "Keto-enol Tautomerization Reverse"
                )
            ):
                return False
        return True

    @property
    def meta_required(self) -> interfaces.MetaKeyPacket:
        return interfaces.MetaKeyPacket(operator_keys={"name"})


@typing.final
@dataclasses.dataclass(frozen=True)
class Allowed_Elements_Filter(metadata.ReactionFilterBase):
    # only allow reactions with specified elements in reactants.
    # does not check hydrogen
    def __call__(self, recipe: interfaces.ReactionExplicit) -> bool:
        if (
            recipe.operator.meta is None
            or recipe.operator.meta["allowed_elements"][0] == "All"
        ):
            return True
        for mol in recipe.reactants:
            if not isinstance(mol.item, interfaces.MolDatRDKit):
                raise NotImplementedError(
                    f"""Filter only implemented for molecule type \
                        MolDatRDKit, not {type(mol.item)}"""
                )
            for atom in mol.item.rdkitmol.GetAtoms():
                if (
                    atom.GetSymbol()
                    not in recipe.operator.meta["allowed_elements"]
                ):
                    return False
        return True

    @property
    def meta_required(self) -> interfaces.MetaKeyPacket:
        return interfaces.MetaKeyPacket(operator_keys={"allowed_elements"})


@typing.final
@dataclasses.dataclass(frozen=True, slots=True)
class Max_Atoms_Filter(metadata.ReactionFilterBase):
    # similar to engine.filter.reaction.max_atoms
    # take a dict, check multiple atoms {6: 12, 8: 2}
    max_atoms_dict: dict

    def __call__(self, recipe: interfaces.ReactionExplicit) -> bool:
        if self.max_atoms_dict is None:
            return True
        for atom_number, max_atoms in self.max_atoms_dict.items():
            for mol in recipe.products:
                if not isinstance(mol.item, interfaces.MolDatRDKit):
                    raise NotImplementedError(
                        f"""Filter only implemented for molecule type \
                            MolDatRDKit, not {type(mol.item)}"""
                    )
                atom_q = rdqueries.AtomNumEqualsQueryAtom(atom_number)
                atom_matches = mol.item.rdkitmol.GetAtomsMatchingQuery(atom_q)
                if len(atom_matches) > max_atoms:
                    return False
        return True

    @property
    def meta_required(self) -> interfaces.MetaKeyPacket:
        return interfaces.MetaKeyPacket()

@typing.final
@dataclasses.dataclass(frozen=True)
class Free_Radical_Polymerization_Filter(metadata.ReactionFilterBase):
    __slots__ = ("frp_types")
    frp_types: tuple

    # only allow products with functional groups that allow for FRP
    def __call__(self, recipe: interfaces.ReactionExplicit) -> bool:
        if self.frp_types is None:
            return True
        
        frp_substructures = {
            "ionic": Chem.MolFromSmarts("[CX3H1](=[O,S])[#6]"),
            "ROMP": Chem.MolFromSmarts("[RC]=[RC]"),
            "vinyl": Chem.MolFromSmarts("[!RC]=[C]"),
            "ROP": Chem.MolFromSmarts("[O,S;X1]=[C,S;X3;R]-[N,O,S;R]"),
        }

        for mol in recipe.products:
            if not isinstance(mol.item, interfaces.MolDatRDKit):
                raise NotImplementedError(
                    f"""Filter only implemented for molecule type \
                        MolDatRDKit, not {type(mol.item)}"""
                )
            
            # loop through FRP types and see if possible
            # return at first instance of available polymerization
            for frp_type in self.frp_types:
                if mol.item.rdkitmol.GetSubstructMatch(frp_substructures[frp_type]):
                    return True
            
            # if no matches
            return False
            
        return True

    @property
    def meta_required(self) -> interfaces.MetaKeyPacket:
        return interfaces.MetaKeyPacket(operator_keys={"frp_types"})

def get_smiles_from_file(file_name):
    def is_valid_smiles(smiles_string):
        try:
            mol = Chem.MolFromSmiles(smiles_string)
            return mol is not None
        except ValueError:
            return False

    if not isinstance(file_name, str):
        return file_name
    if is_valid_smiles(file_name):
        return [file_name]
    with open(
        file_name  # , encoding="utf-8"
    ) as f:
        lines = f.readlines()
    clean_list = list()
    for i in lines:
        if i != "\n":
            clean_list.append(i.strip())
    return clean_list


def generate_network(
    job_name="default_job",
    starters=False,
    helpers=False,
    gen=1,
    direction="forward",
    molecule_thermo_calculator=None,
    max_rxn_thermo_change=15,
    frp_types=None,
    max_atoms=None,  # {"C": 20}
    allow_multiple_reactants="default",  # forward allowed, retro no
    targets=None,  # string or list, set, etc.
):
    if not starters:
        raise Exception("At least one starter is needed to generate a network")

    starters = get_smiles_from_file(starters)
    helpers = get_smiles_from_file(helpers)
    targets = get_smiles_from_file(targets)

    print(f"Job name: {job_name}")
    print(f"Job type: synthetic network expansion {direction}")
    print("Job started on:", datetime.now())
    start_time = time.time()

    engine = dn.create_engine()
    network = engine.new_network()

    if helpers:
        for smiles in helpers:
            network.add_mol(engine.mol.rdkit(smiles))

    my_start_i = -1
    for smiles in starters:
        if my_start_i == -1:
            my_start_i = network.add_mol(engine.mol.rdkit(smiles))
        else:
            network.add_mol(engine.mol.rdkit(smiles))

    if direction == "forward":
        smarts_list = op_smarts
    elif direction == "retro":
        smarts_list = op_retro_smarts

    for smarts in smarts_list:
        if smarts.kekulize_flag is False:
            network.add_op(
                engine.op.rdkit(
                    smarts.smarts,
                    drop_errors=True,
                ),
                meta={
                    "name": smarts.name,
                    "reactants_stoi": smarts.reactants_stoi,
                    "products_stoi": smarts.products_stoi,
                    "enthalpy_correction": smarts.enthalpy_correction,
                    "ring_issue": smarts.ring_issue,
                    "kekulize_flag": smarts.kekulize_flag,
                    "Retro_Not_Aromatic": smarts.Retro_Not_Aromatic,
                    "number_of_steps": smarts.number_of_steps,
                    "allowed_elements": smarts.allowed_elements,
                    "Reaction_type": smarts.reaction_type,
                    "Reaction_direction": direction,
                },
            )
        if smarts.kekulize_flag is True:
            network.add_op(
                engine.op.rdkit(
                    smarts.smarts,
                    kekulize=True,
                    drop_errors=True,
                ),
                meta={
                    "name": smarts.name,
                    "reactants_stoi": smarts.reactants_stoi,
                    "products_stoi": smarts.products_stoi,
                    "enthalpy_correction": smarts.enthalpy_correction,
                    "ring_issue": smarts.ring_issue,
                    "kekulize_flag": smarts.kekulize_flag,
                    "Retro_Not_Aromatic": smarts.Retro_Not_Aromatic,
                    "number_of_steps": smarts.number_of_steps,
                    "allowed_elements": smarts.allowed_elements,
                    "Reaction_type": smarts.reaction_type,
                    "Reaction_direction": direction,
                },
            )

    strat = engine.strat.cartesian(network)

    periodic_table = Chem.GetPeriodicTable()

    if max_atoms is None:
        max_atoms_dict_num = None
    else:
        max_atoms_dict_num = dict()
        for key in max_atoms:
            max_atoms_dict_num[periodic_table.GetAtomicNumber(key)] = max_atoms[
                key
            ]

    if direction == "forward":
        reaction_plan = (
            Max_Atoms_Filter(max_atoms_dict_num)
            >> Ring_Issues_Filter()
            >> Enol_filter_forward()
            >> Check_balance_filter()
            >> Allowed_Elements_Filter()
            >> Chem_Rxn_dH_Calculator(
                "dH", "forward", molecule_thermo_calculator
            )
            >> Rxn_dH_Filter(max_rxn_thermo_change, "dH")
            >> Free_Radical_Polymerization_Filter(frp_types)
        )
        recipe_filter = None

    elif direction == "retro":
        reaction_plan = (
            Max_Atoms_Filter(max_atoms_dict_num)
            >> Ring_Issues_Filter()
            >> Retro_Not_Aromatic_Filter()
            >> Enol_filter_retro()
            >> Allowed_Elements_Filter()
            >> Check_balance_filter()
            >> Chem_Rxn_dH_Calculator("dH", "retro", molecule_thermo_calculator)
            >> Rxn_dH_Filter(max_rxn_thermo_change, "dH")
            >> Free_Radical_Polymerization_Filter("frp_types")
        )
        recipe_filter = Cross_Reaction_Filter(tuple(range(my_start_i)))

    if allow_multiple_reactants != "default":
        if allow_multiple_reactants is True:
            recipe_filter = None
        elif allow_multiple_reactants is False:
            recipe_filter = Cross_Reaction_Filter(tuple(range(my_start_i)))

    bundle_filter = engine.filter.bundle.coreactants(tuple(range(my_start_i)))

    ini_number = len(network.mols)

    strat.expand(
        num_iter=gen,
        reaction_plan=reaction_plan,
        bundle_filter=bundle_filter,
        recipe_filter=recipe_filter,
        save_unreactive=False,
    )

    if targets is not None:
        print("Checking for targets")
        to_check = set()
        if isinstance(targets, str):
            to_check.add(Chem.MolToSmiles(Chem.MolFromSmiles(targets)))
        else:
            for i in targets:
                to_check.add(Chem.MolToSmiles(Chem.MolFromSmiles(i)))

        for mol in network.mols:
            if (
                network.reactivity[network.mols.i(mol.uid)] is True
                and mol.uid in to_check
            ):
                print("Target found for", mol.uid)

    print("Number of generations:", gen)
    print("Number of operators:", len(network.ops))
    print("Number of molecules before expantion:", ini_number)
    print("Number of molecules after expantion:", len(network.mols))
    print("Number of reactions:", len(network.rxns))

    end_time = time.time()
    elapsed_time = (end_time - start_time) / 60
    print(
        "Time used for network generation:",
        "{:.2f}".format(elapsed_time),
        "minutes",
    )
    print()

    network.save_to_file(f"{job_name}_{direction}_saved_network")

    return network
