"""Enzymatic network generation code."""

import collections.abc
import dataclasses
import re
import time
import typing
from datetime import datetime
from pathlib import Path

import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, rdqueries
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
from rdkit.Chem.rdmolops import GetFormalCharge

import doranet as dn
from doranet import interfaces, metadata


def clean_SMILES(smiles):
    mol = Chem.MolFromSmiles(smiles)
    Chem.rdmolops.RemoveStereochemistry(mol)
    cpd_smiles = Chem.MolToSmiles(mol)
    return cpd_smiles


def neutralise_charges(smiles):
    mol = Chem.MolFromSmiles(smiles)
    patts = (
        # Imidazoles
        ("[n+;H]", "n"),
        # Amines
        ("[N+;!H0]", "N"),
        # Carboxylic acids and alcohols
        ("[$([O-]);!$([O-][#7])]", "O"),
        # Thiols
        ("[S-;X1]", "S"),
        # Sulfonamides
        ("[$([N-;X2]S(=O)=O)]", "N"),
        # Enamines
        ("[$([N-;X2][C,N]=C)]", "N"),
        # Tetrazoles
        ("[n-]", "[nH]"),
        # Sulfoxides
        ("[$([S-]=O)]", "S"),
        # Amides
        ("[$([N-]C=O)]", "N"),
    )
    reactions = [
        (AllChem.MolFromSmarts(x), AllChem.MolFromSmiles(y, False))
        for x, y in patts
    ]
    for reactant, product in reactions:
        while mol.HasSubstructMatch(reactant):
            rms = AllChem.ReplaceSubstructs(mol, reactant, product)
            mol = rms[0]
    Chem.SanitizeMol(mol)
    new_smiles = Chem.MolToSmiles(mol)
    return new_smiles


bio_rules_path = Path(__file__).parent / "JN3604IMT_rules.tsv"
cofactors_path = Path(__file__).parent / "all_cofactors.tsv"

bio_rules = pd.read_csv(bio_rules_path, sep="\t")
cofactors = pd.read_csv(cofactors_path, sep="\t")

cofactors_dict = dict()  # {"NAD_Cof": SMILES,}
cofactors_set = set()  # { SMILES,}
for idx, x in enumerate(cofactors["SMILES"]):
    cofactors_dict[cofactors["#ID"][idx]] = Chem.MolToSmiles(
        Chem.MolFromSmiles(x)
    )
    cofactors_set.add(Chem.MolToSmiles(Chem.MolFromSmiles(x)))

cofactors_clean_dict = dict()
cofactors_clean = set()
for idx, x in enumerate(cofactors["SMILES"]):
    cofactors_clean_dict[cofactors["#ID"][idx]] = clean_SMILES(x)
    cofactors_clean.add(clean_SMILES(x))


@typing.final
@dataclasses.dataclass(frozen=True, slots=True)
class SMILESCalculator(metadata.MolPropertyCalc[float]):
    # Calculate SMILES for molecules and save in network
    smiles_key: collections.abc.Hashable

    @property
    def key(self) -> collections.abc.Hashable:
        return self.smiles_key

    @property
    def meta_required(self) -> interfaces.MetaKeyPacket:
        return interfaces.MetaKeyPacket(molecule_keys={self.smiles_key})

    @property
    def resolver(self) -> metadata.MetaDataResolverFunc[float]:
        return metadata.TrivialMetaDataResolverFunc

    def __call__(
        self,
        data: interfaces.DataPacketE[interfaces.MolDatBase],
        prev_value: typing.Optional[float] = None,
    ) -> typing.Optional[float]:
        if prev_value is not None:
            return prev_value
        item = data.item
        if data.meta is not None and self.smiles_key in data.meta:
            return None
        return item.smiles


@typing.final
@dataclasses.dataclass(frozen=True)
class Reaction_Type_Filter(interfaces.RecipeFilter):
    # used for bio rxns, check reactants
    Allow_multi_reactants: bool

    def __call__(self, recipe: interfaces.RecipeExplicit) -> bool:
        if not self.Allow_multi_reactants:
            reas = set()
            for mol in recipe.reactants:
                SMILES = mol.meta["SMILES"]
                if clean_SMILES(SMILES) not in cofactors_clean:
                    reas.add(clean_SMILES(SMILES))
            if len(reas) != 1:
                # only allow A + cofactor or A + A
                return False
        reas_type = recipe.operator.meta["Reactants"].split(";")
        for idx, mol in enumerate(recipe.reactants):
            # check if reactant type is correct
            SMILES = mol.meta["SMILES"]
            if (
                reas_type[idx] == "Any"
                and clean_SMILES(SMILES) in cofactors_clean
            ):
                return False
            if (
                reas_type[idx] != "Any"
                and clean_SMILES(SMILES) != cofactors_clean_dict[reas_type[idx]]
            ):
                return False
        return True

    @property
    def meta_required(self) -> interfaces.MetaKeyPacket:
        return interfaces.MetaKeyPacket(
            operator_keys={"Reactants", "SMARTS"}, molecule_keys={"SMILES"}
        )


@typing.final
@dataclasses.dataclass(frozen=True)
class Product_Filter(metadata.ReactionFilterBase):
    # used for bio rxns, check if product type is correct.
    # Filter out radicals and fragments
    def __call__(self, recipe: interfaces.ReactionExplicit) -> bool:
        pros_type = recipe.operator.meta["Products"].split(";")
        for idx, mol in enumerate(recipe.products):
            if (
                pros_type[idx] != "Any"
                and clean_SMILES(mol.item.smiles)
                != cofactors_clean_dict[pros_type[idx]]
            ):
                return False
            if (
                Descriptors.NumRadicalElectrons(
                    Chem.MolFromSmiles(mol.item.smiles)
                )
                != 0
            ):
                return False
            if "." in mol.item.smiles:
                return False
        return True

    @property
    def meta_required(self) -> interfaces.MetaKeyPacket:
        return interfaces.MetaKeyPacket(operator_keys={"Products", "Name"})


@typing.final
@dataclasses.dataclass(frozen=True)
class Check_balance_filter(metadata.ReactionFilterBase):
    # used for bio rxns, check balance
    def __call__(self, recipe: interfaces.ReactionExplicit) -> bool:
        if True:
            charge_diff = 0
            reactants_dict = dict()
            products_dict = dict()
            pattern = r"([A-Z][a-z]*)(\d*)"
            for mol in recipe.reactants:
                charge_diff += GetFormalCharge(mol.item.rdkitmol)
                smiles = CalcMolFormula(mol.item.rdkitmol)
                matches = re.findall(pattern, smiles)
                for match in matches:
                    element, count = match
                    count = int(count) if count else 1
                    reactants_dict[element] = (
                        reactants_dict.get(element, 0) + count
                    )
            for mol in recipe.products:
                charge_diff -= GetFormalCharge(mol.item.rdkitmol)
                smiles = CalcMolFormula(mol.item.rdkitmol)
                matches = re.findall(pattern, smiles)
                for match in matches:
                    element, count = match
                    count = int(count) if count else 1
                    products_dict[element] = (
                        products_dict.get(element, 0) + count
                    )
            if charge_diff and "H" in reactants_dict:
                reactants_dict["H"] = reactants_dict["H"] - charge_diff
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
@dataclasses.dataclass(frozen=True, slots=True)
class Chem_Rxn_dH_Calculator(metadata.RxnPropertyCalc[float]):
    dH_key: collections.abc.Hashable
    direction: str
    rxn_dH: typing.Callable[[str], float]

    # rxn_dH: user function for dH of reaction
    # input: {reactants:(SMILES,), products:(SMILES,)} output: float
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
        if self.rxn_dH is None:
            return "No_Thermo"
        reas = list()
        pros = list()
        for mol in data.reactants:
            reas.append(mol.item.smiles)
        for mol in data.products:
            pros.append(mol.item.smiles)
        dH = self.rxn_dH({"reactants": reas, "products": pros})
        if dH is None:
            return float("nan")
        if data.operator.meta["enthalpy_correction"] is not None:
            dH = dH + data.operator.meta["enthalpy_correction"]
        # if self.direction == "forward":
        #     dH = dH
        if self.direction == "retro":
            dH = -dH
        return round(dH, 4)

    @property
    def meta_required(self) -> interfaces.MetaKeyPacket:
        return interfaces.MetaKeyPacket(
            operator_keys={
                "enthalpy_correction",
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
        if dH < self.max_dH:
            return True
        return False

    @property
    def meta_required(self) -> interfaces.MetaKeyPacket:
        return interfaces.MetaKeyPacket()


@typing.final
@dataclasses.dataclass(frozen=True, slots=True)
class Max_Atoms_Filter(metadata.ReactionFilterBase):
    # similar to engine.filter.reaction.max_atoms
    # take a dict, check multiple atoms {6: 12, 8: 2}
    # not checking cofactors
    max_atoms_dict: dict
    cofactors: set

    def __call__(self, recipe: interfaces.ReactionExplicit) -> bool:
        if self.max_atoms_dict is None:
            return True
        for atom_number, max_atoms in self.max_atoms_dict.items():
            for mol in recipe.products:
                rdkit_mol = mol.item.rdkitmol
                Chem.rdmolops.RemoveStereochemistry(rdkit_mol)
                if Chem.MolToSmiles(rdkit_mol) not in self.cofactors:
                    atom_query = rdqueries.AtomNumEqualsQueryAtom(atom_number)
                    atom_matches = mol.item.rdkitmol.GetAtomsMatchingQuery(
                        atom_query
                    )
                    if len(atom_matches) > max_atoms:
                        return False
        return True

    @property
    def meta_required(self) -> interfaces.MetaKeyPacket:
        return interfaces.MetaKeyPacket()


def generate_network(
    job_name="default_job",
    starters=False,
    gen=1,
    direction="forward",
    rxn_thermo_calculator=None,
    max_rxn_thermo_change=15,
    max_atoms=None,  # {"C": 20}
    allow_multiple_reactants=False,
    targets=None,  # string or list, set, etc.
    excluded_cofactors=("CARBONYL_CoF", "AMINO_CoF"),
):
    if not starters:
        raise Exception("At least one starter is needed to generate a network")

    print(f"Job Name: {job_name}")
    print("Job Started On:", datetime.now())
    start_time = time.time()

    engine = dn.create_engine()
    network = engine.new_network()

    for key in cofactors_dict:
        if excluded_cofactors is None or key not in excluded_cofactors:
            # add cofactors to network, they're like helpers in chem expansion
            network.add_mol(
                engine.mol.rdkit(cofactors_dict[key]),
                meta={
                    "SMILES": Chem.MolToSmiles(
                        Chem.MolFromSmiles(cofactors_dict[key])
                    )
                },
            )

    my_start_i = -1
    for smiles in starters:
        if my_start_i == -1:
            my_start_i = network.add_mol(
                engine.mol.rdkit(smiles),
                meta={"SMILES": Chem.MolToSmiles(Chem.MolFromSmiles(smiles))},
            )
        else:
            network.add_mol(
                engine.mol.rdkit(smiles),
                meta={"SMILES": Chem.MolToSmiles(Chem.MolFromSmiles(smiles))},
            )

    for idx, x in enumerate(bio_rules["SMARTS"]):
        reas_types = bio_rules["Reactants"][idx].split(";")
        pros_types = bio_rules["Products"][idx].split(";")

        if excluded_cofactors is None or (
            not set(excluded_cofactors) & set(reas_types)
            and not set(excluded_cofactors) & set(pros_types)
        ):
            network.add_op(
                engine.op.rdkit(
                    x,
                    kekulize=False,
                    drop_errors=True,
                ),
                meta={
                    "name": bio_rules["Name"][idx],
                    "Reactants": bio_rules["Reactants"][idx],
                    "Products": bio_rules["Products"][idx],
                    "Comments": bio_rules["Comments"][idx],
                    "SMARTS": x,
                    "reactants_stoi": (1,) * len(x.split(">>")[0].split(".")),
                    "products_stoi": (1,) * len(x.split(">>")[1].split(".")),
                    "enthalpy_correction": 0,
                    "Reaction_type": "Enzymatic",
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

    SMILES_Cal = SMILESCalculator("SMILES")
    Product_check = Product_Filter()
    coreactants_filter = engine.filter.bundle.coreactants(
        tuple(range(my_start_i))
    )

    reaction_plan = (
        Max_Atoms_Filter(max_atoms_dict_num, cofactors_clean)
        >> SMILES_Cal
        >> Product_check
        >> Check_balance_filter()
        >> Chem_Rxn_dH_Calculator("dH", direction, rxn_thermo_calculator)
        >> Rxn_dH_Filter(max_rxn_thermo_change, "dH")
    )

    Type_Filter = Reaction_Type_Filter(allow_multiple_reactants)
    ini_number = len(network.mols)

    strat.expand(
        num_iter=gen,
        reaction_plan=reaction_plan,
        bundle_filter=coreactants_filter,
        recipe_filter=Type_Filter,
        save_unreactive=False,
    )

    if targets is not None:
        print("Checking for targets")
        to_check = set()
        if isinstance(targets, str):
            to_check.add(clean_SMILES(targets))
        else:
            for i in targets:
                to_check.add(clean_SMILES(i))

        for mol in network.mols:
            if (
                network.reactivity[network.mols.i(mol.uid)] is True
                and clean_SMILES(mol.uid) in to_check
            ):
                print("Target found for", mol.uid)

    print("Number of generations:", gen)
    print("Number of operators loaded:", len(network.ops))
    print(
        "Number of molecules before expantion (including cofactors):",
        ini_number,
    )
    print(
        "Number of molecules after expantion (including cofactors):",
        len(network.mols),
    )
    print("Number of reactions:", len(network.rxns))

    end_time = time.time()
    elapsed_time = (end_time - start_time) / 60
    print("time used:", "{:.2f}".format(elapsed_time), "minutes")

    network.save_to_file(f"{job_name}_saved_network")

    return network
