"""Synthetic network generation code."""

import collections.abc
import dataclasses
import re
import time
import typing
from datetime import datetime

from rdkit import Chem
from rdkit.Chem import rdqueries, rdRascalMCES
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
class Regioselectivity_filter(metadata.ReactionFilterBase):
    __slots__ = ("regio_user", "direction")
    regio_user: typing.Optional[typing.Collection]
    direction: str

    def __call__(self, recipe: interfaces.ReactionExplicit) -> bool:
        if self.regio_user is None:  # if user doesn't use it
            return True
        regio_rxn_tup = recipe.operator.meta["regioselectivity"]
        if regio_rxn_tup is None:  # if rxn doesn't have it
            return True
        Markovnikov, Anti_Markovnikov, Baeyer_Villiger, Zaitsev, Hofmann = (
            False,
            False,
            False,
            False,
            False,
        )
        if regio_rxn_tup[0] in self.regio_user:
            if regio_rxn_tup[0] == "Markovnikov":
                Markovnikov = True
            elif regio_rxn_tup[0] == "Anti-Markovnikov":
                Anti_Markovnikov = True
            elif regio_rxn_tup[0] == "Baeyer-Villiger":
                Baeyer_Villiger = True
            elif regio_rxn_tup[0] == "Zaitsev":
                Zaitsev = True
            elif regio_rxn_tup[0] == "Hofmann":
                Hofmann = True
        else:
            return True
        main_rea = recipe.reactants[regio_rxn_tup[1]].item
        if not isinstance(main_rea, interfaces.MolDatRDKit):
            raise NotImplementedError(
                f"""Filter only implemented for molecule type \
                    MolDatRDKit, not {type(main_rea)}"""
            )
        main_pro = recipe.products[regio_rxn_tup[2]].item
        if not isinstance(main_pro, interfaces.MolDatRDKit):
            raise NotImplementedError(
                f"""Filter only implemented for molecule type \
                    MolDatRDKit, not {type(main_pro)}"""
            )
        reactant = main_rea.rdkitmol
        product = main_pro.rdkitmol
        if self.direction == "retro":
            reactant, product = product, reactant
        if Zaitsev or Hofmann:
            reactant, product = product, reactant
        mol1 = reactant  # reactant
        num_atoms1 = mol1.GetNumAtoms()
        # add Hs so terminal double bonds like C=CC, and
        # small mols like CC=CC can be properly matched
        mol1 = Chem.AddHs(mol1)
        mol2 = product
        num_atoms2 = mol2.GetNumAtoms()
        mol2 = Chem.AddHs(mol2)
        # rdRascalMCES opts
        opts = rdRascalMCES.RascalOptions()
        opts.similarityThreshold = 0  # so it works for small mols
        opts.maxBondMatchPairs = 3000  # default 1k.
        opts.ringMatchesRingOnly = True
        # Find MCES
        results = rdRascalMCES.FindMCES(mol1, mol2, opts)
        try:
            res = results[0]
        except IndexError:
            return True  # if molecule too big to match, let it pass
        matching_bonds = res.bondMatches()
        matching_atoms = res.atomMatches()
        mol1_all_bond_indices = [bond.GetIdx() for bond in mol1.GetBonds()]
        mol1_matching_bonds = {t[0] for t in matching_bonds}
        mol1_rxn_bond_index_list = [
            item
            for item in mol1_all_bond_indices
            if item not in mol1_matching_bonds
        ]
        if len(mol1_rxn_bond_index_list) != 1:
            raise ValueError(
                "Number of transformed bonds in the reactant is not 1"
            )
        mol1_rxn_bond_index = mol1_rxn_bond_index_list[0]
        mol1_rxn_bond_obj = mol1.GetBondWithIdx(mol1_rxn_bond_index)
        mol1_atom1_idx = mol1_rxn_bond_obj.GetBeginAtomIdx()
        mol1_atom1 = mol1_rxn_bond_obj.GetBeginAtom()
        mol1_atom2_idx = mol1_rxn_bond_obj.GetEndAtomIdx()
        mol1_atom2 = mol1_rxn_bond_obj.GetEndAtom()

        def find_ketone_carbon(mol, atom1, atom2, index1, index2):
            oxygen_atom_num = 8
            CO_atom = []
            for neighbor in atom1.GetNeighbors():
                # Check if neighbor is Oxygen
                if neighbor.GetAtomicNum() == oxygen_atom_num:
                    bond = mol.GetBondBetweenAtoms(index1, neighbor.GetIdx())
                    if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                        CO_atom.append((atom1, index2))
            for neighbor in atom2.GetNeighbors():
                # Check if neighbor is Oxygen
                if neighbor.GetAtomicNum() == oxygen_atom_num:
                    bond = mol.GetBondBetweenAtoms(index2, neighbor.GetIdx())
                    if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                        CO_atom.append((atom2, index1))
            if len(CO_atom) == 1:
                return CO_atom[0]
            diketones_case = 2
            return len(CO_atom) == diketones_case
            # elif len(CO_atom) == 2: # 1,2-diketones
            #     return True
            # else:    # if no ketone carbon is found
            #     return False

        if Markovnikov or Anti_Markovnikov:
            if num_atoms2 - num_atoms1 == 1:
                mol1_atom1_non_h_neighbors = [
                    nbr.GetIdx()
                    for nbr in mol1_atom1.GetNeighbors()
                    if nbr.GetAtomicNum() != 1
                ]
                num_mol1_atom1_non_h_neighbors = len(mol1_atom1_non_h_neighbors)
                mol1_atom2_non_h_neighbors = [
                    nbr.GetIdx()
                    for nbr in mol1_atom2.GetNeighbors()
                    if nbr.GetAtomicNum() != 1
                ]
                num_mol1_atom2_non_h_neighbors = len(mol1_atom2_non_h_neighbors)

                mol2_atom1_idx = next(
                    second
                    for first, second in matching_atoms
                    if first == mol1_atom1_idx
                )
                mol2_atom2_idx = next(
                    second
                    for first, second in matching_atoms
                    if first == mol1_atom2_idx
                )
                mol2_atom1 = mol2.GetAtomWithIdx(mol2_atom1_idx)
                mol2_atom2 = mol2.GetAtomWithIdx(mol2_atom2_idx)
                mol2_atom1_non_h_neighbors = [
                    nbr.GetIdx()
                    for nbr in mol2_atom1.GetNeighbors()
                    if nbr.GetAtomicNum() != 1
                ]
                num_mol2_atom1_non_h_neighbors = len(mol2_atom1_non_h_neighbors)
                mol2_atom2_non_h_neighbors = [
                    nbr.GetIdx()
                    for nbr in mol2_atom2.GetNeighbors()
                    if nbr.GetAtomicNum() != 1
                ]
                num_mol2_atom2_non_h_neighbors = len(mol2_atom2_non_h_neighbors)
                if (
                    num_mol1_atom1_non_h_neighbors
                    == num_mol1_atom2_non_h_neighbors
                ):
                    return True
                mo1_more_sub = (
                    "atom1"
                    if num_mol1_atom1_non_h_neighbors
                    > num_mol1_atom2_non_h_neighbors
                    else "atom2"
                )
                # check which atom in mol2 gets new group
                if (
                    num_mol2_atom1_non_h_neighbors
                    > num_mol1_atom1_non_h_neighbors
                ):
                    atom_increased = "atom1"
                elif (
                    num_mol2_atom2_non_h_neighbors
                    > num_mol1_atom2_non_h_neighbors
                ):
                    atom_increased = "atom2"
                else:
                    raise ValueError(
                        "Something wrong with the addition reaction"
                    )
                if Markovnikov:
                    return mo1_more_sub == atom_increased
                if Anti_Markovnikov:
                    return mo1_more_sub != atom_increased
            elif num_atoms2 - num_atoms1 > 1:  # addition of alcohols
                return True
            else:
                raise ValueError("Something wrong with the addition reaction")

        elif Baeyer_Villiger:
            carbon_atom_num = 6
            oxygen_atom_num = 8
            mol1_ketone = find_ketone_carbon(
                mol1, mol1_atom1, mol1_atom2, mol1_atom1_idx, mol1_atom2_idx
            )
            if type(mol1_ketone) is tuple:
                mol1_CO_C = mol1_ketone[0]
                mol1_rebond_C_idx = mol1_ketone[1]
                two_neighbors = list()
                for neighbor in mol1_CO_C.GetNeighbors():
                    # 2 neighboring C atoms
                    if neighbor.GetAtomicNum() == carbon_atom_num:
                        two_neighbors.append(neighbor)
                # redefine mol1, mol2 as the 2 neighboring atoms of C=O
                mol1_atom1, mol1_atom2 = two_neighbors
                mol1_atom1_idx = mol1_atom1.GetIdx()
                mol1_atom2_idx = mol1_atom2.GetIdx()
                mol1_unchanged_C_idx = (
                    mol1_atom1_idx
                    if mol1_atom2_idx == mol1_rebond_C_idx
                    else mol1_atom2_idx
                )
                # check if unchanged is a C=O group, if so let it pass
                mol1_unchanged_C = mol1.GetAtomWithIdx(mol1_unchanged_C_idx)
                for neighbor in mol1_unchanged_C.GetNeighbors():
                    # Check if neighbor is O
                    if neighbor.GetAtomicNum() == oxygen_atom_num:
                        bond = mol1.GetBondBetweenAtoms(
                            mol1_unchanged_C_idx, neighbor.GetIdx()
                        )
                        if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                            return True
                mol1_atom1_non_h_neighbors = [
                    nbr.GetIdx()
                    for nbr in mol1_atom1.GetNeighbors()
                    if nbr.GetAtomicNum() == carbon_atom_num
                ]
                num_mol1_atom1_non_h_neighbors = len(mol1_atom1_non_h_neighbors)
                mol1_atom2_non_h_neighbors = [
                    nbr.GetIdx()
                    for nbr in mol1_atom2.GetNeighbors()
                    if nbr.GetAtomicNum() == carbon_atom_num
                ]
                num_mol1_atom2_non_h_neighbors = len(mol1_atom2_non_h_neighbors)
                if (
                    num_mol1_atom1_non_h_neighbors
                    == num_mol1_atom2_non_h_neighbors
                ):
                    return True
                if (
                    num_mol1_atom1_non_h_neighbors
                    > num_mol1_atom2_non_h_neighbors
                ):
                    mo1_more_sub_idx = mol1_atom1_idx
                else:
                    mo1_more_sub_idx = mol1_atom2_idx
                return mo1_more_sub_idx == mol1_rebond_C_idx
            elif mol1_ketone is True:  # 1,2-diketones, let it pass
                return True
            else:  # no ketone
                raise ValueError(
                    "No carbonyl group reacted in the Baeyer Villiger rxn"
                )

        elif Zaitsev or Hofmann:
            carbon_atom_num = 6
            if num_atoms2 - num_atoms1 == 1:
                mol2_atom1_idx = next(
                    second
                    for first, second in matching_atoms
                    if first == mol1_atom1_idx
                )
                mol2_atom2_idx = next(
                    second
                    for first, second in matching_atoms
                    if first == mol1_atom2_idx
                )
                mol2_all_atom_indices = [
                    atom.GetIdx() for atom in mol2.GetAtoms()
                ]
                mol2_matching_atoms = {t[1] for t in matching_atoms}
                mol2_leaving_group_index = next(
                    item
                    for item in mol2_all_atom_indices
                    if (
                        item not in mol2_matching_atoms
                        and mol2.GetAtomWithIdx(item).GetAtomicNum() != 1
                    )
                )
                mol2_base_C = next(
                    atom
                    for atom in mol2.GetAtomWithIdx(
                        mol2_leaving_group_index
                    ).GetNeighbors()
                    if atom.GetAtomicNum() != 1
                )
                mol2_loseH_C_idx = (
                    mol2_atom1_idx
                    if mol2_base_C.GetIdx() == mol2_atom2_idx
                    else mol2_atom2_idx
                )
                base_neighbors = [  # find all neighbor C who has H
                    atom
                    for atom in mol2_base_C.GetNeighbors()
                    if (
                        atom.GetAtomicNum() == carbon_atom_num
                        and any(
                            neighbor.GetAtomicNum() == 1
                            for neighbor in atom.GetNeighbors()
                        )
                    )
                ]
                num_neighbor_list = list()
                for atom in base_neighbors:
                    num_neighbor = len(
                        [
                            nbr.GetIdx()
                            for nbr in atom.GetNeighbors()
                            if nbr.GetAtomicNum() == carbon_atom_num
                        ]
                    )
                    num_neighbor_list.append(num_neighbor)
                if Zaitsev:
                    max_value = max(num_neighbor_list)
                    max_neighbor_atom_indices = {
                        base_neighbors[i].GetIdx()
                        for i, x in enumerate(num_neighbor_list)
                        if x == max_value
                    }
                    return mol2_loseH_C_idx in max_neighbor_atom_indices
                else:  # Hofmann
                    min_value = min(num_neighbor_list)
                    min_neighbor_atom_indices = {
                        base_neighbors[i].GetIdx()
                        for i, x in enumerate(num_neighbor_list)
                        if x == min_value
                    }
                    return mol2_loseH_C_idx in min_neighbor_atom_indices
            else:
                raise ValueError(
                    "The leaving group might be too big for the Zaitsev filter"
                )

    @property
    def meta_required(self) -> interfaces.MetaKeyPacket:
        return interfaces.MetaKeyPacket(
            operator_keys={
                "Reaction_direction",
                "regioselectivity",
            }
        )


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


def validate_regioselectivity(regioselectivity):
    allowed_values = {
        "Markovnikov",
        "Anti-Markovnikov",
        "Zaitsev",
        "Hofmann",
        "Baeyer-Villiger",
    }
    if regioselectivity is None:
        return None
    elif regioselectivity == "all":
        return {
            "Markovnikov",
            "Anti-Markovnikov",
            "Zaitsev",
            "Hofmann",
            "Baeyer-Villiger",
        }
    elif isinstance(regioselectivity, (set, list, tuple)):
        # Check if all elements in the collection are valid
        if not all(item in allowed_values for item in regioselectivity):
            raise ValueError(
                "Invalid regioselectivity values. Allowed values:"
                "Markovnikov, Anti-Markovnikov, Zaitsev, Hofmann,"
                "Baeyer-Villiger"
            )
    else:
        raise ValueError(
            "Invalid regioselectivity values. Allowed values: None, all, or a"
            "collection containing Markovnikov, Anti-Markovnikov, Zaitsev,"
            "Hofmann, Baeyer-Villiger"
        )
    return regioselectivity  # If all checks pass, input is valid


def generate_network(
    job_name="default_job",
    starters=False,
    helpers=False,
    gen=1,
    direction="forward",
    molecule_thermo_calculator=None,
    max_rxn_thermo_change=15,
    max_atoms=None,  # {"C": 20}
    allow_multiple_reactants="default",  # forward allowed, retro no
    targets=None,  # string or list, set, etc.
    regioselectivity=None,  # None, "all", or use a selection:
    # ("Markovnikov","Anti-Markovnikov","Zaitsev","Hofmann","Baeyer-Villiger")
):
    if not starters:
        raise Exception("At least one starter is needed to generate a network")

    user_regio = validate_regioselectivity(regioselectivity)

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
                    "regioselectivity": smarts.regioselectivity,
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
                    "regioselectivity": smarts.regioselectivity,
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
            >> Regioselectivity_filter(user_regio, "forward")
            >> Chem_Rxn_dH_Calculator(
                "dH", "forward", molecule_thermo_calculator
            )
            >> Rxn_dH_Filter(max_rxn_thermo_change, "dH")
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
            >> Regioselectivity_filter(user_regio, "retro")
            >> Chem_Rxn_dH_Calculator("dH", "retro", molecule_thermo_calculator)
            >> Rxn_dH_Filter(max_rxn_thermo_change, "dH")
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
