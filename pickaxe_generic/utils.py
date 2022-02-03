"""
Contains classes which define and implement utility functions.

Classes:

    RxnTracker
      RxnTrackerSingle*
      RxnTrackerDepthFirst*
"""

from abc import ABC, abstractmethod
from itertools import chain, product as iterproduct
from typing import Collection, Dict, Generator, Iterable, List, Optional, Sequence, Set

from containers import ObjectLibrary
from datatypes import Identifier, RxnDatBase


class RxnTracker(ABC):
    """
    Interface representing an object which analyzes rxn network connections.

    Classes implementing this interface are able to create retrosynthetic trees
    based on a precalculated reaction network tree.

    Parameters
    ----------
    target : Identifier
        Unique ID of target molecule.
    reagent_table : Sequence[Identifier] (default: tuple())
        Contains unique IDs of reagents which do not need to be synthesized.
    fail_on_unknown_reagent : bool (default: False)
        If True, do not return paths which require reagents not in
        reagent_table.
    """

    @abstractmethod
    def getParentChains(
        self,
        target: Identifier,
        reagent_table: Sequence[Identifier] = tuple(),
        fail_on_unknown_reagent: bool = False,
    ) -> Iterable[Iterable[Iterable[RxnDatBase]]]:
        """
        Gets parent chains for a particular target molecule.

        Parameters
        ----------
        target : Identifier
            Unique id of target molecule.
        reagent_table : Sequence[Identifier]
            Sequence of reagents which are considered "basic" and which the tree
            search will consider leaf nodes.
        fail_on_unknown_reagent : bool
            If tree requires unlisted reagents, do not return.
        """


class RxnTrackerSingle(RxnTracker):
    """Implements RxnTracker interface; only compatible with reactions
    involving a single reactant and product.  DEVELOPMENT ONLY"""

    _mol_lookup: Dict[Identifier, Identifier]

    def __init__(self, rxn_lib: ObjectLibrary[RxnDatBase]) -> None:
        self._mol_lookup = {}
        self._rxn_lib = rxn_lib
        for rxnid in rxn_lib.ids():
            product_mol = sorted(rxn_lib[rxnid].products)[0]
            if product_mol not in self._mol_lookup:
                self._mol_lookup[product_mol] = []
            self._mol_lookup[product_mol].append(rxnid)

    def _getchains(
        self, cur_mol: Identifier, cur_mols: Optional[Set[Identifier]] = None
    ):
        if cur_mols is None:
            cur_mols = {cur_mol}
        noReactions = True
        if cur_mol not in self._mol_lookup:
            yield list()
        else:
            for rxnid in self._mol_lookup[cur_mol]:
                reactant = sorted(self._rxn_lib[rxnid].reactants)[0]
                if reactant in cur_mols:
                    continue
                noReactions = False
                for rxnpath in self._getchains(reactant, cur_mols.union({reactant})):
                    rxnpath.append(rxnid)
                    yield rxnpath
            if noReactions:
                yield list()

    def getParentChains(
        self,
        target: Identifier,
        reagent_table: Sequence[Identifier] = None,
        fail_on_unknown_reagent: bool = None,
    ) -> Iterable[Iterable[Iterable[RxnDatBase]]]:
        if reagent_table is not None or fail_on_unknown_reagent is not None:
            raise NotImplementedError("Arguments besides target not supported.")
        return ([path] for path in self._getchains(target))


class RxnTrackerDepthFirst(RxnTracker):
    """Implements RxnTracker interface; stores lookups as a hash table within
    the object.  Will eventually deprecate this functionality when the
    ObjectLibrary interface is updated to include native search functionality.
    """

    _mol_lookup: Dict[Identifier, Identifier]
    _rxn_lib: ObjectLibrary[RxnDatBase]

    def __init__(self, rxn_lib: ObjectLibrary[RxnDatBase]) -> None:
        self._mol_lookup = {}
        self._rxn_lib = rxn_lib
        for rxnid in rxn_lib.ids():
            for product_mol in rxn_lib[rxnid].products:
                if product_mol not in self._mol_lookup:
                    self._mol_lookup[product_mol] = []
                self._mol_lookup[product_mol].append(rxnid)

    def _getchains(
        self,
        cur_gen_mols: Collection[Identifier],
        prev_gens_mols: Set[Iterable[Identifier]] = None,
        prev_gens_rxns: Set[Iterable[Identifier]] = None,
        reagent_table: Optional[Iterable[Identifier]] = None,
        fail_on_unknown_reagent: bool = False,
    ) -> Generator[List[Set[Identifier]], None, None]:
        if len(cur_gen_mols) == 0:
            yield []
            return
        if prev_gens_mols is None:
            prev_gens_mols = set()
        if prev_gens_rxns is None:
            prev_gens_rxns = set()
        if reagent_table is None:
            reagent_table = set()
        rxnsets: List[List[RxnDatBase]] = []
        for mol in cur_gen_mols:
            if mol in reagent_table:
                continue
            elif mol not in self._mol_lookup:
                rxnsets.append([])
                continue
            newrxnset: List[RxnDatBase] = [
                rxn
                for rxn in self._mol_lookup[mol]
                if rxn not in prev_gens_rxns
                and all(
                    mol not in prev_gens_mols for mol in self._rxn_lib[rxn].reactants
                )
            ]
            rxnsets.append(newrxnset)
        if not fail_on_unknown_reagent:
            rxnsets = [rxnset for rxnset in rxnsets if len(rxnset) > 0]
        if len(rxnsets) == 0:
            yield []
            return
        tested_combos = set()
        for rxncombo in iterproduct(*rxnsets):
            rxncombo = frozenset(rxncombo)
            if rxncombo in tested_combos:
                continue
            else:
                tested_combos.add(frozenset(rxncombo))
            required_reagents = set(
                mol
                for mol in chain(*(self._rxn_lib[rxn].reactants for rxn in rxncombo))
                if mol not in reagent_table
            )
            if len(required_reagents) == 0:
                yield [rxncombo]
                continue
            for path in self._getchains(
                required_reagents,
                prev_gens_mols.union(cur_gen_mols),
                prev_gens_rxns.union(rxncombo),
                reagent_table,
                fail_on_unknown_reagent,
            ):
                path.append(rxncombo)
                yield path

    def getParentChains(
        self,
        target: Identifier,
        reagent_table: Sequence[Identifier] = tuple(),
        fail_on_unknown_reagent: bool = False,
    ) -> Generator[List[Set[RxnDatBase]], None, None]:
        if fail_on_unknown_reagent and not reagent_table:
            raise ValueError(
                "reagent table must be specified if fail_on_unknown_reagent is True"
            )
        return (
            path
            for path in self._getchains(
                [target],
                reagent_table=reagent_table,
                fail_on_unknown_reagent=fail_on_unknown_reagent,
            )
        )
