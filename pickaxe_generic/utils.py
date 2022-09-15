"""
Contains classes which define and implement utility functions.

Classes:

    RxnTracker
      RxnTrackerSingle*
      RxnTrackerDepthFirst*
"""

from abc import ABC, abstractmethod
from itertools import chain
from itertools import product as iterproduct
from re import X
from types import NoneType
from typing import (
    Callable,
    Collection,
    Dict,
    Generator,
    Generic,
    Iterable,
    Iterator,
    List,
    Optional,
    Sequence,
    Set,
    TypeVar,
    overload,
)

from . import interfaces


class RxnTrackerSingle(interfaces.RxnTracker):
    """Implements RxnTracker interface; only compatible with reactions
    involving a single reactant and product.  DEVELOPMENT ONLY"""

    _mol_lookup: Dict[interfaces.Identifier, list[interfaces.Identifier]]

    def __init__(
        self, rxn_lib: interfaces.ObjectLibrary[interfaces.RxnDatBase]
    ) -> None:
        self._mol_lookup = {}
        self._rxn_lib = rxn_lib
        for rxnid in rxn_lib.ids():
            product_mol = sorted(rxn_lib[rxnid].products)[0]
            if product_mol not in self._mol_lookup:
                self._mol_lookup[product_mol] = []
            self._mol_lookup[product_mol].append(rxnid)

    def _getchains(
        self,
        cur_mol: interfaces.Identifier,
        cur_mols: Optional[Set[interfaces.Identifier]] = None,
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
                for rxnpath in self._getchains(
                    reactant, cur_mols.union({reactant})
                ):
                    rxnpath.append(rxnid)
                    yield rxnpath
            if noReactions:
                yield list()

    def getParentChains(
        self,
        target: interfaces.Identifier,
        reagent_table: Sequence[interfaces.Identifier] = None,
        fail_on_unknown_reagent: bool = None,
        max_depth: Optional[int] = None,
    ) -> Iterable[Iterable[Iterable[interfaces.RxnDatBase]]]:
        if reagent_table is not None or fail_on_unknown_reagent is not None:
            raise NotImplementedError("Arguments besides target not supported.")
        return ([path] for path in self._getchains(target))


class RxnTrackerDepthFirst(interfaces.RxnTracker):
    """Implements RxnTracker interface; stores lookups as a hash table within
    the object.  Will eventually deprecate this functionality when the
    ObjectLibrary interface is updated to include native search functionality.
    """

    _mol_lookup: Dict[interfaces.Identifier, list[interfaces.Identifier]]
    _rxn_lib: interfaces.ObjectLibrary[interfaces.RxnDatBase]

    def __init__(
        self, rxn_lib: interfaces.ObjectLibrary[interfaces.RxnDatBase]
    ) -> None:
        self._mol_lookup = {}
        self._rxn_lib = rxn_lib
        for rxnid in rxn_lib.ids():
            for product_mol in rxn_lib[rxnid].products:
                if product_mol not in self._mol_lookup:
                    self._mol_lookup[product_mol] = []
                self._mol_lookup[product_mol].append(rxnid)

    def _getchains(
        self,
        cur_gen_mols: Collection[interfaces.Identifier],
        prev_gens_mols: Set[interfaces.Identifier] = None,
        prev_gens_rxns: Set[interfaces.Identifier] = None,
        reagent_table: Optional[Iterable[interfaces.Identifier]] = None,
        fail_on_unknown_reagent: bool = False,
        depth: Optional[int] = None,
    ) -> Generator[list[frozenset[interfaces.Identifier]], None, None]:
        if depth is not None:
            if depth <= 0:
                return
            new_depth = depth - 1
        else:
            new_depth = None
        if len(cur_gen_mols) == 0:
            yield []
            return
        if prev_gens_mols is None:
            prev_gens_mols = set()
        if prev_gens_rxns is None:
            prev_gens_rxns = set()
        if reagent_table is None:
            reagent_table = set()
        rxnsets: List[List[interfaces.Identifier]] = []
        for mol in cur_gen_mols:
            if mol in reagent_table:
                continue
            elif mol not in self._mol_lookup:
                rxnsets.append([])
                continue
            newrxnset: List[interfaces.Identifier] = [
                rxn
                for rxn in self._mol_lookup[mol]
                if rxn not in prev_gens_rxns
                and all(
                    mol not in prev_gens_mols
                    for mol in self._rxn_lib[rxn].reactants
                )
            ]
            rxnsets.append(newrxnset)
        if not fail_on_unknown_reagent:
            rxnsets = [rxnset for rxnset in rxnsets if len(rxnset) > 0]
            if len(rxnsets) == 0:  # these two lines
                return  # close the loophole
        if len(rxnsets) == 0:
            yield []
            return
        tested_combos = set()
        for rxntuple in iterproduct(*rxnsets):
            rxncombo: frozenset[interfaces.Identifier] = frozenset(rxntuple)
            if rxncombo in tested_combos:
                continue
            else:
                tested_combos.add(frozenset(rxncombo))
            required_reagents = set(
                mol
                for mol in chain(
                    *(self._rxn_lib[rxn].reactants for rxn in rxncombo)
                )
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
                new_depth,
            ):
                path.append(rxncombo)
                yield path

    def getParentChains(
        self,
        target: interfaces.Identifier,
        reagent_table: Sequence[interfaces.Identifier] = tuple(),
        fail_on_unknown_reagent: bool = False,
        max_depth: Optional[int] = None,
    ) -> Generator[list[frozenset[interfaces.Identifier]], None, None]:
        if fail_on_unknown_reagent and not reagent_table:
            raise ValueError(
                "reagent table must be specified if fail_on_unknown_reagent is "
                "True"
            )
        return (
            path
            for path in self._getchains(
                [target],
                reagent_table=reagent_table,
                fail_on_unknown_reagent=fail_on_unknown_reagent,
                depth=max_depth,
            )
        )


T = TypeVar("T")


def logreduce(function: Callable[[T, T], T], iterable: Iterable[T]) -> T:
    """
    logreduce(function, iterable[, initial]) -> value

    Apply a function of two arguments (satisfying the associative property)
    cumulatively to the items of a sequence or iterable, from left to right, so
    as to reduce the iterable to a single value.  For example, reduce(lambda x,
    y: x+y, [1, 2, 3, 4, 5]) calculates (((1+2)+(3+4))+5).  If initial is
    present, it is placed before the items of the iterable in the calculation,
    and serves as a default when the iterable is empty.

    Memory: maximum of log2(n) objects of iterable stored vs 2 for reduce()
    Speed: log2(n) calls to function() vs n for reduce()
    """
    i = iter(iterable)
    try:
        r_val, stop = _logreduce(function, i, 0)
    except StopIteration:
        raise TypeError("logreduce() of empty iterable with no initial value")
    n = 0
    try:
        while not stop:
            n += 1
            new_val, stop = _logreduce(function, i, n - 1)
            r_val = function(r_val, new_val)
    except StopIteration:
        ...
    return r_val


def _logreduce(
    f: Callable[[T, T], T], i: Iterator[T], n: int
) -> tuple[T, bool]:
    if n == 0:
        x = next(i)
        try:
            y = next(i)
        except StopIteration:
            return x, True
        return f(x, y), False
    x, stop = _logreduce(f, i, n - 1)
    if stop:
        return x, True
    try:
        y, stop = _logreduce(f, i, n - 1)
    except StopIteration:
        return x, True
    return f(x, y), False
