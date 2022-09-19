"""
Contains classes which define and implement utility functions.

Classes:

    RxnTracker
      RxnTrackerSingle*
      RxnTrackerDepthFirst*
"""

import collections.abc
import itertools
import typing

from . import interfaces


class RxnTrackerDepthFirst(interfaces.RxnTracker):
    """Implements RxnTracker interface; stores lookups as a hash table within
    the object.  Will eventually deprecate this functionality when the
    ObjectLibrary interface is updated to include native search functionality.
    """

    _mol_lookup: dict[interfaces.Identifier, list[interfaces.Identifier]]
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
        cur_gen_mols: collections.abc.Collection[interfaces.Identifier],
        prev_gens_mols: set[interfaces.Identifier] = None,
        prev_gens_rxns: set[interfaces.Identifier] = None,
        reagent_table: typing.Optional[
            collections.abc.Iterable[interfaces.Identifier]
        ] = None,
        fail_on_unknown_reagent: bool = False,
        depth: typing.Optional[int] = None,
    ) -> collections.abc.Generator[
        list[frozenset[interfaces.Identifier]], None, None
    ]:
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
        rxnsets: list[list[interfaces.Identifier]] = []
        for mol in cur_gen_mols:
            if mol in reagent_table:
                continue
            elif mol not in self._mol_lookup:
                rxnsets.append([])
                continue
            newrxnset: list[interfaces.Identifier] = [
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
        for rxntuple in itertools.product(*rxnsets):
            rxncombo: frozenset[interfaces.Identifier] = frozenset(rxntuple)
            if rxncombo in tested_combos:
                continue
            else:
                tested_combos.add(frozenset(rxncombo))
            required_reagents = set(
                mol
                for mol in itertools.chain(
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
        reagent_table: collections.abc.Sequence[
            interfaces.Identifier
        ] = tuple(),
        fail_on_unknown_reagent: bool = False,
        max_depth: typing.Optional[int] = None,
    ) -> collections.abc.Generator[
        list[frozenset[interfaces.Identifier]], None, None
    ]:
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


T = typing.TypeVar("T")


def logreduce(
    function: collections.abc.Callable[[T, T], T],
    iterable: collections.abc.Iterable[T],
) -> T:
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
    f: collections.abc.Callable[[T, T], T],
    i: collections.abc.Iterator[T],
    n: int,
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
