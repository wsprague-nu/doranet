"""
Contains classes which implement various filter components.

Classes:

    ReactionFilter (deprecated)
      AlwaysTrueFilter*
      ChainFilter*
      LessThanNElementTypeFilter*
      TanimotoSimilarityFilter*
    UIDPreFilter (deprecated)
      AlwaysTrueUIDPreFilter*
      CoreactantUIDPreFilter*
    MolFilter
      MolFilterMetaVal*
      MolFilterMetaExist*
    RecipeFilter
      CoreactantFilter*
      

"""

import collections.abc
import dataclasses
import typing

import rdkit
import rdkit.Chem

from . import interfaces


class AlwaysTrueFilter(interfaces.ReactionFilter):
    def __call__(self, operator, reactants, products):
        return True


class ChainFilter(interfaces.ReactionFilter):
    def __init__(
        self, filters: collections.abc.Iterable[interfaces.ReactionFilter]
    ):
        self._filters = filters

    def __call__(self, operator, reactants, products):
        return all(
            (filter(operator, reactants, products) for filter in self._filters)
        )


class LessThanNElementTypeFilter(interfaces.ReactionFilter):
    def __init__(self, n: int, proton_number: int):
        self._n = n
        self._p = proton_number
        self._q = rdkit.Chem.rdqueries.AtomNumEqualsQueryAtom(proton_number)

    def __call__(self, operator, reactants, products):
        for mol in products:
            if isinstance(mol, interfaces.MolDatRDKit):
                if len(mol.rdkitmol.GetAtomsMatchingQuery(self._q)) >= self._n:
                    return False
        return True

    def __getstate__(self):
        return (self._n, self._p)

    def __setstate__(self, arg) -> None:
        self._n = arg[0]
        self._p = arg[1]
        self._q = rdkit.Chem.rdqueries.AtomNumEqualsQueryAtom(self._p)


class TanimotoSimilarityFilter(interfaces.ReactionFilter):
    def __init__(self, n: float, smi: str):
        self._n = n
        self._s = smi
        self._tmol = rdkit.Chem.MolFromSmiles(self._s)
        self._tfp = rdkit.Chem.RDKFingerprint(self._tmol)

    def __call__(self, operator, reactants, products):
        for mol in products:
            if isinstance(mol, interfaces.MolDatRDKit):
                mol_fp = rdkit.Chem.RDKFingerprint(mol.rdkitmol)
                similarity = rdkit.DataStructs.TanimotoSimilarity(
                    mol_fp, self._tfp
                )

                if similarity > self._n:
                    return True
        return False


class AlwaysTrueUIDPreFilter(interfaces.UIDPreFilter):
    def __call__(
        self,
        operator: interfaces.Identifier,
        reactants: collections.abc.Sequence[interfaces.Identifier],
    ) -> bool:
        return True


class CoreactantUIDPreFilter(interfaces.UIDPreFilter):
    def __init__(
        self, coreactants: collections.abc.Collection[interfaces.Identifier]
    ):
        self._coreactants = frozenset(coreactants)

    def __call__(
        self,
        operator: interfaces.Identifier,
        reactants: collections.abc.Sequence[interfaces.Identifier],
    ) -> bool:
        for uid in reactants:
            if uid not in self._coreactants:
                return True
        return False


@dataclasses.dataclass(frozen=True)
class MolFilterMetaVal(interfaces.MolFilter):
    __slots__ = ("_key", "_val")
    _key: collections.abc.Hashable
    _val: typing.Any

    def __call__(
        self,
        mol: interfaces.MolDatBase,
        meta: typing.Optional[
            collections.abc.Mapping[collections.abc.Hashable, typing.Any]
        ] = None,
    ) -> bool:
        if meta is None:
            return False
        if self._key not in meta:
            return False
        return meta[self._key] == self._val

    @property
    def meta_required(self) -> interfaces.MetaKeyPacket:
        return interfaces.MetaKeyPacket(frozenset(), frozenset((self._key,)))


@dataclasses.dataclass(frozen=True)
class MolFilterMetaExist(interfaces.MolFilter):
    __slots__ = ("_key", "_val")
    _key: collections.abc.Hashable

    def __call__(
        self,
        mol: interfaces.MolDatBase,
        meta: typing.Optional[
            collections.abc.Mapping[collections.abc.Hashable, typing.Any]
        ] = None,
    ) -> bool:
        if meta is None:
            return False
        if self._key not in meta:
            return False
        return True

    @property
    def meta_required(self) -> interfaces.MetaKeyPacket:
        return interfaces.MetaKeyPacket(frozenset(), frozenset((self._key,)))


@dataclasses.dataclass
class CoreactantFilter(interfaces.RecipeFilter):
    __slots__ = "_coreactants"
    _coreactants: set[interfaces.MolIndex]

    def __init__(
        self, coreactants: typing.Iterable[interfaces.MolIndex]
    ) -> None:
        self._coreactants = set(coreactants)

    def __call__(self, recipe: interfaces.RecipeExplicit) -> bool:
        if all(mol.i in self._coreactants for mol in recipe.reactants):
            return False
        return True


def ReplaceNewValue(
    key: collections.abc.Hashable, old_value: typing.Any, new_value: typing.Any
) -> bool:
    if old_value != new_value:
        return True
    return False


class ReplaceBlacklist:
    def __init__(
        self,
        blacklist_keys: collections.abc.Collection[collections.abc.Hashable],
    ) -> None:
        self._blacklist_keys = frozenset(blacklist_keys)

    def __call__(
        self,
        key: collections.abc.Hashable,
        old_value: typing.Any,
        new_value: typing.Any,
    ) -> bool:
        if key in self._blacklist_keys:
            if old_value is None:
                return True
            elif old_value is False and new_value is True:
                return True
            return False
        return ReplaceNewValue(key, old_value, new_value)


"""
class DefaultMetaDataUpdate:
    def __call__(
        self,
        unit: interfaces.ReactionExplicit,
        network: interfaces.ChemNetwork,
    ) -> Generator[
        tuple[
            Optional[tuple[interfaces.MolIndex, Hashable]],
            Optional[tuple[interfaces.OpIndex, Hashable]],
        ],
        None,
        None,
    ]:
        if unit.operator_meta is not None:
            opIndex = network.ops.i(unit.operator.uid)
            for key, value in unit.operator_meta.items():
                if network.op_meta(opIndex, key) != value:
                    network.op_meta(opIndex, key, value)
                    yield (None, (opIndex, key))
        if unit.reactants_meta is not None:
            for index, reactant in enumerate(unit.reactants):
                molIndex = network.mols.i(reactant.uid)
                for key, value in unit.reactants_meta[index]:
                    if network.mol_meta(molIndex, key) != value:
                        network.mol_meta(molIndex, key, value)
                        yield ((molIndex, key), None)
        if unit.products_meta is not None:
            for index, product in enumerate(unit.products):
                molIndex = network.mols.i(product.uid)
                for key, value in unit.products_meta[index]:
                    if network.mol_meta(molIndex, key) != value:
                        network.mol_meta(molIndex, key, value)
                        yield ((molIndex, key), None)


# def __call__(self, )
"""
