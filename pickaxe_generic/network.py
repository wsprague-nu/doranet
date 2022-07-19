from abc import ABC, abstractmethod
from collections.abc import Iterator, Mapping, Sequence
from dataclasses import dataclass
from typing import Generic, NewType, Optional, TypeVar, Union, overload

from pickaxe_generic.containers import DataUnitGen
from pickaxe_generic.datatypes import (
    DataUnit,
    Identifier,
    MolDatBase,
    OpDatBase,
)

_MolIndex = NewType("_MolIndex", int)
_OpIndex = NewType("_OpIndex", int)
_RxnIndex = NewType("_RxnIndex", int)
_I_T = TypeVar("_I_T", bound=int)
_ID_T = TypeVar("_ID_T", bound=Identifier)


@dataclass(frozen=True, order=True, slots=True)
class Reaction:
    operator: _OpIndex
    reactants: tuple[_MolIndex, ...]
    products: tuple[_MolIndex, ...]


@dataclass(frozen=True)
class __ValueQueryData(Generic[DataUnitGen, _I_T]):
    __slots__ = ("_list", "_map")
    _list: Sequence[DataUnitGen]
    _map: Mapping[Identifier, _I_T]

    def __contains__(self, item: Union[Identifier, DataUnitGen]) -> bool:
        if isinstance(item, DataUnit):
            item = item.uid
        return item in self._map.keys()

    @overload
    def __getitem__(self, item: slice) -> Sequence[DataUnitGen]:
        ...

    @overload
    def __getitem__(self, item: Union[_I_T, Identifier]) -> DataUnitGen:
        ...

    def __getitem__(self, item: Union[slice, _I_T, Identifier]):
        if isinstance(item, slice):
            return self._list[item]
        if isinstance(item, int):
            return self._list[item]
        return self._list[self._map[item]]

    def i(self, uid: Identifier) -> _I_T:
        return self._map[uid]

    def uid(self, i: _I_T) -> Identifier:
        return self._list[i].uid

    def __len__(self) -> int:
        return len(self._list)

    def __iter__(self) -> Iterator[DataUnitGen]:
        return iter(self._list)


@dataclass(frozen=True)
class __ValueQueryAssoc(Generic[_ID_T, _I_T]):
    __slots__ = ("_list", "_map")
    _list: Sequence[_ID_T]
    _map: Mapping[_ID_T, _I_T]

    @overload
    def __getitem__(self, item: slice) -> Sequence[_ID_T]:
        ...

    @overload
    def __getitem__(self, item: _I_T) -> _ID_T:
        ...

    def __getitem__(self, item: Union[slice, _I_T]):
        if isinstance(item, slice):
            return self._list[item]
        return self._list[item]

    def i(self, item: _ID_T) -> _I_T:
        return self._map[item]

    def __len__(self) -> int:
        return len(self._list)

    def __iter__(self) -> Iterator[_ID_T]:
        return iter(self._list)


class ChemNetwork(ABC):
    __slots__ = ()

    @abstractmethod
    def __init__(self) -> None:
        ...


class ChemNetworkBin(ChemNetwork):
    __slots__ = (
        "_mol_list",
        "_op_list",
        "_rxn_list",
        "_mol_map",
        "_op_map",
        "_rxn_map",
        "_mol_query",
        "_op_query",
        "_rxn_query",
    )

    def __init__(self) -> None:
        self._mol_list: Sequence[MolDatBase] = []
        self._op_list: Sequence[OpDatBase] = []
        self._rxn_list: Sequence[Reaction] = []
        self._mol_map: Mapping[Identifier, _MolIndex] = {}
        self._op_map: Mapping[Identifier, _OpIndex] = {}
        self._rxn_map: Mapping[Reaction, _RxnIndex] = {}
        self._mol_query: Optional[
            __ValueQueryData[MolDatBase, _MolIndex]
        ] = None
        self._op_query: Optional[__ValueQueryData[OpDatBase, _OpIndex]] = None
        self._rxn_query: Optional[__ValueQueryAssoc[Reaction, _RxnIndex]] = None

    @property
    def mols(self) -> __ValueQueryData[MolDatBase, _MolIndex]:
        if self._mol_query is None:
            self._mol_query = __ValueQueryData(self._mol_list, self._mol_map)
        return self._mol_query

    @property
    def ops(self) -> __ValueQueryData[OpDatBase, _OpIndex]:
        if self._op_query is None:
            self._op_query = __ValueQueryData(self._op_list, self._op_map)
        return self._op_query

    @property
    def rxn(self) -> __ValueQueryAssoc[Reaction, _RxnIndex]:
        if self._rxn_query is None:
            self._rxn_query = __ValueQueryAssoc(self._rxn_list, self._rxn_map)
        return self._rxn_query
