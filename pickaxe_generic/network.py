from abc import ABC, abstractmethod
from collections.abc import Iterator
from typing import Generic, Optional, Union, overload

from pickaxe_generic.containers import DataUnitGen
from pickaxe_generic.datatypes import (
    Identifier,
    MolDatBase,
    OpDatBase,
    RxnDatBase,
)


class ChemNetwork(ABC):
    __slots__ = ()

    @abstractmethod
    def __init__(self) -> None:
        ...


class __ValueQuery(Generic[DataUnitGen]):
    __slots__ = ("_list", "_map")

    def __init__(
        self, i_list: list[DataUnitGen], i_map: dict[Identifier, int]
    ) -> None:
        self._list = i_list
        self._map = i_map

    @overload
    def __getitem__(self, i: slice) -> list[DataUnitGen]:
        ...

    @overload
    def __getitem__(self, i: Union[int, Identifier]) -> DataUnitGen:
        ...

    def __getitem__(self, i: Union[slice, int, Identifier]):
        if isinstance(i, slice):
            return self._list[i]
        if isinstance(i, int):
            return self._list[i]
        return self._list[self._map[i]]

    def __len__(self) -> int:
        return len(self._list)

    def __iter__(self) -> Iterator[DataUnitGen]:
        return iter(self._list)


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
        self._mol_list: list[MolDatBase] = []
        self._op_list: list[OpDatBase] = []
        self._rxn_list: list[RxnDatBase] = []
        self._mol_map: dict[Identifier, int] = {}
        self._op_map: dict[Identifier, int] = {}
        self._rxn_map: dict[Identifier, int] = {}
        self._mol_query: Optional[__ValueQuery[MolDatBase]] = None
        self._op_query: Optional[__ValueQuery[OpDatBase]] = None
        self._rxn_query: Optional[__ValueQuery[RxnDatBase]] = None

    @property
    def mols(self) -> __ValueQuery[MolDatBase]:
        if self._mol_query is None:
            self._mol_query = __ValueQuery(self._mol_list, self._mol_map)
        return self._mol_query

    @property
    def ops(self) -> __ValueQuery[OpDatBase]:
        if self._op_query is None:
            self._op_query = __ValueQuery(self._op_list, self._op_map)
        return self._op_query

    @property
    def rxns(self) -> __ValueQuery[RxnDatBase]:
        if self._rxn_query is None:
            self._rxn_query = __ValueQuery(self._rxn_list, self._rxn_map)
        return self._rxn_query

    # def rxn_reactants(self, rxn: RxnDatBase)
