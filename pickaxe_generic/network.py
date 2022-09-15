from collections.abc import Collection, Hashable, Iterator, Mapping, Sequence
from dataclasses import dataclass
from functools import cache
from gzip import open as gzopen
from pickle import dump, load
from typing import (
    Any,
    Generic,
    NewType,
    Optional,
    Protocol,
    TypeVar,
    Union,
    overload,
)

from . import interfaces


@dataclass(frozen=True, order=True)
class MolSlot:
    __slots__ = (
        "operator",
        "molecule",
        "argnum",
        "operator_meta",
        "molecule_meta",
    )
    operator: interfaces.OpDatBase
    molecule: interfaces.MolDatBase
    argnum: int
    operator_meta: Optional[Mapping]
    molecule_meta: Optional[Mapping]


def recipe_from_explicit(
    recipe_explicit: interfaces.RecipeExplicit,
) -> interfaces.Recipe:
    return interfaces.Recipe(
        interfaces.OpIndex(recipe_explicit.operator.i),
        tuple(
            interfaces.MolIndex(data.i) for data in recipe_explicit.reactants
        ),
    )


@dataclass(frozen=True)
class _ValueQueryData(Generic[interfaces.T_data, interfaces.T_int]):
    __slots__ = ("_list", "_map")
    _list: Sequence[interfaces.T_data]
    _map: Mapping[interfaces.Identifier, interfaces.T_int]

    def __contains__(
        self, item: Union[interfaces.Identifier, interfaces.T_data]
    ) -> bool:
        if isinstance(item, interfaces.DataUnit):
            item = item.uid
        return item in self._map.keys()

    @overload
    def __getitem__(self, item: slice) -> Sequence[interfaces.T_data]:
        ...

    @overload
    def __getitem__(
        self, item: Union[interfaces.T_int, interfaces.Identifier]
    ) -> interfaces.T_data:
        ...

    def __getitem__(
        self, item: Union[slice, interfaces.T_int, interfaces.Identifier]
    ):
        if isinstance(item, slice):
            return self._list[item]
        if isinstance(item, int):
            return self._list[item]
        return self._list[self._map[item]]

    def i(self, uid: interfaces.Identifier) -> interfaces.T_int:
        return self._map[uid]

    def keys(self) -> Collection[interfaces.Identifier]:
        return self._map.keys()

    def uid(self, i: interfaces.T_int) -> interfaces.Identifier:
        return self._list[i].uid

    def __len__(self) -> int:
        return len(self._list)

    def __iter__(self) -> Iterator[interfaces.T_data]:
        return iter(self._list)


@dataclass(frozen=True)
class _ValueQueryAssoc(Generic[interfaces.T_id, interfaces.T_int]):
    __slots__ = ("_list", "_map")
    _list: Sequence[interfaces.T_id]
    _map: Mapping[interfaces.T_id, interfaces.T_int]

    @overload
    def __getitem__(self, item: slice) -> Sequence[interfaces.T_id]:
        ...

    @overload
    def __getitem__(self, item: interfaces.T_int) -> interfaces.T_id:
        ...

    def __getitem__(self, item: Union[slice, interfaces.T_int]):
        if isinstance(item, slice):
            return self._list[item]
        return self._list[item]

    def i(self, item: interfaces.T_id) -> interfaces.T_int:
        return self._map[item]

    def __len__(self) -> int:
        return len(self._list)

    def __iter__(self) -> Iterator[interfaces.T_id]:
        return iter(self._list)


class ChemNetworkBasic(interfaces.ChemNetwork):
    __slots__ = (
        "_mol_list",
        "_op_list",
        "_rxn_list",
        "_mol_map",
        "_op_map",
        "_rxn_map",
        "_mol_meta",
        "_op_meta",
        "_rxn_meta",
        "_mol_producers",
        "_mol_consumers",
        "_compat_table",
        "_mol_query",
        "_op_query",
        "_rxn_query",
        "_reactive_list",
    )

    def __init__(self) -> None:
        self._mol_list: list[interfaces.MolDatBase] = []
        self._op_list: list[interfaces.OpDatBase] = []
        self._rxn_list: list[interfaces.Reaction] = []

        self._mol_map: dict[interfaces.Identifier, interfaces.MolIndex] = {}
        self._op_map: dict[interfaces.Identifier, interfaces.OpIndex] = {}
        self._rxn_map: dict[interfaces.Reaction, interfaces.RxnIndex] = {}

        self._mol_meta: list[dict] = []
        self._op_meta: list[dict] = []
        self._rxn_meta: list[dict] = []

        self._mol_producers: list[list[interfaces.RxnIndex]] = []
        self._mol_consumers: list[list[interfaces.RxnIndex]] = []

        self._compat_table: list[Sequence[list[interfaces.MolIndex]]] = []

        self._mol_query: Optional[
            _ValueQueryData[interfaces.MolDatBase, interfaces.MolIndex]
        ] = None
        self._op_query: Optional[
            _ValueQueryData[interfaces.OpDatBase, interfaces.OpIndex]
        ] = None
        self._rxn_query: Optional[
            _ValueQueryAssoc[interfaces.Reaction, interfaces.RxnIndex]
        ] = None

        self._reactive_list: list[bool] = []

    @property
    def mols(
        self,
    ) -> _ValueQueryData[interfaces.MolDatBase, interfaces.MolIndex]:
        if self._mol_query is None:
            self._mol_query = _ValueQueryData(self._mol_list, self._mol_map)
        return self._mol_query

    @property
    def ops(self) -> _ValueQueryData[interfaces.OpDatBase, interfaces.OpIndex]:
        if self._op_query is None:
            self._op_query = _ValueQueryData(self._op_list, self._op_map)
        return self._op_query

    @property
    def rxns(
        self,
    ) -> _ValueQueryAssoc[interfaces.Reaction, interfaces.RxnIndex]:
        if self._rxn_query is None:
            self._rxn_query = _ValueQueryAssoc(self._rxn_list, self._rxn_map)
        return self._rxn_query

    def mol_meta(self, index: interfaces.MolIndex, key: Hashable, value=None):
        if value is None:
            if key not in self._mol_meta[index]:
                return None
            return self._mol_meta[index][key]
        self._mol_meta[index][key] = value

    def op_meta(self, index: interfaces.OpIndex, key: Hashable, value=None):
        if value is None:
            if key not in self._op_meta[index]:
                return None
            return self._op_meta[index][key]
        self._op_meta[index][key] = value

    def rxn_meta(self, index: interfaces.RxnIndex, key: Hashable, value=None):
        if value is None:
            if key not in self._rxn_meta[index]:
                return None
            return self._rxn_meta[index][key]
        self._rxn_meta[index][key] = value

    def mol_metas(
        self,
        indices: Optional[Sequence[interfaces.MolIndex]] = None,
        keys: Optional[Collection[Hashable]] = None,
    ) -> Sequence[Mapping[Hashable, Any]]:
        if indices is None:
            if keys is None:
                return self._mol_meta
            return tuple(
                {
                    key: self._mol_meta[i][key]
                    for key in keys
                    if key in self._mol_meta[i]
                }
                for i in range(len(self._mol_list))
            )
        if keys is None:
            return tuple(self._mol_meta[i] for i in indices)
        return tuple(
            {
                key: self._mol_meta[i][key]
                for key in keys
                if key in self._mol_meta[i]
            }
            for i in indices
        )

    def op_metas(
        self,
        indices: Optional[Sequence[interfaces.OpIndex]] = None,
        keys: Optional[Collection[Hashable]] = None,
    ) -> Sequence[Mapping[Hashable, Any]]:
        if indices is None:
            if keys is None:
                return self._op_meta
            return tuple(
                {
                    key: self._op_meta[i][key]
                    for key in keys
                    if key in self._op_meta[i]
                }
                for i in range(len(self._op_list))
            )
        if keys is None:
            return tuple(self._op_meta[i] for i in indices)
        return tuple(
            {
                key: self._op_meta[i][key]
                for key in keys
                if key in self._op_meta[i]
            }
            for i in indices
        )

    def rxn_metas(
        self,
        indices: Optional[Sequence[interfaces.RxnIndex]] = None,
        keys: Optional[Collection[Hashable]] = None,
    ) -> Sequence[Mapping[Hashable, Any]]:
        if indices is None:
            if keys is None:
                return self._rxn_meta
            return tuple(
                {
                    key: self._rxn_meta[i][key]
                    for key in keys
                    if key in self._rxn_meta[i]
                }
                for i in range(len(self._rxn_list))
            )
        if keys is None:
            return tuple(self._rxn_meta[i] for i in indices)
        return tuple(
            {
                key: self._rxn_meta[i][key]
                for key in keys
                if key in self._rxn_meta[i]
            }
            for i in indices
        )

    def compat_table(
        self, index: int
    ) -> Sequence[Sequence[interfaces.MolIndex]]:
        return self._compat_table[index]

    def consumers(
        self, mol: Union[int, interfaces.MolDatBase, interfaces.Identifier]
    ) -> Collection[interfaces.RxnIndex]:
        if isinstance(mol, int):
            return self._mol_consumers[interfaces.MolIndex(mol)]
        elif isinstance(mol, interfaces.MolDatBase):
            return self._mol_consumers[self._mol_map[mol.uid]]
        return self._mol_consumers[self._mol_map[mol]]

    def producers(
        self, mol: Union[int, interfaces.MolDatBase, interfaces.Identifier]
    ) -> Collection[interfaces.RxnIndex]:
        if isinstance(mol, int):
            return self._mol_producers[interfaces.MolIndex(mol)]
        elif isinstance(mol, interfaces.MolDatBase):
            return self._mol_producers[self._mol_map[mol.uid]]
        return self._mol_producers[self._mol_map[mol]]

    def add_mol(
        self,
        mol: interfaces.MolDatBase,
        meta: Optional[Mapping] = None,
        reactive: bool = True,
        custom_compat: Optional[
            Collection[tuple[interfaces.OpIndex, int]]
        ] = None,
    ) -> interfaces.MolIndex:
        # if already in database, return existing index
        mol_uid = mol.uid
        if mol_uid in self._mol_map:
            mol_index = self._mol_map[mol_uid]
            if meta is not None:
                self._mol_meta[mol_index].update(meta)

            if not reactive:
                return mol_index

            # if newly reactive, fill in compat table
            if not self._reactive_list[mol_index]:
                self._reactive_list[mol_index] = True
                if custom_compat is None:
                    for i, op in enumerate(self.ops):
                        for argnum in range(len(op)):
                            if op.compat(mol, argnum):
                                self._compat_table[i][argnum].append(mol_index)
                else:
                    for op_index, argnum in custom_compat:
                        self._compat_table[op_index][argnum].append(mol_index)
            return mol_index

        # add mol to main mol list
        mol_index = interfaces.MolIndex(len(self._mol_list))
        self._mol_list.append(mol)

        # add mol id to UID mapping
        self._mol_map[mol_uid] = mol_index

        # extend consumer/producer table
        self._mol_consumers.append([])
        self._mol_producers.append([])

        # add mol metadata to table
        if meta is None:
            self._mol_meta.append({})
        else:
            self._mol_meta.append(dict(meta))

        self._reactive_list.append(reactive)
        if reactive:
            if custom_compat is None:
                # test operator compatibility and add to table
                for i, op in enumerate(self.ops):
                    for argnum in range(len(op)):
                        if op.compat(mol, argnum):
                            self._compat_table[i][argnum].append(mol_index)
            else:
                for op_index, argnum in custom_compat:
                    self._compat_table[op_index][argnum].append(mol_index)

        return mol_index

    def add_op(
        self, op: interfaces.OpDatBase, meta: Optional[Mapping] = None
    ) -> interfaces.OpIndex:
        # if already in database, return existing index
        op_uid = op.uid
        if op_uid in self._op_map:
            op_index = self._op_map[op_uid]
            if meta is not None:
                self._op_meta[op_index].update(meta)
            return op_index

        # add op to main op list
        op_index = interfaces.OpIndex(len(self._op_list))
        self._op_list.append(op)

        # add op id to UID mapping
        self._op_map[op_uid] = op_index

        # add mol metadata to table
        if meta is None:
            self._op_meta.append({})
        else:
            self._op_meta.append(dict(meta))

        # test operator compatibility and add to table
        self._compat_table.append(
            tuple(
                [
                    [
                        interfaces.MolIndex(mol_index)
                        for mol_index, mol in enumerate(self._mol_list)
                        if self._reactive_list[mol_index]
                        and op.compat(mol, argnum)
                    ]
                    for argnum in range(len(op))
                ]
            )
        )

        return op_index

    def add_rxn(
        self,
        rxn: Optional[interfaces.Reaction] = None,
        op: Optional[interfaces.OpIndex] = None,
        reactants: Optional[Sequence[interfaces.MolIndex]] = None,
        products: Optional[Sequence[interfaces.MolIndex]] = None,
        meta: Optional[Mapping] = None,
    ) -> interfaces.RxnIndex:

        if rxn is None:
            if op is None or reactants is None or products is None:
                raise ValueError(
                    f"op ({op}), reactants ({reactants}), and products ({products}) must all be specified if reaction is None"
                )
            rxn = interfaces.Reaction(op, tuple(reactants), tuple(products))

        # if already in database, return existing index
        if rxn in self._rxn_map:
            rxn_index = self._rxn_map[rxn]
            if meta is not None:
                self._rxn_meta[rxn_index].update(meta)
            return rxn_index

        # sanity check that all reactants and products exist in the network
        if (
            max(max(rxn.reactants), max(rxn.products)) >= len(self._mol_list)
            or min(min(rxn.reactants), min(rxn.products)) < 0
        ):
            raise IndexError(
                f"One of the molecule components for reaction {rxn} is not in the network."
            )
        # sanity check that operator exists in the network
        if rxn.operator >= len(self._op_list):
            raise IndexError(
                f"The operator for reaction {rxn} is not in the network."
            )

        # add rxn to main rxn list
        rxn_index = interfaces.RxnIndex(len(self._rxn_list))
        self._rxn_list.append(rxn)

        # add rxn to index mapping
        self._rxn_map[rxn] = rxn_index

        # add consumption/production mappings
        for i in rxn.reactants:
            self._mol_consumers[i].append(rxn_index)
        for i in rxn.products:
            self._mol_producers[i].append(rxn_index)

        # add rxn metadata to table
        if meta is None:
            self._rxn_meta.append({})
        else:
            self._rxn_meta.append(dict(meta))

        return rxn_index


def dump_network_to_file(
    network: interfaces.ChemNetwork, filepath: str = "network.dat"
) -> None:
    with gzopen(filepath, "wb") as fout:
        dump(network, fout)


def load_network_from_file(
    filepath: str = "network.dat",
) -> interfaces.ChemNetwork:
    with gzopen(filepath, "rb") as fin:
        return load(fin)
