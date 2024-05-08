"""Contains classes for network representation and data storage."""

import collections.abc
import dataclasses
import gzip
import pickle
import typing

from doranet import interfaces


@dataclasses.dataclass(frozen=True, slots=True, order=True)
class MolSlot:
    operator: interfaces.OpDatBase
    molecule: interfaces.MolDatBase
    argnum: int
    operator_meta: typing.Optional[collections.abc.Mapping]
    molecule_meta: typing.Optional[collections.abc.Mapping]


def recipe_from_explicit(
    recipe_explicit: interfaces.RecipeExplicit,
) -> interfaces.Recipe:
    return interfaces.Recipe(
        interfaces.OpIndex(recipe_explicit.operator.i),
        tuple(
            interfaces.MolIndex(data.i) for data in recipe_explicit.reactants
        ),
    )


def index_iterable_guard(
    item: typing.Union[
        interfaces.T_int, collections.abc.Iterable[interfaces.T_int]
    ],
) -> typing.TypeGuard[collections.abc.Iterable[interfaces.T_int]]:
    return hasattr(item, "__iter__")


@dataclasses.dataclass(frozen=True, slots=True)
class _ValueQueryData(
    interfaces.ValueQueryData[interfaces.T_data, interfaces.T_int]
):
    _list: collections.abc.Sequence[interfaces.T_data]
    _map: collections.abc.Mapping[interfaces.Identifier, interfaces.T_int]
    _meta: collections.abc.Sequence[collections.abc.MutableMapping]

    def __contains__(
        self, item: typing.Union[interfaces.Identifier, interfaces.T_data]
    ) -> bool:
        if isinstance(item, interfaces.DataUnit):
            item = item.uid
        return item in self._map

    @typing.overload
    def __getitem__(
        self, item: slice
    ) -> collections.abc.Sequence[interfaces.T_data]: ...

    @typing.overload
    def __getitem__(
        self, item: typing.Union[interfaces.T_int, interfaces.Identifier]
    ) -> interfaces.T_data: ...

    def __getitem__(
        self, item: typing.Union[slice, interfaces.T_int, interfaces.Identifier]
    ):
        if isinstance(item, slice):
            return self._list[item]
        if isinstance(item, int):
            return self._list[item]
        return self._list[self._map[item]]

    def i(self, uid: interfaces.Identifier) -> interfaces.T_int:
        return self._map[uid]

    def keys(self) -> collections.abc.KeysView[interfaces.Identifier]:
        return self._map.keys()

    @typing.overload
    def meta(
        self,
        indices: interfaces.T_int,
        keys: typing.Optional[
            collections.abc.Iterable[collections.abc.Hashable]
        ] = None,
    ) -> collections.abc.Mapping[collections.abc.Hashable, typing.Any]: ...

    @typing.overload
    def meta(
        self,
        indices: typing.Optional[
            collections.abc.Iterable[interfaces.T_int]
        ] = None,
        keys: typing.Optional[
            collections.abc.Iterable[collections.abc.Hashable]
        ] = None,
    ) -> collections.abc.Iterable[
        collections.abc.Mapping[collections.abc.Hashable, typing.Any]
    ]: ...

    def meta(
        self,
        indices: typing.Optional[
            typing.Union[
                interfaces.T_int, collections.abc.Iterable[interfaces.T_int]
            ]
        ] = None,
        keys: typing.Optional[
            collections.abc.Iterable[collections.abc.Hashable]
        ] = None,
    ) -> typing.Union[
        collections.abc.Iterable[
            collections.abc.Mapping[collections.abc.Hashable, typing.Any]
        ],
        collections.abc.Mapping[collections.abc.Hashable, typing.Any],
    ]:
        targets: collections.abc.Iterable[int]
        if indices is None:
            targets = range(len(self._meta))
        elif index_iterable_guard(indices):
            targets = indices
        elif isinstance(indices, int):
            if keys is None:
                return self._meta[indices]
            else:
                m: collections.abc.Mapping = self._meta[indices]
                return {key: m[key] for key in keys if key in m}
        else:
            raise TypeError(
                f"Invalid argument type for `indices`: {type(indices)}"
            )
        if keys is None:
            return tuple(self._meta[i] for i in targets)
        mt = self._meta
        return tuple(
            {key: mt[i][key] for key in keys if key in mt[i]} for i in targets
        )

    def set_meta(
        self, index: interfaces.T_int, values: collections.abc.Mapping
    ) -> None:
        self._meta[index].update(values)

    def uid(self, i: interfaces.T_int) -> interfaces.Identifier:
        return self._list[i].uid

    def __len__(self) -> int:
        return len(self._list)

    def __iter__(self) -> collections.abc.Iterator[interfaces.T_data]:
        return iter(self._list)


@dataclasses.dataclass(frozen=True, slots=True)
class _ValueQueryAssoc(typing.Generic[interfaces.T_id, interfaces.T_int]):
    _list: collections.abc.Sequence[interfaces.T_id]
    _map: collections.abc.Mapping[interfaces.T_id, interfaces.T_int]
    _meta: collections.abc.Sequence[collections.abc.MutableMapping]

    @typing.overload
    def __getitem__(
        self, item: slice
    ) -> collections.abc.Sequence[interfaces.T_id]: ...

    @typing.overload
    def __getitem__(self, item: interfaces.T_int) -> interfaces.T_id: ...

    def __getitem__(self, item: typing.Union[slice, interfaces.T_int]):
        if isinstance(item, slice):
            return self._list[item]
        return self._list[item]

    def i(self, item: interfaces.T_id) -> interfaces.T_int:
        return self._map[item]

    @typing.overload
    def meta(
        self,
        indices: interfaces.T_int,
        keys: typing.Optional[
            collections.abc.Iterable[collections.abc.Hashable]
        ] = None,
    ) -> collections.abc.Mapping[collections.abc.Hashable, typing.Any]: ...

    @typing.overload
    def meta(
        self,
        indices: typing.Optional[
            collections.abc.Iterable[interfaces.T_int]
        ] = None,
        keys: typing.Optional[
            collections.abc.Iterable[collections.abc.Hashable]
        ] = None,
    ) -> collections.abc.Iterable[
        collections.abc.Mapping[collections.abc.Hashable, typing.Any]
    ]: ...

    def meta(
        self,
        indices: typing.Optional[
            typing.Union[
                interfaces.T_int, collections.abc.Iterable[interfaces.T_int]
            ]
        ] = None,
        keys: typing.Optional[
            collections.abc.Iterable[collections.abc.Hashable]
        ] = None,
    ) -> typing.Union[
        collections.abc.Iterable[
            collections.abc.Mapping[collections.abc.Hashable, typing.Any]
        ],
        collections.abc.Mapping[collections.abc.Hashable, typing.Any],
    ]:
        targets: collections.abc.Iterable[int]
        if indices is None:
            targets = range(len(self._meta))
        elif index_iterable_guard(indices):
            targets = indices
        elif isinstance(indices, int):
            if keys is None:
                return self._meta[indices]
            else:
                m: collections.abc.Mapping = self._meta[indices]
                return {key: m[key] for key in keys if key in m}
        else:
            raise TypeError(
                f"Invalid argument type for `indices`: {type(indices)}"
            )
        if keys is None:
            return tuple(self._meta[i] for i in targets)
        mt = self._meta
        return tuple(
            {key: mt[i][key] for key in keys if key in mt[i]} for i in targets
        )

    def set_meta(
        self, index: interfaces.T_int, values: collections.abc.Mapping
    ) -> None:
        self._meta[index].update(values)

    def __len__(self) -> int:
        return len(self._list)

    def __iter__(self) -> collections.abc.Iterator[interfaces.T_id]:
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

        self._compat_table: list[
            collections.abc.Sequence[list[interfaces.MolIndex]]
        ] = []

        self._mol_query: typing.Optional[
            _ValueQueryData[interfaces.MolDatBase, interfaces.MolIndex]
        ] = None
        self._op_query: typing.Optional[
            _ValueQueryData[interfaces.OpDatBase, interfaces.OpIndex]
        ] = None
        self._rxn_query: typing.Optional[
            _ValueQueryAssoc[interfaces.Reaction, interfaces.RxnIndex]
        ] = None

        self._reactive_list: list[bool] = []

    @property
    def mols(
        self,
    ) -> _ValueQueryData[interfaces.MolDatBase, interfaces.MolIndex]:
        if self._mol_query is None:
            self._mol_query = _ValueQueryData(
                self._mol_list, self._mol_map, self._mol_meta
            )
        return self._mol_query

    @property
    def ops(self) -> _ValueQueryData[interfaces.OpDatBase, interfaces.OpIndex]:
        if self._op_query is None:
            self._op_query = _ValueQueryData(
                self._op_list, self._op_map, self._op_meta
            )
        return self._op_query

    @property
    def rxns(
        self,
    ) -> _ValueQueryAssoc[interfaces.Reaction, interfaces.RxnIndex]:
        if self._rxn_query is None:
            self._rxn_query = _ValueQueryAssoc(
                self._rxn_list, self._rxn_map, self._rxn_meta
            )
        return self._rxn_query

    def compat_table(
        self, index: interfaces.OpIndex
    ) -> collections.abc.Sequence[
        collections.abc.Sequence[interfaces.MolIndex]
    ]:
        return self._compat_table[index]

    def consumers(
        self,
        mol: typing.Union[int, interfaces.MolDatBase, interfaces.Identifier],
    ) -> collections.abc.Collection[interfaces.RxnIndex]:
        if isinstance(mol, int):
            return self._mol_consumers[interfaces.MolIndex(mol)]
        elif isinstance(mol, interfaces.MolDatBase):
            return self._mol_consumers[self._mol_map[mol.uid]]
        return self._mol_consumers[self._mol_map[mol]]

    def producers(
        self,
        mol: typing.Union[int, interfaces.MolDatBase, interfaces.Identifier],
    ) -> collections.abc.Collection[interfaces.RxnIndex]:
        if isinstance(mol, int):
            return self._mol_producers[interfaces.MolIndex(mol)]
        elif isinstance(mol, interfaces.MolDatBase):
            return self._mol_producers[self._mol_map[mol.uid]]
        return self._mol_producers[self._mol_map[mol]]

    def add_mol(
        self,
        mol: interfaces.MolDatBase,
        meta: typing.Optional[collections.abc.Mapping] = None,
        reactive: typing.Optional[bool] = None,
        _custom_compat: typing.Optional[
            collections.abc.Collection[tuple[interfaces.OpIndex, int]]
        ] = None,
    ) -> interfaces.MolIndex:
        if _custom_compat is not None and reactive is False:
            raise ValueError(
                "`_custom_compat` cannot be specified when `reactive` is False"
            )
        # if already in database, return existing index
        mol_uid = mol.uid
        if mol_uid in self._mol_map:
            mol_index = self._mol_map[mol_uid]
            if meta is not None:
                self._mol_meta[mol_index].update(meta)

            if reactive is not True:
                return mol_index

            # if newly reactive, fill in compat table
            if not self._reactive_list[mol_index]:
                self._reactive_list[mol_index] = True
                if _custom_compat is None:
                    for i, op in enumerate(self.ops):
                        for argnum in range(len(op)):
                            if op.compat(mol, argnum):
                                self._compat_table[i][argnum].append(mol_index)
                else:
                    for op_index, argnum in _custom_compat:
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

        self._reactive_list.append(reactive is not False)
        if reactive is not False:
            if _custom_compat is None:
                # test operator compatibility and add to table
                for i, op in enumerate(self.ops):
                    for argnum in range(len(op)):
                        if op.compat(mol, argnum):
                            self._compat_table[i][argnum].append(mol_index)
            else:
                for op_index, argnum in _custom_compat:
                    self._compat_table[op_index][argnum].append(mol_index)

        return mol_index

    def add_op(
        self,
        op: interfaces.OpDatBase,
        meta: typing.Optional[collections.abc.Mapping] = None,
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
        operator: typing.Optional[interfaces.OpIndex] = None,
        reactants: typing.Optional[
            collections.abc.Sequence[interfaces.MolIndex]
        ] = None,
        products: typing.Optional[
            collections.abc.Sequence[interfaces.MolIndex]
        ] = None,
        meta: typing.Optional[collections.abc.Mapping] = None,
        rxn: typing.Optional[interfaces.Reaction] = None,
    ) -> interfaces.RxnIndex:
        if rxn is None:
            if operator is None or reactants is None or products is None:
                raise ValueError(
                    f"""op ({operator}), reactants ({reactants}), and products
                        ({products}) must all be specified if reaction is
                        None"""
                )
            rxn = interfaces.Reaction(
                operator, tuple(reactants), tuple(products)
            )

        # if already in database, return existing index
        if rxn in self._rxn_map:
            rxn_index = self._rxn_map[rxn]
            if meta is not None:
                self._rxn_meta[rxn_index].update(meta)
            return rxn_index

        # sanity check that all reactants and products exist in the network
        if (
            max(*rxn.reactants, *rxn.products) >= len(self._mol_list)
            or min(*rxn.reactants, *rxn.products) < 0
        ):
            raise IndexError(
                f"""One of the molecule components for reaction {rxn} is not in
                    the network."""
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

    @property
    def reactivity(self) -> collections.abc.Sequence[bool]:
        return self._reactive_list


@dataclasses.dataclass(frozen=True, slots=True)
class ChemNetworkFacadeMetaTrigger(interfaces.ChemNetwork):
    network: interfaces.ChemNetwork
    metakeys: interfaces.MetaKeyPacket
    updated_mol_meta: set[interfaces.MolIndex] = dataclasses.field(
        default_factory=set
    )
    updated_op_meta: set[interfaces.OpIndex] = dataclasses.field(
        default_factory=set
    )

    @property
    def mols(self):
        return self.network.mols

    @property
    def ops(self):
        return self.network.ops

    @property
    def rxns(self):
        return self.network.rxns

    def compat_table(self, index):
        return self.network.compat_table(index)

    def consumers(self, mol):
        return self.network.consumers(mol)

    def producers(self, mol):
        return self.network.producers(mol)

    def add_mol(self, mol, meta=None, reactive=None, _custom_compat=None):
        if meta is not None and mol in self.network.mols:
            i = self.network.mols.i(mol.uid)
            exist_keys = self.network.mols.meta(i, meta.keys())
            for key in exist_keys:
                if exist_keys[key] != meta[key]:
                    self.updated_mol_meta.add(i)
                    continue
        return self.network.add_mol(mol, meta, reactive, _custom_compat)

    def add_op(self, op, meta=None):
        if meta is not None and op in self.network.ops:
            i = self.network.ops.i(op.uid)
            exist_keys = self.network.ops.meta(i, meta.keys())
            for key in exist_keys:
                if exist_keys[key] != meta[key]:
                    self.updated_op_meta.add(i)
                    continue
        return self.network.add_op(op, meta)

    def add_rxn(
        self, operator=None, reactants=None, products=None, meta=None, rxn=None
    ):
        return self.network.add_rxn(operator, reactants, products, meta, rxn)

    @property
    def reactivity(self):
        return self.network.reactivity

    def save_to_file(self, filename, path="", minimal=False, ext=".pgnet"):
        self.network.save_to_file(filename, path, minimal, ext)


def dump_network_to_file(
    network: interfaces.ChemNetwork, filepath: str = "network.dat"
) -> None:
    with gzip.open(filepath, "wb") as fout:
        pickle.dump(network, fout)


def load_network_from_file(
    filepath: str = "network.dat",
) -> interfaces.ChemNetwork:
    with gzip.open(filepath, "rb") as fin:
        return pickle.load(fin)
