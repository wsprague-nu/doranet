"""Logic and interfaces for metadata processing."""

import abc
import collections.abc
import dataclasses
import itertools
import operator
import typing

from doranet import interfaces, utils


class MetaSink(typing.Protocol):
    __slots__ = ()

    @property
    @abc.abstractmethod
    def meta_required(self) -> interfaces.MetaKeyPacket: ...


_T = typing.TypeVar("_T")
_U = typing.TypeVar("_U")
MetaDataResolverFunc = collections.abc.Callable[[_T, _T], _T]


def TrivialMetaDataResolverFunc(a: _T, b: _T) -> _T:
    return a


class LocalPropertyCalc(abc.ABC, typing.Generic[_T]):
    @property
    @abc.abstractmethod
    def key(self) -> collections.abc.Hashable: ...

    @property
    @abc.abstractmethod
    def meta_required(self) -> interfaces.MetaKeyPacket: ...

    @property
    @abc.abstractmethod
    def resolver(self) -> MetaDataResolverFunc[_T]: ...

    @typing.final
    def __and__(
        self, other: typing.Union["PropertyCompositor", "LocalPropertyCalc"]
    ) -> "PropertyCompositor":
        return MergePropertyCompositor(
            _as_property_compositor(self), _as_property_compositor(other)
        )

    @typing.final
    def __add__(
        self, other: typing.Union["PropertyCompositor", "LocalPropertyCalc"]
    ) -> "PropertyCompositor":
        return FunctionPropertyCompositor(
            operator.add,
            _as_property_compositor(self),
            _as_property_compositor(other),
        )

    @typing.final
    def __sub__(
        self, other: typing.Union["PropertyCompositor", "LocalPropertyCalc"]
    ) -> "PropertyCompositor":
        return FunctionPropertyCompositor(
            operator.sub,
            _as_property_compositor(self),
            _as_property_compositor(other),
        )

    @typing.final
    def __mul__(
        self, other: typing.Union["PropertyCompositor", "LocalPropertyCalc"]
    ) -> "PropertyCompositor":
        return FunctionPropertyCompositor(
            operator.mul,
            _as_property_compositor(self),
            _as_property_compositor(other),
        )

    @typing.final
    def __truediv__(
        self, other: typing.Union["PropertyCompositor", "LocalPropertyCalc"]
    ) -> "PropertyCompositor":
        return FunctionPropertyCompositor(
            operator.truediv,
            _as_property_compositor(self),
            _as_property_compositor(other),
        )

    @typing.final
    def __pow__(
        self, other: typing.Union["PropertyCompositor", "LocalPropertyCalc"]
    ) -> "PropertyCompositor":
        return FunctionPropertyCompositor(
            operator.pow,
            _as_property_compositor(self),
            _as_property_compositor(other),
        )

    @typing.final
    def __rshift__(
        self,
        other: typing.Union[
            "RxnAnalysisStep",
            "PropertyCompositor",
            "ReactionFilterBase",
            "LocalPropertyCalc",
        ],
    ) -> "RxnAnalysisStep":
        return as_rxn_analysis_step(self) >> as_rxn_analysis_step(other)


class MolPropertyCalc(LocalPropertyCalc[_T]):
    @abc.abstractmethod
    def __call__(
        self,
        data: interfaces.DataPacketE[interfaces.MolDatBase],
        prev_value: typing.Optional[_T] = None,
    ) -> typing.Optional[_T]: ...


class MolPropertyFromRxnCalc(LocalPropertyCalc[_T]):
    @abc.abstractmethod
    def __call__(
        self,
        data: interfaces.DataPacketE[interfaces.MolDatBase],
        rxn: interfaces.ReactionExplicit,
        prev_value: typing.Optional[_T] = None,
    ) -> typing.Optional[_T]: ...


class OpPropertyCalc(LocalPropertyCalc[_T]):
    @abc.abstractmethod
    def __call__(
        self,
        data: interfaces.DataPacketE[interfaces.OpDatBase],
        prev_value: typing.Optional[_T] = None,
    ) -> typing.Optional[_T]: ...


class OpPropertyFromRxnCalc(LocalPropertyCalc[_T]):
    @abc.abstractmethod
    def __call__(
        self,
        data: interfaces.DataPacketE[interfaces.OpDatBase],
        rxn: interfaces.ReactionExplicit,
        prev_value: typing.Optional[_T] = None,
    ) -> typing.Optional[_T]: ...


class RxnPropertyCalc(LocalPropertyCalc[_T]):
    @abc.abstractmethod
    def __call__(
        self,
        data: interfaces.ReactionExplicit,
        prev_value: typing.Optional[_T] = None,
    ) -> typing.Optional[_T]: ...


@dataclasses.dataclass(frozen=True)
class KeyOutput:
    mol_keys: frozenset[collections.abc.Hashable]
    op_keys: frozenset[collections.abc.Hashable]
    rxn_keys: frozenset[collections.abc.Hashable]

    def __or__(self, other: "KeyOutput") -> "KeyOutput":
        return KeyOutput(
            self.mol_keys | other.mol_keys,
            self.op_keys | other.op_keys,
            self.rxn_keys | other.rxn_keys,
        )

    def __and__(self, other: "KeyOutput") -> "KeyOutput":
        new_mol_keys = self.mol_keys | other.mol_keys
        if len(new_mol_keys) > len(self.mol_keys) + len(other.mol_keys):
            raise KeyError(
                f"""Conflicting molecule metadata key outputs {self.mol_keys &
                    other.mol_keys}; separate expressions with >> or combine
                    using other operator"""
            )
        new_op_keys = self.op_keys | other.op_keys
        if len(new_op_keys) > len(self.op_keys) + len(other.op_keys):
            raise KeyError(
                f"""Conflicting operator metadata key outputs {self.op_keys &
                    other.op_keys}; separate expressions with >> or combine
                    using other operator"""
            )
        new_rxn_keys = self.rxn_keys | other.rxn_keys
        if len(new_rxn_keys) > len(self.rxn_keys) + len(other.rxn_keys):
            raise KeyError(
                f"""Conflicting reaction metadata key outputs {self.rxn_keys &
                    other.rxn_keys}; separate expressions with >> or combine
                    using other operator"""
            )
        return KeyOutput(new_mol_keys, new_op_keys, new_rxn_keys)


class ReactionFilterBase(abc.ABC):
    __slots__ = ()

    @abc.abstractmethod
    def __call__(self, recipe: interfaces.ReactionExplicit) -> bool: ...

    @property
    def meta_required(self) -> interfaces.MetaKeyPacket:
        return interfaces.MetaKeyPacket()

    @typing.final
    def __and__(self, other: "ReactionFilterBase") -> "ReactionFilterBase":
        return ReactionFilterAnd(self, other)

    @typing.final
    def __invert__(self) -> "ReactionFilterBase":
        return ReactionFilterInv(self)

    @typing.final
    def __or__(self, other: "ReactionFilterBase") -> "ReactionFilterBase":
        return ReactionFilterOr(self, other)

    @typing.final
    def __xor__(self, other: "ReactionFilterBase") -> "ReactionFilterBase":
        return ReactionFilterXor(self, other)

    @typing.final
    def __rshift__(
        self,
        other: typing.Union[
            "RxnAnalysisStep",
            "PropertyCompositor",
            "ReactionFilterBase",
            LocalPropertyCalc,
        ],
    ) -> "RxnAnalysisStep":
        return as_rxn_analysis_step(self) >> as_rxn_analysis_step(other)


@dataclasses.dataclass(frozen=True)
class ReactionFilterAnd(ReactionFilterBase):
    __slots__ = ("_filter1", "_filter2")

    _filter1: ReactionFilterBase
    _filter2: ReactionFilterBase

    def __call__(self, recipe: interfaces.ReactionExplicit) -> bool:
        return self._filter1(recipe) and self._filter2(recipe)

    @property
    def meta_required(self) -> interfaces.MetaKeyPacket:
        return self._filter1.meta_required + self._filter2.meta_required


@dataclasses.dataclass(frozen=True)
class ReactionFilterInv(ReactionFilterBase):
    __slots__ = ("_filter",)
    _filter: ReactionFilterBase

    def __call__(self, recipe: interfaces.ReactionExplicit) -> bool:
        return not self._filter(recipe)

    @property
    def meta_required(self) -> interfaces.MetaKeyPacket:
        return self._filter.meta_required


@dataclasses.dataclass(frozen=True)
class ReactionFilterOr(ReactionFilterBase):
    __slots__ = ("_filter1", "_filter2")
    _filter1: ReactionFilterBase
    _filter2: ReactionFilterBase

    def __call__(self, recipe: interfaces.ReactionExplicit) -> bool:
        return self._filter1(recipe) or self._filter2(recipe)

    @property
    def meta_required(self) -> interfaces.MetaKeyPacket:
        return self._filter1.meta_required + self._filter2.meta_required


@dataclasses.dataclass(frozen=True)
class ReactionFilterXor(ReactionFilterBase):
    __slots__ = ("_filter1", "_filter2")
    _filter1: ReactionFilterBase
    _filter2: ReactionFilterBase

    def __call__(self, recipe: interfaces.ReactionExplicit) -> bool:
        return self._filter1(recipe) != self._filter2(recipe)

    @property
    def meta_required(self) -> interfaces.MetaKeyPacket:
        return self._filter1.meta_required + self._filter2.meta_required


class PropertyCompositor(abc.ABC):
    @abc.abstractmethod
    def __call__(
        self, rxn: interfaces.ReactionExplicit
    ) -> "MetaPropertyState": ...

    @property
    @abc.abstractmethod
    def keys(self) -> KeyOutput: ...

    @typing.final
    def __and__(
        self, other: typing.Union["PropertyCompositor", LocalPropertyCalc]
    ) -> "PropertyCompositor":
        return MergePropertyCompositor(self, _as_property_compositor(other))

    @typing.final
    def __add__(
        self, other: typing.Union["PropertyCompositor", LocalPropertyCalc]
    ) -> "PropertyCompositor":
        return FunctionPropertyCompositor(
            operator.add, self, _as_property_compositor(other)
        )

    @typing.final
    def __sub__(
        self, other: typing.Union["PropertyCompositor", LocalPropertyCalc]
    ) -> "PropertyCompositor":
        return FunctionPropertyCompositor(
            operator.sub, self, _as_property_compositor(other)
        )

    @typing.final
    def __mul__(
        self, other: typing.Union["PropertyCompositor", LocalPropertyCalc]
    ) -> "PropertyCompositor":
        return FunctionPropertyCompositor(
            operator.mul, self, _as_property_compositor(other)
        )

    @typing.final
    def __truediv__(
        self, other: typing.Union["PropertyCompositor", LocalPropertyCalc]
    ) -> "PropertyCompositor":
        return FunctionPropertyCompositor(
            operator.truediv, self, _as_property_compositor(other)
        )

    @typing.final
    def __pow__(
        self, other: typing.Union["PropertyCompositor", LocalPropertyCalc]
    ) -> "PropertyCompositor":
        return FunctionPropertyCompositor(
            operator.pow, self, _as_property_compositor(other)
        )

    @typing.final
    def __rshift__(
        self,
        other: typing.Union[
            "RxnAnalysisStep",
            "PropertyCompositor",
            ReactionFilterBase,
            LocalPropertyCalc,
        ],
    ) -> "RxnAnalysisStep":
        return as_rxn_analysis_step(self) >> as_rxn_analysis_step(other)

    @property
    @abc.abstractmethod
    def meta_required(self) -> interfaces.MetaKeyPacket: ...

    @property
    @abc.abstractmethod
    def resolver(self) -> "MetaUpdateResolver": ...


class MergePropertyCompositor(PropertyCompositor):
    __slots__ = ("_comp1", "_comp2")

    _comp1: PropertyCompositor
    _comp2: PropertyCompositor
    _keys: KeyOutput

    def __init__(
        self, comp1: PropertyCompositor, comp2: PropertyCompositor
    ) -> None:
        self._comp1 = comp1
        self._comp2 = comp2
        self._keys = comp1.keys & comp2.keys

    def __call__(self, rxn: interfaces.ReactionExplicit) -> "MetaPropertyState":
        state1 = self._comp1(rxn)
        state2 = self._comp2(rxn)
        state1.mol_info.update(state2.mol_info)
        state1.op_info.update(state2.op_info)
        state1.rxn_info.update(state2.rxn_info)
        return state1

    @property
    def keys(self) -> KeyOutput:
        return self._keys

    @property
    def meta_required(self) -> interfaces.MetaKeyPacket:
        return self._comp1.meta_required + self._comp2.meta_required

    @property
    def resolver(self) -> "MetaUpdateResolver":
        return self._comp1.resolver | self._comp2.resolver


@dataclasses.dataclass(frozen=True)
class FunctionPropertyCompositor(PropertyCompositor):
    __slots__ = ("_func", "_comp1", "_comp2")

    _func: collections.abc.Callable[[typing.Any, typing.Any], typing.Any]
    _comp1: PropertyCompositor
    _comp2: PropertyCompositor

    def __call__(self, rxn: interfaces.ReactionExplicit) -> "MetaPropertyState":
        state1 = self._comp1(rxn)
        state2 = self._comp2(rxn)

        for key in state1.mol_info.keys() & state2.mol_info.keys():
            propstate1 = state1.mol_info[key]
            propstate2 = state2.mol_info[key]
            propstate_reduced = MetaPropertyStateSingleProp(
                {}, propstate2.resolver
            )
            for mol in propstate1.data.keys() & propstate2.data.keys():
                propstate_reduced.data[mol] = self._func(
                    propstate1.data[mol], propstate2.data[mol]
                )
            propstate1.overwrite(propstate2).overwrite(propstate_reduced)
        for key in state1.op_info.keys() & state2.op_info.keys():
            propstate1 = state1.op_info[key]
            propstate2 = state2.op_info[key]
            propstate_reduced = MetaPropertyStateSingleProp(
                {}, propstate2.resolver
            )
            for op in propstate1.data.keys() & propstate2.data.keys():
                propstate_reduced.data[op] = self._func(
                    propstate1.data[op], propstate2.data[op]
                )
            propstate1.overwrite(propstate2).overwrite(propstate_reduced)
        for key in state1.rxn_info.keys() & state2.rxn_info.keys():
            propstate1 = state1.rxn_info[key]
            propstate2 = state2.rxn_info[key]
            propstate_reduced = MetaPropertyStateSingleProp(
                {}, propstate2.resolver
            )
            for rxn_id in propstate1.data.keys() & propstate2.data.keys():
                propstate_reduced.data[rxn_id] = self._func(
                    propstate1.data[rxn_id], propstate2.data[rxn_id]
                )
            propstate1.overwrite(propstate2).overwrite(propstate_reduced)

        return state1

    @property
    def keys(self) -> KeyOutput:
        return self._comp1.keys | self._comp2.keys

    @property
    def meta_required(self) -> interfaces.MetaKeyPacket:
        return self._comp1.meta_required + self._comp2.meta_required

    @property
    def resolver(self) -> "MetaUpdateResolver":
        return self._comp1.resolver | self._comp2.resolver


@dataclasses.dataclass(frozen=True)
class MolPropertyCompositor(PropertyCompositor, typing.Generic[_T]):
    __slots__ = ("_calc",)

    _calc: MolPropertyCalc[_T]

    def __call__(self, rxn: interfaces.ReactionExplicit) -> "MetaPropertyState":
        mols = itertools.chain(rxn.reactants, rxn.products)
        props = {
            mol.item.uid: calc
            for mol, calc in ((mol, self._calc(mol)) for mol in mols)
            if calc is not None
        }
        single_state = MetaPropertyStateSingleProp(props, self._calc.resolver)
        return MetaPropertyState({self._calc.key: single_state}, {}, {})

    @property
    def keys(self) -> KeyOutput:
        return KeyOutput(frozenset((self._calc.key,)), frozenset(), frozenset())

    @property
    def meta_required(self) -> interfaces.MetaKeyPacket:
        return self._calc.meta_required

    @property
    def resolver(self) -> "MetaUpdateResolver":
        calc = self._calc
        return MetaUpdateResolver({calc.key: calc.resolver}, {}, {})


@dataclasses.dataclass(frozen=True)
class MolRxnPropertyCompositor(PropertyCompositor, typing.Generic[_T]):
    __slots__ = ("_calc",)

    _calc: MolPropertyFromRxnCalc[_T]

    def __call__(self, rxn: interfaces.ReactionExplicit) -> "MetaPropertyState":
        mols = itertools.chain(rxn.reactants, rxn.products)
        props = {
            mol.item.uid: calc
            for mol, calc in ((mol, self._calc(mol, rxn)) for mol in mols)
            if calc is not None
        }
        single_state = MetaPropertyStateSingleProp(props, self._calc.resolver)
        return MetaPropertyState({self._calc.key: single_state}, {}, {})

    @property
    def keys(self) -> KeyOutput:
        return KeyOutput(frozenset((self._calc.key,)), frozenset(), frozenset())

    @property
    def meta_required(self) -> interfaces.MetaKeyPacket:
        return self._calc.meta_required

    @property
    def resolver(self) -> "MetaUpdateResolver":
        calc = self._calc
        return MetaUpdateResolver({calc.key: calc.resolver}, {}, {})


@dataclasses.dataclass(frozen=True)
class OpPropertyCompositor(PropertyCompositor, typing.Generic[_T]):
    __slots__ = ("_calc",)

    _calc: OpPropertyCalc[_T]

    def __call__(self, rxn: interfaces.ReactionExplicit) -> "MetaPropertyState":
        calc = self._calc(rxn.operator)
        if calc is None:
            return MetaPropertyState({}, {}, {})
        props = {rxn.operator.item.uid: calc}
        single_state = MetaPropertyStateSingleProp(props, self._calc.resolver)
        return MetaPropertyState({}, {self._calc.key: single_state}, {})

    @property
    def keys(self) -> KeyOutput:
        return KeyOutput(frozenset(), frozenset((self._calc.key,)), frozenset())

    @property
    def meta_required(self) -> interfaces.MetaKeyPacket:
        return self._calc.meta_required

    @property
    def resolver(self) -> "MetaUpdateResolver":
        calc = self._calc
        return MetaUpdateResolver({}, {calc.key: calc.resolver}, {})


@dataclasses.dataclass(frozen=True)
class OpRxnPropertyCompositor(PropertyCompositor, typing.Generic[_T]):
    __slots__ = ("_calc",)

    _calc: OpPropertyFromRxnCalc[_T]

    def __call__(self, rxn: interfaces.ReactionExplicit) -> "MetaPropertyState":
        calc = self._calc(rxn.operator, rxn)
        if calc is None:
            return MetaPropertyState({}, {}, {})
        props = {rxn.operator.item.uid: calc}
        single_state = MetaPropertyStateSingleProp(props, self._calc.resolver)
        return MetaPropertyState({}, {self._calc.key: single_state}, {})

    @property
    def keys(self) -> KeyOutput:
        return KeyOutput(frozenset(), frozenset((self._calc.key,)), frozenset())

    @property
    def meta_required(self) -> interfaces.MetaKeyPacket:
        return self._calc.meta_required

    @property
    def resolver(self) -> "MetaUpdateResolver":
        calc = self._calc
        return MetaUpdateResolver({}, {calc.key: calc.resolver}, {})


@dataclasses.dataclass(frozen=True)
class RxnPropertyCompositor(PropertyCompositor, typing.Generic[_T]):
    __slots__ = ("_calc",)

    _calc: RxnPropertyCalc[_T]

    def __call__(self, rxn: interfaces.ReactionExplicit) -> "MetaPropertyState":
        calc = self._calc(rxn)
        if calc is None:
            return MetaPropertyState({}, {}, {})
        props = {rxn.uid: calc}
        single_state = MetaPropertyStateSingleProp(
            props,
            self._calc.resolver,  # type: ignore
        )
        return MetaPropertyState({}, {}, {self._calc.key: single_state})

    @property
    def keys(self) -> KeyOutput:
        return KeyOutput(frozenset(), frozenset(), frozenset((self._calc.key,)))

    @property
    def meta_required(self) -> interfaces.MetaKeyPacket:
        return self._calc.meta_required

    @property
    def resolver(self) -> "MetaUpdateResolver":
        calc = self._calc
        return MetaUpdateResolver({}, {}, {calc.key: calc.resolver})


@dataclasses.dataclass
class MetaPropertyStateSingleProp(typing.Generic[_T]):
    # note: this class is intended to be mutable; it will change after merging!
    __slots__ = ("data", "resolver")

    data: dict[interfaces.Identifier, _T]
    resolver: MetaDataResolverFunc[_T]

    def overwrite(
        self, other: "MetaPropertyStateSingleProp[_T]"
    ) -> "MetaPropertyStateSingleProp[_T]":
        self.data.update(other.data)
        self.resolver = other.resolver
        return self

    def __or__(
        self, other: "MetaPropertyStateSingleProp[_T]"
    ) -> "MetaPropertyStateSingleProp[_T]":
        resolved_props: dict[interfaces.Identifier, _T] = {}
        if len(self.data) == 0:
            return other
        if len(other.data) == 0:
            self.resolver = other.resolver
            return self
        common_keys = self.data.keys() & other.data.keys()
        if len(common_keys) == 0:
            other.data.update(self.data)
            return other
        for item_key in common_keys:
            resolved_props[item_key] = self.resolver(
                self.data[item_key], other.data[item_key]
            )
        other.data.update(self.data)
        other.data.update(resolved_props)
        return other


@dataclasses.dataclass
class MetaPropertyState:
    # note: this class is intended to be mutable; it will change after merging!
    __slots__ = ("mol_info", "op_info", "rxn_info")

    mol_info: dict[collections.abc.Hashable, MetaPropertyStateSingleProp]
    op_info: dict[collections.abc.Hashable, MetaPropertyStateSingleProp]
    rxn_info: dict[collections.abc.Hashable, MetaPropertyStateSingleProp]

    def __or__(self, other: "MetaPropertyState") -> "MetaPropertyState":
        if len(self.mol_info) == 0:
            self.mol_info = other.mol_info
        elif len(other.mol_info) != 0:
            common_keys = self.mol_info.keys() & other.mol_info.keys()
            if len(common_keys) == 0:
                self.mol_info.update(other.mol_info)
            else:
                resolved_info: dict[collections.abc.Hashable, typing.Any] = {}
                for prop_key in common_keys:
                    resolved_info[prop_key] = (
                        self.mol_info[prop_key] | other.mol_info[prop_key]
                    )
                self.mol_info.update(other.mol_info)
                self.mol_info.update(resolved_info)
        if len(self.op_info) == 0:
            self.op_info = other.op_info
        elif len(other.op_info) != 0:
            common_keys = self.op_info.keys() & other.op_info.keys()
            if len(common_keys) == 0:
                self.op_info.update(other.op_info)
            else:
                resolved_info = {}
                for prop_key in common_keys:
                    resolved_info[prop_key] = (
                        self.op_info[prop_key] | other.op_info[prop_key]
                    )
                self.op_info.update(other.op_info)
                self.op_info.update(resolved_info)
        if len(self.rxn_info) == 0:
            self.rxn_info = other.rxn_info
        elif len(other.rxn_info) != 0:
            common_keys = self.rxn_info.keys() & other.rxn_info.keys()
            if len(common_keys) == 0:
                self.rxn_info.update(other.rxn_info)
            else:
                resolved_info = {}
                for prop_key in common_keys:
                    resolved_info[prop_key] = (
                        self.rxn_info[prop_key] | other.rxn_info[prop_key]
                    )
                self.rxn_info.update(other.rxn_info)
                self.rxn_info.update(resolved_info)

        return self


class RxnAnalysisStep(abc.ABC):
    @abc.abstractmethod
    def execute(
        self,
        rxns: collections.abc.Iterable[
            tuple[interfaces.ReactionExplicit, bool]
        ],
    ) -> collections.abc.Iterable[tuple[interfaces.ReactionExplicit, bool]]: ...

    @typing.final
    def __rshift__(
        self,
        other: typing.Union[
            "RxnAnalysisStep",
            PropertyCompositor,
            ReactionFilterBase,
            LocalPropertyCalc,
        ],
    ) -> "RxnAnalysisStep":
        if isinstance(self, RxnAnalysisStep) and isinstance(
            other, RxnAnalysisStep
        ):
            return RxnAnalysisStepCompound(self, other)
        return as_rxn_analysis_step(self) >> as_rxn_analysis_step(other)

    @property
    @abc.abstractmethod
    def meta_required(self) -> interfaces.MetaKeyPacket: ...

    @property
    @abc.abstractmethod
    def resolver(self) -> "MetaUpdateResolver": ...


@dataclasses.dataclass(frozen=True)
class RxnAnalysisStepCompound(RxnAnalysisStep):
    __slots__ = ("step1", "step2")

    step1: RxnAnalysisStep
    step2: RxnAnalysisStep

    def execute(
        self,
        rxns: collections.abc.Iterable[
            tuple[interfaces.ReactionExplicit, bool]
        ],
    ) -> collections.abc.Iterable[tuple[interfaces.ReactionExplicit, bool]]:
        return self.step2.execute(
            rxn_out for rxn_out in self.step1.execute(rxns)
        )

    @property
    def meta_required(self) -> interfaces.MetaKeyPacket:
        return self.step1.meta_required + self.step2.meta_required

    @property
    def resolver(self) -> "MetaUpdateResolver":
        return self.step1.resolver | self.step2.resolver


def _mmd(
    i1: typing.Optional[
        collections.abc.Mapping[collections.abc.Hashable, typing.Any]
    ],
    i2: typing.Optional[
        collections.abc.Mapping[collections.abc.Hashable, typing.Any]
    ],
) -> typing.Optional[
    collections.abc.Mapping[collections.abc.Hashable, typing.Any]
]:
    if i1 is None:
        if i2 is None:
            return None
        return i2
    elif i2 is None:
        return i1
    if not isinstance(i1, dict):
        return dict(i1) | i2  # type: ignore [operator]
    return i1 | i2  # type: ignore [operator]


def metalib_to_rxn_meta(
    metalib: MetaPropertyState,
    rxns: collections.abc.Iterable[tuple[interfaces.ReactionExplicit, bool]],
) -> collections.abc.Iterable[tuple[interfaces.ReactionExplicit, bool]]:
    mol_info: dict[
        interfaces.Identifier, dict[collections.abc.Hashable, typing.Any]
    ] = {
        mol.item.uid: dict()
        for mol in itertools.chain(
            *(
                (itertools.chain(rxn[0].reactants, rxn[0].products))
                for rxn in rxns
            )
        )
    }
    for meta_key, mol_dict in metalib.mol_info.items():
        for mol_id, key_val in mol_dict.data.items():
            mol_info[mol_id][meta_key] = key_val
    op_info: dict[
        interfaces.Identifier, dict[collections.abc.Hashable, typing.Any]
    ] = {rxn[0].operator.item.uid: dict() for rxn in rxns}
    for meta_key, op_dict in metalib.op_info.items():
        for op_id, key_val in op_dict.data.items():
            op_info[op_id][meta_key] = key_val
    rxn_info: dict[
        interfaces.Identifier, dict[collections.abc.Hashable, typing.Any]
    ] = {rxn[0].uid: dict() for rxn in rxns}
    for meta_key, rxn_dict in metalib.rxn_info.items():
        for rxn_id, key_val in rxn_dict.data.items():
            rxn_info[rxn_id][meta_key] = key_val
    for rxn, passed_filter in rxns:
        yield (
            interfaces.ReactionExplicit(
                interfaces.DataPacketE(
                    rxn.operator.i,
                    rxn.operator.item,
                    _mmd(rxn.operator.meta, op_info[rxn.operator.item.uid]),
                ),
                tuple(
                    interfaces.DataPacketE(
                        mol.i, mol.item, _mmd(mol.meta, mol_info[mol.item.uid])
                    )
                    for mol in rxn.reactants
                ),
                tuple(
                    interfaces.DataPacketE(
                        mol.i, mol.item, _mmd(mol.meta, mol_info[mol.item.uid])
                    )
                    for mol in rxn.products
                ),
                _mmd(
                    rxn.reaction_meta,
                    rxn_info[rxn.uid],
                ),
            ),
            passed_filter,
        )


@dataclasses.dataclass(frozen=True)
class RxnAnalysisStepProp(RxnAnalysisStep):
    __slots__ = ("_prop",)

    _prop: PropertyCompositor

    def execute(
        self,
        rxns: collections.abc.Iterable[
            tuple[interfaces.ReactionExplicit, bool]
        ],
    ) -> collections.abc.Iterable[tuple[interfaces.ReactionExplicit, bool]]:
        rxn_list = list(rxns)
        meta_lib_generator = (self._prop(rxn[0]) for rxn in rxn_list if rxn[1])
        try:
            prop_map: MetaPropertyState = utils.logreduce(
                operator.or_, meta_lib_generator
            )
        except TypeError:
            prop_map = MetaPropertyState({}, {}, {})
        return metalib_to_rxn_meta(prop_map, rxn_list)

    @property
    def meta_required(self) -> interfaces.MetaKeyPacket:
        return self._prop.meta_required

    @property
    def resolver(self) -> "MetaUpdateResolver":
        return self._prop.resolver


class RxnAnalysisStepFilter(RxnAnalysisStep):
    __slots__ = ("_arg",)

    def __init__(self, arg: ReactionFilterBase) -> None:
        self._arg = arg

    def execute(
        self,
        rxns: collections.abc.Iterable[
            tuple[interfaces.ReactionExplicit, bool]
        ],
    ) -> collections.abc.Iterable[tuple[interfaces.ReactionExplicit, bool]]:
        for rxn, pass_filter in rxns:
            if not pass_filter:
                yield rxn, False
            else:
                yield rxn, self._arg(rxn)

    @property
    def meta_required(self) -> interfaces.MetaKeyPacket:
        return self._arg.meta_required

    @property
    def resolver(self) -> "MetaUpdateResolver":
        return MetaUpdateResolver({}, {}, {})


def _as_property_compositor(
    arg: typing.Union[PropertyCompositor, LocalPropertyCalc],
) -> PropertyCompositor:
    if isinstance(arg, PropertyCompositor):
        return arg
    elif isinstance(arg, MolPropertyCalc):
        return MolPropertyCompositor(arg)
    elif isinstance(arg, MolPropertyFromRxnCalc):
        return MolRxnPropertyCompositor(arg)
    elif isinstance(arg, OpPropertyCalc):
        return OpPropertyCompositor(arg)
    elif isinstance(arg, OpPropertyFromRxnCalc):
        return OpRxnPropertyCompositor(arg)
    elif isinstance(arg, RxnPropertyCalc):
        return RxnPropertyCompositor(arg)
    elif isinstance(arg, LocalPropertyCalc):
        raise TypeError(
            f"""Argument is of type {type(arg)}, relevant PropertyCompositor has
                not been provided"""
        )
    raise TypeError(
        f"""Argument is of type {type(arg)}, must be of type PropertyCompositor
            or LocalPropertyCalc"""
    )


def as_rxn_analysis_step(
    arg: typing.Union[
        RxnAnalysisStep,
        PropertyCompositor,
        ReactionFilterBase,
        LocalPropertyCalc,
    ],
) -> RxnAnalysisStep:
    if isinstance(arg, PropertyCompositor):
        return RxnAnalysisStepProp(arg)
    elif isinstance(arg, ReactionFilterBase):
        return RxnAnalysisStepFilter(arg)
    elif isinstance(arg, LocalPropertyCalc):
        return RxnAnalysisStepProp(_as_property_compositor(arg))
    elif isinstance(arg, RxnAnalysisStep):
        return arg
    raise TypeError(
        f"""Argument is of type {type(arg)}, must be of type RxnAnalysisStep,
            PropertyCompositor, ReactionFilterBase, or LocalPropertyCalc"""
    )


def _compose_property_function(
    func: collections.abc.Callable[[typing.Any, typing.Any], typing.Any],
    comp1: PropertyCompositor,
    comp2: PropertyCompositor,
) -> PropertyCompositor:
    return FunctionPropertyCompositor(
        func, _as_property_compositor(comp1), _as_property_compositor(comp2)
    )


def _merge_metas(
    x: collections.abc.Mapping[collections.abc.Hashable, MetaDataResolverFunc],
    y: collections.abc.Mapping[collections.abc.Hashable, MetaDataResolverFunc],
) -> collections.abc.Mapping[collections.abc.Hashable, MetaDataResolverFunc]:
    if isinstance(x, dict):
        x | y  # type: ignore [operator]
    return {**x, **y}


@typing.final
@dataclasses.dataclass(frozen=True)
class MetaUpdateResolver:
    __slots__ = ("mol_updates", "op_updates", "rxn_updates")

    mol_updates: collections.abc.Mapping[
        collections.abc.Hashable, MetaDataResolverFunc
    ]
    op_updates: collections.abc.Mapping[
        collections.abc.Hashable, MetaDataResolverFunc
    ]
    rxn_updates: collections.abc.Mapping[
        collections.abc.Hashable, MetaDataResolverFunc
    ]

    @typing.final
    def __or__(self, other: "MetaUpdateResolver") -> "MetaUpdateResolver":
        return MetaUpdateResolver(
            _merge_metas(self.mol_updates, other.mol_updates),
            _merge_metas(self.op_updates, other.op_updates),
            _merge_metas(self.rxn_updates, other.rxn_updates),
        )
