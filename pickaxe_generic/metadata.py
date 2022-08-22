import operator
from abc import ABC, abstractmethod
from collections.abc import Collection, Mapping
from dataclasses import dataclass
from email.generator import Generator
from functools import reduce
from itertools import chain
from multiprocessing.sharedctypes import Value
from operator import add, or_
from typing import (
    Any,
    Callable,
    Generic,
    Hashable,
    Iterable,
    Optional,
    Protocol,
    TypeVar,
    Union,
    final,
)

from pickaxe_generic.datatypes import (
    DataPacket,
    DataPacketE,
    DataUnit,
    Identifier,
    MetaKeyPacket,
    MolDatBase,
    OpDatBase,
)
from pickaxe_generic.network import ChemNetwork, ReactionExplicit
from pickaxe_generic.utils import logreduce


class MetaSink(Protocol):
    __slots__ = ()

    @property
    @abstractmethod
    def meta_required(self) -> MetaKeyPacket:
        ...


_T = TypeVar("_T")
_U = TypeVar("_U")
MetaDataResolverFunc = Callable[[_T, _T], _T]


class LocalPropertyCalc(ABC, Generic[_T]):
    @property
    @abstractmethod
    def key(self) -> Hashable:
        ...

    @property
    @abstractmethod
    def meta_required(self) -> MetaKeyPacket:
        ...

    @property
    @abstractmethod
    def resolver(self) -> MetaDataResolverFunc[_T]:
        ...

    @final
    def __and__(
        self, other: Union["PropertyCompositor", "LocalPropertyCalc"]
    ) -> "PropertyCompositor":
        return MergePropertyCompositor(
            _as_property_compositor(self), _as_property_compositor(other)
        )

    @final
    def __add__(
        self, other: Union["PropertyCompositor", "LocalPropertyCalc"]
    ) -> "PropertyCompositor":
        return FunctionPropertyCompositor(
            operator.add,
            _as_property_compositor(self),
            _as_property_compositor(other),
        )

    @final
    def __sub__(
        self, other: Union["PropertyCompositor", "LocalPropertyCalc"]
    ) -> "PropertyCompositor":
        return FunctionPropertyCompositor(
            operator.sub,
            _as_property_compositor(self),
            _as_property_compositor(other),
        )

    @final
    def __mul__(
        self, other: Union["PropertyCompositor", "LocalPropertyCalc"]
    ) -> "PropertyCompositor":
        return FunctionPropertyCompositor(
            operator.mul,
            _as_property_compositor(self),
            _as_property_compositor(other),
        )

    @final
    def __truediv__(
        self, other: Union["PropertyCompositor", "LocalPropertyCalc"]
    ) -> "PropertyCompositor":
        return FunctionPropertyCompositor(
            operator.truediv,
            _as_property_compositor(self),
            _as_property_compositor(other),
        )

    @final
    def __pow__(
        self, other: Union["PropertyCompositor", "LocalPropertyCalc"]
    ) -> "PropertyCompositor":
        return FunctionPropertyCompositor(
            operator.pow,
            _as_property_compositor(self),
            _as_property_compositor(other),
        )

    @final
    def __rshift__(
        self,
        other: Union[
            "RxnAnalysisStep",
            "PropertyCompositor",
            "ReactionFilterBase",
            "LocalPropertyCalc",
        ],
    ) -> "RxnAnalysisStep":
        return as_rxn_analysis_step(self) >> as_rxn_analysis_step(other)


class MolPropertyCalc(LocalPropertyCalc[_T]):
    @abstractmethod
    def __call__(
        self, data: DataPacketE[MolDatBase], prev_value: Optional[_T] = None
    ) -> Optional[_T]:
        ...


class MolPropertyFromRxnCalc(LocalPropertyCalc[_T]):
    @abstractmethod
    def __call__(
        self,
        data: DataPacketE[MolDatBase],
        rxn: ReactionExplicit,
        prev_value: Optional[_T] = None,
    ) -> Optional[_T]:
        ...


class OpPropertyCalc(LocalPropertyCalc[_T]):
    @abstractmethod
    def __call__(
        self, data: DataPacketE[OpDatBase], prev_value: Optional[_T] = None
    ) -> Optional[_T]:
        ...


class OpPropertyFromRxnCalc(LocalPropertyCalc[_T]):
    @abstractmethod
    def __call__(
        self,
        data: DataPacketE[OpDatBase],
        rxn: ReactionExplicit,
        prev_value: Optional[_T] = None,
    ) -> Optional[_T]:
        ...


class RxnPropertyCalc(LocalPropertyCalc[_T]):
    @abstractmethod
    def __call__(
        self, data: ReactionExplicit, prev_value: Optional[_T] = None
    ) -> Optional[_T]:
        ...


@dataclass(frozen=True)
class KeyOutput:
    mol_keys: frozenset[Hashable]
    op_keys: frozenset[Hashable]
    rxn_keys: frozenset[Hashable]

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
                f"Conflicting molecule metadata key outputs {self.mol_keys & other.mol_keys}; separate expressions with >> or combine using other operator"
            )
        new_op_keys = self.op_keys | other.op_keys
        if len(new_op_keys) > len(self.op_keys) + len(other.op_keys):
            raise KeyError(
                f"Conflicting operator metadata key outputs {self.op_keys & other.op_keys}; separate expressions with >> or combine using other operator"
            )
        new_rxn_keys = self.rxn_keys | other.rxn_keys
        if len(new_rxn_keys) > len(self.rxn_keys) + len(other.rxn_keys):
            raise KeyError(
                f"Conflicting reaction metadata key outputs {self.rxn_keys & other.rxn_keys}; separate expressions with >> or combine using other operator"
            )
        return KeyOutput(new_mol_keys, new_op_keys, new_rxn_keys)


class ReactionFilterBase(ABC):
    __slots__ = ()

    @abstractmethod
    def __call__(self, recipe: ReactionExplicit) -> bool:
        ...

    @property
    def meta_required(self) -> MetaKeyPacket:
        return MetaKeyPacket()

    @final
    def __and__(self, other: "ReactionFilterBase") -> "ReactionFilterBase":
        return ReactionFilterAnd(self, other)

    @final
    def __invert__(self) -> "ReactionFilterBase":
        return ReactionFilterInv(self)

    @final
    def __or__(self, other: "ReactionFilterBase") -> "ReactionFilterBase":
        return ReactionFilterOr(self, other)

    @final
    def __xor__(self, other: "ReactionFilterBase") -> "ReactionFilterBase":
        return ReactionFilterXor(self, other)

    @final
    def __rshift__(
        self,
        other: Union[
            "RxnAnalysisStep",
            "PropertyCompositor",
            "ReactionFilterBase",
            LocalPropertyCalc,
        ],
    ) -> "RxnAnalysisStep":
        return as_rxn_analysis_step(self) >> as_rxn_analysis_step(other)


@dataclass(frozen=True)
class ReactionFilterAnd(ReactionFilterBase):
    __slots__ = ("_filter1", "_filter2")

    _filter1: ReactionFilterBase
    _filter2: ReactionFilterBase

    def __call__(self, recipe: ReactionExplicit) -> bool:
        return self._filter1(recipe) and self._filter2(recipe)

    @property
    def meta_required(self) -> MetaKeyPacket:
        return self._filter1.meta_required + self._filter2.meta_required


@dataclass(frozen=True)
class ReactionFilterInv(ReactionFilterBase):
    __slots__ = ("_filter",)
    _filter: ReactionFilterBase

    def __call__(self, recipe: ReactionExplicit) -> bool:
        return not self._filter(recipe)

    @property
    def meta_required(self) -> MetaKeyPacket:
        return self._filter.meta_required


@dataclass(frozen=True)
class ReactionFilterOr(ReactionFilterBase):
    __slots__ = ("_filter1", "_filter2")
    _filter1: ReactionFilterBase
    _filter2: ReactionFilterBase

    def __call__(self, recipe: ReactionExplicit) -> bool:
        return self._filter1(recipe) or self._filter2(recipe)

    @property
    def meta_required(self) -> MetaKeyPacket:
        return self._filter1.meta_required + self._filter2.meta_required


@dataclass(frozen=True)
class ReactionFilterXor(ReactionFilterBase):
    __slots__ = ("_filter1", "_filter2")
    _filter1: ReactionFilterBase
    _filter2: ReactionFilterBase

    def __call__(self, recipe: ReactionExplicit) -> bool:
        return self._filter1(recipe) != self._filter2(recipe)

    @property
    def meta_required(self) -> MetaKeyPacket:
        return self._filter1.meta_required + self._filter2.meta_required


class PropertyCompositor(ABC):
    @abstractmethod
    def __call__(self, rxn: ReactionExplicit) -> "MetaPropertyState":
        ...

    @property
    @abstractmethod
    def keys(self) -> KeyOutput:
        ...

    @final
    def __and__(
        self, other: Union["PropertyCompositor", LocalPropertyCalc]
    ) -> "PropertyCompositor":
        return MergePropertyCompositor(self, _as_property_compositor(other))

    @final
    def __add__(
        self, other: Union["PropertyCompositor", LocalPropertyCalc]
    ) -> "PropertyCompositor":
        return FunctionPropertyCompositor(
            operator.add, self, _as_property_compositor(other)
        )

    @final
    def __sub__(
        self, other: Union["PropertyCompositor", LocalPropertyCalc]
    ) -> "PropertyCompositor":
        return FunctionPropertyCompositor(
            operator.sub, self, _as_property_compositor(other)
        )

    @final
    def __mul__(
        self, other: Union["PropertyCompositor", LocalPropertyCalc]
    ) -> "PropertyCompositor":
        return FunctionPropertyCompositor(
            operator.mul, self, _as_property_compositor(other)
        )

    @final
    def __truediv__(
        self, other: Union["PropertyCompositor", LocalPropertyCalc]
    ) -> "PropertyCompositor":
        return FunctionPropertyCompositor(
            operator.truediv, self, _as_property_compositor(other)
        )

    @final
    def __pow__(
        self, other: Union["PropertyCompositor", LocalPropertyCalc]
    ) -> "PropertyCompositor":
        return FunctionPropertyCompositor(
            operator.pow, self, _as_property_compositor(other)
        )

    @final
    def __rshift__(
        self,
        other: Union[
            "RxnAnalysisStep",
            "PropertyCompositor",
            ReactionFilterBase,
            LocalPropertyCalc,
        ],
    ) -> "RxnAnalysisStep":
        return as_rxn_analysis_step(self) >> as_rxn_analysis_step(other)

    @property
    @abstractmethod
    def meta_required(self) -> MetaKeyPacket:
        ...


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

    def __call__(self, rxn: ReactionExplicit) -> "MetaPropertyState":
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
    def meta_required(self) -> MetaKeyPacket:
        return self._comp1.meta_required + self._comp2.meta_required


@dataclass(frozen=True)
class FunctionPropertyCompositor(PropertyCompositor):
    __slots__ = ("_func", "_comp1", "_comp2")

    _func: Callable[[Any, Any], Any]
    _comp1: PropertyCompositor
    _comp2: PropertyCompositor

    def __call__(self, rxn: ReactionExplicit) -> "MetaPropertyState":
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
    def meta_required(self) -> MetaKeyPacket:
        return self._comp1.meta_required + self._comp2.meta_required


@dataclass(frozen=True)
class MolPropertyCompositor(PropertyCompositor, Generic[_T]):
    __slots__ = ("_calc",)

    _calc: MolPropertyCalc[_T]

    def __call__(self, rxn: ReactionExplicit) -> "MetaPropertyState":
        mols = chain(rxn.reactants, rxn.products)
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
    def meta_required(self) -> MetaKeyPacket:
        return self._calc.meta_required


@dataclass(frozen=True)
class MolRxnPropertyCompositor(PropertyCompositor, Generic[_T]):
    __slots__ = ("_calc",)

    _calc: MolPropertyFromRxnCalc[_T]

    def __call__(self, rxn: ReactionExplicit) -> "MetaPropertyState":
        mols = chain(rxn.reactants, rxn.products)
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
    def meta_required(self) -> MetaKeyPacket:
        return self._calc.meta_required


@dataclass(frozen=True)
class OpPropertyCompositor(PropertyCompositor, Generic[_T]):
    __slots__ = ("_calc",)

    _calc: OpPropertyCalc[_T]

    def __call__(self, rxn: ReactionExplicit) -> "MetaPropertyState":
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
    def meta_required(self) -> MetaKeyPacket:
        return self._calc.meta_required


@dataclass(frozen=True)
class OpRxnPropertyCompositor(PropertyCompositor, Generic[_T]):
    __slots__ = ("_calc",)

    _calc: OpPropertyFromRxnCalc[_T]

    def __call__(self, rxn: ReactionExplicit) -> "MetaPropertyState":
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
    def meta_required(self) -> MetaKeyPacket:
        return self._calc.meta_required


@dataclass(frozen=True)
class RxnPropertyCompositor(PropertyCompositor, Generic[_T]):
    __slots__ = ("_calc",)

    _calc: RxnPropertyCalc[_T]

    def __call__(self, rxn: ReactionExplicit) -> "MetaPropertyState":
        calc = self._calc(rxn)
        if calc is None:
            return MetaPropertyState({}, {}, {})
        props = {rxn.uid: calc}
        single_state = MetaPropertyStateSingleProp(
            props, self._calc.resolver  # type: ignore
        )
        return MetaPropertyState({}, {}, {self._calc.key: single_state})

    @property
    def keys(self) -> KeyOutput:
        return KeyOutput(frozenset(), frozenset(), frozenset((self._calc.key,)))

    @property
    def meta_required(self) -> MetaKeyPacket:
        return self._calc.meta_required


@dataclass
class MetaPropertyStateSingleProp(Generic[_T]):
    # note: this class is intended to be mutable; it will change after merging!
    __slots__ = ("data", "resolver")

    data: dict[Identifier, _T]
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
        resolved_props: dict[Identifier, _T] = {}
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


@dataclass
class MetaPropertyState:
    # note: this class is intended to be mutable; it will change after merging!
    __slots__ = ("mol_info", "op_info", "rxn_info")

    mol_info: dict[Hashable, MetaPropertyStateSingleProp]
    op_info: dict[Hashable, MetaPropertyStateSingleProp]
    rxn_info: dict[Hashable, MetaPropertyStateSingleProp]

    def __or__(self, other: "MetaPropertyState") -> "MetaPropertyState":
        if len(self.mol_info) == 0:
            self.mol_info = other.mol_info
        elif len(other.mol_info) != 0:
            common_keys = self.mol_info.keys() & other.mol_info.keys()
            if len(common_keys) == 0:
                self.mol_info.update(other.mol_info)
            else:
                resolved_info: dict[Hashable, Any] = {}
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


class RxnAnalysisStep(ABC):
    @abstractmethod
    def execute(
        self, rxns: Iterable[tuple[ReactionExplicit, bool]]
    ) -> Iterable[tuple[ReactionExplicit, bool]]:
        ...

    @final
    def __rshift__(
        self,
        other: Union[
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

    @abstractmethod
    @property
    def meta_required(self) -> MetaKeyPacket:
        ...


@dataclass(frozen=True)
class RxnAnalysisStepCompound(RxnAnalysisStep):
    __slots__ = ("step1", "step2")

    step1: RxnAnalysisStep
    step2: RxnAnalysisStep

    def execute(
        self, rxns: Iterable[tuple[ReactionExplicit, bool]]
    ) -> Iterable[tuple[ReactionExplicit, bool]]:
        return self.step2.execute(
            rxn_out for rxn_out in self.step1.execute(rxns)
        )

    @property
    def meta_required(self) -> MetaKeyPacket:
        return self.step1.meta_required + self.step2.meta_required


def _mmd(
    i1: Optional[Mapping[Hashable, Any]], i2: Optional[Mapping[Hashable, Any]]
) -> Optional[Mapping[Hashable, Any]]:
    if i1 is None:
        if i2 is None:
            return None
        return i2
    elif i2 is None:
        return i1
    if not isinstance(i1, dict):
        return dict(i1) | i2
    return i1 | i2


def metalib_to_rxn_meta(
    metalib: MetaPropertyState, rxns: Iterable[tuple[ReactionExplicit, bool]]
) -> Iterable[tuple[ReactionExplicit, bool]]:
    mol_info: dict[Identifier, dict[Hashable, Any]] = {}
    for meta_key, mol_dict in metalib.mol_info.items():
        for mol_id, key_val in mol_dict.data.items():
            mol_info[mol_id][meta_key] = key_val
    op_info: dict[Identifier, dict[Hashable, Any]] = {}
    for meta_key, op_dict in metalib.op_info.items():
        for op_id, key_val in op_dict.data.items():
            op_info[op_id][meta_key] = key_val
    rxn_info: dict[Identifier, dict[Hashable, Any]] = {}
    for meta_key, rxn_dict in metalib.rxn_info.items():
        for rxn_id, key_val in rxn_dict.data.items():
            rxn_info[rxn_id][meta_key] = key_val
    for rxn, passed_filter in rxns:
        yield (
            ReactionExplicit(
                rxn.operator,
                rxn.reactants,
                rxn.products,
                _mmd(
                    rxn.reaction_meta,
                    rxn_info[rxn.uid],
                ),
            ),
            passed_filter,
        )


@dataclass(frozen=True)
class RxnAnalysisStepProp(RxnAnalysisStep):
    __slots__ = ("_prop",)

    def __init__(self, arg: PropertyCompositor) -> None:
        self._prop = arg

    def execute(
        self, rxns: Iterable[tuple[ReactionExplicit, bool]]
    ) -> Iterable[tuple[ReactionExplicit, bool]]:
        rxn_list = list(rxns)
        meta_lib_generator = (self._prop(rxn[0]) for rxn in rxn_list if rxn[1])
        prop_map: MetaPropertyState = logreduce(
            operator.or_, meta_lib_generator
        )
        return metalib_to_rxn_meta(prop_map, rxn_list)

    @property
    def meta_required(self) -> MetaKeyPacket:
        return self._prop.meta_required


class RxnAnalysisStepFilter(RxnAnalysisStep):
    __slots__ = ("_arg",)

    def __init__(self, arg: ReactionFilterBase) -> None:
        self._arg = arg

    def execute(
        self, rxns: Iterable[tuple[ReactionExplicit, bool]]
    ) -> Iterable[tuple[ReactionExplicit, bool]]:
        for rxn, pass_filter in rxns:
            if not pass_filter or not self._arg(rxn):
                yield rxn, False
            else:
                yield rxn, True

    @property
    def meta_required(self) -> MetaKeyPacket:
        return self._arg.meta_required


class MetaUpdate(ABC):
    @abstractmethod
    def __call__(
        self, reaction: ReactionExplicit, network: ChemNetwork
    ) -> Optional[ReactionExplicit]:
        ...

    @final
    def __add__(self, other: "MetaUpdate") -> "MetaUpdate":
        return MetaUpdateMultiple((self, other))


@final
class MetaUpdateMultiple(MetaUpdate):
    __slots__ = "_composed_updates"
    _composed_updates: tuple[MetaUpdate, ...]

    def __init__(self, update_list: Collection[MetaUpdate]) -> None:
        updates: list[MetaUpdate] = []
        for update in update_list:
            if isinstance(update, MetaUpdateMultiple):
                updates.extend(update._composed_updates)
            else:
                updates.append(update)
        self._composed_updates = tuple(updates)

    def __call__(
        self, reaction: ReactionExplicit, network: ChemNetwork
    ) -> ReactionExplicit:
        cur_reaction = reaction
        for update in self._composed_updates:
            new_reaction = update(cur_reaction, network)
            if new_reaction is not None:
                cur_reaction = new_reaction
        return cur_reaction


def _as_property_compositor(
    arg: Union[PropertyCompositor, LocalPropertyCalc]
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
            f"Argument is of type {type(arg)}, relevant PropertyCompositor has not been provided"
        )
    raise TypeError(
        f"Argument is of type {type(arg)}, must be of type PropertyCompositor or LocalPropertyCalc"
    )


def as_rxn_analysis_step(
    arg: Union[
        RxnAnalysisStep,
        PropertyCompositor,
        ReactionFilterBase,
        LocalPropertyCalc,
    ]
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
        f"Argument is of type {type(arg)}, must be of type RxnAnalysisStep, PropertyCompositor, ReactionFilterBase, or LocalPropertyCalc"
    )


def _compose_property_function(
    func: Callable[[Any, Any], Any],
    comp1: PropertyCompositor,
    comp2: PropertyCompositor,
) -> PropertyCompositor:
    return FunctionPropertyCompositor(
        func, _as_property_compositor(comp1), _as_property_compositor(comp2)
    )
