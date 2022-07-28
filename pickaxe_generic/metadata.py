from abc import ABC, abstractmethod
from collections.abc import Collection, Mapping
from functools import reduce
from operator import add, or_
from typing import Optional, Protocol, final

from pickaxe_generic.datatypes import MetaKeyPacket
from pickaxe_generic.network import ChemNetwork, ReactionExplicit


class MetaSink(Protocol):
    __slots__ = ()

    @abstractmethod
    @property
    def meta_required(self) -> MetaKeyPacket:
        ...


class LocalPropertyCalc(MetaSink, ABC):
    __slots__ = ()

    @abstractmethod
    def __call__(self, reaction: ReactionExplicit) -> Mapping:
        ...

    @final
    def __add__(self, other: "LocalPropertyCalc") -> "LocalPropertyCalc":
        return LocalPropertyCalcMultiple((self, other))


@final
class LocalPropertyCalcMultiple(LocalPropertyCalc):
    __slots__ = "_composed_properties"
    _composed_properties: tuple[LocalPropertyCalc, ...]

    def __init__(self, prop_list: Collection[LocalPropertyCalc]) -> None:
        props: list[LocalPropertyCalc] = []
        for prop in prop_list:
            if isinstance(prop, LocalPropertyCalcMultiple):
                props.extend(prop._composed_properties)
            else:
                props.append(prop)
        self._composed_properties = tuple(props)

    def __call__(self, reaction: ReactionExplicit) -> dict:
        return reduce(
            or_,
            (dict(prop(reaction)) for prop in self._composed_properties),
        )

    @property
    def meta_required(self) -> MetaKeyPacket:
        return reduce(
            add, (prop.meta_required for prop in self._composed_properties)
        )


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


@final
class MetaUpdateDefault(MetaUpdate):
    __slots__ = ()

    def __call__(
        self, reaction: ReactionExplicit, network: ChemNetwork
    ) -> None:
        op_meta = reaction.operator_meta
        if op_meta is not None:
            op_index = network.ops.i(reaction.operator.uid)
            for key, value in op_meta:
                network.op_meta(op_index, key, value)
        react_meta = reaction.reactants_meta
        if react_meta is not None:
            ...
