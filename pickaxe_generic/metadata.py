from abc import ABC, abstractmethod
from collections.abc import Collection, Mapping
from functools import reduce
from operator import add, or_
from typing import Protocol, final

from pickaxe_generic.datatypes import MetaKeyPacket
from pickaxe_generic.network import ReactionExplicit


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
