from abc import ABC, abstractmethod

from pickaxe_generic.datatypes import MetaKeyPacket


class MetaSink(ABC):
    __slots__ = ()

    @abstractmethod
    @property
    def meta_required(self) -> MetaKeyPacket:
        ...
