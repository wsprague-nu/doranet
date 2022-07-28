from abc import ABC, abstractmethod

from pickaxe_generic.filters import MetaKeyPacket


class MetaSink(ABC):
    @abstractmethod
    def meta_required(self) -> MetaKeyPacket:
        ...
