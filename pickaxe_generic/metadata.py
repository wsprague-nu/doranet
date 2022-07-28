from abc import ABC, abstractmethod
from dataclasses import dataclass

from pickaxe_generic.datatypes import MetaKeyPacket


class MetaSink(ABC):
    @abstractmethod
    def meta_required(self) -> MetaKeyPacket:
        ...
