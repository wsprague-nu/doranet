from abc import ABC, abstractmethod
from dataclasses import dataclass


@dataclass(frozen=True)
class MetaKeyPacket:
    operator_keys: frozenset = frozenset()
    molecule_keys: frozenset = frozenset()

    def __add__(self, other: "MetaKeyPacket") -> "MetaKeyPacket":
        return MetaKeyPacket(
            self.operator_keys.union(other.operator_keys),
            self.molecule_keys.union(other.molecule_keys),
        )


class MetaSink(ABC):
    @abstractmethod
    def meta_required(self) -> MetaKeyPacket:
        ...
