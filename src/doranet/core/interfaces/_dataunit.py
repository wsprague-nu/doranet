import abc
import typing

from ._identifier import Identifier

if typing.TYPE_CHECKING:
    from doranet.core.interfaces import NetworkEngine


T_ci = typing.TypeVar("T_ci", contravariant=False)


class DataUnit(abc.ABC):
    """
    Basic data storage object.

    Object which provides a unique, hashable identifier, and can serve up a
    binary form of the object.

    Attributes
    ----------
    blob : bytes
        Binary representation of object.
    uid : pickaxe_generic.interfaces.Identifier
        Unique identifier of object.
    """

    __slots__ = ()

    @property
    @abc.abstractmethod
    def blob(self) -> bytes:
        """
        Binary representation of object.

        Must be able to initialize object when passed to __setstate__ method of
        any sibling class with same immediate parent (viz. initialize a
        MolDatBasicV2, even if obtained from a MolDatBasicV1).
        """

    @classmethod
    @abc.abstractmethod
    def from_bytes(
        cls: type[T_ci], data: bytes, engine: "NetworkEngine"
    ) -> T_ci:
        """
        Produce object from bytestring.

        Bytestring should be derived from the .blob property.
        """

    @property
    @abc.abstractmethod
    def uid(self) -> Identifier:
        """
        Return unique identifier of object.

        Must be hashable in order to facilitate lookup tables utilizing hashes.
        """
