import abc
import typing
import warnings

from ._identifier import Identifier

if typing.TYPE_CHECKING:
    from doranet.core.interfaces import NetworkEngine


T_ci = typing.TypeVar("T_ci", contravariant=False)


class _DataUnitDeprecated(abc.ABC):
    """
    Basic data storage object (DEPRECATED).

    Object which provides a unique, hashable identifier, and can serve up a
    binary form of the object.

    .. deprecated:: 0.6.0
        `_DataUnitDeprecated` features will be removed in DORAnet 0.7.0, they
        are replaced by `DataUnit` native features.

    Attributes
    ----------
    blob : bytes
        Binary representation of object.
    uid : doranet.interfaces.Identifier
        Unique identifier of object.
    """

    @property
    @abc.abstractmethod
    def blob(self) -> bytes:
        """
        Binary representation of object.

        Must be able to initialize object when passed to __setstate__ method of
        any sibling class with same immediate parent (viz. initialize a
        MolDatBasicV2, even if obtained from a MolDatBasicV1).

        .. deprecated:: 0.6.0
          `blob` property will be removed in DORAnet 0.7.0, it is replaced by
          `to_bytes_fast()` to distinguish it from `from_bytes_small`.
        """

    @classmethod
    @abc.abstractmethod
    def from_bytes(
        cls: type[T_ci], data: bytes, engine: "NetworkEngine"
    ) -> T_ci:
        """
        Produce object from bytestring.

        Bytestring should be derived from the .blob property.

        .. deprecated:: 0.6.0
          `from_bytes` will be removed in DORAnet 0.7.0, it is replaced by
          `from_bytes_fast` to distinguish it from `from_bytes_small`.
        """

    @property
    @abc.abstractmethod
    def uid(self) -> Identifier:
        """
        Return unique identifier of object.

        Must be hashable in order to facilitate lookup tables utilizing hashes.

        .. deprecated:: 0.6.0
          `uid` property will be removed in DORAnet 0.7.0, it is replaced by
          `get_uid()` to clarify that it is a method, not an attribute.
        """


class DataUnit(_DataUnitDeprecated):
    """
    Basic data storage object.

    An implementation of DataUnit should provide methods which enable rapid
    storage and transmission of the contained data type. It should also be
    fully immutable (e.g., unable to be modified by outside processes),
    otherwise the assumptions of DORAnet do not hold.
    """

    @abc.abstractmethod
    def to_bytes_fast(self) -> bytes:
        """
        Serialize `self` to non-unique binary form, optimized for speed.

        This method is used for pickling/unpickling and is meant for local
        interprocess communication purposes (CPU bottleneck).

        The returned value should satisfy
        `self.__class__.from_bytes_fast(self.to_bytes_fast()) == self` and may
        require overriding `from_bytes_fast()` as well.
        """

    @classmethod
    @abc.abstractmethod
    def from_bytes_fast(cls, bstring: bytes) -> typing.Self:
        """
        Deserialize `Self` from speed-optimized binary form.

        This method is used for transmitting the `DataUnit` for I/O bottlenecked
        processes and size optimization, for example database communication.

        The returned value should satisfy
        `self.__class__.from_bytes_fast(self.to_bytes_fast()) == self` and may
        require overriding `to_bytes_fast()` as well.

        Parameters
        ----------
        bstring : bytes
            Binary string containing speed-optimized binary representation.
        """

    def to_bytes_small(self) -> bytes:
        """
        Serialize `self` to non-unique binary form, optimized for size.

        This method is used for transmitting the `DataUnit` for I/O bottlenecked
        processes and size optimization, for example database communication.

        The returned value should satisfy
        `self.__class__.from_bytes_small(self.to_bytes_small()) == self` and may
        require overriding `from_bytes_small()` as well.
        """
        return self.to_bytes_fast()

    @classmethod
    def from_bytes_small(cls, bstring: bytes) -> typing.Self:
        """
        Deserialize `Self` from size-optimized binary form.

        This method is used for transmitting the `DataUnit` for I/O bottlenecked
        processes and size optimization, for example database communication.

        The returned value should satisfy
        `self.__class__.from_bytes_small(self.to_bytes_small()) == self` and may
        require overriding `to_bytes_small()` as well.

        Parameters
        ----------
        bstring : bytes
            Binary string containing size-optimized binary representation.
        """
        return cls.from_bytes_fast(bstring)

    @abc.abstractmethod
    def get_uid(self) -> str:
        """
        Return unique identifier of object.

        Does not have to define object, only be unique with respect to other
        objects of the same class.
        """

    def __hash__(self) -> int:
        return hash(self.get_uid())

    def __repr__(self) -> str:
        return f"{type(self)}.from_bytes_fast({self.to_bytes_fast()!r})"

    def __str__(self) -> str:
        return self.get_uid()

    @typing.final
    def __eq__(self, other: object) -> bool:
        return (
            type(other) is type(self)
            and isinstance(other, type(self))  # probably overkill
            and self.get_uid() == other.get_uid()
        )

    @typing.final
    def __lt__(self, other: object) -> bool:
        if not isinstance(other, DataUnit):
            raise TypeError(
                f"'<' not supported between instances of '{type(self)}' and "
                f"'{type(other)}'"
            )
        return (type(self).__name__, self.get_uid()) < (
            type(other).__name__,
            other.get_uid(),
        )

    @typing.final
    @property
    def blob(self) -> bytes:
        warnings.warn(
            "DataUnit.blob is a deprecated function, please use"
            "self.to_bytes_fast() instead",
            DeprecationWarning,
            stacklevel=1,
        )
        return self.to_bytes_fast()

    @typing.final
    @property
    def uid(self) -> str:
        warnings.warn(
            "DataUnit.uid is a deprecated function, please use"
            "self.get_uid() instead",
            DeprecationWarning,
            stacklevel=1,
        )
        return self.get_uid()

    @typing.final
    @classmethod
    def from_bytes(cls, data: bytes, engine: "NetworkEngine") -> typing.Self:
        warnings.warn(
            "DataUnit.from_bytes() is a deprecated function, please use"
            "DataUnit.from_bytes_fast() instead",
            DeprecationWarning,
            stacklevel=1,
        )
        return cls.from_bytes_fast(data)
