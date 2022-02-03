"""
Contains classes which define and implement relevant smart containers.

Classes:

    ObjectLibrary
      ObjectLibraryBasic*
      ObjectLibraryKeyVal*
"""

from abc import ABC, abstractmethod
from typing import (
    Callable,
    Dict,
    final,
    Generator,
    Generic,
    Iterable,
    Iterator,
    Optional,
    TypeVar,
    Union,
)

from pickaxe_generic.datatypes import DataUnit, Identifier

DataUnitGen = TypeVar("DataUnitGen", bound=DataUnit)


class ObjectLibrary(ABC, Generic[DataUnitGen]):
    """
    Interface representing library of data.

    Classes implementing this interface manage multiple instances of a hashable
    object, and may have responsibility for synchronization with external
    databases which may also manage this information (be volatile).  Contained
    objects must have a "uid" attribute which contains a hashable unique id.

    Current implementations assume that this library will never shrink or remove
    entries.
    """

    __slots__ = ()

    @abstractmethod
    def add(self, obj: Union[Iterable[DataUnitGen], DataUnitGen]) -> None:
        """
        Adds an object or multiple objects to the library (as long as it is
        unique to the other objects in the library).

        Parameters
        ----------
        obj : Union[Iterable[DataUnit], DataUnit]
            Object(s) to be added.
        """

    @abstractmethod
    def ids(self) -> Iterable[Identifier]:
        """
        Returns a set of keys used in the library.

        Returns
        -------
        Iterable[Identifier]
            Ids of objects in the library.
        """
        pass

    @abstractmethod
    def __contains__(self, item: DataUnitGen) -> bool:
        pass

    @abstractmethod
    def __getitem__(self, id: Identifier) -> DataUnitGen:
        pass

    @abstractmethod
    def __iter__(self) -> Iterator[DataUnitGen]:
        pass

    @abstractmethod
    def __len__(self) -> int:
        pass


@final
class ObjectLibraryBasic(ObjectLibrary, Generic[DataUnitGen]):
    """
    Minimal class implementing the ObjectLibrary interface.

    This class wraps a dict.

    Parameters
    ----------
    objects : Optional[Iterable[DataUnit]]
        Objects which are to be included in the object library.
    """

    __slots__ = ("_lookup",)

    _lookup: Dict[Identifier, DataUnitGen]

    def __init__(
        self,
        objects: Optional[Union[Iterable[DataUnitGen], DataUnitGen]] = None,
    ) -> None:
        self._lookup = {}
        if isinstance(objects, Iterable):
            for item in objects:
                self._lookup[item.uid] = item
        elif objects is not None:
            self._lookup[objects.uid] = objects

    def add(self, obj: Union[Iterable[DataUnitGen], DataUnitGen]) -> None:
        if isinstance(obj, Iterable):
            for item in obj:
                self._lookup[item.uid] = item
        else:
            self._lookup[obj.uid] = obj

    def ids(self) -> Generator[Identifier, None, None]:
        return (key for key in self._lookup.keys())

    def __contains__(self, item: DataUnitGen) -> bool:
        return item.uid in self._lookup

    def __getitem__(self, id: Identifier) -> DataUnitGen:
        return self._lookup[id]

    def __iter__(self) -> Iterator[DataUnitGen]:
        return iter(self._lookup.values())

    def __len__(self) -> int:
        return len(self._lookup)


@final
class ObjectLibraryKeyVal(ObjectLibrary, Generic[DataUnitGen]):
    """
    Minimal class implementing the ObjectLibrary interface; stores binary
    representation and initializes using passed function.

    This class wraps a dict.

    Parameters
    ----------
    objects : Optional[Iterable[DataUnitGen]]
        Objects which are to be included in the object library.
    initializer : Callable[[bytes],DataUnitGen]
        Initializer which can convert bytes to relevant data type.
    """

    __slots__ = (
        "_initializer",
        "_lookup",
    )

    _lookup: Dict[Identifier, bytes]
    _initializer: Callable[[bytes], DataUnitGen]

    def __init__(
        self,
        objects: Optional[Union[Iterable[DataUnitGen], DataUnitGen]] = None,
        initializer: Optional[Callable[[bytes], DataUnitGen]] = None,
    ) -> None:
        if initializer is None:
            raise ValueError("initializer must be specified")
        self._lookup = {}
        self._initializer = initializer
        if isinstance(objects, Iterable):
            for item in objects:
                self._lookup[item.uid] = item.blob
        elif objects is not None:
            self._lookup[objects.uid] = objects.blob

    def add(self, obj: Union[Iterable[DataUnitGen], DataUnitGen]) -> None:
        if isinstance(obj, Iterable):
            for item in obj:
                self._lookup[item.uid] = item.blob
        else:
            self._lookup[obj.uid] = obj.blob

    def ids(self) -> Generator[Identifier, None, None]:
        return (key for key in self._lookup.keys())

    def __contains__(self, item: DataUnitGen) -> bool:
        return item.uid in self._lookup

    def __getitem__(self, id: Identifier) -> DataUnitGen:
        return self._initializer(self._lookup[id])

    def __iter__(self) -> Iterator[DataUnitGen]:
        return (self._initializer(value) for value in self._lookup.values())

    def __len__(self) -> int:
        return len(self._lookup)
