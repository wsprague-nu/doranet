"""
Contains classes which define and implement relevant smart containers.

Classes:

    interfaces.ObjectLibrary (deprecated)
      ObjectLibraryBasic*
      ObjectLibraryKeyVal*
"""

import collections.abc
import typing

from . import interfaces


@typing.final
class ObjectLibraryBasic(
    interfaces.ObjectLibrary, typing.Generic[interfaces.T_data]
):
    """
    Minimal class implementing the ObjectLibrary interface.

    This class wraps a dict.

    Parameters
    ----------
    objects : Optional[Iterable[DataUnit]]
        Objects which are to be included in the object library.
    """

    __slots__ = ("_lookup",)

    _lookup: dict[interfaces.Identifier, interfaces.T_data]

    def __init__(
        self,
        objects: typing.Optional[
            typing.Union[
                collections.abc.Iterable[interfaces.T_data], interfaces.T_data
            ]
        ] = None,
    ) -> None:
        self._lookup = {}
        if isinstance(objects, collections.abc.Iterable):
            for item in objects:
                self._lookup[item.uid] = item
        elif objects is not None:
            self._lookup[objects.uid] = objects

    def add(
        self,
        obj: typing.Union[
            collections.abc.Iterable[interfaces.T_data], interfaces.T_data
        ],
    ) -> None:
        if isinstance(obj, collections.abc.Iterable):
            for item in obj:
                self._lookup[item.uid] = item
        else:
            self._lookup[obj.uid] = obj

    def ids(
        self,
    ) -> collections.abc.Generator[interfaces.Identifier, None, None]:
        return (key for key in self._lookup)

    def __contains__(self, item: interfaces.T_data) -> bool:
        return item.uid in self._lookup

    def __getitem__(self, item: interfaces.Identifier) -> interfaces.T_data:
        return self._lookup[item]

    def __iter__(self) -> collections.abc.Iterator[interfaces.T_data]:
        return iter(self._lookup.values())

    def __len__(self) -> int:
        return len(self._lookup)


@typing.final
class ObjectLibraryKeyVal(
    interfaces.ObjectLibrary, typing.Generic[interfaces.T_data]
):
    """
    Minimal class implementing the ObjectLibrary interface.

    This class stores binary representation and initializes using passed
    function.

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

    _lookup: dict[interfaces.Identifier, bytes]

    def __init__(
        self,
        objects: typing.Optional[
            typing.Union[
                collections.abc.Iterable[interfaces.T_data], interfaces.T_data
            ]
        ] = None,
        initializer: typing.Optional[
            collections.abc.Callable[[bytes], interfaces.T_data]
        ] = None,
    ) -> None:
        if initializer is None:
            raise ValueError("initializer must be specified")
        self._lookup = {}
        self._initializer = initializer
        if isinstance(objects, collections.abc.Iterable):
            for item in objects:
                self._lookup[item.uid] = item.blob
        elif objects is not None:
            self._lookup[objects.uid] = objects.blob

    def add(
        self,
        obj: typing.Union[
            collections.abc.Iterable[interfaces.T_data], interfaces.T_data
        ],
    ) -> None:
        if isinstance(obj, collections.abc.Iterable):
            for item in obj:
                self._lookup[item.uid] = item.blob
        else:
            self._lookup[obj.uid] = obj.blob

    def ids(
        self,
    ) -> collections.abc.Generator[interfaces.Identifier, None, None]:
        return (key for key in self._lookup)

    def __contains__(self, item: interfaces.T_data) -> bool:
        return item.uid in self._lookup

    def __getitem__(self, item: interfaces.Identifier) -> interfaces.T_data:
        return self._initializer(self._lookup[item])

    def __iter__(self) -> collections.abc.Iterator[interfaces.T_data]:
        return (self._initializer(value) for value in self._lookup.values())

    def __len__(self) -> int:
        return len(self._lookup)
