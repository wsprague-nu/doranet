"""Interfaces for unique identifiers."""

import collections.abc
import typing


# Identifier = str
class Identifier(collections.abc.Hashable, typing.Protocol):
    """
    Unique identifier/name of object.

    Should be orderable, hashable, and immutable.

    Methods
    -------
    __hash__:
        Hash function for object.
    __eq__:
        Equality function for object.
    __lt__:
        Ordering function for object.
    """

    def __hash__(self) -> int:
        """
        Hashes object to integer.

        Returns
        -------
        int
            Integer representing hashed value of object.  Should be
            (effectively) unique.
        """
        ...

    def __eq__(self, value: object, /) -> bool:
        """
        Compare object to others of similar type.

        If x == y, then x is equivalent to y.  This enables hashtables.

        Parameters
        ----------
        value
            Object to be compared.

        Returns
        -------
        bool
            True if object is equivalent to value, False otherwise.
        """
        ...

    def __lt__(self, other, /) -> bool:
        """
        Compare object to others of similar type.

        If canonical, this functionality allows canonical sorting.

        Parameters
        ----------
        other
            Object to be compared.

        Returns
        -------
        bool
            True if other is after self when ordered, False otherwise.
        """
        ...
