"""Essential Operator object type."""

import abc
import typing
from collections.abc import Iterable

from ._dataunit import DataUnit


class Operator(DataUnit):
    """
    Interface representing an operator data object.

    Classes implementing this interface manage information about a single
    operator which acts on a `DataUnit` and can generate iterables of
    product `DataUnit`s.

    Attributes
    ----------
    blob : bytes
        Binary representation of operator.
    uid : doranet.interfaces.Identifier
        Unique identifier of operator.

    Methods
    -------
    __call__:
        Generate reactions from molecules.
    __len__:
        Number of arguments in operator.
    compat:
        Check compatibility of molecules with operator argument.
    """

    @abc.abstractmethod
    def __call__(self, *reactants: DataUnit) -> Iterable[Iterable[DataUnit]]:
        """
        React a sequence of `DataUnit` objects using internal operator.

        For every combination of reaction sites which is possible, there is a
        set of product molecules. This method returns an iterator over each of
        these sets, which may not be unique.

        Parameters
        ----------
        *reactants : DataUnit
            Reactants which match the arguments in the operator.

        Returns
        -------
        Iterable[Iterable[DataUnit]]
            Iterable of reaction product sets.
        """

    @abc.abstractmethod
    def num_args(self) -> int:
        """
        Return number of arguments in operator.

        Returns
        -------
        int
            Number of arguments in operator.
        """

    @typing.final
    def __len__(self) -> int:
        """
        Return number of arguments in operator.

        Returns
        -------
        int
            Number of arguments in operator.
        """
        return self.num_args()
