"""
Contains classes which implement global hook functions.
"""

import dataclasses
import typing

from . import interfaces


@dataclasses.dataclass(slots=True)
class NumberIterCondition:
    _num_iter: int

    def __call__(
        self, network: interfaces.ChemNetwork
    ) -> typing.Optional[bool]:
        self._num_iter -= 1
        if self._num_iter <= 0:
            return False
        return None


@dataclasses.dataclass(frozen=True, slots=True)
class MaxMoleculesCondition:
    _max_mols: int

    def __call__(
        self, network: interfaces.ChemNetwork
    ) -> typing.Optional[bool]:
        if len(network.mols) > self._max_mols:
            return False
        return None
