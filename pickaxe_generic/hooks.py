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
