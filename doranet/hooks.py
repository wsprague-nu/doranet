"""Contains classes which implement global hook functions."""

import dataclasses

from doranet import interfaces


@dataclasses.dataclass(slots=True)
class NumberIterCondition(interfaces.GlobalUpdateHook):
    _num_iter: int

    def __call__(
        self, network: interfaces.ChemNetwork
    ) -> interfaces.GlobalHookReturnValue:
        self._num_iter -= 1
        if self._num_iter <= 0:
            return interfaces.GlobalHookReturnValue.STOP
        return interfaces.GlobalHookReturnValue.CONTINUE


@dataclasses.dataclass(frozen=True, slots=True)
class MaxMoleculesCondition(interfaces.GlobalUpdateHook):
    _max_mols: int

    def __call__(
        self, network: interfaces.ChemNetwork
    ) -> interfaces.GlobalHookReturnValue:
        if len(network.mols) > self._max_mols:
            return interfaces.GlobalHookReturnValue.STOP
        return interfaces.GlobalHookReturnValue.CONTINUE


@dataclasses.dataclass(frozen=True, slots=True)
class TargetMoleculeCondition(interfaces.GlobalUpdateHook):
    _target_mol: interfaces.MolDatBase

    def __call__(
        self, network: interfaces.ChemNetwork
    ) -> interfaces.GlobalHookReturnValue:
        if (
            self._target_mol in network.mols
            and network.reactivity[network.mols.i(self._target_mol.uid)]
        ):
            return interfaces.GlobalHookReturnValue.STOP
        return interfaces.GlobalHookReturnValue.CONTINUE
