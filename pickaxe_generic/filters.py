from abc import ABC, abstractmethod
from typing import Collection, Iterable, Sequence

from rdkit.Chem.rdqueries import AtomNumEqualsQueryAtom

from pickaxe_generic.datatypes import (
    Identifier,
    MolDatBase,
    MolDatRDKit,
    OpDatBase,
)


class ReactionFilter(ABC):
    @abstractmethod
    def __call__(
        self,
        operator: OpDatBase,
        reactants: Sequence[MolDatBase],
        products: Sequence[MolDatBase],
    ) -> bool:
        pass


class AlwaysTrueFilter(ReactionFilter):
    def __call__(self, operator, reactants, products):
        return True


class ChainFilter(ReactionFilter):
    def __init__(self, filters: Iterable[ReactionFilter]):
        self._filters = filters

    def __call__(self, operator, reactants, products):
        return all(
            (filter(operator, reactants, products) for filter in self._filters)
        )


class LessThanNElementTypeFilter(ReactionFilter):
    def __init__(self, n: int, proton_number: int):
        self._n = n
        self._p = proton_number
        self._q = AtomNumEqualsQueryAtom(proton_number)

    def __call__(self, operator, reactants, products):
        for mol in products:
            if isinstance(mol, MolDatRDKit):
                if len(mol.rdkitmol.GetAtomsMatchingQuery(self._q)) >= self._n:
                    return False
        return True

    def __getstate__(self):
        return (self._n, self._p)

    def __setstate__(self, arg) -> None:
        self._n = arg[0]
        self._p = arg[1]
        self._q = AtomNumEqualsQueryAtom(self._p)


class UIDPreFilter(ABC):
    @abstractmethod
    def __call__(
        self,
        operator: Identifier,
        reactants: Sequence[Identifier],
    ) -> bool:
        pass


class AlwaysTrueUIDPreFilter(UIDPreFilter):
    def __call__(
        self, operator: Identifier, reactants: Sequence[Identifier]
    ) -> bool:
        return True


class CoreactantUIDPreFilter(UIDPreFilter):
    def __init__(self, coreactants: Collection[Identifier]):
        self._coreactants = frozenset(coreactants)

    def __call__(
        self, operator: Identifier, reactants: Sequence[Identifier]
    ) -> bool:
        for uid in reactants:
            if uid not in self._coreactants:
                return True
        return False

class TanimotoSimilarityFilter(ReactionFilter):
    def __init__(self, n: float, smi: str):
        self._n = n
        self._s = smi
        self._tmol = Chem.MolFromSmiles(self._s)
        self._tfp = Chem.RDKFingerprint(self._tmol)

    def __call__(self, operator, reactants, products):
        for mol in products:
            if isinstance(mol, MolDatRDKit):
                mol_fp = Chem.RDKFingerprint(mol.rdkitmol)
                similarity = DataStructs.TanimotoSimilarity(mol_fp, self._tfp)
                
                if similarity > self._n:
                    return True
        return False
