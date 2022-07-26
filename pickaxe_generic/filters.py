from abc import ABC, abstractmethod
from collections.abc import Collection, Iterable, Mapping, Sequence
from dataclasses import dataclass
from typing import Optional, Protocol, final

from rdkit.Chem import MolFromSmiles, RDKFingerprint
from rdkit.Chem.rdqueries import AtomNumEqualsQueryAtom
from rdkit.DataStructs import TanimotoSimilarity

from pickaxe_generic.datatypes import (
    Identifier,
    MolDatBase,
    MolDatRDKit,
    OpDatBase,
)
from pickaxe_generic.network import (
    ChemNetwork,
    ReactionExplicit,
    RecipeExplicit,
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
        self._tmol = MolFromSmiles(self._s)
        self._tfp = RDKFingerprint(self._tmol)

    def __call__(self, operator, reactants, products):
        for mol in products:
            if isinstance(mol, MolDatRDKit):
                mol_fp = RDKFingerprint(mol.rdkitmol)
                similarity = TanimotoSimilarity(mol_fp, self._tfp)

                if similarity > self._n:
                    return True
        return False


@dataclass(frozen=True)
class MetaKeyPacket:
    operator_keys: frozenset = frozenset()
    reactant_keys: frozenset = frozenset()

    def __add__(self, other: 'MetaKeyPacket') -> 'MetaKeyPacket':
        return MetaKeyPacket(self.operator_keys.union(other.operator_keys),self.reactant_keys.union(other.reactant_keys))


class RecipeFilter(ABC):
    __slots__ = ()

    @abstractmethod
    def __call__(self, recipe: RecipeExplicit) -> bool:
        ...

    @property
    def meta_required(self) -> MetaKeyPacket:
        return MetaKeyPacket()

    @final
    def __and__(self, other: "RecipeFilter") -> "RecipeFilter":
        return RecipeFilterAnd(self, other)

    @final
    def __inv__(self, other: "RecipeFilter") -> "RecipeFilter":
        return RecipeFilterInv(self, other)

    @final
    def __or__(self, other: "RecipeFilter") -> "RecipeFilter":
        return RecipeFilterOr(self, other)

    @final
    def __xor__(self, other: "RecipeFilter") -> "RecipeFilter":
        return RecipeFilterXor(self, other)


@dataclass(frozen=True)  # type: ignore
def RecipeFilterAnd(RecipeFilter):
    __slots__ = ("_filter1", "_filter2")

    _filter1: RecipeFilter
    _filter2: RecipeFilter

    def __call__(self, recipe: RecipeExplicit) -> bool:
        return self._filter1(recipe) and self._filter2(recipe)

    @property
    def meta_required(self) -> MetaKeyPacket:
        return self._filter1.meta_required + self._filter2.meta_required


@dataclass(frozen=True)  # type: ignore
def RecipeFilterInv(RecipeFilter):
    __slots__ = ("_filter",)
    _filter: RecipeFilter

    def __call__(self, recipe: RecipeExplicit) -> bool:
        return not self._filter(recipe)

    @property
    def meta_required(self) -> MetaKeyPacket:
        return self._filter.meta_required


@dataclass(frozen=True)  # type: ignore
def RecipeFilterOr(RecipeFilter):
    __slots__ = ("_filter1", "_filter2")
    _filter1: RecipeFilter
    _filter2: RecipeFilter

    def __call__(self, recipe: RecipeExplicit) -> bool:
        return self._filter1(recipe) or self._filter2(recipe)

    @property
    def meta_required(self) -> MetaKeyPacket:
        return self._filter1.meta_required + self._filter2.meta_required


@dataclass(frozen=True)  # type: ignore
def RecipeFilterXor(RecipeFilter):
    __slots__ = ("_filter1", "_filter2")
    _filter1: RecipeFilter
    _filter2: RecipeFilter

    def __call__(self, recipe: RecipeExplicit) -> bool:
        return self._filter1(recipe) != self._filter2(recipe)

    @property
    def meta_required(self) -> MetaKeyPacket:
        return self._filter1.meta_required + self._filter2.meta_required


class RankValue(Protocol):
    @abstractmethod
    def __lt__(self, other: "RankValue") -> bool:
        ...


class RecipeRanker(Protocol):
    __slots__ = ()

    @abstractmethod
    def __call__(self, recipe: RecipeExplicit) -> Optional[RankValue]:
        ...


class ReactionFilterBase(ABC):
    __slots__ = ()

    @abstractmethod
    def __call__(self, recipe: ReactionExplicit) -> bool:
        ...

    @property
    def meta_required(self) -> MetaKeyPacket:
        return MetaKeyPacket()

    @final
    def __and__(self, other: "ReactionFilterBase") -> "ReactionFilterBase":
        return ReactionFilterAnd(self, other)

    @final
    def __inv__(self, other: "ReactionFilterBase") -> "ReactionFilterBase":
        return ReactionFilterInv(self, other)

    @final
    def __or__(self, other: "ReactionFilterBase") -> "ReactionFilterBase":
        return ReactionFilterOr(self, other)

    @final
    def __xor__(self, other: "ReactionFilterBase") -> "ReactionFilterBase":
        return ReactionFilterXor(self, other)


@dataclass(frozen=True)  # type: ignore
def ReactionFilterAnd(ReactionFilterBase):
    __slots__ = ("_filter1", "_filter2")

    _filter1: ReactionFilterBase
    _filter2: ReactionFilterBase

    def __call__(self, recipe: ReactionExplicit) -> bool:
        return self._filter1(recipe) and self._filter2(recipe)

    @property
    def meta_required(self) -> MetaKeyPacket:
        return self._filter1.meta_required + self._filter2.meta_required


@dataclass(frozen=True)  # type: ignore
def ReactionFilterInv(ReactionFilterBase):
    __slots__ = ("_filter",)
    _filter: ReactionFilterBase

    def __call__(self, recipe: ReactionExplicit) -> bool:
        return not self._filter(recipe)

    @property
    def meta_required(self) -> MetaKeyPacket:
        return self._filter.meta_required


@dataclass(frozen=True)  # type: ignore
def ReactionFilterOr(ReactionFilterBase):
    __slots__ = ("_filter1", "_filter2")
    _filter1: ReactionFilterBase
    _filter2: ReactionFilterBase

    def __call__(self, recipe: ReactionExplicit) -> bool:
        return self._filter1(recipe) or self._filter2(recipe)

    @property
    def meta_required(self) -> MetaKeyPacket:
        return self._filter1.meta_required + self._filter2.meta_required


@dataclass(frozen=True)  # type: ignore
def ReactionFilterXor(ReactionFilterBase):
    __slots__ = ("_filter1", "_filter2")
    _filter1: ReactionFilterBase
    _filter2: ReactionFilterBase

    def __call__(self, recipe: ReactionExplicit) -> bool:
        return self._filter1(recipe) != self._filter2(recipe)

    @property
    def meta_required(self) -> MetaKeyPacket:
        return self._filter1.meta_required + self._filter2.meta_required


def LocalMetaDataCalculator(ABC):
    __slots__ = ()

    # @abstractmethod


# def __call__(self, )
