from abc import ABC, abstractmethod
from collections.abc import Collection, Iterable, Mapping, Sequence
from dataclasses import dataclass
from typing import Any, Generator, Hashable, Optional, Protocol, Union, final

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
    _MolIndex,
    _OpIndex,
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
    molecule_keys: frozenset = frozenset()

    def __add__(self, other: "MetaKeyPacket") -> "MetaKeyPacket":
        return MetaKeyPacket(
            self.operator_keys.union(other.operator_keys),
            self.molecule_keys.union(other.molecule_keys),
        )


class MolFilter(ABC):
    __slots__ = ()

    @abstractmethod
    def __call__(
        self, mol: MolDatBase, meta: Optional[Mapping[Hashable, Any]] = None
    ) -> bool:
        ...

    @property
    def meta_required(self) -> MetaKeyPacket:
        return MetaKeyPacket()

    @final
    def __and__(self, other: "MolFilter") -> "MolFilter":
        return MolFilterAnd(self, other)

    @final
    def __inv__(self) -> "MolFilter":
        return MolFilterInv(self)

    @final
    def __or__(self, other: "MolFilter") -> "MolFilter":
        return MolFilterOr(self, other)

    @final
    def __xor__(self, other: "MolFilter") -> "MolFilter":
        return MolFilterXor(self, other)


@dataclass(frozen=True)
class MolFilterAnd(MolFilter):
    __slots__ = ("_filter1", "_filter2")

    _filter1: MolFilter
    _filter2: MolFilter

    def __call__(
        self, mol: MolDatBase, meta: Optional[Mapping[Hashable, Any]] = None
    ) -> bool:
        return self._filter1(mol, meta) and self._filter2(mol, meta)

    @property
    def meta_required(self) -> MetaKeyPacket:
        return self._filter1.meta_required + self._filter2.meta_required


@dataclass(frozen=True)
class MolFilterInv(MolFilter):
    __slots__ = ("_filter",)
    _filter: MolFilter

    def __call__(
        self, mol: MolDatBase, meta: Optional[Mapping[Hashable, Any]] = None
    ) -> bool:
        return not self._filter(mol, meta)

    @property
    def meta_required(self) -> MetaKeyPacket:
        return self._filter.meta_required


@dataclass(frozen=True)
class MolFilterOr(MolFilter):
    __slots__ = ("_filter1", "_filter2")
    _filter1: MolFilter
    _filter2: MolFilter

    def __call__(
        self, mol: MolDatBase, meta: Optional[Mapping[Hashable, Any]] = None
    ) -> bool:
        return self._filter1(mol, meta) or self._filter2(mol, meta)

    @property
    def meta_required(self) -> MetaKeyPacket:
        return self._filter1.meta_required + self._filter2.meta_required


@dataclass(frozen=True)
class MolFilterXor(MolFilter):
    __slots__ = ("_filter1", "_filter2")
    _filter1: MolFilter
    _filter2: MolFilter

    def __call__(
        self, mol: MolDatBase, meta: Optional[Mapping[Hashable, Any]] = None
    ) -> bool:
        return self._filter1(mol, meta) != self._filter2(mol, meta)

    @property
    def meta_required(self) -> MetaKeyPacket:
        return self._filter1.meta_required + self._filter2.meta_required


@dataclass(frozen=True)
class MolFilterMetaVal(MolFilter):
    __slots__ = ("_key", "_val")
    _key: Hashable
    _val: Any

    def __call__(
        self, mol: MolDatBase, meta: Optional[Mapping[Hashable, Any]] = None
    ) -> bool:
        if meta is None:
            return False
        if self._key not in meta:
            return False
        return meta[self._key] == self._val

    @property
    def meta_required(self) -> MetaKeyPacket:
        return MetaKeyPacket(frozenset(), frozenset((self._key,)))


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
    def __inv__(self) -> "RecipeFilter":
        return RecipeFilterInv(self)

    @final
    def __or__(self, other: "RecipeFilter") -> "RecipeFilter":
        return RecipeFilterOr(self, other)

    @final
    def __xor__(self, other: "RecipeFilter") -> "RecipeFilter":
        return RecipeFilterXor(self, other)


@dataclass(frozen=True)
class RecipeFilterAnd(RecipeFilter):
    __slots__ = ("_filter1", "_filter2")

    _filter1: RecipeFilter
    _filter2: RecipeFilter

    def __call__(self, recipe: RecipeExplicit) -> bool:
        return self._filter1(recipe) and self._filter2(recipe)

    @property
    def meta_required(self) -> MetaKeyPacket:
        return self._filter1.meta_required + self._filter2.meta_required


@dataclass(frozen=True)
class RecipeFilterInv(RecipeFilter):
    __slots__ = ("_filter",)
    _filter: RecipeFilter

    def __call__(self, recipe: RecipeExplicit) -> bool:
        return not self._filter(recipe)

    @property
    def meta_required(self) -> MetaKeyPacket:
        return self._filter.meta_required


@dataclass(frozen=True)
class RecipeFilterOr(RecipeFilter):
    __slots__ = ("_filter1", "_filter2")
    _filter1: RecipeFilter
    _filter2: RecipeFilter

    def __call__(self, recipe: RecipeExplicit) -> bool:
        return self._filter1(recipe) or self._filter2(recipe)

    @property
    def meta_required(self) -> MetaKeyPacket:
        return self._filter1.meta_required + self._filter2.meta_required


@dataclass(frozen=True)
class RecipeFilterXor(RecipeFilter):
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
    @abstractmethod
    def __call__(self, recipe: RecipeExplicit) -> Optional[RankValue]:
        ...

    @property
    def meta_required(self) -> MetaKeyPacket:
        return MetaKeyPacket()


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
    def __inv__(self) -> "ReactionFilterBase":
        return ReactionFilterInv(self)

    @final
    def __or__(self, other: "ReactionFilterBase") -> "ReactionFilterBase":
        return ReactionFilterOr(self, other)

    @final
    def __xor__(self, other: "ReactionFilterBase") -> "ReactionFilterBase":
        return ReactionFilterXor(self, other)


@dataclass(frozen=True)
class ReactionFilterAnd(ReactionFilterBase):
    __slots__ = ("_filter1", "_filter2")

    _filter1: ReactionFilterBase
    _filter2: ReactionFilterBase

    def __call__(self, recipe: ReactionExplicit) -> bool:
        return self._filter1(recipe) and self._filter2(recipe)

    @property
    def meta_required(self) -> MetaKeyPacket:
        return self._filter1.meta_required + self._filter2.meta_required


@dataclass(frozen=True)
class ReactionFilterInv(ReactionFilterBase):
    __slots__ = ("_filter",)
    _filter: ReactionFilterBase

    def __call__(self, recipe: ReactionExplicit) -> bool:
        return not self._filter(recipe)

    @property
    def meta_required(self) -> MetaKeyPacket:
        return self._filter.meta_required


@dataclass(frozen=True)
class ReactionFilterOr(ReactionFilterBase):
    __slots__ = ("_filter1", "_filter2")
    _filter1: ReactionFilterBase
    _filter2: ReactionFilterBase

    def __call__(self, recipe: ReactionExplicit) -> bool:
        return self._filter1(recipe) or self._filter2(recipe)

    @property
    def meta_required(self) -> MetaKeyPacket:
        return self._filter1.meta_required + self._filter2.meta_required


@dataclass(frozen=True)
class ReactionFilterXor(ReactionFilterBase):
    __slots__ = ("_filter1", "_filter2")
    _filter1: ReactionFilterBase
    _filter2: ReactionFilterBase

    def __call__(self, recipe: ReactionExplicit) -> bool:
        return self._filter1(recipe) != self._filter2(recipe)

    @property
    def meta_required(self) -> MetaKeyPacket:
        return self._filter1.meta_required + self._filter2.meta_required


class MetaDataCalculatorLocal(Protocol):
    @abstractmethod
    def __call__(self, unit: Union[ReactionExplicit, RecipeExplicit]) -> None:
        ...

    @property
    @abstractmethod
    def meta_required(self) -> MetaKeyPacket:
        ...


class MetaDataUpdate(Protocol):
    @abstractmethod
    def __call__(
        self, unit: ReactionExplicit, network: ChemNetwork
    ) -> Generator[tuple[Optional[_MolIndex], Optional[_OpIndex]], None, None]:
        ...


class DefaultMetaDataUpdate:
    def __call__(
        self, unit: ReactionExplicit, network: ChemNetwork
    ) -> Generator[tuple[Optional[_MolIndex], Optional[_OpIndex]], None, None]:
        if unit.operator_meta is not None:
            opIndex = network.ops.i(unit.operator.uid)
            for key, value in unit.operator_meta.items():
                network.op_meta(opIndex, key, value)
            yield (None, opIndex)
        if unit.reactants_meta is not None:
            for index, reactant in enumerate(unit.reactants):
                molIndex = network.mols.i(reactant.uid)
                for key, value in unit.reactants_meta[index]:
                    network.mol_meta(molIndex, key, value)
                yield (molIndex, None)
        if unit.products_meta is not None:
            for index, product in enumerate(unit.products):
                molIndex = network.mols.i(product.uid)
                for key, value in unit.products_meta[index]:
                    network.mol_meta(molIndex, key, value)
                yield (molIndex, None)


# def __call__(self, )
