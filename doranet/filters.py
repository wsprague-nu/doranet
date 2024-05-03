"""Contains classes which implement various filter components."""

import collections.abc
import dataclasses
import typing

import rdkit
import rdkit.Chem
import rdkit.Chem.rdqueries

from doranet import interfaces, metadata

# class AlwaysTrueFilter(interfaces.ReactionFilter):
#     def __call__(self, operator, reactants, products):
#         return True


# class ChainFilter(interfaces.ReactionFilter):
#     def __init__(
#         self, filters: collections.abc.Iterable[interfaces.ReactionFilter]
#     ):
#         self._filters = filters

#     def __call__(self, operator, reactants, products):
#         return all(
#             (
#                 filter(operator, reactants, products)
#                 for filter in self._filters
#             )  # yo
#         )


# class LessThanNElementTypeFilter(interfaces.ReactionFilter):
#     def __init__(self, n: int, proton_number: int):
#         self._n = n
#         self._p = proton_number
#         self._q = rdkit.Chem.rdqueries.AtomNumEqualsQueryAtom(proton_number)

#     def __call__(self, operator, reactants, products):
#         for mol in products:
#             if (
#                 isinstance(mol, interfaces.MolDatRDKit)
#                 and len(mol.rdkitmol.GetAtomsMatchingQuery(self._q)) >=
#                     self._n
#             ):
#                 return False
#         return True

#     def __getstate__(self):
#         return (self._n, self._p)

#     def __setstate__(self, arg) -> None:
#         self._n = arg[0]
#         self._p = arg[1]
#         self._q = rdkit.Chem.rdqueries.AtomNumEqualsQueryAtom(self._p)


# class TanimotoSimilarityFilter(interfaces.ReactionFilter):
#     def __init__(self, n: float, smi: str):
#         self._n = n
#         self._s = smi
#         self._tmol = rdkit.Chem.MolFromSmiles(self._s)
#         self._tfp = rdkit.Chem.RDKFingerprint(self._tmol)

#     def __call__(self, operator, reactants, products):
#         for mol in products:
#             if isinstance(mol, interfaces.MolDatRDKit):
#                 mol_fp = rdkit.Chem.RDKFingerprint(mol.rdkitmol)
#                 similarity = rdkit.DataStructs.TanimotoSimilarity(
#                     mol_fp, self._tfp
#                 )

#                 if similarity > self._n:
#                     return True
#         return False


# class AlwaysTrueUIDPreFilter(interfaces.UIDPreFilter):
#     def __call__(
#         self,
#         operator: interfaces.Identifier,
#         reactants: collections.abc.Sequence[interfaces.Identifier],
#     ) -> bool:
#         return True


# @dataclasses.dataclass(frozen=True)
# class CoreactantUIDPreFilter(interfaces.UIDPreFilter):
#     coreactants: collections.abc.Container[interfaces.Identifier]

#     def __call__(
#         self,
#         operator: interfaces.Identifier,
#         reactants: collections.abc.Sequence[interfaces.Identifier],
#     ) -> bool:
#         return any(uid not in self.coreactants for uid in reactants)


@typing.final
@dataclasses.dataclass(frozen=True)
class MolFilterMetaVal(interfaces.MolFilter):
    __slots__ = ("key", "val")
    key: collections.abc.Hashable
    val: typing.Any

    def __call__(
        self,
        mol: interfaces.DataPacket[interfaces.MolDatBase],
        op: typing.Optional[interfaces.DataPacket[interfaces.OpDatBase]],
        arg_num: typing.Optional[int],
    ) -> bool:
        if mol.meta is None:
            return False
        if self.key not in mol.meta:
            return False
        return mol.meta[self.key] == self.val

    @property
    def meta_required(self) -> interfaces.MetaKeyPacket:
        return interfaces.MetaKeyPacket(frozenset(), frozenset((self.key,)))


@typing.final
@dataclasses.dataclass(frozen=True)
class MolFilterMetaExist(interfaces.MolFilter):
    __slots__ = ("key",)
    key: collections.abc.Hashable

    def __call__(
        self,
        mol: interfaces.DataPacket[interfaces.MolDatBase],
        op: typing.Optional[interfaces.DataPacket[interfaces.OpDatBase]],
        arg_num: typing.Optional[int],
    ) -> bool:
        if mol.meta is None or self.key not in mol.meta:
            return False
        return True

    @property
    def meta_required(self) -> interfaces.MetaKeyPacket:
        return interfaces.MetaKeyPacket(frozenset(), frozenset((self.key,)))


@typing.final
@dataclasses.dataclass(frozen=True, slots=True)
class MolFilterMetaFunc(interfaces.MolFilter):
    key: collections.abc.Hashable
    pred: collections.abc.Callable[[typing.Any], bool]
    unknown_pass: bool = False

    def __call__(
        self,
        mol: interfaces.DataPacket[interfaces.MolDatBase],
        op: typing.Optional[interfaces.DataPacket[interfaces.OpDatBase]],
        arg_num: typing.Optional[int],
    ) -> bool:
        meta = mol.meta
        key = self.key
        if meta is None or key not in meta:
            return self.unknown_pass
        return self.pred(meta[key])

    @property
    def meta_required(self) -> interfaces.MetaKeyPacket:
        return interfaces.MetaKeyPacket(frozenset(), frozenset((self.key,)))


@typing.final
@dataclasses.dataclass(frozen=True, slots=True)
class MolFilterIndex(interfaces.MolFilter):
    indices: collections.abc.Container[interfaces.MolIndex]
    whitelist: bool = False

    def __call__(
        self,
        mol: interfaces.DataPacket[interfaces.MolDatBase],
        op: typing.Optional[interfaces.DataPacket[interfaces.OpDatBase]],
        arg_num: typing.Optional[int],
    ) -> bool:
        if mol.i in self.indices:
            return self.whitelist
        return not self.whitelist

    @property
    def meta_required(self) -> interfaces.MetaKeyPacket:
        return interfaces.MetaKeyPacket()


@typing.final
@dataclasses.dataclass(frozen=True, slots=True)
class BundleFilterCoreactants(interfaces.BundleFilter):
    coreagents: collections.abc.Container[interfaces.MolIndex]

    def __call__(
        self, bundle: interfaces.RecipeBundle
    ) -> collections.abc.Iterable[interfaces.RecipeBundle]:
        mol_args = bundle.args
        n_args = len(mol_args)
        coreagents = self.coreagents
        coreagents_args = tuple(
            tuple(mol for mol in mols if mol.i in coreagents)
            for mols in mol_args
        )
        for i in range(n_args):
            yield interfaces.RecipeBundle(
                bundle.operator,
                coreagents_args[:i]
                + (
                    tuple(
                        mol for mol in mol_args[i] if mol.i not in coreagents
                    ),
                )
                + mol_args[i + 1 :],
            )

    @property
    def meta_required(self) -> interfaces.MetaKeyPacket:
        return interfaces.MetaKeyPacket()


@typing.final
@dataclasses.dataclass(frozen=True)
class CoreactantFilter(interfaces.RecipeFilter):
    __slots__ = ("coreactants",)
    coreactants: collections.abc.Container[interfaces.MolIndex]

    def __call__(self, recipe: interfaces.RecipeExplicit) -> bool:
        if all(mol.i in self.coreactants for mol in recipe.reactants):
            return False
        return True

    @property
    def meta_required(self) -> interfaces.MetaKeyPacket:
        return interfaces.MetaKeyPacket()


@typing.final
@dataclasses.dataclass(frozen=True)
class GenerationFilter(metadata.ReactionFilterBase):
    __slots__ = ("max_gens", "gen_key")

    max_gens: int
    gen_key: collections.abc.Hashable

    def __call__(self, recipe: interfaces.ReactionExplicit) -> bool:
        if all(
            mol.meta is not None
            and self.gen_key in mol.meta
            and mol.meta[self.gen_key] + 1 < self.max_gens
            for mol in recipe.reactants
        ):
            return True
        return False

    @property
    def meta_required(self) -> interfaces.MetaKeyPacket:
        return interfaces.MetaKeyPacket(frozenset((self.gen_key,)))


@dataclasses.dataclass(frozen=True, slots=True)
class ReactionFilterMaxAtoms(metadata.ReactionFilterBase):
    max_atoms: int
    query: typing.Any

    @classmethod
    def from_num(
        cls, max_atoms: int, proton_number: typing.Optional[int] = None
    ) -> "ReactionFilterMaxAtoms":
        if proton_number is None:
            return ReactionFilterMaxAtoms(max_atoms, None)
        return ReactionFilterMaxAtoms(
            max_atoms,
            rdkit.Chem.rdqueries.AtomNumEqualsQueryAtom(proton_number),
        )

    def __call__(self, recipe: interfaces.ReactionExplicit):
        for mol in recipe.products:
            if not isinstance(mol.item, interfaces.MolDatRDKit):
                raise NotImplementedError(
                    f"""Counting # of atoms in non-RDKit molecules is not yet
                        supported (found {repr(mol)})"""
                )
            if (
                len(mol.item.rdkitmol.GetAtomsMatchingQuery(self.query))
                > self.max_atoms
            ):
                return False
        return True


def ReplaceNewValue(
    key: collections.abc.Hashable, old_value: typing.Any, new_value: typing.Any
) -> bool:
    if old_value != new_value:
        return True
    return False


class ReplaceBlacklist:
    def __init__(
        self,
        blacklist_keys: collections.abc.Collection[collections.abc.Hashable],
    ) -> None:
        self._blacklist_keys = frozenset(blacklist_keys)

    def __call__(
        self,
        key: collections.abc.Hashable,
        old_value: typing.Any,
        new_value: typing.Any,
    ) -> bool:
        if key in self._blacklist_keys:
            return (old_value is None) or (
                old_value is False and new_value is True
            )
        return ReplaceNewValue(key, old_value, new_value)


"""
class DefaultMetaDataUpdate:
    def __call__(
        self,
        unit: interfaces.ReactionExplicit,
        network: interfaces.ChemNetwork,
    ) -> Generator[
        tuple[
            Optional[tuple[interfaces.MolIndex, Hashable]],
            Optional[tuple[interfaces.OpIndex, Hashable]],
        ],
        None,
        None,
    ]:
        if unit.operator_meta is not None:
            opIndex = network.ops.i(unit.operator.uid)
            for key, value in unit.operator_meta.items():
                if network.op_meta(opIndex, key) != value:
                    network.op_meta(opIndex, key, value)
                    yield (None, (opIndex, key))
        if unit.reactants_meta is not None:
            for index, reactant in enumerate(unit.reactants):
                molIndex = network.mols.i(reactant.uid)
                for key, value in unit.reactants_meta[index]:
                    if network.mol_meta(molIndex, key) != value:
                        network.mol_meta(molIndex, key, value)
                        yield ((molIndex, key), None)
        if unit.products_meta is not None:
            for index, product in enumerate(unit.products):
                molIndex = network.mols.i(product.uid)
                for key, value in unit.products_meta[index]:
                    if network.mol_meta(molIndex, key) != value:
                        network.mol_meta(molIndex, key, value)
                        yield ((molIndex, key), None)


# def __call__(self, )
"""
