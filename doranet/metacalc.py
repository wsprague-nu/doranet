"""Contains metadata calculator classes."""

import collections.abc
import dataclasses
import typing

import rdkit.Chem.rdMolDescriptors

from doranet import interfaces, metadata


@typing.final
@dataclasses.dataclass(frozen=True)
class GenerationCalculator(metadata.MolPropertyFromRxnCalc[int]):
    __slots__ = ("gen_key",)

    gen_key: collections.abc.Hashable

    @property
    def key(self) -> collections.abc.Hashable:
        return self.gen_key

    @property
    def meta_required(self) -> interfaces.MetaKeyPacket:
        return interfaces.MetaKeyPacket(
            molecule_keys=frozenset((self.gen_key,))
        )

    @property
    def resolver(self) -> metadata.MetaDataResolverFunc[int]:
        return min

    def __call__(
        self,
        data: interfaces.DataPacketE[interfaces.MolDatBase],
        rxn: interfaces.ReactionExplicit,
        prev_value: typing.Optional[int] = None,
    ) -> typing.Optional[int]:
        if data in rxn.reactants:
            return None
        cur_gen = None
        for reactant in rxn.reactants:
            if reactant.meta is None or self.gen_key not in reactant.meta:
                return None
            if cur_gen is None:
                cur_gen = reactant.meta[self.gen_key] + 1
            cur_gen = max(cur_gen, reactant.meta[self.gen_key] + 1)
        if cur_gen is None:
            return None
        if prev_value is not None and prev_value < cur_gen:
            return None
        if (
            data.meta is not None
            and self.gen_key in data.meta
            and data.meta[self.gen_key] <= cur_gen
        ):
            return None
        return cur_gen


@typing.final
@dataclasses.dataclass(frozen=True, slots=True)
class MassWasteCalculator(metadata.MolPropertyFromRxnCalc[float]):
    masswaste_key: collections.abc.Hashable
    mw_key: collections.abc.Hashable

    @property
    def key(self) -> collections.abc.Hashable:
        return self.masswaste_key

    @property
    def meta_required(self) -> interfaces.MetaKeyPacket:
        return interfaces.MetaKeyPacket(
            molecule_keys=frozenset((self.masswaste_key, self.mw_key))
        )

    @property
    def resolver(self) -> metadata.MetaDataResolverFunc[float]:
        return min

    def __call__(
        self,
        data: interfaces.DataPacketE[interfaces.MolDatBase],
        rxn: interfaces.ReactionExplicit,
        prev_value: typing.Optional[float] = None,
    ) -> typing.Optional[float]:
        if (
            sum(1 for mol in rxn.products if mol.item == data.item)
            - sum(1 for mol in rxn.reactants if mol.item == data.item)
            <= 0
        ):
            return None
        if data.meta is None:
            return None
        masseff_key = self.masswaste_key
        mw_key = self.mw_key
        if any(
            mol.meta is None or masseff_key not in mol.meta
            for mol in rxn.reactants
        ):
            return None
        if any(
            mol.meta is None or mw_key not in mol.meta for mol in rxn.products
        ):
            return None
        reactants_sum: float = sum(
            mol.meta[masseff_key]  # type: ignore
            for mol in rxn.reactants
            if mol.item != data.item
        )
        products_sum: float = sum(
            mol.meta[mw_key]  # type: ignore
            for mol in rxn.products
            if mol.item != data.item  # type: ignore
        )
        final_result = reactants_sum + products_sum
        return final_result


@typing.final
@dataclasses.dataclass(frozen=True, slots=True)
class MolWeightCalculator(metadata.MolPropertyCalc[int]):
    weight_key: collections.abc.Hashable

    @property
    def key(self) -> collections.abc.Hashable:
        return self.weight_key

    @property
    def meta_required(self) -> interfaces.MetaKeyPacket:
        return interfaces.MetaKeyPacket()

    @property
    def resolver(self) -> metadata.MetaDataResolverFunc[int]:
        return metadata.TrivialMetaDataResolverFunc

    def __call__(
        self,
        data: interfaces.DataPacketE[interfaces.MolDatBase],
        prev_value: typing.Optional[int] = None,
    ) -> typing.Optional[int]:
        if prev_value is not None:
            return prev_value
        item = data.item
        if not isinstance(item, interfaces.MolDatRDKit):
            raise NotImplementedError(
                f"""MolWeightCalculator has not been implemented for molecule
                    type {type(item)}"""
            )
        return rdkit.Chem.rdMolDescriptors.CalcExactMolWt(item.rdkitmol)
