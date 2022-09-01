from dataclasses import dataclass
from typing import Optional

from pickaxe_generic.datatypes import DataPacketE, MetaKeyPacket, MolDatBase
from pickaxe_generic.engine import create_engine
from pickaxe_generic.filters import CoreactantFilter
from pickaxe_generic.metadata import (
    LocalPropertyCalc,
    MetaDataResolverFunc,
    MetaUpdateResolver,
    MolPropertyCalc,
    MolPropertyFromRxnCalc,
    ReactionFilterBase,
)
from pickaxe_generic.network import ChemNetworkBasic, ReactionExplicit
from pickaxe_generic.strategies import PriorityQueueStrategyBasic

network = ChemNetworkBasic()
engine = create_engine()

main_reagent = "C1(C)OC(=O)C=C(O)C=1"  # TAL

mol_smiles = (
    "CO",
    "O",
    "[H][H]",
    "N",
    "OO",
    "C=C",
    "CC=O",
)

op_kekule_smarts = (
    # hydrogenation of alkene
    "[C+0:1]=[C+0:2].[H][H]>>[*:1][*:2]",
    # hydrogenation of ketone
    "[C+0:1]=[O+0:2].[H][H]>>[*:1][*:2]",
    # diels-alder, something isn't working here
    "[C+0:1]=[C+0:2][C+0:3]=[C+0:4].[C+0:5]=[C+0:6]>>[*:1]1[*:2]=[*:3][*:4][*:5][*:6]1",
    # hydrolysis of ether
    "[*+0:1][O+0:2]!@[*+0:3].[O+0H2:4]>>[*:1][*:2].[*:3][*:4]",
    "[*+0:1][O+0:2]@[*+0:3].[O+0H2:4]>>([*:1][*:2].[*:3][*:4])",
    # keto-enol tautomerization
    "[C+0H:1][C+0:2]=[O+0:3]>>[*:1]=[*:2][*:3]",
    "[C+0:1]=[C+0:2][O+0H:3]>>[*:1][*:2]=[*:3]",
    # epoxidation of alkene
    "[C+0:1]=[C+0:2].[O+0H:3][O+0H:4]>>[*:1]1[*:2][*:3]1.[*:4]",
)

network.add_mol(engine.Mol(main_reagent), meta={"gen": 0})
coreagents = []
for smiles in mol_smiles:
    coreagents.append(network.add_mol(engine.Mol(smiles), meta={"gen": 0}))

for smarts in op_kekule_smarts:
    network.add_op(engine.Op(smarts, kekulize=True))

strategy = PriorityQueueStrategyBasic(network)


@dataclass(frozen=True)
class GenerationFilter(ReactionFilterBase):
    __slots__ = ("gen_key", "final_gen")

    gen_key: str
    final_gen: int

    def __call__(self, recipe: ReactionExplicit) -> bool:
        if all(
            mol.meta is None
            or self.gen_key not in mol.meta
            or mol.meta[self.gen_key] + 1 < self.final_gen
            for mol in recipe.reactants
        ):
            return True
        return False

    @property
    def meta_required(self) -> MetaKeyPacket:
        return MetaKeyPacket(frozenset((self.gen_key,)))


@dataclass(frozen=True)
class GenerationCalculator(MolPropertyFromRxnCalc[int]):
    __slots__ = ("gen_key",)

    gen_key: str

    @property
    def key(self) -> str:
        return self.gen_key

    @property
    def meta_required(self) -> MetaKeyPacket:
        return MetaKeyPacket(molecule_keys=frozenset((self.gen_key,)))

    @property
    def resolver(self) -> MetaDataResolverFunc[int]:
        return lambda x, y: min(x, y)

    def __call__(
        self,
        data: DataPacketE[MolDatBase],
        rxn: ReactionExplicit,
        prev_value: Optional[int] = None,
    ) -> Optional[int]:
        if data in rxn.reactants:
            return None
        cur_gen = 1
        for reactant in rxn.reactants:
            if reactant.meta is not None and self.gen_key in reactant.meta:
                cur_gen = max(cur_gen, reactant.meta[self.gen_key] + 1)
        if prev_value is not None and prev_value < cur_gen:
            return None
        if (
            data.meta is not None
            and self.gen_key in data.meta
            and data.meta[self.gen_key] <= cur_gen
        ):
            return None
        return cur_gen


strategy.expand(
    max_recipes=20,
    heap_size=20,
    recipe_filter=CoreactantFilter(coreagents),
    mc_local=GenerationCalculator("gen") >> GenerationFilter("gen", 2),
    mc_update=MetaUpdateResolver({"gen": lambda x, y: min(x, y)}, {}, {}),
)

print(network)
