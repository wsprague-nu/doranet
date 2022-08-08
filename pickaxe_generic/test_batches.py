from pickaxe_generic.engine import create_engine
from pickaxe_generic.filters import CoreactantFilter
from pickaxe_generic.network import ChemNetworkBasic
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

network.add_mol(engine.Mol(main_reagent))
coreagents = []
for smiles in mol_smiles:
    coreagents.append(network.add_mol(engine.Mol(smiles)))

for smarts in op_kekule_smarts:
    network.add_op(engine.Op(smarts, kekulize=True))

strategy = PriorityQueueStrategyBasic(network)

strategy.expand(
    max_recipes=10, heap_size=10, recipe_filter=CoreactantFilter(coreagents)
)

print(network)
