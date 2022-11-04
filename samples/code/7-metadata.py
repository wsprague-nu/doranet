import pickaxe_generic as pg

engine = pg.create_engine()

network = engine.new_network()

reagents = [
    "[H][H]",  # hydrogen
    "O",  # water
    "CO",  # methanol
    "CCO",  # ethanol
    "CC(O)=O",  # acetic acid
]

operator_smarts = {
    "ester_hydrolysis_nonring": "[O&+0:1]=[C&+0:2]-&!@[O&+0&H0:3].[O&+0&H2:4]>>[*:1]=[*:2]-[*:4].[*:3]",
    "ester_hydrolysis_ring": "[O&+0:1]=[C&+0:2]-&@[O&+0&H0:3].[O&+0&H2:4]>>([*:1]=[*:2]-[*:4].[*:3])",
    "esterification": "[O&+0:1]=[C&+0:2]-[O&+0&H1:3].[O&+0&H1:4]>>[*:1]=[*:2]-[*:4].[*:3]",
    "esterification_intra": "([O&+0:1]=[C&+0:2]-[O&+0&H1:3].[O&+0&H1:4])>>[*:1]=[*:2]-[*:4].[*:3]",
    "hydrogenation of carbonyl": "[C+0:1]=[O+0:2].[H][H]>>[*:1][*:2]",
}

for smiles in reagents:
    network.add_mol(engine.mol.rdkit(smiles), meta={"gen": 0, "waste": 0})

for name, smarts in operator_smarts.items():
    network.add_op(engine.op.rdkit(smarts), meta={"name": name})

network.save_to_file("7-metadata")

from pprint import pprint

network = engine.network_from_file("7-metadata")

strat = engine.strat.cartesian(network)

gen_calc = engine.meta.generation("gen")
mw_calc = engine.meta.mw("mw")
mass_efficiency_calc = engine.meta.masswaste("waste", "mw")
gen_filter = engine.filter.reaction.generation(3, "gen")
reaction_plan = (gen_calc & mw_calc) >> gen_filter >> mass_efficiency_calc

strat.expand(num_iter=5, reaction_plan=reaction_plan)

pprint(
    [
        (i, v[0].smiles, v[1])
        for i, v in enumerate(
            zip(network.mols, network.mols.meta(keys=["gen", "waste"]))
        )
    ]
)

pprint(list(enumerate(network.rxns)))
