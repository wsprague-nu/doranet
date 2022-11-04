import pprint

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

network.save_to_file("8-global-hooks")

network_no_hook = engine.network_from_file("8-global-hooks")
network_with_hook = engine.network_from_file("8-global-hooks")

strat_no_hook = engine.strat.cartesian(network_no_hook)
strat_with_hook = engine.strat.cartesian(network_with_hook)

mol_limit_hook = engine.hook.max_mols(10)
gen_calc = engine.meta.generation("gen")

strat_no_hook.expand(num_iter=3, reaction_plan=gen_calc)
strat_with_hook.expand(
    num_iter=3, global_hooks=[mol_limit_hook], reaction_plan=gen_calc
)

pprint.pprint(
    [
        (i, v[0].smiles, v[1])
        for i, v in enumerate(
            zip(network_no_hook.mols, network_no_hook.mols.meta(keys=["gen"]))
        )
    ]
)

print()

pprint.pprint(
    [
        (i, v[0].smiles, v[1])
        for i, v in enumerate(
            zip(
                network_with_hook.mols,
                network_with_hook.mols.meta(keys=["gen"]),
            )
        )
    ]
)

network = engine.network_from_file("8-global-hooks")

strat = engine.strat.cartesian(network)

target_hook = engine.hook.target(engine.mol.rdkit("CC(O)O"))
gen_calc = engine.meta.generation("gen")

strat.expand(num_iter=3, global_hooks=[target_hook], reaction_plan=gen_calc)

pprint.pprint(
    [
        (i, v[0].smiles, v[1])
        for i, v in enumerate(
            zip(network.mols, network.mols.meta(keys=["gen"]))
        )
    ]
)
