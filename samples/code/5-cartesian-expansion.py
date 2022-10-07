import pprint

import pickaxe_generic as pg

engine = pg.create_engine()

initial_reactant_smiles = [
    "[H][H]",  # hydrogen
    "CC=O",  # acetaldehyde
    "CC(C)=O",  # acetone
    "CCCO",  # propanol
    "C=CC=C",  # butadiene
]

operator_smarts = {
    "hydrogenation of alkene/carbonyl": "[C,O;+0:1]=[C&+0:2].[#1][#1]>>[*:1]-[*:2]"
}

network = engine.new_network()
for smiles in initial_reactant_smiles:
    network.add_mol(engine.mol.rdkit(smiles))
for name, smarts in operator_smarts.items():
    network.add_op(engine.op.rdkit(smarts), meta={"name": name})

for i in range(len(network.mols)):
    network.mols.set_meta(i, {"generation": 0})
network.mols.set_meta(0, {"generation": -1})
network.mols.set_meta(4, {"generation": -1})
strat = engine.strat.cartesian(network, gen_key="generation")
strat.expand(num_iter=2)
pprint.pprint(list(enumerate(network.mols)))
pprint.pprint(list(network.rxns))
pprint.pprint(list(network.mols.meta()))

operator_smarts_big = {
    "hydrogenation of alkene/carbonyl": "[C,O;+0:1]=[C&+0:2].[#1][#1]>>[*:1]-[*:2]",
    "hydrogenolysis/hydrodeoxygenation of C-O bond (non-ring)": "[C&+0:1]-&!@[O&+0:2].[#1][#1]>>[*:1].[*:2]",
    "hydrogenolysis of C-O bond (ring)": "[C&+0:1]-&@[O&+0:2].[#1][#1]>>([*:1].[*:2])",
    "hydrolysis of ester (non-ring)": "[O&+0:1]=[C&+0:2]-&!@[O&+0&H0:3].[O&+0&H2:4]>>[*:1]=[*:2]-[*:4].[*:3]",
    "hydrolysis of ester (ring)": "[O&+0:1]=[C&+0:2]-&@[O&+0&H0:3].[O&+0&H2:4]>>([*:1]=[*:2]-[*:4].[*:3])",
    "hydrolysis of ether (non-ring)": "[C&+0:1]-[O&+0&H0:2]-&!@[C&+0:3].[O&+0&H2:4]>>[*:1]-[*:2].[*:3]-[*:4]",
    "hydrolysis of ether (ring)": "[C&+0:1]-[O&+0&H0:2]-&@[C&+0:3].[O&+0&H2:4]>>([*:1]-[*:2].[*:3]-[*:4])",
}
