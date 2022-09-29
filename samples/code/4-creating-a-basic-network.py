import pickaxe_generic as pg

engine = pg.create_engine()

water = engine.mol.rdkit("O")
ethanol = engine.mol.rdkit("CCO")
acetone = engine.mol.rdkit("CC(C)=O")
butanone = engine.mol.rdkit("CCC(C)=O")
methyl_butanoate = engine.mol.rdkit("CCCC(=O)OC")
delta_valerolactone = engine.mol.rdkit("O=C1CCCCO1")
hydroxyvaleric_acid = engine.mol.rdkit("O=C(O)CCCCO")

aldol_condensation = engine.op.rdkit(
    "[O&+0:1]=[C&+0:2]-[C&+0;H2,H3:3].[C&+0:4]=[O&+0:5]>>[*:1]=[*:2]-[*:3]=[*:4].[*:5]"
)
ester_hydrolysis_nonring = engine.op.rdkit(
    "[O&+0:1]=[C&+0:2]-&!@[O&+0&H0:3].[O&+0&H2:4]>>[*:1]=[*:2]-[*:4].[*:3]"
)
ester_hydrolysis_ring = engine.op.rdkit(
    "[O&+0:1]=[C&+0:2]-&@[O&+0&H0:3].[O&+0&H2:4]>>([*:1]=[*:2]-[*:4].[*:3])"
)
esterification = engine.op.rdkit(
    "[O&+0:1]=[C&+0:2]-[O&+0&H1:3].[O&+0&H1:4]>>[*:1]=[*:2]-[*:4].[*:3]"
)
esterification_intra = engine.op.rdkit(
    "([O&+0:1]=[C&+0:2]-[O&+0&H1:3].[O&+0&H1:4])>>[*:1]=[*:2]-[*:4].[*:3]"
)

network = engine.new_network()

water_i = network.add_mol(water)

ester_hydrolysis_ring = engine.op.rdkit(
    "[O&+0:1]=[C&+0:2]-&@[O&+0&H0:3].[O&+0&H2:4]>>([*:1]=[*:2]-[*:4].[*:3])"
)
network.add_op(ester_hydrolysis_ring)

network.add_mol(delta_valerolactone)

hydroxyvaleric_acid = ester_hydrolysis_ring(delta_valerolactone, water)[0][0]

water_i = network.mols.i(water.uid)
delta_valerolactone_i = network.mols.i(delta_valerolactone.uid)
ester_hydrolysis_ring_i = network.ops.i(ester_hydrolysis_ring.uid)

hydroxyvaleric_acid = ester_hydrolysis_ring(delta_valerolactone, water)[0][0]
hydroxyvaleric_acid_i = network.add_mol(hydroxyvaleric_acid)

reaction_i = network.add_rxn(
    operator=ester_hydrolysis_ring_i,
    reactants=(delta_valerolactone_i, water_i),
    products=(hydroxyvaleric_acid_i,),
)

print(reaction_i)

print(network.rxns[reaction_i])

print(network.consumers(0))
print(network.producers(0))
print(network.producers(2))
