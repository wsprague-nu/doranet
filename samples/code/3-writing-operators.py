import pickaxe_generic as pg

engine = pg.create_engine()

water = engine.mol.rdkit("O")
acetone = engine.mol.rdkit("CC(C)=O")
butanone = engine.mol.rdkit("CCC(C)=O")
succinaldehyde = engine.mol.rdkit("O=CCCC=O")
tetrahydrofuran = engine.mol.rdkit("OC1CCC(O)O1")

aldol_condensation = engine.op.rdkit(
    "[O&+0:1]=[C&+0:2]-[C&+0;H2,H3:3].[C&+0:4]=[O&+0:5]>>[*:1]=[*:2]-[*:3]=[*:4].[*:5]"
)

print(aldol_condensation.smarts)
