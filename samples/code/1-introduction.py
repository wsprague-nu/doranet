import pickaxe_generic as pg

engine = pg.create_engine()

acetone = engine.mol.rdkit("CC(=O)C")
print(acetone)

cond_aldol = engine.op.rdkit(
    "[O&+0:1]=[C&+0:2]-[C&+0;H2,H3:3].[C&+0:4]=[O&+0:5]>>[*:1]=[*:2]-[*:3]=[*:4].[*:5]"
)

products = cond_aldol(acetone, acetone)
print(products)
