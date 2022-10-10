import pickaxe_generic as pg

engine = pg.create_engine()

ozone = engine.mol.rdkit("[O-][O+]=O")
alkyne = engine.mol.rdkit("CCC#C")
fluorine = engine.mol.rdkit("F")

op = engine.op.rdkit(
    "[C&+0:1]#[C&+0:2].[F&H1&+0:3].[F&H1&+0:4]>>[*:1]-[*:2](-[*:3])-[*:4]"
)

print(op(alkyne, fluorine, fluorine))
