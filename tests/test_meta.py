"""Test metadata updates."""

import doranet as dn


def test_meta_update():
    engine = dn.create_engine()
    network = engine.new_network()
    network.add_mol(engine.mol.rdkit("C#C"), meta={"gen": 0})
    network.add_op(engine.op.rdkit("[C:1]#[C:2]>>[*:1]=[*:2]"))
    network.add_op(engine.op.rdkit("[C:1]#[C:2]>>[*:1]-[*:2]"))
    network.add_op(engine.op.rdkit("[C:1]=[C:2]>>[*:1]-[*:2]"))
    network.add_op(engine.op.rdkit("[C:1]-[C:2]>>[*:1].[*:2]"))
    network.add_op(engine.op.rdkit("[C:1]-[C:2]>>[*:1]=[*:2]"))

    strat = engine.strat.cartesian(network)

    strat.expand(reaction_plan=dn.metacalc.GenerationCalculator("gen"))

    assert network.mols.meta(dn.interfaces.MolIndex(0), ("gen",))["gen"] == 0
    assert network.mols.meta(dn.interfaces.MolIndex(1), ("gen",))["gen"] == 1
    assert network.mols.meta(dn.interfaces.MolIndex(2), ("gen",))["gen"] == 1
    assert network.mols.meta(dn.interfaces.MolIndex(3), ("gen",))["gen"] == 2  # noqa: PLR2004
