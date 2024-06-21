"""Test tutorial exercises."""

import os
import tempfile

import doranet as dn


def test_tutorial_2():
    engine = dn.create_engine()
    acetone = engine.mol.rdkit("CC(=O)C")

    assert repr(acetone) == "MolDatBasic('CC(C)=O')"
    assert acetone.uid == "CC(C)=O"
    assert acetone.smiles == "CC(C)=O"

    aldol_condensation = engine.op.rdkit(
        "[O&+0:1]=[C&+0:2]-[C&+0;H2,H3:3].[C&+0:4]=[O&+0:5]>>[*:1]=[*:2]-[*:3]=[*:4].[*:5]"
    )

    assert (
        repr(aldol_condensation(acetone, acetone))
        == "((MolDatBasic('CC(=O)C=C(C)C'), MolDatBasic('O')), (MolDatBasic('CC(=O)C=C(C)C'), MolDatBasic('O')))"
    )


def test_tutorial_3():
    engine = dn.create_engine()

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

    assert (
        repr(aldol_condensation(acetone, acetone))
        == "((MolDatBasic('CC(=O)C=C(C)C'), MolDatBasic('O')), (MolDatBasic('CC(=O)C=C(C)C'), MolDatBasic('O')))"
    )

    assert (
        repr(aldol_condensation(butanone, acetone))
        == "((MolDatBasic('CC(=O)C(C)=C(C)C'), MolDatBasic('O')), (MolDatBasic('CCC(=O)C=C(C)C'), MolDatBasic('O')))"
    )

    assert aldol_condensation.compat(acetone, 0) is True

    assert aldol_condensation.compat(water, 1) is False

    assert aldol_condensation(acetone, water) == tuple()

    assert (
        aldol_condensation.smarts
        == "[O&+0:1]=[C&+0:2]-[C&+0;H2,H3:3].[C&+0:4]=[O&+0:5]>>[*:1]=[*:2]-[*:3]=[*:4].[*:5]"
    )

    assert len(aldol_condensation) == 2  # noqa: PLR2004

    try:
        aldol_condensation(acetone)
        raise AssertionError()
    except RuntimeError:
        assert True

    esterification = engine.op.rdkit(
        "[O&+0:1]=[C&+0:2]-[O&+0&H1:3].[O&+0&H1:4]>>[*:1]=[*:2]-[*:4].[*:3]"
    )

    assert (
        repr(esterification(hydroxyvaleric_acid, ethanol))
        == "((MolDatBasic('CCOC(=O)CCCCO'), MolDatBasic('O')),)"
    )

    esterification_intra = engine.op.rdkit(
        "([O&+0:1]=[C&+0:2]-[O&+0&H1:3].[O&+0&H1:4])>>[*:1]=[*:2]-[*:4].[*:3]"
    )

    assert (
        repr(esterification_intra(hydroxyvaleric_acid))
        == "((MolDatBasic('O=C1CCCCO1'), MolDatBasic('O')),)"
    )

    try:
        esterification_intra(hydroxyvaleric_acid, ethanol)
        raise AssertionError()
    except RuntimeError:
        assert True

    ester_hydrolysis_incorrect = engine.op.rdkit(
        "[O&+0:1]=[C&+0:2]-[O&+0&H0:3].[O&+0&H2:4]>>[*:1]=[*:2]-[*:4].[*:3]"
    )

    assert (
        repr(ester_hydrolysis_incorrect(methyl_butanoate, water))
        == "((MolDatBasic('CCCC(=O)O'), MolDatBasic('CO')),)"
    )

    assert (
        repr(ester_hydrolysis_incorrect(delta_valerolactone, water))
        == "((MolDatBasic('CCCCC(=O)O'), MolDatBasic('CCCCO')),)"
    )

    ester_hydrolysis_nonring = engine.op.rdkit(
        "[O&+0:1]=[C&+0:2]-&!@[O&+0&H0:3].[O&+0&H2:4]>>[*:1]=[*:2]-[*:4].[*:3]"
    )
    ester_hydrolysis_ring = engine.op.rdkit(
        "[O&+0:1]=[C&+0:2]-&@[O&+0&H0:3].[O&+0&H2:4]>>([*:1]=[*:2]-[*:4].[*:3])"
    )

    assert (
        repr(ester_hydrolysis_nonring(methyl_butanoate, water))
        == "((MolDatBasic('CCCC(=O)O'), MolDatBasic('CO')),)"
    )

    assert ester_hydrolysis_nonring(delta_valerolactone, water) == tuple()

    assert ester_hydrolysis_ring(methyl_butanoate, water) == tuple()

    assert (
        repr(ester_hydrolysis_ring(delta_valerolactone, water))
        == "((MolDatBasic('O=C(O)CCCCO'),),)"
    )


def test_tutorial_4():
    engine = dn.create_engine()

    water = engine.mol.rdkit("O")
    ethanol = engine.mol.rdkit("CCO")
    butanone = engine.mol.rdkit("CCC(C)=O")
    delta_valerolactone = engine.mol.rdkit("O=C1CCCCO1")
    hydroxyvaleric_acid = engine.mol.rdkit("O=C(O)CCCCO")

    ester_hydrolysis_ring = engine.op.rdkit(
        "[O&+0:1]=[C&+0:2]-&@[O&+0&H0:3].[O&+0&H2:4]>>([*:1]=[*:2]-[*:4].[*:3])"
    )

    network = engine.new_network()

    water_i = network.add_mol(water)

    assert water_i == 0

    assert repr(network.mols[water_i]) == "MolDatBasic('O')"

    assert water in network.mols
    assert water.uid in network.mols
    assert butanone not in network.mols

    assert network.mols.i(water.uid) == 0

    assert repr(network.mols[water_i]) == "MolDatBasic('O')"

    assert network.add_mol(water) == 0

    assert len(network.mols) == 1

    ester_hydrolysis_ring_i = network.add_op(ester_hydrolysis_ring)

    assert ester_hydrolysis_ring_i == 0

    network.add_mol(delta_valerolactone)

    hydroxyvaleric_acid = ester_hydrolysis_ring(delta_valerolactone, water)[0][
        0
    ]

    assert repr(hydroxyvaleric_acid) == "MolDatBasic('O=C(O)CCCCO')"

    water_i = network.mols.i(water.uid)
    delta_valerolactone_i = network.mols.i(delta_valerolactone.uid)
    ester_hydrolysis_ring_i = network.ops.i(ester_hydrolysis_ring.uid)

    hydroxyvaleric_acid = ester_hydrolysis_ring(delta_valerolactone, water)[0][
        0
    ]
    hydroxyvaleric_acid_i = network.add_mol(hydroxyvaleric_acid)

    reaction_i = network.add_rxn(
        operator=ester_hydrolysis_ring_i,
        reactants=(delta_valerolactone_i, water_i),
        products=(hydroxyvaleric_acid_i,),
    )

    assert reaction_i == 0
    assert network.rxns[0] == dn.interfaces.Reaction(0, (1, 0), (2,))

    assert network.consumers(0) == [0]
    assert network.producers(0) == []
    assert network.producers(2) == [0]

    with tempfile.TemporaryDirectory() as temp_dir:
        network.save_to_file("saved_network", temp_dir)
        assert os.listdir(temp_dir) == ["saved_network.pgnet"]

        network_loaded = engine.network_from_file("saved_network", temp_dir)
        assert (
            repr(list(network_loaded.mols))
            == "[MolDatBasic('O'), MolDatBasic('O=C1CCCCO1'), MolDatBasic('O=C(O)CCCCO')]"
        )

    assert network.compat_table(ester_hydrolysis_ring_i) == ([1], [0])

    network.add_mol(water, reactive=False)

    ethanol_i = network.add_mol(ethanol, meta={"is_alcohol": True})
    water_i = network.add_mol(
        water, meta={"is_alcohol": False, "solvent": True}
    )

    network.mols.set_meta(water_i, {"is_water": True, "is_toxic": False})

    assert network.mols.meta(water_i) == {
        "is_alcohol": False,
        "solvent": True,
        "is_water": True,
        "is_toxic": False,
    }
    assert network.mols.meta(water_i, ["is_alcohol"]) == {"is_alcohol": False}
    assert network.mols.meta([ethanol_i], ["is_alcohol"]) == (
        {"is_alcohol": True},
    )
    assert network.mols.meta([water_i, ethanol_i], keys=["is_alcohol"]) == (
        {"is_alcohol": False},
        {"is_alcohol": True},
    )
    assert network.mols.meta([ethanol_i, water_i], keys=["solvent"]) == (
        {},
        {"solvent": True},
    )
    assert network.mols.meta(keys=["is_alcohol"]) == (
        {"is_alcohol": False},
        {},
        {},
        {"is_alcohol": True},
    )
    assert network.mols.meta() == (
        {
            "is_alcohol": False,
            "solvent": True,
            "is_water": True,
            "is_toxic": False,
        },
        {},
        {},
        {"is_alcohol": True},
    )
