"""Unit tests for network module."""

from typing import Iterable, Tuple

from network import (
    CartesianStrategy,
    Identifier,
    MolDatBasic,
    ObjectLibraryBasic,
    OpDatBasic,
    RxnDatBasic,
)
from pytest import fixture, mark
from rdkit.Chem import MolFromSmiles
from rdkit.Chem.rdchem import Mol as RDKitMol
from rdkit.Chem.rdChemReactions import ChemicalReaction as RDKitRxn
from rdkit.Chem.rdChemReactions import ReactionFromSmarts


class TestMolDatBasic:
    @mark.parametrize(
        "init",
        [
            MolFromSmiles("c1c(N(=O)=O)cccc1", sanitize=False),
            "c1c(N(=O)=O)cccc1",
            MolDatBasic("c1c(N(=O)=O)cccc1").blob,
        ],
    )
    def test_init(self, init):
        """Test initialization from different formats."""
        mol = MolDatBasic(init)
        rdkitmol = MolFromSmiles("O=[N+]([O-])c1ccccc1")
        assert hasattr(mol, "blob")
        assert mol.inchikey == "LQNUZADURLCDLV-UHFFFAOYSA-N"
        assert isinstance(mol.rdkitmol, RDKitMol)
        assert mol.smiles == "O=[N+]([O-])c1ccccc1"
        assert mol.uid == "O=[N+]([O-])c1ccccc1"

    def test_ordering(self):
        """Test comparison between molecule objects."""
        mol1 = MolDatBasic("O=[N+]([O-])c1ccccc1")
        mol2 = MolDatBasic("CCC(=O)[O-]")
        assert mol2 < mol1
        assert not mol1 < mol2

    def test_hash(self):
        """Test hash function."""
        mol = MolDatBasic("O=[N+]([O-])c1ccccc1")
        assert hash(mol.uid) == hash("O=[N+]([O-])c1ccccc1")

    @mark.parametrize(
        "init",
        [
            MolFromSmiles("c1c(N(=O)=O)cccc1", sanitize=False),
            "c1c(N(=O)=O)cccc1",
            MolDatBasic("c1c(N(=O)=O)cccc1").blob,
        ],
    )
    def test_equality(self, init):
        """Test equality of molecule objects."""
        mol = MolDatBasic(init)
        molequal = MolDatBasic("O=[N+]([O-])c1ccccc1")
        molunequal1 = MolDatBasic("c1c(N(=O)=O)cccc1", sanitize=False)
        molunequal2 = MolDatBasic("CC")
        assert mol == molequal
        assert mol != molunequal1
        assert mol != molunequal2


class TestOpDatBasic:
    @mark.parametrize(
        "init",
        [
            ReactionFromSmarts(
                "[C:1](=[O:2])-[OD1].[N!H0:3]>>[C:1](=[O:2])[N:3]"
            ),
            ReactionFromSmarts(
                "[N!H0:3].[C:1](=[O:2])-[OD1]>>[C:1](=[O:2])[N:3]"
            ),
            "[C:1](=[O:2])-[OD1].[N!H0:3]>>[C:1](=[O:2])[N:3]",
            "[N!H0:3].[C:1](=[O:2])-[OD1]>>[C:1](=[O:2])[N:3]",
            OpDatBasic("[C:1](=[O:2])-[OD1].[N!H0:3]>>[C:1](=[O:2])[N:3]").blob,
            OpDatBasic("[N!H0:3].[C:1](=[O:2])-[OD1]>>[C:1](=[O:2])[N:3]").blob,
        ],
    )
    def test_init(self, init):
        """Test initialization from different formats."""
        op = OpDatBasic(init)
        assert hasattr(op, "blob")
        assert isinstance(op.rdkitrxn, RDKitRxn)
        assert (
            op.smarts == "[C:1](=[O:2])-[O&D1].[N&!H0:3]>>[C:1](=[O:2])[N:3]"
            or op.smarts == "[N&!H0:3].[C:1](=[O:2])-[O&D1]>>[C:1](=[O:2])[N:3]"
        )
        assert op.uid == (
            ("O[CH:1]=[O:2]", "[NH3:3]"),
            ("[CH:1](=[O:2])[NH2:3]",),
        )

    def test_ordering(self):
        """Test comparison between operator objects."""
        op1 = OpDatBasic("[C:1](=[O:2])O.[N:3]>>[C:1](=[O:2])[N:3]")
        op2 = OpDatBasic("[CH1:1][OH:2].[OH][C:3]=[O:4]>>[C:1][O:2][C:3]=[O:4]")
        assert op1 < op2
        assert not op2 < op1
        assert not op1 == op2

    @mark.parametrize(
        "init",
        [
            ReactionFromSmarts(
                "[C:1](=[O:2])-[OD1].[N!H0:3]>>[C:1](=[O:2])[N:3]"
            ),
            "[C:1](=[O:2])-[OD1].[N!H0:3]>>[C:1](=[O:2])[N:3]",
            OpDatBasic("[C:1](=[O:2])-[OD1].[N!H0:3]>>[C:1](=[O:2])[N:3]").blob,
        ],
    )
    def test_equality(self, init):
        """Test equality of operator objects."""
        op = OpDatBasic(init)
        assert op == OpDatBasic(
            "[C:1](=[O:2])-[OD1].[N!H0:3]>>[C:1](=[O:2])[N:3]"
        )
        assert op == OpDatBasic(
            "[CH1:1](=[O:2])[OD1].[N:3]>>[C:1](=[O:2])[N:3]"
        )
        assert op != OpDatBasic("[C:1]=[O,N:2]>>[C:1][*:2]")
        assert op != OpDatBasic(
            "[CH1:1][OH:2].[OH][C:3]=[O:4]>>[C:1][O:2][C:3]=[O:4]"
        )

    def test_hash(self):
        """Test hash method."""
        assert isinstance(
            hash(OpDatBasic("[C:1]=[O,N:2]>>[C:1][*:2]").uid), int
        )

    def test_compat(self):
        """Test compat method."""
        op = OpDatBasic("[C:1](=[O:2])-[OD1].[N!H0:3]>>[C:1](=[O:2])[N:3]")
        mol1 = MolDatBasic("C(COC(=O)O)C(=O)O")
        mol2 = MolDatBasic("NC")
        mol3 = MolDatBasic("CC")
        mol4 = MolDatBasic("")
        assert op.compat(mol1, 0)
        assert not op.compat(mol1, 1)
        assert op.compat(mol2, 1)
        assert not op.compat(mol2, 0)
        assert not op.compat(mol3, 0)
        assert not op.compat(mol3, 1)
        assert not op.compat(mol4, 0)
        assert not op.compat(mol4, 1)

    def test_react_1(self):
        """Test react method."""
        op = OpDatBasic("[C:1](=[O:2])-[OD1].[N!H0:3]>>[C:1](=[O:2])[N:3]")
        mol1 = MolDatBasic("C(COC(=O)O)C(=O)O")
        mol2 = MolDatBasic("NC")
        rxn1, rxn2 = op((mol1, mol2))
        rxn1_expected = ("CNC(=O)OCCC(=O)O",)
        rxn2_expected = ("CNC(=O)CCOC(=O)O",)
        assert frozenset(
            (tuple(mol.uid for mol in rxn1), tuple(mol.uid for mol in rxn2))
        ) == frozenset((rxn1_expected, rxn2_expected))

    def test_react_2(self):
        op = OpDatBasic("[C:1][C:2][N:3]>>[C:1].[C:2][N:3]")
        mol = MolDatBasic("NCCC(C)CN")
        prod1_expected = frozenset(("CC(C)CN", "CN"))
        prod2_expected = frozenset(("CCCCN", "CN"))
        prods = op((mol,))
        assert frozenset(
            (
                frozenset(tuple(mol.uid for mol in prods[0])),
                frozenset(tuple(mol.uid for mol in prods[1])),
            )
        ) == frozenset((prod1_expected, prod2_expected))


class TestRxnDatBasic:
    @fixture(scope="class")
    def rxnpacket_gen(
        self,
    ) -> Tuple[Identifier, Tuple[Identifier, Identifier]]:
        op = OpDatBasic("[C:1](=[O:2])-[OD1].[N!H0:3]>>[C:1](=[O:2])[N:3]")
        rmol1 = MolDatBasic("C(COC(=O)O)C(=O)O")
        rmol2 = MolDatBasic("NC")
        return (op.uid, (rmol1.uid, rmol2.uid))

    @fixture(scope="class")
    def rxn1_gen(self, rxnpacket_gen) -> RxnDatBasic:
        pmol = MolDatBasic("CNC(=O)OCCC(=O)O")
        return RxnDatBasic(rxnpacket_gen[0], rxnpacket_gen[1], (pmol.uid,))

    @fixture(scope="class")
    def rxn1_reversed_gen(self, rxnpacket_gen) -> RxnDatBasic:
        op = OpDatBasic("[N!H0:3].[C:1](=[O:2])-[OD1]>>[C:1](=[O:2])[N:3]")
        pmol = MolDatBasic("CNC(=O)OCCC(=O)O")
        return RxnDatBasic(op.uid, rxnpacket_gen[1], (pmol.uid,))

    @fixture(scope="class")
    def rxn2_gen(self, rxnpacket_gen) -> RxnDatBasic:
        pmol = MolDatBasic("CNC(=O)CCOC(=O)O")
        return RxnDatBasic(rxnpacket_gen[0], rxnpacket_gen[1], (pmol.uid,))

    def test_init(self, rxn1_gen):
        """Test initialization of reaction from bytes."""
        rxn2 = RxnDatBasic(reaction=rxn1_gen.blob)
        assert rxn1_gen == rxn2

    def test_ordering_1(self, rxn1_gen, rxn2_gen):
        """Test comparison between operator objects."""
        assert rxn2_gen < rxn1_gen

    def test_ordering_2(self, rxn1_gen, rxn2_gen):
        assert not rxn1_gen < rxn2_gen

    def test_equality_1(self, rxn1_gen):
        """Test equality between operator objects."""
        assert rxn1_gen == rxn1_gen

    def test_equality_2(self, rxn2_gen):
        assert rxn2_gen == rxn2_gen

    def test_equality_3(self, rxn1_gen, rxn2_gen):
        assert not rxn1_gen == rxn2_gen

    def test_equality_4(self, rxn1_gen, rxn1_reversed_gen):
        assert rxn1_gen == rxn1_reversed_gen

    def test_hash(self, rxn1_gen):
        """Test hash method."""
        assert isinstance(hash(rxn1_gen.uid), int)


@mark.parametrize(
    "init",
    [
        [
            MolDatBasic(smiles)
            for smiles in (
                "CC",  # ethane
                "NC",  # hydrogen cyanide
                "CCC",  # propane
                "C=C",  # ethene
                "CBr",  # bromomethane
                "C#N",  # hydrocyanic acid
                "[Na+].[Cl-]",  # sodium chloride
                "CC(O)C",  # 2-propanol
                "CC(=O)C",  # 2-propanone
                "CC(CC)C",  # 2-methylbutane
                "CC(C)CC(=O)",  # 2-methylbutanal
                "c1c(N(=O)=O)cccc1",  # nitrobenzene
                "CC(C)(C)CC",  # 2,2-dimethylbutane
                "C=1CCCCC1",  # cyclohexene
                "c1ccccc1",  # benzene
                "C1OC1CC",  # ethyloxirane
                "c1cc2ccccc2cc1",  # naphthalene
                "CCC(=O)[O-]",  # ionized propanoic acid
                "c1cccc[n+]1CC(=O)O",  # 1-carboxylmethyl pyridinium
                "c1ccccc1N",
                "[nH]1ccc2c1cc(N)cc2",
            )
        ],
        (MolDatBasic("CC"), MolDatBasic("NC")),
        MolDatBasic("CC"),
        [],
        None,
    ],
    scope="class",
)
class TestObjectLibraryBasic_MolDatBasic:
    @fixture(scope="function")
    def mol_lib_gen(self, init):
        return ObjectLibraryBasic(init)

    @fixture(scope="class")
    def mol_add(self):
        return MolDatBasic("CCCC")

    def test_init(self, init):
        """Test initialization of object library."""
        lib = ObjectLibraryBasic(init)
        assert ObjectLibraryBasic(init) is not None
        if init is None:
            return
        elif isinstance(init, Iterable):
            for item in init:
                assert item in lib
        else:
            assert init in lib

    def test_add_1(self, mol_lib_gen, mol_add):
        """Test adding an object to the library."""
        assert mol_add not in mol_lib_gen
        mol_lib_gen.add(mol_add)
        assert mol_add in mol_lib_gen

    def test_add_2(self, mol_lib_gen, mol_add):
        mol2 = MolDatBasic("C1=CC=C2[C-]=CC=CC2=C1.[K+]")
        assert mol_add not in mol_lib_gen
        assert mol2 not in mol_lib_gen
        mol_lib_gen.add((mol_add, mol2))
        assert mol_add in mol_lib_gen
        assert mol2 in mol_lib_gen

    def test_keys_1(self, mol_lib_gen, mol_add):
        """Test adding an object and checking for its UID in keys."""
        assert mol_add.uid not in mol_lib_gen.ids()
        mol_lib_gen.add(mol_add)
        assert mol_add.uid in mol_lib_gen.ids()

    def test_keys_2(self, mol_lib_gen, mol_add):
        mol2 = MolDatBasic("C1=CC=C2[C-]=CC=CC2=C1.[K+]")
        assert mol_add.uid not in mol_lib_gen.ids()
        assert mol2.uid not in mol_lib_gen.ids()
        mol_lib_gen.add((mol_add, mol2))
        assert mol_add.uid in mol_lib_gen.ids()
        assert mol2.uid in mol_lib_gen.ids()

    def test_iter(self, mol_lib_gen):
        """Test iteration through the library."""
        for item in mol_lib_gen:
            assert isinstance(item, MolDatBasic)

    def test_len(self, init, mol_lib_gen):
        """Test length of library."""
        if isinstance(init, Iterable):
            assert len(mol_lib_gen) == len(init)
        elif init is not None:
            assert len(mol_lib_gen) == 1
        else:
            assert len(mol_lib_gen) == 0


@mark.parametrize(
    "init",
    [
        [
            OpDatBasic(smarts)
            for smarts in (
                "[C:1][C:1]>>[C:1]",
                "[C:1]=[O,N:2]>>[C:1][*:2]",
                "[CH1:1][OH:2].[OH][C:3]=[O:4]>>[C:1][O:2][C:3]=[O:4]",
            )
        ],
        (
            OpDatBasic("[C:1][C:1]>>[C:1]"),
            OpDatBasic("[C:1]=[O,N:2]>>[C:1][*:2]"),
        ),
        OpDatBasic("[C:1][C:1]>>[C:1]"),
        [],
        None,
    ],
    scope="class",
)
class TestObjectLibraryBasic_OpDatBasic:
    @fixture(scope="function")
    def op_lib_gen(self, init):
        return ObjectLibraryBasic(init)

    @fixture(scope="class")
    def op_add(self):
        return OpDatBasic("[*:1][Nh2:2]>>[*:1][Nh0:2](~[OD1])~[OD1]")

    def test_init(self, init):
        """Test initialization of object library."""
        lib = ObjectLibraryBasic(init)
        assert ObjectLibraryBasic(init) is not None
        if init is None:
            return
        elif isinstance(init, Iterable):
            for item in init:
                assert item in lib
        else:
            assert init in lib

    def test_add_1(self, op_lib_gen, op_add):
        """Test adding an object to the library."""
        assert op_add not in op_lib_gen
        op_lib_gen.add(op_add)
        assert op_add in op_lib_gen

    def test_add_2(self, op_lib_gen, op_add):
        op2 = OpDatBasic("[C:1](=[O:2])-[OD1].[N!H0:3]>>[C:1](=[O:2])[N:3]")
        assert op_add not in op_lib_gen
        assert op2 not in op_lib_gen
        op_lib_gen.add((op_add, op2))
        assert op_add in op_lib_gen
        assert op2 in op_lib_gen

    def test_keys_1(self, op_lib_gen, op_add):
        """Test adding an object and checking for its UID in keys."""
        assert op_add.uid not in op_lib_gen.ids()
        op_lib_gen.add(op_add)
        assert op_add.uid in op_lib_gen.ids()

    def test_keys_2(self, op_lib_gen, op_add):
        op2 = OpDatBasic("[C:1](=[O:2])-[OD1].[N!H0:3]>>[C:1](=[O:2])[N:3]")
        assert op_add.uid not in op_lib_gen.ids()
        assert op2.uid not in op_lib_gen.ids()
        op_lib_gen.add((op_add, op2))
        assert op_add.uid in op_lib_gen.ids()
        assert op2.uid in op_lib_gen.ids()

    def test_iter(self, op_lib_gen):
        """Test iteration through the library."""
        for item in op_lib_gen:
            assert isinstance(item, OpDatBasic)

    def test_len(self, init, op_lib_gen):
        """Test length of library."""
        if isinstance(init, Iterable):
            assert len(op_lib_gen) == len(init)
        elif init is not None:
            assert len(op_lib_gen) == 1
        else:
            assert len(op_lib_gen) == 0

    def test_call(self, init, op_lib_gen):
        """Test referencing objects by UID."""
        if isinstance(init, Iterable):
            for item in init:
                assert item == op_lib_gen[item.uid]
        elif init is not None:
            assert init == op_lib_gen[init.uid]


class TestCartesianStrategy:
    @mark.parametrize(
        "mols,ops",
        [
            (ObjectLibraryBasic(), ObjectLibraryBasic()),
            (ObjectLibraryBasic(MolDatBasic("CC")), ObjectLibraryBasic()),
            (
                ObjectLibraryBasic(),
                ObjectLibraryBasic(
                    OpDatBasic(
                        "[C:1](=[O:2])-[OD1].[N!H0:3]>>[C:1](=[O:2])[N:3]"
                    )
                ),
            ),
            (
                ObjectLibraryBasic(MolDatBasic("CC")),
                ObjectLibraryBasic(
                    OpDatBasic(
                        "[C:1](=[O:2])-[OD1].[N!H0:3]>>[C:1](=[O:2])[N:3]"
                    )
                ),
            ),
            (
                ObjectLibraryBasic((MolDatBasic("CC"), MolDatBasic("CCC"))),
                ObjectLibraryBasic(
                    (
                        OpDatBasic("[C:1][C:1]>>[C:1]"),
                        OpDatBasic(
                            "[C:1](=[O:2])-[OD1].[N!H0:3]>>[C:1](=[O:2])[N:3]"
                        ),
                    ),
                ),
            ),
        ],
    )
    def test_init(self, mols, ops):
        assert CartesianStrategy(mols, ops, ObjectLibraryBasic()) is not None

    @mark.parametrize(
        "mols_gen,ops_gen,num_mols,num_rxns",
        [
            (
                ObjectLibraryBasic(MolDatBasic("CC")),
                ObjectLibraryBasic(
                    OpDatBasic(
                        "[C:1](=[O:2])-[OD1].[N!H0:3]>>[C:1](=[O:2])[N:3]"
                    )
                ),
                1,
                0,
            ),
            (
                ObjectLibraryBasic((MolDatBasic("CC"), MolDatBasic("CCC"))),
                ObjectLibraryBasic(
                    (
                        OpDatBasic("[C:1][C:2]>>[C:1]"),
                        OpDatBasic(
                            "[C:1](=[O:2])-[OD1].[N!H0:3]>>[C:1](=[O:2])[N:3]"
                        ),
                    ),
                ),
                3,
                3,
            ),
            (
                ObjectLibraryBasic((MolDatBasic("CCC"), MolDatBasic("CCCCC"))),
                ObjectLibraryBasic(
                    (
                        OpDatBasic("[C:1][C:2][C:3]>>[C:2]"),
                        OpDatBasic(
                            "[C:1](=[O:2])-[OD1].[N!H0:3]>>[C:1](=[O:2])[N:3]"
                        ),
                    ),
                ),
                3,
                2,
            ),
            (
                ObjectLibraryBasic((MolDatBasic("CCC"), MolDatBasic("CCCCC"))),
                ObjectLibraryBasic(
                    (
                        OpDatBasic("[CH3:1][C:2][C:3]>>[CH3:1][C:2].[C:3]"),
                        OpDatBasic(
                            "[C:1](=[O:2])-[OD1].[N!H0:3]>>[C:1](=[O:2])[N:3]"
                        ),
                    ),
                ),
                4,
                2,
            ),
            (
                ObjectLibraryBasic((MolDatBasic("CCN"), MolDatBasic("CCCCC"))),
                ObjectLibraryBasic(
                    (OpDatBasic("[C:1][C:2]>>[C:1]"),),
                ),
                7,
                12,
            ),
            (
                ObjectLibraryBasic((MolDatBasic("CCN"), MolDatBasic("CCCCC"))),
                ObjectLibraryBasic(
                    (OpDatBasic("[C:1][C:2]>>[C:1].[C:2]"),),
                ),
                7,
                7,
            ),
            (
                ObjectLibraryBasic((MolDatBasic("CCCC(O)C"),)),
                ObjectLibraryBasic(
                    (OpDatBasic("[CH2:3]-[CH1:2]-[OH]>>[*:3]=[*:2]"),),
                ),
                2,
                1,
            ),
            (
                ObjectLibraryBasic(
                    (
                        MolDatBasic("C(COC(=O)O)C(=O)O"),
                        MolDatBasic("CC(=O)O"),
                        MolDatBasic("NC"),
                    )
                ),
                ObjectLibraryBasic(
                    (
                        OpDatBasic(
                            "[C:1](=[O:2])-[OD1].[N!H0:3][CH3]>>[C:1](=[O:2])[N:3]"
                        ),
                    ),
                ),
                7,
                5,
            ),
            (
                ObjectLibraryBasic(
                    (
                        MolDatBasic("C(COC(=O)O)C(=O)O"),
                        MolDatBasic("CC(=O)O"),
                        MolDatBasic("NC"),
                    )
                ),
                ObjectLibraryBasic(
                    (
                        OpDatBasic(
                            "[C:1](=[O:2])-[OD1].[N!H0:3][CH3]>>[C:1](=[O:2])[N:3]"
                        ),
                        OpDatBasic("[C:1][C:2]>>[C:1].[C:2]"),
                    ),
                ),
                14,
                20,
            ),
        ],
    )
    def test_expand_1(self, mols_gen, ops_gen, num_mols, num_rxns):
        rxn_lib = ObjectLibraryBasic()
        strat = CartesianStrategy(mols_gen, ops_gen, rxn_lib)
        strat.expand()
        print(mols_gen)
        assert len(mols_gen) == num_mols
        assert len(rxn_lib) == num_rxns

    @mark.parametrize(
        "mols_gen,ops_gen,num_mols",
        [
            (
                ObjectLibraryBasic(MolDatBasic("CC")),
                ObjectLibraryBasic(
                    OpDatBasic(
                        "[C:1](=[O:2])-[OD1].[N!H0:3]>>[C:1](=[O:2])[N:3]"
                    )
                ),
                1,
            ),
            (
                ObjectLibraryBasic((MolDatBasic("CC"), MolDatBasic("CCC"))),
                ObjectLibraryBasic(
                    (
                        OpDatBasic("[C:1][C:2]>>[C:1]"),
                        OpDatBasic(
                            "[C:1](=[O:2])-[OD1].[N!H0:3]>>[C:1](=[O:2])[N:3]"
                        ),
                    ),
                ),
                3,
            ),
            (
                ObjectLibraryBasic((MolDatBasic("CCC"), MolDatBasic("CCCCC"))),
                ObjectLibraryBasic(
                    (
                        OpDatBasic("[C:1][C:2][C:3]>>[C:2]"),
                        OpDatBasic(
                            "[C:1](=[O:2])-[OD1].[N!H0:3]>>[C:1](=[O:2])[N:3]"
                        ),
                    ),
                ),
                3,
            ),
            (
                ObjectLibraryBasic((MolDatBasic("CCC"), MolDatBasic("CCCCC"))),
                ObjectLibraryBasic(
                    (
                        OpDatBasic("[CH3:1][C:2][C:3]>>[CH3:1][C:2].[C:3]"),
                        OpDatBasic(
                            "[C:1](=[O:2])-[OD1].[N!H0:3]>>[C:1](=[O:2])[N:3]"
                        ),
                    ),
                ),
                4,
            ),
            (
                ObjectLibraryBasic((MolDatBasic("CCN"), MolDatBasic("CCCCC"))),
                ObjectLibraryBasic(
                    (OpDatBasic("[C:1][C:2]>>[C:1]"),),
                ),
                7,
            ),
            (
                ObjectLibraryBasic((MolDatBasic("CCN"), MolDatBasic("CCCCC"))),
                ObjectLibraryBasic(
                    (OpDatBasic("[C:1][C:2]>>[C:1].[C:2]"),),
                ),
                7,
            ),
            (
                ObjectLibraryBasic((MolDatBasic("CCCC(O)C"),)),
                ObjectLibraryBasic(
                    (OpDatBasic("[CH2:3]-[CH1:2]-[OH]>>[*:3]=[*:2]"),),
                ),
                2,
            ),
            (
                ObjectLibraryBasic(
                    (
                        MolDatBasic("C(COC(=O)O)C(=O)O"),
                        MolDatBasic("CC(=O)O"),
                        MolDatBasic("NC"),
                    )
                ),
                ObjectLibraryBasic(
                    (
                        OpDatBasic(
                            "[C:1](=[O:2])-[OD1].[N!H0:3][CH3]>>[C:1](=[O:2])[N:3]"
                        ),
                    ),
                ),
                7,
            ),
            (
                ObjectLibraryBasic(
                    (
                        MolDatBasic("C(COC(=O)O)C(=O)O"),
                        MolDatBasic("CC(=O)O"),
                        MolDatBasic("NC"),
                    )
                ),
                ObjectLibraryBasic(
                    (
                        OpDatBasic(
                            "[C:1](=[O:2])-[OD1].[N!H0:3][CH3]>>[C:1](=[O:2])[N:3]"
                        ),
                        OpDatBasic("[C:1][C:2]>>[C:1].[C:2]"),
                    ),
                ),
                14,
            ),
        ],
    )
    def test_expand_2(self, mols_gen, ops_gen, num_mols):
        rxn_lib = ObjectLibraryBasic()
        strat = CartesianStrategy(mols_gen, ops_gen, rxn_lib)
        strat.expand(max_mols=20)
        print(mols_gen)
        assert len(mols_gen) == 20 or len(mols_gen) == num_mols

    @mark.parametrize(
        "mols_gen,ops_gen,num_rxns",
        [
            (
                ObjectLibraryBasic(MolDatBasic("CC")),
                ObjectLibraryBasic(
                    OpDatBasic(
                        "[C:1](=[O:2])-[OD1].[N!H0:3]>>[C:1](=[O:2])[N:3]"
                    )
                ),
                0,
            ),
            (
                ObjectLibraryBasic((MolDatBasic("CC"), MolDatBasic("CCC"))),
                ObjectLibraryBasic(
                    (
                        OpDatBasic("[C:1][C:2]>>[C:1]"),
                        OpDatBasic(
                            "[C:1](=[O:2])-[OD1].[N!H0:3]>>[C:1](=[O:2])[N:3]"
                        ),
                    ),
                ),
                3,
            ),
            (
                ObjectLibraryBasic((MolDatBasic("CCC"), MolDatBasic("CCCCC"))),
                ObjectLibraryBasic(
                    (
                        OpDatBasic("[C:1][C:2][C:3]>>[C:2]"),
                        OpDatBasic(
                            "[C:1](=[O:2])-[OD1].[N!H0:3]>>[C:1](=[O:2])[N:3]"
                        ),
                    ),
                ),
                2,
            ),
            (
                ObjectLibraryBasic((MolDatBasic("CCC"), MolDatBasic("CCCCC"))),
                ObjectLibraryBasic(
                    (
                        OpDatBasic("[CH3:1][C:2][C:3]>>[CH3:1][C:2].[C:3]"),
                        OpDatBasic(
                            "[C:1](=[O:2])-[OD1].[N!H0:3]>>[C:1](=[O:2])[N:3]"
                        ),
                    ),
                ),
                2,
            ),
            (
                ObjectLibraryBasic((MolDatBasic("CCN"), MolDatBasic("CCCCC"))),
                ObjectLibraryBasic(
                    (OpDatBasic("[C:1][C:2]>>[C:1]"),),
                ),
                12,
            ),
            (
                ObjectLibraryBasic((MolDatBasic("CCN"), MolDatBasic("CCCCC"))),
                ObjectLibraryBasic(
                    (OpDatBasic("[C:1][C:2]>>[C:1].[C:2]"),),
                ),
                7,
            ),
            (
                ObjectLibraryBasic((MolDatBasic("CCCC(O)C"),)),
                ObjectLibraryBasic(
                    (OpDatBasic("[CH2:3]-[CH1:2]-[OH]>>[*:3]=[*:2]"),),
                ),
                1,
            ),
            (
                ObjectLibraryBasic(
                    (
                        MolDatBasic("C(COC(=O)O)C(=O)O"),
                        MolDatBasic("CC(=O)O"),
                        MolDatBasic("NC"),
                    )
                ),
                ObjectLibraryBasic(
                    (
                        OpDatBasic(
                            "[C:1](=[O:2])-[OD1].[N!H0:3][CH3]>>[C:1](=[O:2])[N:3]"
                        ),
                    ),
                ),
                5,
            ),
            (
                ObjectLibraryBasic(
                    (
                        MolDatBasic("C(COC(=O)O)C(=O)O"),
                        MolDatBasic("CC(=O)O"),
                        MolDatBasic("NC"),
                    )
                ),
                ObjectLibraryBasic(
                    (
                        OpDatBasic(
                            "[C:1](=[O:2])-[OD1].[N!H0:3][CH3]>>[C:1](=[O:2])[N:3]"
                        ),
                        OpDatBasic("[C:1][C:2]>>[C:1].[C:2]"),
                    ),
                ),
                20,
            ),
        ],
    )
    def test_expand_3(self, mols_gen, ops_gen, num_rxns):
        rxn_lib = ObjectLibraryBasic()
        strat = CartesianStrategy(mols_gen, ops_gen, rxn_lib)
        strat.expand(max_rxns=5)
        print(mols_gen)
        assert len(rxn_lib) == 5 or len(rxn_lib) == num_rxns

    def test_expand_4(self):
        mol_lib = ObjectLibraryBasic((MolDatBasic("CCCC"),))
        op_lib = ObjectLibraryBasic((OpDatBasic("[C:1][C:2]>>[C:1]"),))
        rxn_lib = ObjectLibraryBasic()
        strat = CartesianStrategy(mol_lib, op_lib, rxn_lib)
        assert len(mol_lib) == 1 and len(rxn_lib) == 0
        strat.refresh()
        strat.expand(num_gens=1)
        assert len(mol_lib) == 4 and len(rxn_lib) == 3
        strat.refresh()
        strat.expand(num_gens=1)
        assert len(mol_lib) == 4 and len(rxn_lib) == 6

    def test_refresh(self):
        mol_lib = ObjectLibraryBasic((MolDatBasic("CCC"),))
        op_lib = ObjectLibraryBasic((OpDatBasic("[C:1][C:2]>>[C:1]"),))
        rxn_lib = ObjectLibraryBasic()
        strat = CartesianStrategy(mol_lib, op_lib, rxn_lib)
        strat.expand()
        assert len(mol_lib) == 3
        assert len(rxn_lib) == 3
        mol_lib.add(MolDatBasic("CCCCC"))
        strat.refresh()
        strat.expand()
        assert len(mol_lib) == 5
        assert len(rxn_lib) == 10
