"""Contains classes which define and implement molecule-operator data units."""

import builtins
import collections.abc
import io
import pickle
import typing

import rdkit
import rdkit.Chem
import rdkit.Chem.rdchem
import rdkit.Chem.rdChemReactions
import rdkit.Chem.rdinchi
import rdkit.Chem.rdmolfiles
import rdkit.Chem.rdmolops

from doranet import interfaces

# some code to make loads more safe to arbitrary code execution
# necessary since external data from a database may be input
# if you are having issues, add relevant classes to _safe_%module%_classes
# or if modules not in builtins required, add other if clauses
_safe_builtins_classes: frozenset[str] = frozenset(
    {
        "frozenset",
        "tuple",
    }
)


class __SafeUnpickler(pickle.Unpickler):
    def find_class(self, module: str, name: str) -> typing.Any:
        if module == "builtins" and name in _safe_builtins_classes:
            return getattr(builtins, name)
        raise pickle.UnpicklingError(f"global '{module}.{name}' is forbidden")


def _loads(string_in: bytes) -> typing.Any:
    return __SafeUnpickler(io.BytesIO(string_in)).load()

    # @final
    # def __getstate__(self) -> bytes:
    #     """
    #     Serializes object based on blob property.
    #     """
    #     return self.blob

    # @abstractmethod
    # def __setstate__(self, data: bytes) -> None:
    #     """
    #     Deserializes object from blob.
    #     """


@typing.final
class MolDatBasicV1(interfaces.MolDatRDKit):
    """
    Version of MolDatRDKit which caches only SMILES and blob.

    Speeds: 5
    """

    __slots__ = ("_blob", "_smiles")
    _blob: bytes
    _smiles: str

    def __init__(
        self,
        molecule: rdkit.Chem.rdchem.Mol | str | bytes,
        sanitize: bool = True,
        neutralize: bool = False,
    ) -> None:
        rdkitmol = self._processinput(molecule, sanitize, neutralize)
        self._buildfrommol(rdkitmol)

    def _buildfrommol(self, in_val: rdkit.Chem.rdchem.Mol) -> None:
        self._blob = in_val.ToBinary()
        self._smiles = rdkit.Chem.rdmolfiles.MolToSmiles(in_val)

    @property
    def blob(self) -> bytes:
        return self._blob

    @property
    def inchikey(self) -> str:
        return rdkit.Chem.rdinchi.MolToInchiKey(self.rdkitmol)

    @property
    def rdkitmol(self) -> rdkit.Chem.rdchem.Mol:
        return rdkit.Chem.rdchem.Mol(self._blob)

    @property
    def smiles(self) -> str:
        return self._smiles

    @property
    def uid(self) -> interfaces.Identifier:
        return self._smiles

    def __repr__(self) -> str:
        return f"MolDatBasic('{self.smiles}')"


@typing.final
class MolDatBasicV2(interfaces.MolDatRDKit):
    """
    Version of MolDatRDKit which caches all values.

    Speeds: 1,3,4,6
    """

    __slots__ = ("_blob", "_inchikey", "_rdkitmol", "_smiles")
    _blob: typing.Optional[bytes]
    _inchikey: typing.Optional[str]
    _rdkitmol: rdkit.Chem.rdchem.Mol
    _smiles: str

    def __init__(
        self,
        molecule: rdkit.Chem.rdchem.Mol | str | bytes,
        sanitize: bool = True,
        neutralize: bool = False,
    ) -> None:
        self._blob = None
        self._inchikey = None
        rdkitmol = self._processinput(molecule, sanitize, neutralize)
        self._buildfrommol(rdkitmol)

    def _buildfrommol(self, in_val: rdkit.Chem.rdchem.Mol) -> None:
        self._rdkitmol = in_val
        self._smiles = rdkit.Chem.rdmolfiles.MolToSmiles(in_val)

    @property
    def blob(self) -> bytes:
        if self._blob is None:
            self._blob = self._rdkitmol.ToBinary()
        return self._blob

    @property
    def inchikey(self) -> str:
        if self._inchikey is None:
            self._inchikey = rdkit.Chem.rdinchi.MolToInchiKey(self._rdkitmol)
        return self._inchikey

    @property
    def rdkitmol(self) -> rdkit.Chem.rdchem.Mol:
        return self._rdkitmol

    @property
    def smiles(self) -> str:
        return self._smiles

    @property
    def uid(self) -> str:
        return self._smiles

    def __repr__(self) -> str:
        return f'MolDatBasic("{self.smiles}")'


@typing.final
class OpDatBasic(interfaces.OpDatRDKit):
    """
    Minimal class implementing the OpDatRDKit interface.

    Classes implementing this interface manage information about a single
    rdkit-compatible SMARTS operator.

    Parameters
    ----------
    operator : rdkit.Chem.rdChemReactions.ChemicalReaction | str | bytes
        SMARTS string which is used to generate operator data, otherwise some
        encoding of relevant data.

    Attributes
    ----------
    blob : bytes
        Binary representation of operator.
    rdkitrxn : rdkit.Chem.rdChemReactions.ChemicalReaction
        RDKit reaction object.
    smarts : str
        SMARTS string representing operator.
    uid : Tuple[Tuple[str, ...], Tuple[str, ...]]
        Unique identifier of object, in this case based on the SMARTS string.
    """

    __slots__ = (
        "_blob",
        "_engine",
        "_kekulize",
        "_drop_errors",
        "_rdkitrxn",
        "_smarts",
        "_templates",
        "_uid",
    )

    _rdkitrxn: rdkit.Chem.rdChemReactions.ChemicalReaction
    _templates: typing.Optional[tuple[rdkit.Chem.rdchem.Mol, ...]]
    _engine: interfaces.NetworkEngine

    _blob: typing.Optional[bytes]
    _smarts: typing.Optional[str]
    _uid: typing.Optional[tuple[tuple[str, ...], tuple[str, ...]]]
    _kekulize: bool
    _drop_errors: bool

    def __init__(
        self,
        operator: rdkit.Chem.rdchem.Mol | str | bytes,
        engine: interfaces.NetworkEngine,
        kekulize: bool = False,
        drop_errors: bool = False,
    ) -> None:
        if isinstance(operator, rdkit.Chem.rdChemReactions.ChemicalReaction):
            self._rdkitrxn = operator
            self._kekulize = kekulize
            self._drop_errors = drop_errors
        elif isinstance(operator, str):
            self._rdkitrxn = rdkit.Chem.rdChemReactions.ReactionFromSmarts(
                operator
            )
            self._kekulize = kekulize
            self._drop_errors = drop_errors
            # SanitizeRxn(self._rdkitrxn)
        elif isinstance(operator, bytes):
            self._rdkitrxn, self._kekulize, self._drop_errors = _loads(operator)
            self._blob = operator
        else:
            raise NotImplementedError("Invalid operator type")
        self._templates = None
        self._engine = engine
        self._blob = None
        self._smarts = None
        self._uid = None

    @property
    def blob(self) -> bytes:
        if self._blob is None:
            self._blob = pickle.dumps(
                (self.rdkitrxn, self._kekulize, self._drop_errors)
            )
        return self._blob

    @property
    def rdkitrxn(self) -> rdkit.Chem.rdChemReactions.ChemicalReaction:
        return self._rdkitrxn

    @property
    def smarts(self) -> str:
        if self._smarts is None:
            self._smarts = rdkit.Chem.rdChemReactions.ReactionToSmarts(
                self._rdkitrxn
            )
        return self._smarts

    @property
    def uid(self) -> str:
        if self._smarts is None:
            self._smarts = rdkit.Chem.rdChemReactions.ReactionToSmarts(
                self._rdkitrxn
            )
        return self._smarts

    def compat(self, mol: interfaces.MolDatBase, arg: int) -> bool:
        if self._templates is None:
            self._templates = self._build_templates()
        if isinstance(mol, interfaces.MolDatRDKit):
            tempmol = mol.rdkitmol
            if self._kekulize:
                tempmol = rdkit.Chem.rdchem.Mol(tempmol)
                rdkit.Chem.rdmolops.Kekulize(tempmol, clearAromaticFlags=True)
            return tempmol.HasSubstructMatch(
                self._templates[arg], useChirality=True
            )
        else:
            return False

    def _build_templates(self) -> tuple[rdkit.Chem.rdchem.Mol, ...]:
        return tuple(self._rdkitrxn.GetReactants())

    def _attempt_reaction(
        self, mols: collections.abc.Iterable[rdkit.Chem.rdchem.Mol]
    ) -> collections.abc.Iterable[rdkit.Chem.rdchem.Mol]:
        try:
            return self._rdkitrxn.RunReactants(mols, maxProducts=0)
        except Exception as err:
            print(type(err))
            raise err

    def _process_new_mol(self, mol) -> interfaces.MolDatRDKit | None:
        if self._drop_errors:
            try:
                return self._engine.mol.rdkit(mol)
            except Exception:
                return None
        return self._engine.mol.rdkit(mol)

    def __call__(
        self, *reactants: interfaces.MolDatBase
    ) -> tuple[tuple[interfaces.MolDatBase, ...], ...]:
        rdkitmols: list[rdkit.Chem.rdchem.Mol] = [
            reactant.rdkitmol
            for reactant in reactants
            if isinstance(reactant, interfaces.MolDatRDKit)
        ]
        if self._kekulize:
            rdkitmols = [
                rdkit.Chem.rdchem.Mol(rdkitmol) for rdkitmol in rdkitmols
            ]
            for rdkitmol in rdkitmols:
                rdkit.Chem.rdmolops.Kekulize(rdkitmol, clearAromaticFlags=True)
        try:
            return tuple(
                product_set  # type: ignore [misc]
                for product_set in (
                    tuple(
                        self._process_new_mol(product) for product in products
                    )
                    for products in self._rdkitrxn.RunReactants(
                        rdkitmols, maxProducts=0
                    )
                )
                if None not in product_set
            )
        except Exception as err:
            raise RuntimeError(
                f"Error occurred when using operator {self} on {reactants}"
            ) from err

    def __len__(self) -> int:
        return (
            self._rdkitrxn.GetNumReactantTemplates()
            + self._rdkitrxn.GetNumAgentTemplates()
        )

    def __lt__(self, other: object) -> bool:
        if isinstance(other, OpDatBasic):
            return self.uid < other.uid
        raise NotImplementedError("Comparison not implemented")

    def __repr__(self) -> str:
        sval = f"OpDatBasic({repr(self.uid)})"
        if self._kekulize:
            sval = sval + ",kekulize"
        return sval
