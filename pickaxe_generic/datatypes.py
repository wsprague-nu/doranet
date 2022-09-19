"""
Contains classes which define and implement molecule-operator data units.

Classes:

    Identifier(Protocol)
    DataUnit
      MolDatBase
        MolDatRDKit
          MolDatBasicV1*
          MolDatBasicV2*
      OpDatBase
        OpDatRDKit
          OpDatBasic*
      RxnDatBase
        RxnDatBasic*
"""

import builtins
import collections.abc
import io
import pickle
import typing

from rdkit.Chem import Mol as BuildMol
from rdkit.Chem import MolToSmiles
from rdkit.Chem.inchi import MolToInchiKey
from rdkit.Chem.rdchem import AtomValenceException, KekulizeException
from rdkit.Chem.rdchem import Mol as RDKitMol
from rdkit.Chem.rdChemReactions import ChemicalReaction as RDKitRxn
from rdkit.Chem.rdChemReactions import ReactionFromSmarts, ReactionToSmarts
from rdkit.Chem.rdmolops import Kekulize

from . import interfaces

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


def loads(string_in: bytes) -> typing.Any:
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
        molecule: RDKitMol | str | bytes,
        sanitize: bool = True,
        neutralize: bool = False,
    ) -> None:
        rdkitmol = self._processinput(molecule, sanitize, neutralize)
        self._buildfrommol(rdkitmol)

    def _buildfrommol(self, in_val: RDKitMol) -> None:
        self._blob = in_val.ToBinary()
        self._smiles = MolToSmiles(in_val)

    @property
    def blob(self) -> bytes:
        return self._blob

    @property
    def inchikey(self) -> str:
        return MolToInchiKey(self.rdkitmol)

    @property
    def rdkitmol(self) -> RDKitMol:
        return BuildMol(self._blob)

    @property
    def smiles(self) -> str:
        return self._smiles

    @property
    def uid(self) -> interfaces.Identifier:
        return self._smiles

    def __repr__(self) -> str:
        return f"MolDatBasic('{self.smiles}')"

    def __setstate__(self, arg: bytes) -> None:
        self._blob = arg
        self._smiles = MolToSmiles(BuildMol(arg))


@typing.final
class MolDatBasicV2(interfaces.MolDatRDKit):
    """
    Version of MolDatRDKit which caches all values.

    Speeds: 1,3,4,6
    """

    __slots__ = ("_blob", "_inchikey", "_rdkitmol", "_smiles")
    _blob: typing.Optional[bytes]
    _inchikey: typing.Optional[str]
    _rdkitmol: RDKitMol
    _smiles: str

    def __init__(
        self,
        molecule: RDKitMol | str | bytes,
        sanitize: bool = True,
        neutralize: bool = False,
    ) -> None:
        self._blob = None
        self._inchikey = None
        rdkitmol = self._processinput(molecule, sanitize, neutralize)
        self._buildfrommol(rdkitmol)

    def _buildfrommol(self, in_val: RDKitMol) -> None:
        self._rdkitmol = in_val
        self._smiles = MolToSmiles(in_val)

    @property
    def blob(self) -> bytes:
        if self._blob is None:
            self._blob = self._rdkitmol.ToBinary()
        return self._blob

    @property
    def inchikey(self) -> str:
        if self._inchikey is None:
            self._inchikey = MolToInchiKey(self._rdkitmol)
        return self._inchikey

    @property
    def rdkitmol(self) -> RDKitMol:
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
    operator : RDKitMol | str | bytes
        SMARTS string which is used to generate operator data, otherwise some
        encoding of relevant data.

    Attributes
    ----------
    blob : bytes
        Binary representation of operator.
    rdkitrxn : RDKitRxn
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
        "_rdkitrxn",
        "_smarts",
        "_templates",
        "_uid",
    )

    _rdkitrxn: RDKitRxn
    _templates: typing.Optional[tuple[RDKitMol, ...]]
    _engine: interfaces.NetworkEngine

    _blob: typing.Optional[bytes]
    _smarts: typing.Optional[str]
    _uid: typing.Optional[tuple[tuple[str, ...], tuple[str, ...]]]
    _kekulize: bool

    def __init__(
        self,
        operator: RDKitMol | str | bytes,
        engine: interfaces.NetworkEngine,
        kekulize_before_operation: bool = False,
    ) -> None:
        if isinstance(operator, RDKitRxn):
            self._rdkitrxn = operator
            self._kekulize = kekulize_before_operation
        elif isinstance(operator, str):
            self._rdkitrxn = ReactionFromSmarts(operator)
            self._kekulize = kekulize_before_operation
            # SanitizeRxn(self._rdkitrxn)
        elif isinstance(operator, bytes):
            self._rdkitrxn, self._kekulize = loads(operator)
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
            self._blob = pickle.dumps((self.rdkitrxn, self._kekulize))
        return self._blob

    @property
    def rdkitrxn(self) -> RDKitRxn:
        return self._rdkitrxn

    @property
    def smarts(self) -> str:
        if self._smarts is None:
            self._smarts = ReactionToSmarts(self._rdkitrxn)
        return self._smarts

    @property
    def uid(self) -> str:
        if self._smarts is None:
            self._smarts = ReactionToSmarts(self._rdkitrxn)
        return self._smarts

    def compat(self, mol: interfaces.MolDatBase, arg: int) -> bool:
        if self._templates is None:
            self._templates = self._build_templates()
        if isinstance(mol, interfaces.MolDatRDKit):
            tempmol = mol.rdkitmol
            if self._kekulize:
                tempmol = BuildMol(tempmol)
                Kekulize(tempmol, clearAromaticFlags=True)
            return tempmol.HasSubstructMatch(
                self._templates[arg], useChirality=True
            )
        else:
            return False

    def _build_templates(self) -> tuple[RDKitMol, ...]:
        return tuple(self._rdkitrxn.GetReactants())

    def _attempt_reaction(
        self, mols: collections.abc.Iterable[RDKitMol]
    ) -> collections.abc.Iterable[RDKitMol]:
        try:
            return self._rdkitrxn.RunReactants(mols, maxProducts=0)
        except Exception as err:
            print(type(err))
            raise err

    def __call__(
        self, reactants: collections.abc.Sequence[interfaces.MolDatBase]
    ) -> tuple[tuple[interfaces.MolDatBase, ...], ...]:
        rdkitmols: list[RDKitMol] = [
            reactant.rdkitmol
            for reactant in reactants
            if isinstance(reactant, interfaces.MolDatRDKit)
        ]
        if self._kekulize:
            rdkitmols = [BuildMol(rdkitmol) for rdkitmol in rdkitmols]
            for rdkitmol in rdkitmols:
                Kekulize(rdkitmol, clearAromaticFlags=True)
        try:
            return tuple(
                tuple(self._engine.Mol(product) for product in products)
                for products in self._rdkitrxn.RunReactants(
                    rdkitmols, maxProducts=0
                )
            )
        except AtomValenceException as err:
            raise ValueError(
                f"Error occurred when using operator {self} on {reactants}"
            ) from err
        except KekulizeException as err:
            raise ValueError(
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


@typing.final
class RxnDatBasic(interfaces.RxnDatBase):
    """
    Minimal class implementing the RxnDatBase interface.

    This class manages a simple combination of operator, reactant, and product
    object without any additional metadata; essentially a dataclass but with the
    added benefit of being archivable in an ObjLib due to inheriting DataUnit.

    Parameters
    ----------
    operator : Identifier
        Operator object ID.
    products : Iterable[Identifier]
        Products of reaction IDs.
    reactants : Iterable[Identifier]
        Reactants involved in reaction IDs.

    Attributes
    ----------
    blob : bytes
        Binary representation of reaction.
    operator : Identifier
        Operator object ID.
    products : Iterable[Identifier]
        Products of reaction IDs.
    reactants : Iterable[Identifier]
        Reactants involved in reaction IDs.
    uid : Tuple[Identifier, Tuple[Identifier,...], Tuple[Identifier,...]]
        Unique identifier of object.
    """

    __slots__ = ("_blob", "_operator", "_products", "_reactants", "_uid")

    _operator: interfaces.Identifier
    _products: frozenset[interfaces.Identifier]
    _reactants: frozenset[interfaces.Identifier]

    _blob: typing.Optional[bytes]
    _uid: typing.Optional[
        tuple[
            interfaces.Identifier,
            tuple[interfaces.Identifier, ...],
            tuple[interfaces.Identifier, ...],
        ]
    ]

    def __init__(
        self,
        operator: typing.Optional[interfaces.Identifier] = None,
        reactants: typing.Optional[
            collections.abc.Iterable[interfaces.Identifier]
        ] = None,
        products: typing.Optional[
            collections.abc.Iterable[interfaces.Identifier]
        ] = None,
        reaction: typing.Optional[bytes] = None,
    ) -> None:
        if reaction is not None:
            data: tuple[
                interfaces.Identifier,
                tuple[interfaces.Identifier, ...],
                tuple[interfaces.Identifier, ...],
            ] = loads(reaction)
            self._operator = data[0]
            self._products = frozenset(data[1])
            self._reactants = frozenset(data[2])
        elif (
            operator is not None
            and reactants is not None
            and products is not None
        ):
            self._operator = operator
            self._reactants = frozenset(reactants)
            self._products = frozenset(products)
        else:
            raise TypeError("Insufficient arguments provided")
        self._blob = None
        self._uid = None

    @property
    def blob(self) -> bytes:
        if self._blob is None:
            self._blob = pickle.dumps(
                tuple((self._operator, self._products, self._reactants))
            )
        return self._blob

    @property
    def operator(self) -> interfaces.Identifier:
        return self._operator

    @property
    def products(self) -> frozenset[interfaces.Identifier]:
        return self._products

    @property
    def reactants(self) -> frozenset[interfaces.Identifier]:
        return self._reactants

    @property
    def uid(
        self,
    ) -> tuple[
        interfaces.Identifier,
        tuple[interfaces.Identifier, ...],
        tuple[interfaces.Identifier, ...],
    ]:
        if self._uid is None:
            self._uid = (
                self._operator,
                tuple(sorted(self._products)),
                tuple(sorted(self._reactants)),
            )
        return self._uid

    def __lt__(self, other: object) -> bool:
        if isinstance(other, RxnDatBasic):
            return self.uid < other.uid
        raise NotImplementedError("Comparison not implemented")

    def __repr__(self) -> str:
        return (
            f"RxnDatBasic(operator={repr(self.uid[0])}, "
            f"reactants={repr(self.uid[2])}, products={repr(self.uid[1])})"
        )
