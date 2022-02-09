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

from __future__ import annotations

import builtins
from abc import ABC, abstractmethod
from io import BytesIO
from pickle import Unpickler, UnpicklingError, dumps
from typing import (
    TYPE_CHECKING,
    Any,
    FrozenSet,
    Iterable,
    List,
    Optional,
    Protocol,
    Sequence,
    Tuple,
    Union,
    final,
)

from rdkit.Chem import Mol as BuildMol
from rdkit.Chem import MolFromSmiles, MolToSmiles
from rdkit.Chem.inchi import MolToInchiKey
from rdkit.Chem.rdchem import AtomValenceException, KekulizeException
from rdkit.Chem.rdchem import Mol as RDKitMol
from rdkit.Chem.rdChemReactions import ChemicalReaction as RDKitRxn
from rdkit.Chem.rdChemReactions import ReactionFromSmarts, ReactionToSmarts
from rdkit.Chem.rdmolops import AssignStereochemistry
from rdkit.Chem.rdmolops import SanitizeMol

if TYPE_CHECKING:
    from pickaxe_generic.engine import NetworkEngine

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


class __SafeUnpickler(Unpickler):
    def find_class(self, module: str, name: str) -> Any:
        if module == "builtins" and name in _safe_builtins_classes:
            return getattr(builtins, name)
        raise UnpicklingError(f"global '{module}.{name}' is forbidden")


def loads(s: bytes) -> Any:
    return __SafeUnpickler(BytesIO(s)).load()


class Identifier(Protocol):
    """
    Orderable, hashable object used as unique identifier.  Ideally immutable
    using public methods.

    Methods
    -------
    __hash__
    __eq__
    __lt__
    """

    def __hash__(self) -> int:
        """
        Hashes object to integer.

        Returns
        -------
        int
            Integer representing hashed value of object.  Should be
            (effectively) unique.
        """

    def __eq__(self, other: object) -> bool:
        """
        Compares object to others of similar type.  Enables hashtables.

        Arguments
        ---------
        other : object
            Object to be compared.

        Returns
        -------
        bool
            True if object is equivalent to other, False otherwise.
        """

    def __lt__(self, other: object) -> bool:
        """
        Compares object to others of similar type.  Allows sorting.

        Arguments
        ---------
        other : object
            Object to be compared.

        Returns
        -------
        bool
            True if object is after self when ordered, False otherwise.
        """


class DataUnit(ABC):
    """
    Object which provides a unique, hashable identifier, a method of ordering,
    and can serve up a binary form of the object.

    Attributes
    ----------
    blob : bytes
        Binary representation of object.
    uid : Identifier
        Unique identifier of object.

    Methods
    -------
    __eq__
    __lt__
    """

    __slots__ = ()

    @property
    @abstractmethod
    def blob(self) -> bytes:
        """
        Binary representation of object.  Must be able to initialize object when
        passed to __init__ method of any subclass of same type (viz. initialize
        a MolDatBasicV2, even if obtained from a MolDatBasicV1).
        """

    @property
    @abstractmethod
    def uid(self) -> Identifier:
        """
        Unique identifier of object.  Must be hashable in order to facilitate
        lookup tables utilizing hashes.
        """

    def __eq__(self, other: object) -> bool:
        """
        Compares object to others of similar type.  Enables hashtables.

        Arguments
        ---------
        other : object
            Object to be compared.

        Returns
        -------
        bool
            True if object is equivalent to other, False otherwise.
        """
        if isinstance(other, DataUnit):
            return self.uid == other.uid
        raise NotImplementedError("Cannot compare these objects")

    @abstractmethod
    def __lt__(self, other: object) -> bool:
        """
        Compares object to others of similar type.  Allows sorting.

        Arguments
        ---------
        other : object
            Object to be compared.

        Returns
        -------
        bool
            True if object is after self when ordered, False otherwise.
        """

    def __getstate__(self):
        return self.blob

    def __setstate__(self, data):
        self.__init__(data)


class MolDatBase(DataUnit):
    """
    Interface representing molecule data.

    Classes implementing this interface manage information about a single
    molecule, allowing for memory management and lumped molecule frameworks.

    Attributes
    ----------
    WORK IN PROGRESS
    """

    __slots__ = ()


class MolDatRDKit(MolDatBase):
    """
    Interface representing an RDKit molecule data object.

    Classes implementing this interface manage information about a single
    rdkit-compatible molecule.  Defines the constructor for this type of
    MolDat, and thus must be subclassed only by implementations with @final.
    Default behavior for __lt__ is sorting by UID.

    Parameters
    ----------
    molecule : Union[RDKitMol, str, bytes]
        Sufficient information to generate molecule in the form of an RDKitMol,
        a SMILES string, or the pickled form of an RDKitMol.
    sanitize : bool (default: True)
        Should be True when using input from non-sanitized sources.  Is not
        applied when initializing from blob.
    neutralize : bool (default: False)
        Should be True if you want hydrogens to be added/subtracted to
        neutralize a molecule and input is a SMILES string or non-neutralized
        molecule.  Is not applied when initializing from blob.

    Attributes
    ----------
    blob : bytes
        Binary representation of molecule.
    inchikey : str
        InChIKey of object.
    rdkitmol : RDKitMol
        RDKit molecule object.
    smiles : str
        SMILES string of object.
    uid : Identifier
        Unique identifier of object.
    """

    __slots__ = ()

    @abstractmethod
    def __init__(
        self,
        molecule: Union[RDKitMol, str, bytes],
        sanitize: bool = True,
        neutralize: bool = False,
    ) -> None:
        pass

    def __lt__(self, other: object) -> bool:
        if isinstance(other, MolDatRDKit):
            return self.uid < other.uid
        else:
            raise TypeError(
                f"Invalid comparison between objects of type {type(self)} and type {type(other)}"
            )

    @property
    @abstractmethod
    def blob(self) -> bytes:
        """
        RDKit-generated bytestring enabling quick regeneration of RDKitMol
        object.
        """

    @property
    @abstractmethod
    def inchikey(self) -> str:
        """
        InChIKey hash of molecule.
        """

    @property
    @abstractmethod
    def rdkitmol(self) -> RDKitMol:
        """
        Reference to RDKit molecule object containing basic properties.
        """

    @property
    @abstractmethod
    def smiles(self) -> str:
        """
        Canonical SMILES string of molecule object.
        """

    def _processinput(
        self,
        molecule: Union[RDKitMol, str, bytes],
        sanitize: bool = True,
        neutralize: bool = False,
    ) -> RDKitMol:
        if isinstance(molecule, bytes):
            rdkitmol = BuildMol(molecule)
            # if sanitize:
            #    SanitizeMol(rdkitmol)
            #    AssignStereochemistry(rdkitmol, True, True, True)
            # if neutralize:
            #    raise NotImplementedError("No neutralize function coded")
        elif isinstance(molecule, RDKitMol):
            rdkitmol = molecule
            # print(MolToSmiles(rdkitmol))
            if sanitize:
                SanitizeMol(rdkitmol)
                AssignStereochemistry(rdkitmol, True, True, True)
            if neutralize:
                raise NotImplementedError("No neutralize function coded")
        elif isinstance(molecule, str):
            if sanitize:
                rdkitmol = MolFromSmiles(molecule, sanitize=True)
            else:
                rdkitmol = MolFromSmiles(molecule, sanitize=False)
            if neutralize:
                raise NotImplementedError("No neutralize function coded")
        else:
            raise TypeError("Invalid molecule type")
        if rdkitmol is None:
            raise TypeError("Invalid molecule information")
        return rdkitmol


@final
class MolDatBasicV1(MolDatRDKit):
    """
    Version of MolDatRDKit which caches only SMILES and blob.

    Speeds: 5
    """

    __slots__ = ("_blob", "_smiles")
    _blob: bytes
    _smiles: str

    def __init__(
        self,
        molecule: Union[RDKitMol, str, bytes],
        sanitize: bool = True,
        neutralize: bool = False,
    ) -> None:
        rdkitmol = self._processinput(molecule, sanitize, neutralize)
        self._blob = rdkitmol.ToBinary()
        self._smiles = MolToSmiles(rdkitmol)

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
    def uid(self) -> Identifier:
        return self._smiles


@final
class MolDatBasicV2(MolDatRDKit):
    """
    Version of MolDatRDKit which caches all values.

    Speeds: 1,3,4,6
    """

    __slots__ = ("_blob", "_inchikey", "_rdkitmol", "_smiles")
    _blob: Optional[bytes]
    _inchikey: Optional[str]
    _rdkitmol: RDKitMol
    _smiles: str

    def __init__(
        self,
        molecule: Union[RDKitMol, str, bytes],
        sanitize: bool = True,
        neutralize: bool = False,
    ) -> None:
        self._blob = None
        self._inchikey = None
        self._rdkitmol = self._processinput(molecule, sanitize, neutralize)
        self._smiles = MolToSmiles(self._rdkitmol)

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


class OpDatBase(DataUnit):
    """
    Interface representing operator data.

    Classes implementing this interface manage information about a single
    operator which acts on MolDatBase and can generate RxnDatBase objects.

    Methods
    ----------
    compat
    __call__
    __len__
    """

    __slots__ = ()

    @abstractmethod
    def compat(self, mol: MolDatBase, arg: int) -> bool:
        """
        Determines compatibility of MolDat object with operator argument.

        Parameters
        ----------
        mol : MolDatBase
            MolDat object which is to be compared.
        arg : int
            Index of argument which is to be compared.

        Returns
        -------
        bool
            True if MolDat may be passed as argument arg to operator.
        """

    @abstractmethod
    def __call__(
        self, reactants: Sequence[MolDatBase]
    ) -> Iterable[Iterable[MolDatBase]]:
        """
        Reacts a sequence of MolDat objects using internal operator and returns
        a sequence of RxnDat objects which contain metadata about potential
        results.

        Parameters
        ----------
        reactants : Sequence[MolDatBase]
            Reactants which match the arguments in the operator.

        Returns
        -------
        Iterable[RxnDatBase]
            Iterable of reactions which are produced by applying the operator.
        """

    @abstractmethod
    def __len__(self) -> int:
        """
        Defined as number of arguments in operator.
        """


class OpDatRDKit(OpDatBase):
    """
    Interface representing an RDKit SMARTS operator.  Agents are treated as
    arguments following reagent arguments.

    Classes implementing this interface manage information about a single
    rdkit-compatible SMARTS operator.

    Attributes
    ----------
    smarts : str
        SMARTS string representing operator.
    rdkitrxn : RDKitRxn
        RDKit reaction object.
    """

    __slots__ = ()

    @abstractmethod
    def __init__(self, operator: Union[RDKitRxn, str, bytes]):
        pass

    @property
    @abstractmethod
    def smarts(self) -> str:
        """
        SMARTS string encoding operator information.
        """

    @property
    @abstractmethod
    def rdkitrxn(self) -> RDKitRxn:
        """
        Reference to RDKit reaction object.
        """


@final
class OpDatBasic(OpDatRDKit):
    """
    Minimal class implementing the OpDatRDKit interface.

    Classes implementing this interface manage information about a single
    rdkit-compatible SMARTS operator.

    Parameters
    ----------
    operator : Union[RDKitRxn, str, bytes]
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
        "_rdkitrxn",
        "_smarts",
        "_templates",
        "_uid",
    )

    _rdkitrxn: RDKitRxn
    _templates: Optional[Tuple[RDKitMol, ...]]
    _engine: NetworkEngine

    _blob: Optional[bytes]
    _smarts: Optional[str]
    _uid: Optional[Tuple[Tuple[str, ...], Tuple[str, ...]]]

    def __init__(
        self, operator: Union[RDKitRxn, str, bytes], engine: NetworkEngine
    ) -> None:
        if isinstance(operator, RDKitRxn):
            self._rdkitrxn = operator
        elif isinstance(operator, str):
            self._rdkitrxn = ReactionFromSmarts(operator)
            # SanitizeRxn(self._rdkitrxn)
        elif isinstance(operator, bytes):
            self._rdkitrxn = loads(operator)
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
            self._blob = dumps(self.rdkitrxn)
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

    def compat(self, mol: MolDatBase, arg: int) -> bool:
        if self._templates is None:
            self._templates = self._build_templates()
        if isinstance(mol, MolDatRDKit):
            return mol.rdkitmol.HasSubstructMatch(
                self._templates[arg], useChirality=True
            )
        else:
            return False

    def _build_templates(self) -> Tuple[RDKitMol, ...]:
        return tuple(self._rdkitrxn.GetReactants())

    def _attempt_reaction(self, mols: Iterable[RDKitMol]) -> Iterable[RDKitMol]:
        try:
            return self._rdkitrxn.RunReactants(mols, maxProducts=0)
        except Exception as e:
            print(type(e))
            raise e

    def __call__(
        self, reactants: Sequence[MolDatBase]
    ) -> Tuple[Tuple[MolDatBase, ...], ...]:
        rdkitmols: List[RDKitMol] = [
            reactant.rdkitmol
            for reactant in reactants
            if isinstance(reactant, MolDatRDKit)
        ]
        try:
            return tuple(
                tuple(self._engine.Mol(product) for product in products)
                for products in self._rdkitrxn.RunReactants(rdkitmols, maxProducts=0)
            )
        except AtomValenceException as e:
            raise ValueError(
                f"Error occurred when using operator {self} on {reactants}"
            ) from e
        except KekulizeException as e:
            raise ValueError(
                f"Error occurred when using operator {self} on {reactants}"
            ) from e

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
        return f"OpDatBasic({repr(self.uid)})"


class RxnDatBase(DataUnit):
    """
    Interface representing reaction data.

    Class implementing this interface manage information about a single reaction
    between several molecules to produce several molecules, with an associated
    operator.

    Attributes
    ----------
    operator : Identifier
        Operator object ID.
    products : Iterable[Identifier]
        Products of reaction IDs.
    reactants : Iterable[Identifier]
        Reactants involved in reaction IDs.
    """

    __slots__ = ()

    @abstractmethod
    def __init__(
        self,
        operator: Optional[Identifier] = None,
        reactants: Optional[Iterable[Identifier]] = None,
        products: Optional[Iterable[Identifier]] = None,
        reaction: Optional[bytes] = None,
    ) -> None:
        pass

    @property
    @abstractmethod
    def operator(self) -> Identifier:
        """
        Operator ID involved in reaction.
        """

    @property
    @abstractmethod
    def products(self) -> Iterable[Identifier]:
        """
        Products involved in reaction IDs.
        """

    @property
    @abstractmethod
    def reactants(self) -> Iterable[Identifier]:
        """
        Reactants involved in reaction IDs.
        """


@final
class RxnDatBasic(RxnDatBase):
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

    _operator: Identifier
    _products: FrozenSet[Identifier]
    _reactants: FrozenSet[Identifier]

    _blob: Optional[bytes]
    _uid: Optional[Tuple[Identifier, Tuple[Identifier, ...], Tuple[Identifier, ...]]]

    def __init__(
        self,
        operator: Optional[Identifier] = None,
        reactants: Optional[Iterable[Identifier]] = None,
        products: Optional[Iterable[Identifier]] = None,
        reaction: Optional[bytes] = None,
    ) -> None:
        if reaction is not None:
            data: Tuple[
                Identifier, Tuple[Identifier, ...], Tuple[Identifier, ...]
            ] = loads(reaction)
            self._operator = data[0]
            self._products = frozenset(data[1])
            self._reactants = frozenset(data[2])
        elif operator is not None and reactants is not None and products is not None:
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
            self._blob = dumps(tuple((self._operator, self._products, self._reactants)))
        return self._blob

    @property
    def operator(self) -> Identifier:
        return self._operator

    @property
    def products(self) -> FrozenSet[Identifier]:
        return self._products

    @property
    def reactants(self) -> FrozenSet[Identifier]:
        return self._reactants

    @property
    def uid(
        self,
    ) -> Tuple[Identifier, Tuple[Identifier, ...], Tuple[Identifier, ...]]:
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
        return f"RxnDatBasic(operator={repr(self.uid[0])}, reactants={repr(self.uid[2])}, products={repr(self.uid[1])})"
