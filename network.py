"""
Contains classes which define and implement molecule-operator networks for use
in reaction network generation.

Classes:

    Identifier(Protocol)
    DataUnit
      MolDatBase(DataUnit)
        MolDatRDKit(MolDatBase)
          MolDatBasic(MolDatRDKit)
      OpDatBase(DataUnit)
        OpDatRDKit(OpDatBase)
          OpDatBasic(OpDatRDKit)
      RxnDatBase(DataUnit)
        RxnDatBasic(RxnDatBase)
    ObjectLibrary
      ObjectLibraryBasic(ObjectLibrary)
    ExpansionStrategy
      CartesianStrategy(ExpansionStrategy)

"""

from abc import ABC, abstractmethod
import builtins
from io import BytesIO
from itertools import chain, product as iterproduct
from pickle import dumps, Unpickler, UnpicklingError
from typing import (
    Any,
    Callable,
    Dict,
    FrozenSet,
    Generator,
    Generic,
    Iterable,
    Iterator,
    List,
    Optional,
    Protocol,
    Sequence,
    Set,
    Tuple,
    TypeVar,
    Union,
    final,
)
from sys import exc_info

import rdkit.Chem.rdChemReactions
from rdkit.Chem import Mol as BuildMol, MolFromSmiles, MolToSmiles
from rdkit.Chem.inchi import MolToInchiKey
from rdkit.Chem.rdchem import (
    AtomValenceException,
    KekulizeException,
    Mol as RDKitMol,
)
from rdkit.Chem.rdChemReactions import ChemicalReaction as RDKitRxn
from rdkit.Chem.rdChemReactions import (
    ReactionFromSmarts,
    ReactionToSmarts,
    SanitizeRxn,
)
from rdkit.Chem.rdmolops import AssignStereochemistry, SanitizeMol


# some code to make loads more safe to arbitrary code execution
# necessary since external data from a database may be input
# if you are having issues, add relevant classes to _safe_%module%_classes
# or if modules not in builtins required, add other if clauses
_safe_builtins_classes = frozenset(
    {
        "frozenset",
        "tuple",
    }
)
_safe_rdkit_chem_rdchemreactions_classes = frozenset(
    {
        "ChemicalReaction",
    }
)


class __SafeUnpickler(Unpickler):
    def find_class(self, module, name):
        if module == "builtins" and name in _safe_builtins_classes:
            return getattr(builtins, name)
        elif (
            module == "rdkit.Chem.rdChemReactions"
            and name in _safe_rdkit_chem_rdchemreactions_classes
        ):
            return getattr(rdkit.Chem.rdChemReactions, name)
        raise UnpicklingError(f"global '{module}.{name}' is forbidden")


def loads(s):
    return __SafeUnpickler(BytesIO(s)).load()


def create_engine(speed: int = 5):
    return NetworkEngineBasic(speed=speed)


def copy_set_values(src: Set, dst: Set) -> Set:
    """
    Copies references from src set to dst set if available, reduces memory
    redundancy.

    Parameters
    ----------
    src : Set
        Set containing source references.
    dst : Set
        Set containing references to be overwritten.

    Returns
    -------
    Set
        Updated version of dst with references from src.
    """
    return (dst & src) | dst


class Identifier(Protocol):
    """
    Orderable, hashable object used as unique identifier.

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
            Integer representing hashed value of object.  Should be nearly
            unique.
        """
        pass

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
        pass

    def __lt__(self, other) -> bool:
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
        pass


class DataUnit(ABC):
    """
    Object which provides a unique, hashable identifier, a method of ordering,
    and can serve up a NamedTuple containing database form of the object.

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
        passed to __init__ method of compatible subclass.
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


DataUnitGen = TypeVar("DataUnitGen", bound=DataUnit)


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

    def __lt__(self, other: object):
        if isinstance(other, MolDatRDKit):
            return self.uid < other.uid
        raise NotImplementedError("Comparison not implemented for this object")

    def __repr__(self):
        return f"MolDatRDKit({self.smiles})"


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
    _engine: "NetworkEngine"

    _blob: Optional[bytes]
    _smarts: Optional[str]
    _uid: Optional[Tuple[Tuple[str, ...], Tuple[str, ...]]]

    def __init__(
        self, operator: Union[RDKitRxn, str, bytes], engine: "NetworkEngine"
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
                for products in self._rdkitrxn.RunReactants(
                    rdkitmols, maxProducts=0
                )
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
    _uid: Optional[
        Tuple[Identifier, Tuple[Identifier, ...], Tuple[Identifier, ...]]
    ]

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
            self._blob = dumps(
                tuple((self._operator, self._products, self._reactants))
            )
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


class ObjectLibrary(ABC, Generic[DataUnitGen]):
    """
    Interface representing library of data.

    Classes implementing this interface manage multiple instances of a hashable
    object, and may have responsibility for synchronization with external
    databases which may also manage this information (be volatile).  Contained
    objects must have a "uid" attribute which contains a hashable unique id.

    Current implementation assumes that this library will never shrink or remove
    entries.
    """

    __slots__ = ()

    @abstractmethod
    def add(self, obj: Union[Iterable[DataUnitGen], DataUnitGen]) -> None:
        """
        Adds an object or multiple objects to the library (as long as it is
        unique to the other objects in the library).

        Parameters
        ----------
        obj : Union[Iterable[DataUnit], DataUnit]
            Object(s) to be added.
        """

    @abstractmethod
    def ids(self) -> Iterable[Identifier]:
        """
        Returns a set of keys used in the library.

        Returns
        -------
        Iterable[Identifier]
            Ids of objects in the library.
        """
        pass

    @abstractmethod
    def __contains__(self, item: DataUnitGen) -> bool:
        pass

    @abstractmethod
    def __getitem__(self, id: Identifier) -> DataUnitGen:
        pass

    @abstractmethod
    def __iter__(self) -> Iterator[DataUnitGen]:
        pass

    @abstractmethod
    def __len__(self) -> int:
        pass


@final
class ObjectLibraryBasic(ObjectLibrary, Generic[DataUnitGen]):
    """
    Minimal class implementing the ObjectLibrary interface.

    This class wraps a dict.

    Parameters
    ----------
    objects : Optional[Iterable[DataUnit]]
        Objects which are to be included in the object library.
    """

    __slots__ = ("_lookup",)

    _lookup: Dict[Identifier, DataUnitGen]

    def __init__(
        self,
        objects: Optional[Union[Iterable[DataUnitGen], DataUnitGen]] = None,
    ) -> None:
        self._lookup = {}
        if isinstance(objects, Iterable):
            for item in objects:
                self._lookup[item.uid] = item
        elif objects is not None:
            self._lookup[objects.uid] = objects

    def add(self, obj: Union[Iterable[DataUnitGen], DataUnitGen]) -> None:
        if isinstance(obj, Iterable):
            for item in obj:
                self._lookup[item.uid] = item
        else:
            self._lookup[obj.uid] = obj

    def ids(self) -> Generator[Identifier, None, None]:
        return (key for key in self._lookup.keys())

    def __contains__(self, item: DataUnitGen) -> bool:
        return item.uid in self._lookup

    def __getitem__(self, id: Identifier) -> DataUnitGen:
        return self._lookup[id]

    def __iter__(self) -> Iterator[DataUnitGen]:
        return iter(self._lookup.values())

    def __len__(self) -> int:
        return len(self._lookup)


@final
class ObjectLibraryKeyVal(ObjectLibrary, Generic[DataUnitGen]):
    """
    Minimal class implementing the ObjectLibrary interface; stores binary
    representation and initializes using passed function.

    This class wraps a dict.

    Parameters
    ----------
    objects : Optional[Iterable[DataUnitGen]]
        Objects which are to be included in the object library.
    initializer : Callable[[bytes],DataUnitGen]
        Initializer which can convert bytes to relevant data type.
    """

    __slots__ = (
        "_initializer",
        "_lookup",
    )

    _lookup: Dict[Identifier, bytes]
    _initializer: Callable[[bytes], DataUnitGen]

    def __init__(
        self,
        objects: Optional[Union[Iterable[DataUnitGen], DataUnitGen]] = None,
        initializer: Callable[[bytes], DataUnitGen] = None,
    ) -> None:
        if initializer is None:
            raise ValueError("initializer must be specified")
        self._lookup = {}
        self._initializer = initializer
        if isinstance(objects, Iterable):
            for item in objects:
                self._lookup[item.uid] = item.blob
        elif objects is not None:
            self._lookup[objects.uid] = objects.blob

    def add(self, obj: Union[Iterable[DataUnitGen], DataUnitGen]) -> None:
        if isinstance(obj, Iterable):
            for item in obj:
                self._lookup[item.uid] = item.blob
        else:
            self._lookup[obj.uid] = obj.blob

    def ids(self) -> Generator[Identifier, None, None]:
        return (key for key in self._lookup.keys())

    def __contains__(self, item: DataUnitGen) -> bool:
        return item.uid in self._lookup

    def __getitem__(self, id: Identifier) -> DataUnitGen:
        return self._initializer(self._lookup[id])

    def __iter__(self) -> Iterator[DataUnitGen]:
        return (self._initializer(value) for value in self._lookup.values())

    def __len__(self) -> int:
        return len(self._lookup)


class ExpansionStrategy(ABC):
    """
    Interface representing a network expansion strategy.

    Classes implementing this interface use information from a molecule and
    operator library to generate new reactions, which are then output to a
    reaction library.
    """

    __slots__ = ()

    @abstractmethod
    def __init__(
        self,
        mol_lib: ObjectLibrary[MolDatBase],
        op_lib: ObjectLibrary[OpDatBase],
        rxn_lib: ObjectLibrary[RxnDatBase],
    ) -> None:
        pass

    @abstractmethod
    def expand(
        self,
        max_rxns: Optional[int] = None,
        max_mols: Optional[int] = None,
        num_gens: Optional[int] = None,
    ) -> None:
        """
        Expand molecule library.

        Parameters
        ----------
        max_rxns : Optional[int] (default: None)
            Limit of new reactions to add.  If None, no limit.
        max_mols : Optional[int] (default: None)
            Limit of new molecules to add.  If None, no limit.
        num_gens : Optional[int] (defauilt: None)
            Maximum generations of reactions to enumerate.  If None, no limit.
        """

    @abstractmethod
    def refresh(self) -> None:
        """
        Refreshes cache of strategy.  Use when mol_lib, op_lib, and rxn_lib are
        volatile and an update is known to have occurred.
        """


@final
class CartesianStrategy(ExpansionStrategy):
    """
    Implements ExpansionStrategy interface via Cartesian product of molecules
    and operators.  The outermost loop is through all operators, and each
    subsequent inner loop iterates through the molecule library for each
    argument.  The strategy then refreshes itself with new molecules and
    proceeds again.

    Parameters
    ----------
    mol_lib : ObjectLibrary
        Library containing MolDatBase objects.
    op_lib : ObjectLibrary
        Library containing OpDatBase objects.
    rxn_lib : ObjectLibrary
        Library containing RxnDatBase objects.
    """

    __slots__ = (
        "_compat_table",
        "_engine",
        "_mol_cache",
        "_op_cache",
        "_recipe_cache",
        "_mol_lib",
        "_op_lib",
        "_rxn_lib",
    )

    # list of compatible molecule uids stored for each argument/operator
    # combination as dict[op][argnum]
    _compat_table: Dict[Identifier, List[List[Identifier]]]

    # dict of molecules whose compatibility has been tested
    _mol_cache: Dict[Identifier, MolDatBase]

    # dict of operators whose compatibility has been tested
    _op_cache: Dict[Identifier, OpDatBase]

    # set of reactions which have already been tried
    _recipe_cache: Set[Tuple[Identifier, FrozenSet[Identifier]]]

    _mol_lib: ObjectLibrary[MolDatBase]
    _op_lib: ObjectLibrary[OpDatBase]
    _rxn_lib: ObjectLibrary[RxnDatBase]

    """
    # set of molecules which can be used as the first argument
    # via stride/offset
    _mol_init: Set[MolDatBase]

    # stride between first argument molecules (used for parallelism)
    _stride: int

    # offset for first argument molecule (used for parallelism)
    _offset: int
    """

    def __init__(
        self,
        mol_lib: ObjectLibrary[MolDatBase],
        op_lib: ObjectLibrary[OpDatBase],
        rxn_lib: ObjectLibrary[RxnDatBase],
        engine: "NetworkEngine",
    ) -> None:
        self._engine = engine
        self._mol_lib = mol_lib
        self._mol_cache = {}
        self._op_lib = op_lib
        self._op_cache = {}
        self._rxn_lib = rxn_lib
        self._recipe_cache = set()

        # initialize compat_table
        self._compat_table = {
            op.uid: [[] for i in range(len(op))] for op in self._op_lib
        }

        # initialize _op_cache
        for op in self._op_lib:
            self._op_cache[op.uid] = op

        # fill _compat_table
        for mol in self._mol_lib:
            self._mol_cache[mol.uid] = mol
            self._add_mol_to_compat(mol)

    def expand(
        self,
        max_rxns: Optional[int] = None,
        max_mols: Optional[int] = None,
        num_gens: Optional[int] = None,
        filter: Callable[[MolDatBase],bool] = lambda _: True,
    ) -> None:
        exhausted: bool = False
        num_mols: int = 0
        num_rxns: int = 0
        gen: int = 0
        while not exhausted:
            if num_gens is not None and gen >= num_gens:
                return
            exhausted = True
            for op_uid in self._op_cache:
                react_uids: Tuple[Identifier, ...]
                for react_uids in iterproduct(*self._compat_table[op_uid]):
                    recipe = (op_uid, frozenset(react_uids))
                    if recipe in self._recipe_cache:
                        continue
                    op = self._op_cache[op_uid]
                    reactants = tuple(
                        self._mol_cache[uid] for uid in react_uids
                    )
                    for productset in op(reactants):
                        prod_uids = frozenset(mol.uid for mol in productset)
                        # print(prod_uids)
                        rxn = self._engine.Rxn(op_uid, react_uids, prod_uids)
                        if rxn in self._rxn_lib:
                            continue
                        temp_mols: List[MolDatBase] = []
                        for product in productset:
                            if (
                                product.uid not in self._mol_cache
                                and product not in self._mol_lib
                                and filter(product)
                            ):
                                temp_mols.append(product)
                                num_mols += 1
                                if max_mols is not None and num_mols > max_mols:
                                    return
                        for mol in temp_mols:
                            self._mol_lib.add(mol)
                        """if self._engine.speed == 4:
                            rxn = self._engine.Rxn(
                                op_uid,
                                react_uids,
                                copy_set_values(
                                    frozenset(self._mol_lib.ids()), rxn.products
                                ),
                            )"""
                        self._rxn_lib.add(rxn)
                        exhausted = False
                        num_rxns += 1
                        if max_rxns is not None and num_rxns >= max_rxns:
                            return
                    self._recipe_cache.add(recipe)
            gen += 1
            self.refresh()

    def refresh(self) -> None:
        # check for molecules in mol_lib which are not in _mol_cache and add them
        if len(self._mol_lib) > len(self._mol_cache):
            for mol_uid in self._mol_lib.ids():
                if mol_uid not in self._mol_cache:
                    mol = self._mol_lib[mol_uid]
                    self._mol_cache[mol_uid] = mol
                    self._add_mol_to_compat(mol)

        # check for molecules in op_lib which are not in _op_cache and add them
        if len(self._op_lib) > len(self._compat_table):
            for op_uid in self._op_lib.ids():
                if op_uid not in self._op_cache:
                    op = self._op_lib[op_uid]
                    self._op_cache[op_uid] = op
                    self._add_op_to_compat(op)

    def _add_mol_to_compat(self, mol: MolDatBase) -> None:
        """
        Add entries to compat_table for new molecule.

        Parameters
        ----------
        mol : MolDatBase
            Molecule object to be added to compat_table.
        """
        for op in self._op_cache.values():
            for arg in range(len(self._compat_table[op.uid])):
                if op.compat(mol, arg):
                    self._compat_table[op.uid][arg].append(mol.uid)
        # if (index % self._stride == 0):
        #    self._mol_init.add(mol)

    def _add_op_to_compat(self, op: OpDatBase) -> None:
        """
        Add entries to compat_table for new operator.

        Parameters
        ----------
        op : OpDatBase
            Operator object to be added to compat_table.
        """
        optable: List[List[Identifier]] = [[] for _ in range(len(op))]
        for arg in range(len(optable)):
            for mol in self._mol_cache.values():
                if op.compat(mol, arg):
                    optable[arg].append(mol.uid)
        self._compat_table[op.uid] = optable


class HybridExpansionStrategy(ExpansionStrategy):
    """
    Interface representing a hybrid network expansion strategy.

    Classes implementing this interface use information from a molecule library
    and multiple operator libraries to generate new reactions, which are then
    output to a reaction library.
    """

    __slots__ = ()

    @abstractmethod
    def __init__(
        self,
        mol_lib: ObjectLibrary[MolDatBase],
        op_libs: Sequence[ObjectLibrary[OpDatBase]],
        rxn_lib: ObjectLibrary[RxnDatBase],
    ) -> None:
        pass


@final
class OrderedCartesianHybridExpansionStrategy(HybridExpansionStrategy):
    """
    Implements HybridExpansionStrategy interface via Cartesian product of
    molecules and operators.  The outermost loop is through all operators, and
    each subsequent inner loop iterates through the molecule library for each
    argument.  The strategy then refreshes itself with new molecules and
    proceeds again.  The ordered condition requires that molecules produced by
    operators from a library further down the sequence cannot be operated on by
    operators from a library toward the beginning of the sequence.

    Parameters
    ----------
    mol_lib : ObjectLibrary[MolDatBase]
        Library containing MolDatBase objects.
    op_libs : Sequence[ObjectLibrary[OpDatBase]]
        Container of ObjectLibraries containing OpDatBase objects.
    rxn_lib : ObjectLibrary[RxnDatBase]
        Library containing RxnDatBase objects.
    """

    __slots__ = (
        "_compat_table",
        "_mol_cache",
        "_op_cache",
        "_recipe_cache",
        "_mol_lib",
        "_op_libs",
        "_rxn_lib",
        "_engine",
    )

    # list of compatible molecule uids stored for each argument/operator
    # combination as list[op_lib_index][op][argnum]
    _compat_table: List[Dict[Identifier, List[List[Identifier]]]]

    # dict of molecules whose compatibility has been tested
    _mol_cache: Dict[Identifier, Tuple[MolDatBase, int]]

    # dict of operators whose compatibility has been tested
    _op_cache: List[Dict[Identifier, OpDatBase]]

    # set of reactions which have already been tried
    _recipe_cache: Set[Tuple[Identifier, FrozenSet[Identifier]]]

    _mol_lib: ObjectLibrary[MolDatBase]
    _op_libs: Sequence[ObjectLibrary[OpDatBase]]
    _rxn_lib: ObjectLibrary[RxnDatBase]

    def __init__(
        self,
        mol_lib: ObjectLibrary[MolDatBase],
        op_libs: Sequence[ObjectLibrary[OpDatBase]],
        rxn_lib: ObjectLibrary[RxnDatBase],
        engine: "NetworkEngine",
    ) -> None:
        self._engine = engine
        self._mol_lib = mol_lib
        self._mol_cache = {}
        self._op_libs = op_libs
        self._op_cache = []
        self._rxn_lib = rxn_lib
        self._recipe_cache = set()

        # initialize compat_table
        self._compat_table = [
            {op.uid: [[] for i in range(len(op))] for op in op_lib}
            for op_lib in self._op_libs
        ]

        # initialize _op_cache
        for i in range(len(self._op_libs)):
            op_lib_cache = {}
            for op in self._op_libs[i]:
                op_lib_cache[op.uid] = op
            self._op_cache.append(op_lib_cache)

        # fill _compat_table
        for mol in self._mol_lib:
            self._mol_cache[mol.uid] = (mol, -1)
            self._add_mol_to_compat(mol, -1)

    def _add_mol_to_compat(self, mol: MolDatBase, step: int) -> None:
        """
        Add entries to compat_table for new molecule.

        Parameters
        ----------
        mol : MolDatBase
            Molecule object to be added to compat_table.
        step : int
            Index of earliest operator library which generated the molecule.
        """
        for i in range(len(self._op_cache)):
            for op in self._op_cache[i]:
                for arg in range(len(self._compat_table[i][op])):
                    if self._op_libs[i][op].compat(mol, arg):
                        self._compat_table[i][op][arg].append(mol.uid)

    def _add_op_to_compat(self, op: OpDatBase, step: int) -> None:
        """
        Add entries to compat_table for new operator.

        Parameters
        ----------
        op : OpDatBase
            Operator object to be added to compat_table.
        step : int
            Index of op_lib operator originates from.
        """
        optable: List[List[Identifier]] = [[] for _ in range(len(op))]
        for arg in range(len(optable)):
            for mol, _ in self._mol_cache.values():
                if op.compat(mol, arg):
                    optable[arg].append(mol.uid)
        self._compat_table[step][op.uid] = optable

    def _compat_table_generator(
        self,
    ) -> Generator[
        Tuple[Tuple[OpDatBase, FrozenSet[RDKitMol]], int], None, None
    ]:
        for lib_index in range(len(self._op_cache)):
            for op_uid in self._op_cache[lib_index]:
                for react_uids in (
                    reactantset
                    for reactantset in iterproduct(
                        *(self._compat_table[lib_index][op_uid])
                    )
                    if all(
                        self._mol_cache[r][1] <= lib_index for r in reactantset
                    )
                ):
                    yield (op_uid, frozenset(react_uids)), lib_index

    def expand(
        self,
        max_rxns: Optional[int] = None,
        max_mols: Optional[int] = None,
        num_gens: Optional[int] = None,
        filter: Callable[[MolDatBase],bool] = lambda _: True,
    ) -> None:
        exhausted: bool = False
        num_mols: int = 0
        num_rxns: int = 0
        gen: int = 0

        while not exhausted:
            if num_gens is not None and gen >= num_gens:
                return
            exhausted = True
            for recipe, lib_index in self._compat_table_generator():
                if recipe in self._recipe_cache:
                    continue
                op = self._op_cache[lib_index][recipe[0]]
                reactants = tuple(self._mol_cache[uid][0] for uid in recipe[1])
                for productset in op(reactants):
                    prod_uids = frozenset(mol.uid for mol in productset)
                    rxn = self._engine.Rxn(recipe[0], recipe[1], prod_uids)
                    if rxn in self._rxn_lib:
                        continue
                    temp_mols: List[MolDatBase] = []
                    for product in productset:
                        if product.uid in self._mol_cache:
                            cur_index = self._mol_cache[product.uid][1]
                            self._mol_cache[product.uid] = (
                                self._mol_cache[product.uid][0],
                                min(lib_index, cur_index),
                            )
                        elif (
                            product.uid not in self._mol_cache
                            and product not in self._mol_lib
                            and filter(product)
                        ):
                            temp_mols.append(product)
                            num_mols += 1
                            if max_mols is not None and num_mols > max_mols:
                                return
                    for mol in temp_mols:
                        self._mol_lib.add(mol)
                        self._mol_cache[mol.uid] = (mol,lib_index)
                        self._add_mol_to_compat(mol,lib_index)
                    self._rxn_lib.add(rxn)
                    exhausted = False
                    num_rxns += 1
                    if max_rxns is not None and num_rxns >= max_rxns:
                        return
                self._recipe_cache.add(recipe)
        gen += 1
        self.refresh()

    def refresh(self) -> None:
        if len(self._mol_lib) > len(self._mol_cache):
            for mol_uid in self._mol_lib.ids():
                if mol_uid not in self._mol_cache:
                    mol = self._mol_lib[mol_uid]
                    self._mol_cache[mol_uid] = mol, -1
                    self._add_mol_to_compat(mol, -1)
        n = 0
        for op_lib, op_cache, compat_table in zip(
            self._op_libs, self._op_cache, self._compat_table
        ):
            if len(op_lib) > len(compat_table):
                for op_uid in op_lib.ids():
                    if op_uid not in op_cache:
                        op = op_lib[op_uid]
                        op_cache[op_uid] = op
                        self._add_op_to_compat(op, n)
            n += 1


class NetworkEngine(ABC):
    """
    Interface representing an object which serves up other objects based on
    configuration parameters.

    Classes implementing this interface determine which type of network objects
    are constructed based on configuration options.
    """

    @property
    @abstractmethod
    def speed(self) -> int:
        """
        Defined speed of engine configuration.
        """

    @abstractmethod
    def Mol(
        self,
        molecule: Union[RDKitMol, str, bytes],
        sanitize: bool = True,
        neutralize: bool = False,
    ) -> MolDatRDKit:
        """
        Initializes a MolDatRDKit object of relevant type.
        """

    @abstractmethod
    def Op(self, operator: Union[RDKitRxn, str, bytes]) -> OpDatRDKit:
        """
        Initializes an OpDatRDKit object of relevant type.
        """

    @abstractmethod
    def Rxn(
        self,
        operator: Optional[Identifier] = None,
        reactants: Optional[Iterable[Identifier]] = None,
        products: Optional[Iterable[Identifier]] = None,
        reaction: Optional[bytes] = None,
    ) -> RxnDatBase:
        """
        Initializes a RxnDatBase object of relevant type.
        """

    @abstractmethod
    def Libs(
        self,
    ) -> Tuple[
        ObjectLibrary[MolDatBase],
        ObjectLibrary[OpDatBase],
        ObjectLibrary[RxnDatBase],
    ]:
        """
        Initializes the three basic ObjectLibraries necessary to run a Strategy.
        """

    @abstractmethod
    def CartesianStrategy(
        self,
        mol_lib: ObjectLibrary[MolDatBase],
        op_lib: ObjectLibrary[OpDatBase],
        rxn_lib: ObjectLibrary[RxnDatBase],
    ) -> CartesianStrategy:
        """
        Initializes a CartesianStrategy of relevant type.
        """


class NetworkEngineBasic(NetworkEngine):
    """
    Implements NetworkEngine class for different speed efficiencies.  Default
    for module.

    Parameters
    ----------
    speed : int (default: 5)
        Integer between 1 and 6 determining speed/memory tradeoff of bulk data.
            6: Maximum caching in RAM, no disk use.
            5: Most data in RAM, no disk use.
            4: Bare minimum data in RAM, no disk use.
            3: Fast primary keys in RAM, disk caches values.
            2: Smallest possible primary keys in RAM, disk caches values.
            1: Fast primary keys and values both stored on disk.
    """

    __slots__ = (
        "_Mol",
        "_Op",
        "_Rxn",
        "_Mol_Lib",
        "_Op_Lib",
        "_Rxn_Lib",
        "_CartesianStrategy",
        "_speed",
    )

    def __init__(self, speed: int = 5):
        if speed == 1:
            # type: ignore
            self._Mol = MolDatBasicV2
            raise NotImplementedError("Speed not yet implemented")
        elif speed == 2:
            self._Mol = MolDatBasicV2
            raise NotImplementedError("Speed not yet implemented")
        elif speed == 3:
            self._Mol = MolDatBasicV2
            raise NotImplementedError("Speed not yet implemented")
        elif speed == 4:
            self._Mol = MolDatBasicV2
            self._Mol_Lib = lambda: ObjectLibraryKeyVal(initializer=self.Mol)
            self._Op_Lib = lambda: ObjectLibraryKeyVal(initializer=self.Op)
            # self._Rxn_Lib = lambda: ObjectLibraryKeyVal(initializer=self.Rxn)
            self._Rxn_Lib = ObjectLibraryBasic
        elif speed == 5:
            self._Mol = MolDatBasicV1
            self._Mol_Lib = ObjectLibraryBasic
            self._Op_Lib = ObjectLibraryBasic
            self._Rxn_Lib = ObjectLibraryBasic
        elif speed == 6:
            self._Mol = MolDatBasicV2
            self._Mol_Lib = ObjectLibraryBasic
            self._Op_Lib = ObjectLibraryBasic
            self._Rxn_Lib = ObjectLibraryBasic
        else:
            raise ValueError(f"speed = {speed} is invalid")
        self._Op = OpDatBasic
        self._Rxn = RxnDatBasic
        self._CartesianStrategy = CartesianStrategy
        self._speed = speed

    @property
    def speed(self) -> int:
        return self._speed

    def Mol(
        self,
        molecule: Union[RDKitMol, str, bytes],
        sanitize: bool = True,
        neutralize: bool = False,
    ) -> MolDatRDKit:
        return self._Mol(
            molecule=molecule, sanitize=sanitize, neutralize=neutralize
        )

    def Op(self, operator: Union[RDKitRxn, str, bytes]) -> OpDatBasic:
        return self._Op(operator=operator, engine=self)

    def Rxn(
        self,
        operator: Optional[Identifier] = None,
        reactants: Optional[Iterable[Identifier]] = None,
        products: Optional[Iterable[Identifier]] = None,
        reaction: Optional[bytes] = None,
    ) -> RxnDatBasic:
        return self._Rxn(
            operator=operator,
            reactants=reactants,
            products=products,
            reaction=reaction,
        )

    def Libs(
        self,
    ) -> Tuple[
        ObjectLibrary[MolDatBase],
        ObjectLibrary[OpDatBase],
        ObjectLibrary[RxnDatBase],
    ]:
        mol_lib: ObjectLibrary[MolDatBase] = self._Mol_Lib()
        op_lib: ObjectLibrary[OpDatBase] = self._Op_Lib()
        rxn_lib: ObjectLibrary[RxnDatBase] = self._Rxn_Lib()
        return (mol_lib, op_lib, rxn_lib)

    def CartesianStrategy(
        self,
        mol_lib: ObjectLibrary[MolDatBase],
        op_lib: ObjectLibrary[OpDatBase],
        rxn_lib: ObjectLibrary[RxnDatBase],
    ) -> CartesianStrategy:
        return self._CartesianStrategy(
            mol_lib=mol_lib, op_lib=op_lib, rxn_lib=rxn_lib, engine=self
        )


class RxnTracker(ABC):
    """
    Interface representing an object which analyzes rxn network connections.

    Classes implementing this interface are able to create retrosynthetic trees
    based on a precalculated reaction network tree.

    Parameters
    ----------
    target : Identifier
        Unique ID of target molecule.
    reagent_table : Sequence[Identifier] (default: tuple())
        Contains unique IDs of reagents which do not need to be synthesized.
    fail_on_unknown_reagent : bool (default: False)
        If True, do not return paths which require reagents not in
        reagent_table.
    """

    def getParentChains(
        self,
        target: Identifier,
        reagent_table: Sequence[Identifier] = tuple(),
        fail_on_unknown_reagent: bool = False,
    ) -> Iterable[Iterable[Iterable[RxnDatBase]]]:
        """
        Gets parent chains for a particular target molecule.

        Parameters
        ----------
        target : Identifier
            Unique id of target molecule.
        reagent_table : Sequence[Identifier]
            Sequence of reagents which are considered "basic" and which the tree
            search will consider leaf nodes.
        fail_on_unknown_reagent : bool
            If tree requires unlisted reagents, do not return.
        """


class RxnTrackerSingle(RxnTracker):
    """Implements RxnTracker interface; only compatible with reactions
    involving a single reactant and product.  DEVELOPMENT ONLY"""

    _mol_lookup: Dict[Identifier, Identifier]

    def __init__(self, rxn_lib: ObjectLibrary[RxnDatBase]) -> None:
        self._mol_lookup = {}
        self._rxn_lib = rxn_lib
        for rxnid in rxn_lib.ids():
            product_mol = sorted(rxn_lib[rxnid].products)[0]
            if product_mol not in self._mol_lookup:
                self._mol_lookup[product_mol] = []
            self._mol_lookup[product_mol].append(rxnid)

    def _getchains(
        self, cur_mol: Identifier, cur_mols: Optional[Set[Identifier]] = None
    ):
        if cur_mols is None:
            cur_mols = {cur_mol}
        noReactions = True
        if cur_mol not in self._mol_lookup:
            yield list()
        else:
            for rxnid in self._mol_lookup[cur_mol]:
                reactant = sorted(self._rxn_lib[rxnid].reactants)[0]
                if reactant in cur_mols:
                    continue
                noReactions = False
                for rxnpath in self._getchains(
                    reactant, cur_mols.union({reactant})
                ):
                    rxnpath.append(rxnid)
                    yield rxnpath
            if noReactions:
                yield list()

    def getParentChains(
        self,
        target: Identifier,
        #        reagent_table: Sequence[Identifier] = tuple(),
        #        fail_on_unknown_reagent: bool = False,
    ) -> Iterable[Iterable[RxnDatBase]]:
        return (path for path in self._getchains(target))


class RxnTrackerDepthFirst(RxnTracker):
    """Implements RxnTracker interface; stores lookups as a hash table within
    the object.  Will eventually deprecate this functionality when the
    ObjectLibrary interface is updated to include native search functionality.
    """

    _mol_lookup: Dict[Identifier, Identifier]
    _rxn_lib: ObjectLibrary[RxnDatBase]

    def __init__(self, rxn_lib: ObjectLibrary[RxnDatBase]) -> None:
        self._mol_lookup = {}
        self._rxn_lib = rxn_lib
        for rxnid in rxn_lib.ids():
            for product_mol in rxn_lib[rxnid].products:
                if product_mol not in self._mol_lookup:
                    self._mol_lookup[product_mol] = []
                self._mol_lookup[product_mol].append(rxnid)

    def _getchains(
        self,
        cur_gen_mols: Iterable[Identifier],
        prev_gens_mols: Optional[Iterable[Identifier]] = None,
        prev_gens_rxns: Optional[Iterable[Identifier]] = None,
        reagent_table: Optional[Iterable[Identifier]] = None,
        fail_on_unknown_reagent: bool = False,
    ) -> Generator[List[Set[Identifier]], None, None]:
        if len(cur_gen_mols) == 0:
            yield []
            return
        if prev_gens_mols is None:
            prev_gens_mols = set()
        if prev_gens_rxns is None:
            prev_gens_rxns = set()
        if reagent_table is None:
            reagent_table = set()
        rxnsets = []
        for mol in cur_gen_mols:
            if mol in reagent_table:
                continue
            elif mol not in self._mol_lookup:
                rxnsets.append([])
                continue
            newrxnset = [
                rxn
                for rxn in self._mol_lookup[mol]
                if rxn not in prev_gens_rxns
                and all(
                    mol not in prev_gens_mols
                    for mol in self._rxn_lib[rxn].reactants
                )
            ]
            rxnsets.append(newrxnset)
        if not fail_on_unknown_reagent:
            rxnsets = [rxnset for rxnset in rxnsets if len(rxnset) > 0]
        if len(rxnsets) == 0:
            yield []
            return
        tested_combos = set()
        for rxncombo in iterproduct(*rxnsets):
            rxncombo = frozenset(rxncombo)
            if rxncombo in tested_combos:
                continue
            else:
                tested_combos.add(frozenset(rxncombo))
            required_reagents = set(
                mol
                for mol in chain(
                    *(self._rxn_lib[rxn].reactants for rxn in rxncombo)
                )
                if mol not in reagent_table
            )
            if len(required_reagents) == 0:
                yield [rxncombo]
                continue
            for path in self._getchains(
                required_reagents,
                prev_gens_mols.union(cur_gen_mols),
                prev_gens_rxns.union(rxncombo),
                reagent_table,
                fail_on_unknown_reagent,
            ):
                path.append(rxncombo)
                yield path

    def getParentChains(
        self,
        target: Identifier,
        reagent_table: Sequence[Identifier] = tuple(),
        fail_on_unknown_reagent: bool = False,
    ) -> Generator[List[Set[RxnDatBase]], None, None]:
        if fail_on_unknown_reagent and not reagent_table:
            raise ValueError(
                "reagent table must be specified if fail_on_unknown_reagent is True"
            )
        return (
            path
            for path in self._getchains(
                [target],
                reagent_table=reagent_table,
                fail_on_unknown_reagent=fail_on_unknown_reagent,
            )
        )
