"""
Contains interfaces for major datatypes in Pickaxe-Generic.

Classes:


"""

import abc
import collections.abc
import typing

import rdkit

T = typing.TypeVar("T")
T_ci = typing.TypeVar("T_ci", contravariant=False)
T_data = typing.TypeVar("T_data", bound="DataUnit")


class Identifier(collections.abc.Hashable, typing.Protocol):
    """
    Orderable, hashable object used as unique identifier.

    This value is ideally immutable using public methods.

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
        Compare object to others of similar type.  Enables hashtables.

        Arguments
        ---------
        other : object
            Object to be compared.

        Returns
        -------
        bool
            True if object is equivalent to other, False otherwise.
        """

    def __lt__(self, other) -> bool:
        """
        Compare object to others of similar type.  Allows sorting.

        Arguments
        ---------
        other
            Object to be compared.

        Returns
        -------
        bool
            True if object is after self when ordered, False otherwise.
        """


class DataUnit(abc.ABC):
    """
    Basic data storage object.

    Object which provides a unique, hashable identifier, a method of ordering,
    and can serve up a binary form of the object.

    Attributes
    ----------
    blob : bytes
        Binary representation of object.
    uid : Identifier
        Unique identifier of object.
    """

    __slots__ = ()

    @property
    @abc.abstractmethod
    def blob(self) -> bytes:
        """
        Binary representation of object.

        Must be able to initialize object when passed to __setstate__ method of
        any subclass of same type (viz. initialize a MolDatBasicV2, even if
        obtained from a MolDatBasicV1).
        """

    @classmethod
    @abc.abstractmethod
    def from_bytes(
        cls: type[T_ci], data: bytes, engine: typing.Optional["NetworkEngine"]
    ) -> T_ci:
        """
        Produce object from bytestring.

        Bytestring should be derived from the .blob property, and should be
        able to initialize any object directly derived from the parent class.
        """

    @property
    @abc.abstractmethod
    def uid(self) -> Identifier:
        """
        Unique identifier of object.

        Must be hashable in order to facilitate lookup tables utilizing hashes.
        """

    @typing.final
    def __getstate__(self) -> bytes:
        return self.blob


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
    molecule : RDKitMol | str
        Sufficient information to generate molecule in the form of an RDKitMol
        or a SMILES string.
    sanitize : bool (default: True)
        Should be True when using input from non-sanitized sources.
    neutralize : bool (default: False)
        Should be True if you want hydrogens to be added/subtracted to
        neutralize a molecule and input is a SMILES string or non-neutralized
        molecule.

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

    @abc.abstractmethod
    def __init__(
        self,
        molecule: rdkit.Chem.rdchem.Mol | str,
        sanitize: bool = True,
        neutralize: bool = False,
    ) -> None:
        pass

    def __lt__(self, other: object) -> bool:
        if isinstance(other, MolDatRDKit):
            return self.uid < other.uid
        else:
            raise TypeError(
                f"Invalid comparison between objects of type {type(self)} and"
                f"type {type(other)}"
            )

    def __setstate__(self, in_val: bytes) -> None:
        rdkitmol: typing.Optional[rdkit.Chem.rdchem.Mol] = rdkit.Chem.Mol(
            in_val
        )
        if rdkitmol is None:
            raise ValueError("Invalid molecule bytestring")
        self._buildfrommol(rdkitmol)

    @abc.abstractmethod
    def _buildfrommol(self, in_val: rdkit.Chem.rdchem.Mol) -> None:
        pass

    @property
    @abc.abstractmethod
    def blob(self) -> bytes:
        """
        RDKit-generated bytestring.

        Enables quick regeneration of RDKitMol object.
        """

    @property
    @abc.abstractmethod
    def inchikey(self) -> str:
        """Return InChIKey hash of molecule."""

    @property
    @abc.abstractmethod
    def rdkitmol(self) -> rdkit.Chem.rdchem.Mol:
        """Return RDKit molecule object containing basic properties."""

    @property
    @abc.abstractmethod
    def smiles(self) -> str:
        """Return canonical SMILES string of molecule object."""

    def _processinput(
        self,
        molecule: rdkit.Chem.rdchem.Mol | str | bytes,
        sanitize: bool = True,
        neutralize: bool = False,
    ) -> rdkit.Chem.rdchem.Mol:
        if isinstance(molecule, bytes):
            rdkitmol = rdkit.Chem.Mol(molecule)
            # if sanitize:
            #    SanitizeMol(rdkitmol)
            #    AssignStereochemistry(rdkitmol, True, True, True)
            # if neutralize:
            #    raise NotImplementedError("No neutralize function coded")
        elif isinstance(molecule, rdkit.Chem.rdchem.Mol):
            rdkitmol = molecule
            # print(MolToSmiles(rdkitmol))
            if sanitize:
                rdkit.Chem.rdmolops.SanitizeMol(rdkitmol)
                rdkit.Chem.rdmolops.AssignStereochemistry(
                    rdkitmol, True, True, True
                )
            if neutralize:
                raise NotImplementedError("No neutralize function coded")
        elif isinstance(molecule, str):
            if sanitize:
                rdkitmol = rdkit.Chem.MolFromSmiles(molecule, sanitize=True)
            else:
                rdkitmol = rdkit.Chem.MolFromSmiles(molecule, sanitize=False)
            if neutralize:
                raise NotImplementedError("No neutralize function coded")
        else:
            raise TypeError("Invalid molecule type")
        if rdkitmol is None:
            raise TypeError("Invalid molecule information")
        return rdkitmol


class OpDatBase(DataUnit):
    """
    Interface representing operator data.

    Classes implementing this interface manage information about a single
    operator which acts on MolDatBase and can generate RxnDatBase objects.

    Methods
    -------
    compat
    __call__
    __len__
    """

    __slots__ = ()

    @abc.abstractmethod
    def compat(self, mol: MolDatBase, arg: int) -> bool:
        """
        Determine compatibility of MolDat object with operator argument.

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

    @abc.abstractmethod
    def __call__(
        self, reactants: collections.abc.Sequence[MolDatBase]
    ) -> collections.abc.Iterable[collections.abc.Iterable[MolDatBase]]:
        """
        React a sequence of MolDat objects using internal operator.

        Return a sequence of RxnDat objects which contain metadata about
        potential results.

        Parameters
        ----------
        reactants : Sequence[MolDatBase]
            Reactants which match the arguments in the operator.

        Returns
        -------
        Iterable[RxnDatBase]
            Iterable of reactions which are produced by applying the operator.
        """

    @abc.abstractmethod
    def __len__(self) -> int:
        """Return number of arguments in operator."""


class OpDatRDKit(OpDatBase):
    """
    Interface representing an RDKit SMARTS operator.

    Agents are treated as arguments following reagent arguments.  Classes
    implementing this interface manage information about a single
    rdkit-compatible SMARTS operator.

    Attributes
    ----------
    smarts : str
        SMARTS string representing operator.
    rdkitrxn : RDKitRxn
        RDKit reaction object.
    """

    __slots__ = ()

    @abc.abstractmethod
    def __init__(self, operator: rdkit.Chem.rdchem.Mol | str | bytes):
        pass

    @property
    @abc.abstractmethod
    def smarts(self) -> str:
        """Return SMARTS string encoding operator information."""

    @property
    @abc.abstractmethod
    def rdkitrxn(self) -> rdkit.Chem.rdChemReactions.ChemicalReaction:
        """Return RDKit reaction object."""


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

    @abc.abstractmethod
    def __init__(
        self,
        operator: typing.Optional[Identifier] = None,
        reactants: typing.Optional[collections.abc.Iterable[Identifier]] = None,
        products: typing.Optional[collections.abc.Iterable[Identifier]] = None,
        reaction: typing.Optional[bytes] = None,
    ) -> None:
        pass

    @property
    @abc.abstractmethod
    def operator(self) -> Identifier:
        """Return ID of operator involved in reaction."""

    @property
    @abc.abstractmethod
    def products(self) -> collections.abc.Iterable[Identifier]:
        """Return IDs of products involved in reaction."""

    @property
    @abc.abstractmethod
    def reactants(self) -> collections.abc.Iterable[Identifier]:
        """Return IDs of reactants involved in reaction."""


class ObjectLibrary(abc.ABC, typing.Generic[T_data]):
    """
    Interface representing library of data.

    Classes implementing this interface manage multiple instances of a hashable
    object, and may have responsibility for synchronization with external
    databases which may also manage this information (be volatile).  Contained
    objects must have a "uid" attribute which contains a hashable unique id.

    Current implementations assume that this library will never shrink or remove
    entries.
    """

    __slots__ = ()

    @abc.abstractmethod
    def add(
        self, obj: typing.Union[collections.abc.Iterable[T_data], T_data]
    ) -> None:
        """
        Add an object or multiple objects to the library.

        This function does not add the new item if it has the same UID as an
        item already in the library.

        Parameters
        ----------
        obj : Union[Iterable[DataUnit], DataUnit]
            Object(s) to be added.
        """

    @abc.abstractmethod
    def ids(self) -> collections.abc.Iterable[Identifier]:
        """
        Return a set of keys used in the library.

        Returns
        -------
        Iterable[Identifier]
            Ids of objects in the library.
        """

    @abc.abstractmethod
    def __contains__(self, item: T_data) -> bool:
        """
        Check if ObjectLibrary contains an object where object.uid == item.uid.

        Parameters
        ----------
        item : DataUnit
            Item to be checked against internal object list.

        Returns
        -------
        bool
            True if ObjectLibrary contains object with same UID.
        """

    @abc.abstractmethod
    def __getitem__(self, item: Identifier) -> T_data:
        """
        Return object where object.uid == item returns True.

        Parameters
        ----------
        item : Identifier
            Item to be checked against internal object list.

        Returns
        -------
        DataUnit
            Object where uid attribute is equal to item.
        """

    @abc.abstractmethod
    def __iter__(self) -> collections.abc.Iterator[T_data]:
        """
        Return an iterator over the objects contained in the ObjectLibrary.

        Returns
        -------
        Iterator[DataUnitGen]
            Iterator over objects contained in the ObjectLibrary.
        """

    @abc.abstractmethod
    def __len__(self) -> int:
        """
        Return the number of items contained in the ObjectLibrary.

        Returns
        -------
        int
            Number of items in ObjectLibrary.
        """


class ExpansionStrategy(abc.ABC):
    """
    Interface representing a network expansion strategy.

    Classes implementing this interface use information from a molecule and
    operator library to generate new reactions, which are then output to a
    reaction library.
    """

    __slots__ = ()

    @abc.abstractmethod
    def __init__(
        self,
        mol_lib: ObjectLibrary[MolDatBase],
        op_lib: ObjectLibrary[OpDatBase],
        rxn_lib: ObjectLibrary[RxnDatBase],
    ) -> None:
        pass

    @abc.abstractmethod
    def expand(
        self,
        max_rxns: typing.Optional[int] = None,
        max_mols: typing.Optional[int] = None,
        num_gens: typing.Optional[int] = None,
        custom_filter: typing.Optional[
            collections.abc.Callable[
                [
                    OpDatBase,
                    collections.abc.Sequence[MolDatBase],
                    collections.abc.Sequence[MolDatBase],
                ],
                bool,
            ]
        ] = None,
        custom_uid_prefilter: typing.Optional[
            collections.abc.Callable[
                [Identifier, collections.abc.Sequence[Identifier]], bool
            ]
        ] = None,
    ) -> None:
        """
        Expand molecule library.

        Parameters
        ----------
        max_rxns : Optional[int] (default: None)
            Limit of new reactions to add.  If None, no limit.
        max_mols : Optional[int] (default: None)
            Limit of new molecules to add.  If None, no limit.
        num_gens : Optional[int] (default: None)
            Maximum generations of reactions to enumerate.  If None, no limit.
        custom_filter: Optional[Callable[[OpDatBase, Sequence[MolDatBase],
                       Sequence[MolDatBase]], bool]] (default: None)
            Filter which selects which reactions to retain.
        custom_uid_prefilter: Optional[Callable[[Identifier,
                              Sequence[Identifier]], bool]]
            Filter which selects which operator UID and reactant UIDs to retain.
        """

    @abc.abstractmethod
    def refresh(self) -> None:
        """
        Refresh active molecules and operators from attached libraries.
        """


class NetworkEngine(abc.ABC):
    """
    Interface representing an object which serves up other objects based on
    configuration parameters.

    Classes implementing this interface determine which type of network objects
    are constructed based on configuration options.
    """

    @property
    @abc.abstractmethod
    def speed(self) -> int:
        """
        Defined speed of engine configuration.
        """

    @abc.abstractmethod
    def Mol(
        self,
        molecule: typing.Union[rdkit.Chem.rdchem.Mol, str, bytes],
        sanitize: bool = True,
        neutralize: bool = False,
    ) -> MolDatRDKit:
        """
        Initializes a MolDatRDKit object of relevant type.
        """

    @abc.abstractmethod
    def Op(
        self,
        operator: typing.Union[
            rdkit.Chem.rdChemReactions.ChemicalReaction, str, bytes
        ],
    ) -> OpDatRDKit:
        """
        Initializes an OpDatRDKit object of relevant type.
        """

    @abc.abstractmethod
    def Rxn(
        self,
        operator: typing.Optional[Identifier] = None,
        reactants: typing.Optional[collections.abc.Iterable[Identifier]] = None,
        products: typing.Optional[collections.abc.Iterable[Identifier]] = None,
        reaction: typing.Optional[bytes] = None,
    ) -> RxnDatBase:
        """
        Initializes a RxnDatBase object of relevant type.
        """

    @abc.abstractmethod
    def Libs(
        self,
    ) -> tuple[
        ObjectLibrary[MolDatBase],
        ObjectLibrary[OpDatBase],
        ObjectLibrary[RxnDatBase],
    ]:
        """
        Initializes the three basic ObjectLibraries necessary to run a Strategy.
        """

    @abc.abstractmethod
    def CartesianStrategy(
        self,
        mol_lib: ObjectLibrary[MolDatBase],
        op_lib: ObjectLibrary[OpDatBase],
        rxn_lib: ObjectLibrary[RxnDatBase],
    ):
        """
        Initializes a CartesianStrategy of relevant type.
        """
