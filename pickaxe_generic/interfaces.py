"""
Contains interfaces for major datatypes in Pickaxe-Generic.

Classes:


"""

import abc
import collections.abc
import typing

import rdkit


T = typing.TypeVar('T')
T_ci = typing.TypeVar('T_ci',contravariant=False)
T_data = typing.TypeVar('T_data',bound='DataUnit')


class Identifier(collections.abc.Hashable,typing.Protocol):
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
    def from_bytes(cls: type[T_ci], data: bytes, engine: typing.Optional[NetworkEngine]) -> T_ci:
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

    @typing.abstractmethod
    def __init__(
        self,
        molecule: RDKitMol | str,
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
        rdkitmol: Optional[RDKitMol] = BuildMol(in_val)
        if rdkitmol is None:
            raise ValueError("Invalid molecule bytestring")
        self._buildfrommol(rdkitmol)

    @abstractmethod
    def _buildfrommol(self, in_val: RDKitMol) -> None:
        pass

    @property
    @abstractmethod
    def blob(self) -> bytes:
        """
        RDKit-generated bytestring.

        Enables quick regeneration of RDKitMol object.
        """

    @property
    @abstractmethod
    def inchikey(self) -> str:
        """Return InChIKey hash of molecule."""

    @property
    @abstractmethod
    def rdkitmol(self) -> RDKitMol:
        """Return RDKit molecule object containing basic properties."""

    @property
    @abstractmethod
    def smiles(self) -> str:
        """Return canonical SMILES string of molecule object."""

    def _processinput(
        self,
        molecule: RDKitMol | str | bytes,
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

