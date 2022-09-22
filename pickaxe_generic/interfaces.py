"""
Contains interfaces for major datatypes in Pickaxe-Generic.

Classes:

    Identifier(Protocol)
    DataUnit
      MolDatBase
        MolDatRDKit
      OpDatBase
        OpDatRDKit
      RxnDatBase (deprecated)
    MetaKeyPacket
    DataPacket
    DataPacketE
    MolFilter
      MolFilterAnd*
      MolFilterInv*
      MolFilterOr*
      MolFilterXor*
    Reaction
    ReactionExplicit
    Recipe
    RecipeExplicit
    RecipeFilter
      RecipeFilterAnd*
      RecipeFilterInv*
      RecipeFilterOr*
      RecipeFilterXor*
    NetworkEngine
    ValueQueryData
    ValueQueryAssoc
    ChemNetwork
    RankValue
    RecipeRanker
    PriorityQueueStrategy
    RxnTracker
    MetaDataCalculatorLocal
    MetaDataUpdate
    ObjectLibrary (deprecated)
    ExpansionStrategy (deprecated)
    ReactionFilter (deprecated)
    UIDPreFilter (deprecated)
"""

import abc
import collections.abc
import dataclasses
import typing

import rdkit
import rdkit.Chem
import rdkit.Chem.rdChemReactions

T = typing.TypeVar("T")
T_ci = typing.TypeVar("T_ci", contravariant=False)
T_ci_co = typing.TypeVar("T_ci_co", covariant=True)
T_data = typing.TypeVar("T_data", bound="DataUnit")
T_id = typing.TypeVar("T_id", bound="Identifier")
T_int = typing.TypeVar("T_int", bound=int)


MolIndex = typing.NewType("MolIndex", int)
OpIndex = typing.NewType("OpIndex", int)
RxnIndex = typing.NewType("RxnIndex", int)


class Identifier(collections.abc.Hashable, typing.Protocol):
    """
    Unique identifier/name of object.

    Should be orderable, hashable, and immutable.

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

    def __eq__(self, other) -> bool:
        """
        Compare object to others of similar type.

        If x == y, then x is equivalent to y.  This enables hashtables.

        Arguments
        ---------
        other
            Object to be compared.

        Returns
        -------
        bool
            True if object is equivalent to other, False otherwise.
        """

    def __lt__(self, other) -> bool:
        """
        Compare object to others of similar type.

        If canonical, this functionality allows canonical sorting.

        Arguments
        ---------
        other
            Object to be compared.

        Returns
        -------
        bool
            True if other is after self when ordered, False otherwise.
        """


class DataUnit(abc.ABC):
    """
    Basic data storage object.

    Object which provides a unique, hashable identifier, and can serve up a
    binary form of the object.

    Attributes
    ----------
    blob : bytes
        Binary representation of object.
    uid : pickaxe_generic.interfaces.Identifier
        Unique identifier of object.
    """

    __slots__ = ()

    @property
    @abc.abstractmethod
    def blob(self) -> bytes:
        """
        Binary representation of object.

        Must be able to initialize object when passed to __setstate__ method of
        any sibling class with same immediate parent (viz. initialize a
        MolDatBasicV2, even if obtained from a MolDatBasicV1).
        """

    @classmethod
    @abc.abstractmethod
    def from_bytes(
        cls: type[T_ci], data: bytes, engine: "NetworkEngine"
    ) -> T_ci:
        """
        Produce object from bytestring.

        Bytestring should be derived from the .blob property.
        """

    @property
    @abc.abstractmethod
    def uid(self) -> Identifier:
        """
        Unique identifier of object.

        Must be hashable in order to facilitate lookup tables utilizing hashes.
        """


class MolDatBase(DataUnit):
    """
    Empty interface signalling molecule status of DataUnit class.

    Classes implementing this interface manage information about a single
    molecule, allowing for memory management and lumped molecule frameworks.

    Attributes
    ----------
    blob : bytes
        Binary representation of object.
    uid : pickaxe_generic.interfaces.Identifier
        Unique identifier of object.
    """

    __slots__ = ()


class MolDatRDKit(MolDatBase):
    """
    Interface representing an RDKit molecule data object.

    Classes implementing this interface manage information about a single
    rdkit-compatible molecule.  Defines the constructor for this type of
    MolDat, and thus must be subclassed only by implementations with
    @typing.final.

    Parameters
    ----------
    molecule : rdkit.Chem.rdchem.Mol | str
        Sufficient information to generate molecule in the form of an RDKit Mol
        or a SMILES string.
    sanitize : bool (default: True)
        Run through RDKit sanitizer first.  Should be True when using input
        from non-sanitized sources.  Can be disabled for speed improvements but
        canonicity is not guaranteed.

    Attributes
    ----------
    blob : bytes
    inchikey : str
    rdkitmol : rdkit.Chem.rdchem.Mol
    smiles : str
    uid : pickaxe_generic.interfaces.Identifier

    Other Parameters
    ----------------
    neutralize : bool (default: False)
        Should be True if you want hydrogens to be added/subtracted to
        neutralize a molecule and input is a SMILES string or non-neutralized
        molecule.
    """

    __slots__ = ()

    @abc.abstractmethod
    def __init__(
        self,
        molecule: typing.Union[rdkit.Chem.rdchem.Mol, str],
        sanitize: bool = True,
        neutralize: bool = False,
    ) -> None:
        ...

    @property
    @abc.abstractmethod
    def blob(self) -> bytes:
        """
        RDKit-generated bytestring.

        Enables quick regeneration of RDKitMol object.

        Returns
        -------
        bytes
            Binary representation of RDKit molecule, generated from native
            RDKit methods.
        """

    @classmethod
    def from_bytes(
        cls: type["MolDatRDKit"], data: bytes, engine: "NetworkEngine"
    ) -> "MolDatRDKit":
        """
        Generate new RDKit molecule from bytestring, according to engine
        configuration.

        Parameters
        ----------
        data : bytes
            Bytestring containing sufficient binary information to initialize
            molecule.
        engine : pickaxe_generic.interfaces.NetworkEngine
            Engine containing settings for molecule initialization.

        Returns
        -------
        pickaxe_generic.interfaces.MolDatRDKit
            Molecule returned from processing bytestring.
        """
        return engine.Mol(data)

    @property
    @abc.abstractmethod
    def inchikey(self) -> str:
        """
        InChIKey hash of molecule.

        Returns
        -------
        str
            InChIKey hash in string form, generated from RDKit.
        """

    @property
    @abc.abstractmethod
    def rdkitmol(self) -> rdkit.Chem.rdchem.Mol:
        """
        RDKit molecule object.

        Returns
        -------
        rdkit.Chem.rdchem.Mol
            RDKit molecule, generated from RDKit.
        """

    @property
    @abc.abstractmethod
    def smiles(self) -> str:
        """
        SMILES string.

        If the molecule has been initialized with sanitize=True, then this is
        the canonical RDKit SMILES string.

        Returns
        -------
        str
            SMILES string.
        """

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
    Interface representing an RDKit operator (ChemicalReaction) data object.

    Classes implementing this interface manage information about a single
    operator which acts on a molecule and can generate tuples of product
    molecules.

    Attributes
    ----------
    blob : bytes
        Binary representation of operator.
    uid : pickaxe_generic.interfaces.Identifier
        Unique identifier of operator.

    Methods
    -------
    __call__
    __len__
    compat
    """

    __slots__ = ()

    @abc.abstractmethod
    def __call__(
        self, reactants: collections.abc.Sequence[MolDatBase]
    ) -> collections.abc.Iterable[collections.abc.Iterable[MolDatBase]]:
        """
        React a sequence of MolDatBase objects using internal operator.

        For every combination of reaction sites which is possible, there is a
        set of product molecules.  This method returns an iterator over each of
        these sets, which may not be unique.

        Parameters
        ----------
        reactants : collections.abc.Sequence[pickaxe_generic.interfaces.MolDatBase]
            Reactants which match the arguments in the operator.

        Returns
        -------
        collections.abc.Iterable[collections.abc.Iterable[pickaxe_generic.interfaces.RxnDatBase]]
            Iterable of reaction product sets.
        """

    @abc.abstractmethod
    def __len__(self) -> int:
        """
        Return number of arguments in operator.

        Returns
        -------
        int
            Number of arguments in operator.
        """

    @abc.abstractmethod
    def compat(self, mol: MolDatBase, arg: int) -> bool:
        """
        Determine compatibility of molecule object with operator argument.

        This feature enables caching of molecule-operator compatibility for more
        efficient network expansion, avoiding the need to test for the presence
        of reactive sites on each molecule each time the operator is called.

        Parameters
        ----------
        mol : pickaxe_generic.interfaces.MolDatBase
            MolDat object which is to be compared.
        arg : int
            Index of argument which is to be compared.

        Returns
        -------
        bool
            True if molecule is compatible with specified operator argument.
        """


class OpDatRDKit(OpDatBase):
    """
    Interface representing an RDKit SMARTS operator.

    Agents are treated as arguments following reagent arguments.  Classes
    implementing this interface manage information about a single
    rdkit-compatible SMARTS operator.

    Attributes
    ----------
    blob : bytes
        Binary representation of operator.
    rdkitrxn : rdkit.Chem.rdChemReactions.ChemicalReaction
        RDKit reaction object.
    smarts : str
        SMARTS string representing operator.
    uid : pickaxe_generic.interfaces.Identifier
        Unique identifier of operator.
    """

    __slots__ = ()

    @abc.abstractmethod
    def __init__(self, operator: rdkit.Chem.rdchem.Mol | str | bytes) -> None:
        ...

    @classmethod
    def from_bytes(
        cls,
        data: bytes,
        engine: "NetworkEngine",
    ) -> "OpDatRDKit":
        """
        Generate new RDKit operator from bytestring, according to engine
        configuration.

        Parameters
        ----------
        data : bytes
            Bytestring containing sufficient binary information to initialize
            operator.
        engine : pickaxe_generic.interfaces.NetworkEngine
            Engine containing settings for operator initialization.

        Returns
        -------
        pickaxe_generic.interfaces.OpDatRDKit
            Operator returned from processing bytestring.
        """
        return engine.Op(data)

    @property
    @abc.abstractmethod
    def rdkitrxn(self) -> rdkit.Chem.rdChemReactions.ChemicalReaction:
        """
        RDKit operator object.

        Returns
        -------
        rdkit.Chem.rdChemReactions.ChemicalReaction
            RDKit operator, generated from RDKit.
        """

    @property
    @abc.abstractmethod
    def smarts(self) -> str:
        """
        SMARTS reaction string corresponding to operator.

        Returns
        -------
        str
            SMARTS string.
        """


class RxnDatBase(DataUnit):
    """
    Interface representing reaction data.

    .. deprecated:: 0.3.0
        Reactions are no longer represented as DataUnits with the advent of
        the ChemNetwork object.  Reactions are instead associations of an
        operator, ordered reactants, ordered products, and metadata.  If a
        single struct is required, suggest using
        pickaxe_generic.interfaces.Reaction if network is available or
        pickaxe_generic.interfaces.ReactionExplicit if not.

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

    @typing.final
    @classmethod
    def from_bytes(
        cls,
        data: bytes,
        engine: "NetworkEngine",
    ) -> "RxnDatBase":
        return engine.Rxn(data)

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


@typing.final
@dataclasses.dataclass(frozen=True)
class MetaKeyPacket:
    """
    Dataclass containing information about metadata keys necessary for filters
    and process-local metadata calculators to operate.

    Attributes
    ----------
    operator_keys : collections.abc.Set[collections.abc.Hashable] (default: frozenset())
        Collection of operator keys.
    molecule_keys : collections.abc.Set[collections.abc.Hashable] (default: frozenset())
        Collection of molecule keys.
    live_operator : bool (default: False)
        Whether initialized operators are required (True) or only index + UID (False).
    live_molecule : bool (default: False)
        Whether initialized molecules are required (True) or only index + UID (False).

    Examples
    --------
    >>> MetaKeyPacket(molecule_keys={"enthalpy","generation"})
    MetaKeyPacket(operator_keys=frozenset(), molecule_keys={'enthalpy', 'generation'}, live_operator=False, live_molecule=False)

    Combining two MetaKeyPackets.

    >>> MetaKeyPacket(molecule_keys={"enthalpy","generation"}) + MetaKeyPacket(molecule_keys={"cost","toxicity"}, live_operator=True)
    MetaKeyPacket(operator_keys=frozenset(), molecule_keys={'enthalpy', 'cost', 'generation', 'toxicity'}, live_operator=True, live_molecule=False)
    """

    operator_keys: collections.abc.Set[collections.abc.Hashable] = frozenset()
    molecule_keys: collections.abc.Set[collections.abc.Hashable] = frozenset()
    live_operator: bool = False
    live_molecule: bool = False

    def __add__(self, other: "MetaKeyPacket") -> "MetaKeyPacket":
        """
        Combine requirements from two MetaKeyPacket objects.

        The new MetaKeyPacket contains the union of the keysets and the result
        of the "or" operator on the boolean values.

        Returns
        -------
        MetaKeyPacket
            The combined requirements of both original MetaKeyPackets.
        """
        return MetaKeyPacket(
            self.operator_keys | other.operator_keys,
            self.molecule_keys | other.molecule_keys,
            self.live_operator or other.live_operator,
            self.live_molecule or other.live_molecule,
        )


@dataclasses.dataclass(frozen=True)
class DataPacket(typing.Generic[T_data]):
    """
    Dataclass containing information about a particular DataUnit.

    Attributes
    ----------
    i : int
        Index of the DataUnit in some ChemNetwork.
    item : typing.Optional[DataUnit]
        The value of the DataUnit itself.
    meta : typing.Optional[collections.abc.Mapping]
        Metadata associated with the DataUnit.
    """

    __slots__ = ("i", "item", "meta")
    i: int
    item: typing.Optional[T_data]
    meta: typing.Optional[collections.abc.Mapping]


@dataclasses.dataclass(frozen=True)
class DataPacketE(DataPacket, typing.Generic[T_data]):
    """mostly unused type???"""

    __slots__ = ("item",)
    item: T_data


class MolFilter(abc.ABC):
    """
    Interface representing a molecule filter.

    Classes which implement this interface provide a method for evaluating a
    molecule's information and returning True or False.

    Attributes
    ----------
    meta_required : MetaKeyPacket
        Metadata required in order for filter to evaluate molecule.

    Notes
    -----
    For developing new MolFilters: if metadata is required in order to run the
    filter, then meta_required must be provided or that metadata is not
    guaranteed to be included.
    """

    __slots__ = ()

    @abc.abstractmethod
    def __call__(self, mol: DataPacket[MolDatBase]) -> bool:
        """
        Evaluate a DataPacket using filter function.

        Parameters
        ----------
        mol : DataPacket[MolDatBase]
            DataPacket containing information about a molecule.

        Returns
        -------
        bool
            Whether or not the DataPacket representing a molecule passes the
            filter.
        """

    @property
    def meta_required(self) -> MetaKeyPacket:
        """
        Specifier for information required by filter function.

        Returns
        -------
        MetaKeyPacket
            MetaKeyPacket containing information on which molecule metadata is
            necessary to run filter, and if live molecule is required.
        """

    @typing.final
    def __and__(self, other: "MolFilter") -> "MolFilterAnd":
        """
        Compose the intersection of two MolFilters.

        Returns
        -------
        MolFilterAnd
            MolFilter which returns the result of `self(mol) and other(mol)`.
        """
        return MolFilterAnd(self, other)

    @typing.final
    def __invert__(self) -> "MolFilterInv":
        """
        Invert the value of the MolFilter.

        Returns
        -------
        MolFilterInv
            MolFilter which returns the result of `not self(mol)`.
        """
        return MolFilterInv(self)

    @typing.final
    def __or__(self, other: "MolFilter") -> "MolFilterOr":
        """
        Compose the union of two MolFilters.

        Returns
        -------
        MolFilterOr
            MolFilter which returns the result of `self(mol) or other(mol)`.
        """
        return MolFilterOr(self, other)

    @typing.final
    def __xor__(self, other: "MolFilter") -> "MolFilterXor":
        """
        Compose the symmetric difference of two MolFilters.

        Returns
        -------
        MolFilterXor
            MolFilter which returns the result of `self(mol) != other(mol)`
        """
        return MolFilterXor(self, other)


@typing.final
@dataclasses.dataclass(frozen=True)
class MolFilterAnd(MolFilter):
    """
    Class which composes the intersection of two filters.

    Notes
    -----
    If initialized with filter1 and filter2, when called will return the result
    of `filter1(mol) and filter2(mol)`.
    """

    __slots__ = ("_filter1", "_filter2")

    _filter1: MolFilter
    _filter2: MolFilter

    def __call__(self, mol: DataPacket[MolDatBase]) -> bool:
        """
        Calculate the intersection result of the composed filters.

        Parameters
        ----------
        mol : DataPacket[MolDatBase]
            DataPacket containing information about a molecule.

        Returns
        -------
        bool
            If initialized with filter1 and filter2, returns the result of
            `filter1(mol) and filter2(mol)`.
        """
        return self._filter1(mol) and self._filter2(mol)

    @property
    def meta_required(self) -> MetaKeyPacket:
        """
        Return the union of the specifiers for the composed filter functions.

        Returns
        -------
        MetaKeyPacket
            The combined metadata key packet for both composed filter functions.
        """
        return self._filter1.meta_required + self._filter2.meta_required


@typing.final
@dataclasses.dataclass(frozen=True)
class MolFilterInv(MolFilter):
    __slots__ = ("_filter",)
    _filter: MolFilter

    def __call__(self, mol: DataPacket[MolDatBase]) -> bool:
        """
        Calculate the inverse result of the composed filter.

        Parameters
        ----------
        mol : DataPacket[MolDatBase]
            DataPacket containing information about a molecule.

        Returns
        -------
        bool
            If initialized with filter1, returns the result of
            `not filter1(mol)`.
        """
        return not self._filter(mol)

    @property
    def meta_required(self) -> MetaKeyPacket:
        """
        Return the specifiers of the composed filter function.

        Returns
        -------
        MetaKeyPacket
            The metadata key packet of the composed filter function.
        """
        return self._filter.meta_required


@typing.final
@dataclasses.dataclass(frozen=True)
class MolFilterOr(MolFilter):
    """
    Class which composes the union of two filters.

    Notes
    -----
    If initialized with filter1 and filter2, when called will return the result
    of `filter1(mol) or filter2(mol)`.
    """

    __slots__ = ("_filter1", "_filter2")
    _filter1: MolFilter
    _filter2: MolFilter

    def __call__(self, mol: DataPacket[MolDatBase]) -> bool:
        """
        Calculate the union result of the composed filters.

        Parameters
        ----------
        mol : DataPacket[MolDatBase]
            DataPacket containing information about a molecule.

        Returns
        -------
        bool
            If initialized with filter1 and filter2, returns the result of
            `filter1(mol) or filter2(mol)`.
        """
        return self._filter1(mol) or self._filter2(mol)

    @property
    def meta_required(self) -> MetaKeyPacket:
        """
        Return the union of the specifiers for the composed filter functions.

        Returns
        -------
        MetaKeyPacket
            The combined metadata key packet for both composed filter functions.
        """
        return self._filter1.meta_required + self._filter2.meta_required


@typing.final
@dataclasses.dataclass(frozen=True)
class MolFilterXor(MolFilter):
    """
    Class which composes the symmetric difference of two filters.

    Notes
    -----
    If initialized with filter1 and filter2, when called will return the result
    of `filter1(mol) != filter2(mol)`.
    """

    __slots__ = ("_filter1", "_filter2")
    _filter1: MolFilter
    _filter2: MolFilter

    def __call__(self, mol: DataPacket[MolDatBase]) -> bool:
        """
        Calculate the symmetric difference result of the composed filters.

        Parameters
        ----------
        mol : DataPacket[MolDatBase]
            DataPacket containing information about a molecule.

        Returns
        -------
        bool
            If initialized with filter1 and filter2, returns the result of
            `filter1(mol) != filter2(mol)`.
        """
        return self._filter1(mol) != self._filter2(mol)

    @property
    def meta_required(self) -> MetaKeyPacket:
        """
        Return the union of the specifiers for the composed filter functions.

        Returns
        -------
        MetaKeyPacket
            The combined metadata key packet for both composed filter functions.
        """
        return self._filter1.meta_required + self._filter2.meta_required


@dataclasses.dataclass(frozen=True, order=True)
class Reaction:
    __slots__ = ("operator", "reactants", "products")
    operator: OpIndex
    reactants: tuple[MolIndex, ...]
    products: tuple[MolIndex, ...]


@dataclasses.dataclass(frozen=True, order=True)
class ReactionExplicit:
    __slots__ = (
        "operator",
        "reactants",
        "products",
        "reaction_meta",
    )
    operator: DataPacketE[OpDatBase]
    reactants: tuple[DataPacketE[MolDatBase], ...]
    products: tuple[DataPacketE[MolDatBase], ...]
    reaction_meta: typing.Optional[collections.abc.Mapping]

    @property
    def uid(
        self,
    ) -> tuple[Identifier, tuple[Identifier, ...], tuple[Identifier, ...]]:
        return (
            self.operator.item.uid,
            tuple(mol.item.uid for mol in self.reactants),
            tuple(mol.item.uid for mol in self.products),
        )


@dataclasses.dataclass(frozen=True)
class Recipe:
    __slots__ = ("operator", "reactants")
    operator: OpIndex
    reactants: tuple[MolIndex, ...]

    def __eq__(self, other: object) -> bool:
        if (
            isinstance(other, Recipe)
            and self.operator == other.operator
            and self.reactants == other.reactants
        ):
            return True
        return False

    def __lt__(self, other: "Recipe") -> bool:
        self_order = sorted(self.reactants, reverse=True)
        other_order = sorted(other.reactants, reverse=True)
        for val_self, val_other in zip(self_order, other_order):
            if val_self < val_other:
                return False
            elif val_other < val_self:
                return True
        if len(self.reactants) < len(other.reactants):
            return False
        elif len(other.reactants) < len(self.reactants):
            return True
        if self.operator < other.operator:
            return False
        elif other.operator < self.operator:
            return True
        return other.reactants < self.reactants


@dataclasses.dataclass(frozen=True)
class RecipeExplicit:
    __slots__ = ("operator", "reactants", "operator_meta", "reactants_meta")
    operator: DataPacket[OpDatBase]
    reactants: tuple[DataPacket[MolDatBase], ...]


class RecipeFilter(abc.ABC):
    __slots__ = ()

    @abc.abstractmethod
    def __call__(self, recipe: RecipeExplicit) -> bool:
        ...

    @property
    def meta_required(self) -> MetaKeyPacket:
        return MetaKeyPacket()

    @typing.final
    def __and__(self, other: "RecipeFilter") -> "RecipeFilter":
        return RecipeFilterAnd(self, other)

    @typing.final
    def __invert__(self) -> "RecipeFilter":
        return RecipeFilterInv(self)

    @typing.final
    def __or__(self, other: "RecipeFilter") -> "RecipeFilter":
        return RecipeFilterOr(self, other)

    @typing.final
    def __xor__(self, other: "RecipeFilter") -> "RecipeFilter":
        return RecipeFilterXor(self, other)


@dataclasses.dataclass(frozen=True)
class RecipeFilterAnd(RecipeFilter):
    __slots__ = ("_filter1", "_filter2")

    _filter1: RecipeFilter
    _filter2: RecipeFilter

    def __call__(self, recipe: RecipeExplicit) -> bool:
        return self._filter1(recipe) and self._filter2(recipe)

    @property
    def meta_required(self) -> MetaKeyPacket:
        return self._filter1.meta_required + self._filter2.meta_required


@dataclasses.dataclass(frozen=True)
class RecipeFilterInv(RecipeFilter):
    __slots__ = ("_filter",)
    _filter: RecipeFilter

    def __call__(self, recipe: RecipeExplicit) -> bool:
        return not self._filter(recipe)

    @property
    def meta_required(self) -> MetaKeyPacket:
        return self._filter.meta_required


@dataclasses.dataclass(frozen=True)
class RecipeFilterOr(RecipeFilter):
    __slots__ = ("_filter1", "_filter2")
    _filter1: RecipeFilter
    _filter2: RecipeFilter

    def __call__(self, recipe: RecipeExplicit) -> bool:
        return self._filter1(recipe) or self._filter2(recipe)

    @property
    def meta_required(self) -> MetaKeyPacket:
        return self._filter1.meta_required + self._filter2.meta_required


@dataclasses.dataclass(frozen=True)
class RecipeFilterXor(RecipeFilter):
    __slots__ = ("_filter1", "_filter2")
    _filter1: RecipeFilter
    _filter2: RecipeFilter

    def __call__(self, recipe: RecipeExplicit) -> bool:
        return self._filter1(recipe) != self._filter2(recipe)

    @property
    def meta_required(self) -> MetaKeyPacket:
        return self._filter1.meta_required + self._filter2.meta_required


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
        kekulize: bool = False,
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
    ) -> "RxnDatBase":
        """
        Initializes a RxnDatBase object of relevant type.
        """

    @abc.abstractmethod
    def Libs(
        self,
    ) -> tuple[
        "ObjectLibrary[MolDatBase]",
        "ObjectLibrary[OpDatBase]",
        "ObjectLibrary[RxnDatBase]",
    ]:
        """
        Initializes the three basic ObjectLibraries necessary to run a Strategy.
        """

    @abc.abstractmethod
    def CartesianStrategy(
        self,
        mol_lib: "ObjectLibrary[MolDatBase]",
        op_lib: "ObjectLibrary[OpDatBase]",
        rxn_lib: "ObjectLibrary[RxnDatBase]",
    ):
        """
        Initializes a CartesianStrategy of relevant type.
        """

    @abc.abstractmethod
    def from_bytes(
        self,
        cls: type[T_ci],
        data: bytes,
    ) -> T_ci:
        ...


class ValueQueryData(typing.Protocol[T_data, T_int]):
    @abc.abstractmethod
    def __contains__(self, item: typing.Union[Identifier, T_data]) -> bool:
        ...

    @typing.overload
    @abc.abstractmethod
    def __getitem__(self, item: slice) -> collections.abc.Sequence[T_data]:
        ...

    @typing.overload
    @abc.abstractmethod
    def __getitem__(self, item: typing.Union[T_int, Identifier]) -> T_data:
        ...

    @abc.abstractmethod
    def __getitem__(self, item: typing.Union[slice, T_int, Identifier]):
        ...

    @abc.abstractmethod
    def i(self, uid: Identifier) -> T_int:
        ...

    @abc.abstractmethod
    def keys(self) -> collections.abc.Collection[Identifier]:
        ...

    @abc.abstractmethod
    def uid(self, i: T_int) -> Identifier:
        ...

    @abc.abstractmethod
    def __len__(self) -> int:
        ...

    @abc.abstractmethod
    def __iter__(self) -> collections.abc.Iterator[T_data]:
        ...


class ValueQueryAssoc(typing.Protocol[T_id, T_int]):
    @typing.overload
    @abc.abstractmethod
    def __getitem__(self, item: slice) -> collections.abc.Sequence[T_id]:
        ...

    @typing.overload
    @abc.abstractmethod
    def __getitem__(self, item: T_int) -> T_id:
        ...

    @abc.abstractmethod
    def __getitem__(self, item: typing.Union[slice, T_int]):
        ...

    @abc.abstractmethod
    def i(self, item: T_id) -> T_int:
        ...

    @abc.abstractmethod
    def __len__(self) -> int:
        ...

    @abc.abstractmethod
    def __iter__(self) -> collections.abc.Iterator[T_id]:
        ...


class ChemNetwork(abc.ABC):
    __slots__ = ()

    @abc.abstractmethod
    def __init__(self) -> None:
        ...

    @property
    @abc.abstractmethod
    def mols(self) -> ValueQueryData[MolDatBase, MolIndex]:
        ...

    @property
    @abc.abstractmethod
    def ops(self) -> ValueQueryData[OpDatBase, OpIndex]:
        ...

    @property
    @abc.abstractmethod
    def rxns(self) -> ValueQueryAssoc[Reaction, RxnIndex]:
        ...

    @abc.abstractmethod
    def mol_meta(
        self, index: MolIndex, key: collections.abc.Hashable, value=None
    ):
        ...

    @abc.abstractmethod
    def op_meta(
        self, index: OpIndex, key: collections.abc.Hashable, value=None
    ):
        ...

    @abc.abstractmethod
    def rxn_meta(
        self, index: RxnIndex, key: collections.abc.Hashable, value=None
    ):
        ...

    @abc.abstractmethod
    def mol_metas(
        self,
        indices: typing.Optional[collections.abc.Sequence[MolIndex]] = None,
        keys: typing.Optional[
            collections.abc.Collection[collections.abc.Hashable]
        ] = None,
    ) -> collections.abc.Sequence[
        collections.abc.Mapping[collections.abc.Hashable, typing.Any]
    ]:
        ...

    @abc.abstractmethod
    def op_metas(
        self,
        indices: typing.Optional[collections.abc.Sequence[OpIndex]] = None,
        keys: typing.Optional[
            collections.abc.Collection[collections.abc.Hashable]
        ] = None,
    ) -> collections.abc.Sequence[
        collections.abc.Mapping[collections.abc.Hashable, typing.Any]
    ]:
        ...

    @abc.abstractmethod
    def rxn_metas(
        self,
        indices: typing.Optional[collections.abc.Sequence[RxnIndex]] = None,
        keys: typing.Optional[
            collections.abc.Collection[collections.abc.Hashable]
        ] = None,
    ) -> collections.abc.Sequence[
        collections.abc.Mapping[collections.abc.Hashable, typing.Any]
    ]:
        ...

    @abc.abstractmethod
    def compat_table(
        self, index: int
    ) -> collections.abc.Sequence[collections.abc.Sequence[MolIndex]]:
        ...

    @abc.abstractmethod
    def consumers(
        self, mol: typing.Union[int, MolDatBase, Identifier]
    ) -> collections.abc.Collection[RxnIndex]:
        ...

    @abc.abstractmethod
    def producers(
        self, mol: typing.Union[int, MolDatBase, Identifier]
    ) -> collections.abc.Collection[RxnIndex]:
        ...

    @abc.abstractmethod
    def add_mol(
        self,
        mol: MolDatBase,
        meta: typing.Optional[collections.abc.Mapping] = None,
        reactive: bool = True,
        custom_compat: typing.Optional[
            collections.abc.Collection[tuple[OpIndex, int]]
        ] = None,
    ) -> MolIndex:
        ...

    @abc.abstractmethod
    def add_op(
        self,
        mol: OpDatBase,
        meta: typing.Optional[collections.abc.Mapping] = None,
    ) -> OpIndex:
        ...

    @abc.abstractmethod
    def add_rxn(
        self,
        rxn: typing.Optional[Reaction] = None,
        op: typing.Optional[OpIndex] = None,
        reactants: typing.Optional[collections.abc.Sequence[MolIndex]] = None,
        products: typing.Optional[collections.abc.Sequence[MolIndex]] = None,
        meta: typing.Optional[collections.abc.Mapping] = None,
    ) -> RxnIndex:
        ...


class RankValue(typing.Protocol):
    @abc.abstractmethod
    def __lt__(self, other: "RankValue") -> bool:
        ...


class RecipeRanker(typing.Protocol):
    @abc.abstractmethod
    def __call__(
        self,
        recipe: RecipeExplicit,
        min_rank: typing.Optional[RankValue] = None,
    ) -> typing.Optional[RankValue]:
        ...

    @property
    def meta_required(self) -> MetaKeyPacket:
        return MetaKeyPacket()


class PriorityQueueStrategy(abc.ABC):
    __slots__ = ()

    @abc.abstractmethod
    def __init__(
        self,
        network: ChemNetwork,
        num_procs: typing.Optional[int] = None,
        blacklist_key: typing.Optional[str] = None,
    ) -> None:
        ...

    @abc.abstractmethod
    def expand(
        self,
        max_recipes: typing.Optional[int] = None,
        heap_size: int = 1,
        batch_size: typing.Optional[int] = None,
        beam_size: typing.Optional[int] = 1,
        # mol_filter_local: Optional[MolFilter] = None,
        # mol_filter: Optional[MolFilter] = None,
        recipe_filter: typing.Optional[RecipeFilter] = None,
        recipe_ranker: typing.Optional[RecipeRanker] = None,
        # mc_local: Optional[MetaDataCalculatorLocal] = None,
        # mc_update: Optional[MetaDataUpdate] = DefaultMetaDataUpdate(),
    ) -> None:
        ...


class RxnTracker(abc.ABC):
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

    @abc.abstractmethod
    def getParentChains(
        self,
        target: Identifier,
        reagent_table: collections.abc.Sequence[Identifier] = tuple(),
        fail_on_unknown_reagent: bool = False,
    ) -> collections.abc.Iterable[
        collections.abc.Iterable[collections.abc.Iterable[Identifier]]
    ]:
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


class MetaDataCalculatorLocal(typing.Protocol):
    @abc.abstractmethod
    def __call__(
        self, unit: typing.Union[ReactionExplicit, RecipeExplicit]
    ) -> None:
        ...

    @property
    @abc.abstractmethod
    def meta_required(self) -> MetaKeyPacket:
        ...


class MetaDataUpdate(typing.Protocol):
    @abc.abstractmethod
    def __call__(
        self, unit: ReactionExplicit, network: ChemNetwork
    ) -> collections.abc.Generator[
        tuple[
            typing.Optional[tuple[MolIndex, collections.abc.Hashable]],
            typing.Optional[tuple[OpIndex, collections.abc.Hashable]],
        ],
        None,
        None,
    ]:
        ...


class ObjectLibrary(abc.ABC, typing.Generic[T_data]):
    """
    Interface representing library of data.

    .. deprecated:: 0.3.0
        `ObjectLibrary` will be removed in pickaxe_generic 0.4.0.  The pattern
        of several ObjectLibraries representing a network has been replaced by
        the ChemNetwork object representing all nodes.  However,
        ObjectLibraries generated through a NetworkEngine will temporarily be
        available via a facade to ChemNetwork.

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
        Iterator[T_data]
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

    .. deprecated:: 0.3.0
        `ExpansionStrategy` will be removed in pickaxe_generic 0.4.0.
        This definition was determined to be too expansive to capture the wide
        number and variety of possible expansion operations.  Recommend use of
        engine.strats.pq() or engine.strats.cartesian().

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


class ReactionFilter(abc.ABC):
    """
    Interface representing an old-style reaction filter.

    .. deprecated:: 0.3.0
        `ReactionFilter` will be removed in pickaxe_generic 0.4.0.
        ReactionFilter is not general enough to enable effective reaction
        filtering; a new struct (ReactionFilterBase) completely replaces it
        paradigm by permitting metadata and integration with analysis flows.
    """

    __slots__ = ()

    @abc.abstractmethod
    def __call__(
        self,
        operator: OpDatBase,
        reactants: collections.abc.Sequence[MolDatBase],
        products: collections.abc.Sequence[MolDatBase],
    ) -> bool:
        pass


class UIDPreFilter(abc.ABC):
    """
    Interface representing an old-style UID prefilter.

    .. deprecated:: 0.3.0
        `UIDPreFilter` will be removed in pickaxe_generic 0.4.0.  UIDPreFilter
        is not general enough to enable effective reaction filtering; a new
        struct (RecipeFilter) completely replaces it paradigm by permitting
        metadata and more general filtering options.
    """

    __slots__ = ()

    @abc.abstractmethod
    def __call__(
        self,
        operator: Identifier,
        reactants: collections.abc.Sequence[Identifier],
    ) -> bool:
        pass
