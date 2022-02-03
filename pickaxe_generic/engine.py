"""
Contains classes which define and implement dependency-injection engines.
Also contains master configuration function create_engine.

Functions:

    create_engine

Classes:

    NetworkEngine
      NetworkEngineBasic*
"""


from abc import ABC, abstractmethod
from msilib.schema import Class
from typing import Any, Iterable, Optional, Tuple, Union

from rdkit.Chem.rdchem import Mol as RDKitMol
from rdkit.Chem.rdChemReactions import ChemicalReaction as RDKitRxn

from pickaxe_generic.containers import (
    ObjectLibrary,
    ObjectLibraryBasic,
    ObjectLibraryKeyVal,
)
from pickaxe_generic.datatypes import (
    Identifier,
    MolDatBase,
    MolDatBasicV1,
    MolDatBasicV2,
    MolDatRDKit,
    OpDatBase,
    OpDatBasic,
    OpDatRDKit,
    RxnDatBase,
    RxnDatBasic,
)
from pickaxe_generic.strategies import CartesianStrategy


def create_engine(speed: int = 5):
    return NetworkEngineBasic(speed=speed)


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
    _Mol: Any
    _Mol_Lib: Any
    _Op_Lib: Any
    _Rxn_Lib: Any

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
        return self._Mol(molecule=molecule, sanitize=sanitize, neutralize=neutralize)

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
