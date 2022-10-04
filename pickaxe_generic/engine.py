"""
Contains classes which define and implement dependency-injection engines.
"""


import base64
import collections.abc
import dataclasses
import gzip
import pickle
import typing
import xml.dom.minidom

import rdkit.Chem
import rdkit.Chem.rdChemReactions

from . import (
    containers,
    datatypes,
    filters,
    interfaces,
    metacalc,
    network,
    strategies,
)


def create_engine(speed: int = 5, np: int = 1) -> interfaces.NetworkEngine:
    """
    Initialize and return a NetworkEngine based on configuration parameters.

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
    np : int (default: 1)
        Number of processes to be used by strategies.
    """
    if np < 1:
        raise ValueError("Must have number of processes greater than 0.")
    return NetworkEngineBasic(speed=speed, np=np)


class _rdkit_op_init(typing.NamedTuple):
    engine: interfaces.NetworkEngine
    optype: type[interfaces.OpDatRDKit]

    def __call__(
        self,
        operator: typing.Union[
            rdkit.Chem.rdChemReactions.ChemicalReaction, str, bytes
        ],
        kekulize: bool = False,
    ) -> interfaces.OpDatRDKit:
        """
        Creates an object which manages an RDKit SMARTS operator.

        Agents are treated as arguments following reagent arguments.  Classes
        implementing this interface manage information about a single
        rdkit-compatible SMARTS operator.

        Parameters
        ----------
        operator : typing.Union[rdkit.Chem.rdChemReactions.ChemicalReaction,
                                str, bytes]
            SMARTS string which is used to generate operator data, otherwise an
            RDKit operator itself.
        kekulize : bool (default: False)
            Whether to kekulize reactants before reaction.
        """
        return self.optype(
            operator=operator, engine=self.engine, kekulize=kekulize
        )


@dataclasses.dataclass(frozen=True, slots=True)
class _cartesian_op_init:
    _engine: interfaces.NetworkEngine

    def __call__(
        self,
        network: interfaces.ChemNetwork,
        gen_key: collections.abc.Hashable = "generation",
    ) -> strategies.CartesianStrategyUpdated:
        return strategies.CartesianStrategyUpdated(
            network, self._engine, gen_key
        )


class NetworkEngineBasic(interfaces.NetworkEngine):
    """
    Implements NetworkEngine class for different speed efficiencies.

    Default for module.

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
    np : int (default: 1)
        Number of processes to be used by strategies.
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
        "_np",
    )
    _Mol: typing.Any
    _Mol_Lib: typing.Any
    _Op_Lib: typing.Any
    _Rxn_Lib: typing.Any

    def __init__(self, speed: int = 5, np: int = 1):
        if speed == 1:
            # type: ignore
            self._Mol = datatypes.MolDatBasicV2
            raise NotImplementedError("Speed not yet implemented")
        elif speed == 2:
            self._Mol = datatypes.MolDatBasicV2
            raise NotImplementedError("Speed not yet implemented")
        elif speed == 3:
            self._Mol = datatypes.MolDatBasicV2
            raise NotImplementedError("Speed not yet implemented")
        elif speed == 4:
            self._Mol = datatypes.MolDatBasicV2
            self._Mol_Lib = lambda: containers.ObjectLibraryKeyVal(
                initializer=self.Mol
            )
            self._Op_Lib = lambda: containers.ObjectLibraryKeyVal(
                initializer=self.Op
            )
            # self._Rxn_Lib = lambda: ObjectLibraryKeyVal(initializer=self.Rxn)
            self._Rxn_Lib = containers.ObjectLibraryBasic
        elif speed == 5:
            self._Mol = datatypes.MolDatBasicV1
            self._Mol_Lib = containers.ObjectLibraryBasic
            self._Op_Lib = containers.ObjectLibraryBasic
            self._Rxn_Lib = containers.ObjectLibraryBasic
        elif speed == 6:
            self._Mol = datatypes.MolDatBasicV2
            self._Mol_Lib = containers.ObjectLibraryBasic
            self._Op_Lib = containers.ObjectLibraryBasic
            self._Rxn_Lib = containers.ObjectLibraryBasic
        else:
            raise ValueError(f"speed = {speed} is invalid")
        self._Op = datatypes.OpDatBasic
        self._Rxn = datatypes.RxnDatBasic
        self._CartesianStrategy = strategies.CartesianStrategy
        self._speed = speed
        if np != 1:
            raise NotImplementedError(
                "Multiprocessing has not yet been implemented"
            )
        self._np = np

    @property
    def speed(self) -> int:
        return self._speed

    @property
    def np(self) -> int:
        return self._np

    def Mol(
        self,
        molecule: typing.Union[rdkit.Chem.rdchem.Mol, str, bytes],
        sanitize: bool = True,
        neutralize: bool = False,
    ) -> interfaces.MolDatRDKit:
        return self._Mol(
            molecule=molecule, sanitize=sanitize, neutralize=neutralize
        )

    def Op(
        self,
        operator: typing.Union[
            rdkit.Chem.rdChemReactions.ChemicalReaction, str, bytes
        ],
        kekulize: bool = False,
    ) -> datatypes.OpDatBasic:
        return self._Op(operator=operator, engine=self, kekulize=kekulize)

    def Rxn(
        self,
        operator: typing.Optional[interfaces.Identifier] = None,
        reactants: typing.Optional[
            collections.abc.Iterable[interfaces.Identifier]
        ] = None,
        products: typing.Optional[
            collections.abc.Iterable[interfaces.Identifier]
        ] = None,
        reaction: typing.Optional[bytes] = None,
    ) -> datatypes.RxnDatBasic:
        return self._Rxn(
            operator=operator,
            reactants=reactants,
            products=products,
            reaction=reaction,
        )

    def Libs(
        self,
    ) -> tuple[
        interfaces.ObjectLibrary[interfaces.MolDatBase],
        interfaces.ObjectLibrary[interfaces.OpDatBase],
        interfaces.ObjectLibrary[interfaces.RxnDatBase],
    ]:
        mol_lib: interfaces.ObjectLibrary[
            interfaces.MolDatBase
        ] = self._Mol_Lib()
        op_lib: interfaces.ObjectLibrary[interfaces.OpDatBase] = self._Op_Lib()
        rxn_lib: interfaces.ObjectLibrary[
            interfaces.RxnDatBase
        ] = self._Rxn_Lib()
        return (mol_lib, op_lib, rxn_lib)

    def CartesianStrategy(
        self,
        mol_lib: interfaces.ObjectLibrary[interfaces.MolDatBase],
        op_lib: interfaces.ObjectLibrary[interfaces.OpDatBase],
        rxn_lib: interfaces.ObjectLibrary[interfaces.RxnDatBase],
    ) -> strategies.CartesianStrategy:
        return self._CartesianStrategy(
            mol_lib=mol_lib, op_lib=op_lib, rxn_lib=rxn_lib, engine=self
        )

    @property
    def mol(self):
        return interfaces.MoleculeTypes(self._Mol)

    @property
    def op(self):
        return interfaces.OperatorTypes(_rdkit_op_init(self, self._Op))

    @property
    def strat(self):
        return interfaces.StrategyTypes(
            _cartesian_op_init(self),
            strategies.PriorityQueueStrategyBasic,
        )

    @property
    def filter(self):
        return interfaces.FilterTypes(
            interfaces.MolFilterTypes(
                filters.MolFilterMetaVal,
                filters.MolFilterMetaExist,
                filters.MolFilterMetaFunc,
            ),
            interfaces.BundleFilterTypes(),
            interfaces.RecipeFilterTypes(filters.CoreactantFilter),
            interfaces.ReactionFilterTypes(filters.GenerationFilter),
        )

    @property
    def meta(self):
        return interfaces.MetaCalcTypes(metacalc.GenerationCalculator)

    def new_network(self):
        return network.ChemNetworkBasic()

    def network_from_file(
        self,
        filename: str,
        path: str = "./",
        ext: str = ".pgnet",
    ) -> interfaces.ChemNetwork:
        filepath = path + filename + ext
        with gzip.open(filepath, "r") as fin:
            md = xml.dom.minidom
            data = md.parse(fin).documentElement
            version = int(data.getAttribute("version"))
            if version == 0:
                subversion = int(data.getAttribute("subversion"))
                if subversion == 0:
                    bvals = base64.urlsafe_b64decode(data.firstChild.data)
                    network: interfaces.ChemNetwork = pickle.loads(bvals)
                    return network
                else:
                    raise NotImplementedError(
                        f"File at {filepath} is incompatible with this version of Pickaxe_Generic, please update (file version={version}.{subversion}, max supported=0.0)"
                    )
            else:
                raise NotImplementedError(
                    f"File at {filepath} is incompatible with this version of Pickaxe_Generic, please update (file version={version}, max supported=0)"
                )
