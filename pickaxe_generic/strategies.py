"""
Contains classes which define and implement network expansion strategies.

Classes:

    ExpansionStrategy
      CartesianStrategy*
      HybridExpansionStrategy
        OrderedCartesianHybridExpansionStrategy*
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from itertools import product as iterproduct
from typing import (
    TYPE_CHECKING,
    Callable,
    Dict,
    FrozenSet,
    Generator,
    List,
    Optional,
    Sequence,
    Set,
    Tuple,
    final,
)

from rdkit.Chem.rdchem import Mol as RDKitMol

from pickaxe_generic.containers import ObjectLibrary
from pickaxe_generic.datatypes import (
    Identifier,
    MolDatBase,
    OpDatBase,
    RxnDatBase,
)
from pickaxe_generic.filters import AlwaysTrueFilter

if TYPE_CHECKING:
    from pickaxe_generic.engine import NetworkEngine


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
        num_gens : Optional[int] (default: None)
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
        filter: Callable[
            [OpDatBase, Sequence[MolDatBase], Sequence[MolDatBase]], bool
        ] = AlwaysTrueFilter(),
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
                        if not filter(op, reactants, tuple(productset)):
                            continue
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
        # check for molecules in mol_lib which are not in _mol_cache and add
        # them
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
        Tuple[Tuple[Identifier, FrozenSet[RDKitMol]], int], None, None
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
        filter: Callable[[MolDatBase], bool] = lambda _: True,
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
                        self._mol_cache[mol.uid] = (mol, lib_index)
                        self._add_mol_to_compat(mol, lib_index)
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
