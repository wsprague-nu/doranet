"""
Contains classes which define and implement network expansion strategies.

Classes:

    ExpansionStrategy
      CartesianStrategy*
      CartesianStrategyParallel*
"""

from abc import ABC, abstractmethod
from concurrent.futures import ProcessPoolExecutor
from itertools import chain, islice
from itertools import product as iterproduct
from typing import (
    Callable,
    Collection,
    Generator,
    Hashable,
    Iterable,
    Optional,
    Protocol,
    Sequence,
    TypeVar,
    Union,
    final,
)

from pickaxe_generic.containers import ObjectLibrary
from pickaxe_generic.datatypes import (
    Identifier,
    MolDatBase,
    OpDatBase,
    RxnDatBase,
)
from pickaxe_generic.filters import (
    DefaultMetaDataUpdate,
    MetaDataCalculatorLocal,
    MetaDataUpdate,
    MetaKeyPacket,
    MolFilter,
    MolFilterMetaVal,
    RankValue,
    RecipeFilter,
    RecipeRanker,
)
from pickaxe_generic.network import ChemNetwork, Recipe, _MolIndex


class _ReactionProvider(Protocol):
    def Rxn(
        self,
        operator: Identifier = None,
        reactants: Collection[Identifier] = None,
        products: Collection[Identifier] = None,
    ) -> RxnDatBase:
        ...


def _generate_recipes_from_compat_table(
    compat: dict[Identifier, list[list[Identifier]]],
    known_cache: Optional[
        set[tuple[Identifier, tuple[Identifier, ...]]]
    ] = None,
    uid_prefilter: Optional[
        Callable[[Identifier, Sequence[Identifier]], bool]
    ] = None,
    add_to_cache: bool = True,
) -> Generator[tuple[Identifier, tuple[Identifier, ...]], None, None]:
    gen_cache = set()
    for op_uid, compat_list in compat.items():
        reactants: tuple[Identifier, ...]
        for reactants in iterproduct(*compat_list):
            recipe = (op_uid, reactants)
            if known_cache is not None:
                if recipe in known_cache:
                    continue
                else:
                    if add_to_cache:
                        known_cache.add(recipe)
            if uid_prefilter is not None and not uid_prefilter(
                op_uid, reactants
            ):
                continue
            recipe = (op_uid, reactants)
            if recipe not in gen_cache:
                gen_cache.add(recipe)
                yield (op_uid, reactants)
            else:
                continue


def _evaluate_reaction(
    operator: OpDatBase,
    reactants: Sequence[MolDatBase],
    engine: _ReactionProvider,
    rxn_filter: Optional[
        Callable[[OpDatBase, Sequence[MolDatBase], Sequence[MolDatBase]], bool]
    ] = None,
    return_rejects: bool = False,
) -> tuple[
    tuple[tuple[RxnDatBase, tuple[MolDatBase, ...]], ...],
    tuple[tuple[RxnDatBase, tuple[MolDatBase, ...]], ...],
]:
    resultslist: list[tuple[RxnDatBase, tuple[MolDatBase, ...]]] = []
    rejectslist: list[tuple[RxnDatBase, tuple[MolDatBase, ...]]] = []
    for productset in operator(reactants):
        productset = tuple(productset)
        if rxn_filter is not None and not rxn_filter(
            operator, reactants, productset
        ):
            if not return_rejects:
                continue
            reactant_uids = tuple((mol.uid for mol in reactants))
            product_uids = tuple((mol.uid for mol in productset))
            reaction = engine.Rxn(operator.uid, reactant_uids, product_uids)
            rejectslist.append((reaction, productset))
            continue
        reactant_uids = tuple((mol.uid for mol in reactants))
        product_uids = tuple((mol.uid for mol in productset))
        reaction = engine.Rxn(operator.uid, reactant_uids, product_uids)
        resultslist.append((reaction, productset))
    return tuple(resultslist), tuple(rejectslist)


def _evaluate_reaction_unpack(
    job: tuple[
        OpDatBase,
        Sequence[MolDatBase],
        _ReactionProvider,
        Optional[
            Callable[
                [OpDatBase, Sequence[MolDatBase], Sequence[MolDatBase]], bool
            ]
        ],
    ]
) -> tuple[
    tuple[tuple[RxnDatBase, tuple[MolDatBase, ...]], ...],
    tuple[tuple[RxnDatBase, tuple[MolDatBase, ...]], ...],
]:

    return _evaluate_reaction(*job)


T = TypeVar("T")


def _chunk_generator(
    n: int, iterable: Iterable[T]
) -> Generator[Iterable[T], None, None]:
    it = iter(iterable)
    while True:
        chunk_it = islice(it, n)
        try:
            first_el = next(chunk_it)
        except StopIteration:
            return
        yield chain((first_el,), chunk_it)


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
        custom_filter: Optional[
            Callable[
                [OpDatBase, Sequence[MolDatBase], Sequence[MolDatBase]], bool
            ]
        ] = None,
        custom_uid_prefilter: Optional[
            Callable[[Identifier, Sequence[Identifier]], bool]
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

    @abstractmethod
    def refresh(self) -> None:
        """
        Refresh active molecules and operators from attached libraries.
        """


@final
class CartesianStrategy(ExpansionStrategy):
    """
    Implements ExpansionStrategy interface via Cartesian product of molecules
    and operators.
    """

    def __init__(
        self,
        mol_lib: ObjectLibrary[MolDatBase],
        op_lib: ObjectLibrary[OpDatBase],
        rxn_lib: ObjectLibrary[RxnDatBase],
        engine: _ReactionProvider,
        blacklist: Optional[set[Identifier]] = None,
    ) -> None:
        """Initialize strategy by attaching libraries.

        Parameters
        ----------
        mol_lib : ObjectLibrary[MolDatBase]
            Library containing molecules.
        op_lib : ObjectLibrary[OpDatBase]
            Library containing operators.
        rxn_lib : ObjectLibrary[RxnDatBase]
            Library containing reactions.
        engine : _ReactionProvider
            Object used to initialize reactions.
        """
        self._mol_lib = mol_lib
        self._op_lib = op_lib
        self._rxn_lib = rxn_lib
        self._engine = engine
        if blacklist is None:
            blacklist = set()
        self._blacklist = blacklist

        # stores active molecule IDs
        self._active_mols: set[Identifier] = set(
            (uid for uid in self._mol_lib.ids())
        )
        # stores active operator IDs
        self._active_ops: set[Identifier] = set(
            (uid for uid in self._op_lib.ids())
        )
        # stores combinations of reagents and operators which have been tried
        self._recipe_cache: set[
            tuple[Identifier, tuple[Identifier, ...]]
        ] = set()

        # create operator-molecule compatibility table
        self._compat_table: dict[Identifier, list[list[Identifier]]] = {}
        for op_uid in self._active_ops:
            self._add_op_to_compat(op_uid)

        # fill operator-molecule compatibility table
        for mol_uid in self._active_mols:
            self._add_mol_to_compat(mol_uid)

    def _add_mol_to_compat(self, mol: Union[Identifier, MolDatBase]) -> None:
        """
        Add entries to compat_table for new molecule.

        Parameters
        ----------
        mol : Identifier | MolDatBase
            Molecule to be added to compat_table.
        """
        if isinstance(mol, MolDatBase):
            mol_uid = mol.uid
        else:
            mol_uid = mol
            mol = self._mol_lib[mol]

        for op_uid in self._active_ops:
            op = self._op_lib[op_uid]
            for arg in range(len(self._compat_table[op_uid])):
                if op.compat(mol, arg):
                    self._compat_table[op_uid][arg].append(mol_uid)

    def _add_op_to_compat(self, op: Union[Identifier, OpDatBase]) -> None:
        """
        Add entries to compat_table for new operator.

        Parameters
        ----------
        op : Identifier | OpDatBase
            Operator to be added to compat_table.
        """
        if isinstance(op, OpDatBase):
            op_uid = op.uid
        else:
            op_uid = op
            op = self._op_lib[op]

        optable: list[list[Identifier]] = [[] for _ in range(len(op))]
        for arg in range(len(optable)):
            for mol_uid in self._active_mols:
                mol = self._mol_lib[mol_uid]
                if op.compat(mol, arg):
                    optable[arg].append(mol_uid)
        self._compat_table[op_uid] = optable

    @property
    def blacklist(self) -> set[Identifier]:
        return self._blacklist

    @blacklist.setter
    def blacklist(self, value: set[Identifier]) -> None:
        self._blacklist = value

    def reset_recipe_cache(self) -> None:
        self._recipe_cache = set()

    def refresh(self) -> None:
        if len(self._mol_lib) > len(self._active_mols):
            for mol_uid in self._mol_lib.ids():
                if (
                    mol_uid not in self._active_mols
                    and mol_uid not in self._blacklist
                ):
                    self._active_mols.add(mol_uid)
                    self._add_mol_to_compat(mol_uid)

        if len(self._op_lib) > len(self._active_ops):
            for op_uid in self._op_lib.ids():
                if op_uid not in self._active_ops:
                    self._active_ops.add(op_uid)
                    self._add_op_to_compat(op_uid)

    def expand(
        self,
        max_rxns: Optional[int] = None,
        max_mols: Optional[int] = None,
        num_gens: Optional[int] = None,
        custom_filter: Optional[
            Callable[
                [OpDatBase, Sequence[MolDatBase], Sequence[MolDatBase]], bool
            ]
        ] = None,
        custom_uid_prefilter: Optional[
            Callable[[Identifier, Sequence[Identifier]], bool]
        ] = None,
        retain_products_to_blacklist: bool = False,
    ) -> None:
        # value used to tell if any new reactions have occurred in a generation
        exhausted: bool = False
        num_mols: int = 0
        num_rxns: int = 0
        gen: int = 0

        while not exhausted:
            if num_gens is not None and gen >= num_gens:
                return
            exhausted = True

            # iterate through possible reactant combinations
            for recipe in _generate_recipes_from_compat_table(
                self._compat_table, self._recipe_cache, custom_uid_prefilter
            ):

                operator = self._op_lib[recipe[0]]
                reactants = tuple(
                    (self._mol_lib[mol_uid] for mol_uid in recipe[1])
                )

                # iterate through evaluated reactions
                results, rejects = _evaluate_reaction(
                    operator,
                    reactants,
                    self._engine,
                    custom_filter,
                    retain_products_to_blacklist,
                )
                for reaction, products in rejects:
                    num_mols += len(products)
                    num_rxns += 1
                    if max_mols is not None and num_mols > max_mols:
                        return
                    if max_rxns is not None and num_rxns > max_rxns:
                        return
                    for mol in products:
                        self._blacklist.add(mol.uid)
                        self._mol_lib.add(mol)
                    self._rxn_lib.add(reaction)
                    exhausted = False
                for reaction, products in results:
                    num_mols += len(products)
                    num_rxns += 1
                    if max_mols is not None and num_mols > max_mols:
                        return
                    if max_rxns is not None and num_rxns > max_rxns:
                        return
                    for mol in products:
                        self._blacklist.discard(mol.uid)
                        self._mol_lib.add(mol)
                    self._rxn_lib.add(reaction)
                    exhausted = False
            gen += 1
            self.refresh()


@final
class CartesianStrategyParallel(ExpansionStrategy):
    """
    Implements ExpansionStrategy interface via Cartesian product of molecules
    and operators.
    """

    def __init__(
        self,
        mol_lib: ObjectLibrary[MolDatBase],
        op_lib: ObjectLibrary[OpDatBase],
        rxn_lib: ObjectLibrary[RxnDatBase],
        engine: _ReactionProvider,
        num_procs: int,
    ) -> None:
        """Initialize strategy by attaching libraries.

        Parameters
        ----------
        mol_lib : ObjectLibrary[MolDatBase]
            Library containing molecules.
        op_lib : ObjectLibrary[OpDatBase]
            Library containing operators.
        rxn_lib : ObjectLibrary[RxnDatBase]
            Library containing reactions.
        engine : _ReactionProvider
            Object used to initialize reactions.
        """
        self._mol_lib = mol_lib
        self._op_lib = op_lib
        self._rxn_lib = rxn_lib
        self._engine = engine
        self._num_procs = num_procs

        # stores active molecule IDs
        self._active_mols: set[Identifier] = set(
            (uid for uid in self._mol_lib.ids())
        )
        # stores active operator IDs
        self._active_ops: set[Identifier] = set(
            (uid for uid in self._op_lib.ids())
        )
        # stores combinations of reagents and operators which have been tried
        self._recipe_cache: set[
            tuple[Identifier, tuple[Identifier, ...]]
        ] = set()

        # create operator-molecule compatibility table
        self._compat_table: dict[Identifier, list[list[Identifier]]] = {}
        for op_uid in self._active_ops:
            self._add_op_to_compat(op_uid)

        # fill operator-molecule compatibility table
        for mol_uid in self._active_mols:
            self._add_mol_to_compat(mol_uid)

    def _add_mol_to_compat(self, mol: Union[Identifier, MolDatBase]) -> None:
        """
        Add entries to compat_table for new molecule.

        Parameters
        ----------
        mol : Identifier | MolDatBase
            Molecule to be added to compat_table.
        """
        if isinstance(mol, MolDatBase):
            mol_uid = mol.uid
        else:
            mol_uid = mol
            mol = self._mol_lib[mol]

        for op_uid in self._active_ops:
            op = self._op_lib[op_uid]
            for arg in range(len(self._compat_table[op_uid])):
                if op.compat(mol, arg):
                    self._compat_table[op_uid][arg].append(mol_uid)

    def _add_op_to_compat(self, op: Union[Identifier, OpDatBase]) -> None:
        """
        Add entries to compat_table for new operator.

        Parameters
        ----------
        op : Identifier | OpDatBase
            Operator to be added to compat_table.
        """
        if isinstance(op, OpDatBase):
            op_uid = op.uid
        else:
            op_uid = op
            op = self._op_lib[op]

        optable: list[list[Identifier]] = [[] for _ in range(len(op))]
        for arg in range(len(optable)):
            for mol_uid in self._active_mols:
                mol = self._mol_lib[mol_uid]
                if op.compat(mol, arg):
                    optable[arg].append(mol_uid)
        self._compat_table[op_uid] = optable

    def reset_recipe_cache(self) -> None:
        self._recipe_cache = set()

    def refresh(self) -> None:
        if len(self._mol_lib) > len(self._active_mols):
            for mol_uid in self._mol_lib.ids():
                if mol_uid not in self._active_mols:
                    self._active_mols.add(mol_uid)
                    self._add_mol_to_compat(mol_uid)

        if len(self._op_lib) > len(self._active_ops):
            for op_uid in self._op_lib.ids():
                if op_uid not in self._active_ops:
                    self._active_ops.add(op_uid)
                    self._add_op_to_compat(op_uid)

    def expand(
        self,
        max_rxns: Optional[int] = None,
        max_mols: Optional[int] = None,
        num_gens: Optional[int] = None,
        custom_filter: Optional[
            Callable[
                [OpDatBase, Sequence[MolDatBase], Sequence[MolDatBase]], bool
            ]
        ] = None,
        custom_uid_prefilter: Optional[
            Callable[[Identifier, Sequence[Identifier]], bool]
        ] = None,
    ) -> None:
        # value used to tell if any new reactions have occurred in a generation
        exhausted: bool = False
        num_mols: int = 0
        num_rxns: int = 0
        gen: int = 0

        while not exhausted:
            if num_gens is not None and gen >= num_gens:
                return
            exhausted = True

            job_generator = _chunk_generator(
                self._num_procs * 100,
                (
                    (
                        self._op_lib[recipe[0]],
                        [self._mol_lib[mol_uid] for mol_uid in recipe[1]],
                        self._engine,
                        custom_filter,
                    )
                    for recipe in _generate_recipes_from_compat_table(
                        self._compat_table,
                        self._recipe_cache,
                        custom_uid_prefilter,
                    )
                ),
            )

            # iterate through possible reactant combinations
            with ProcessPoolExecutor(max_workers=self._num_procs) as executor:
                for jobchunk in job_generator:
                    jobchunk = tuple(jobchunk)
                    for results, _ in executor.map(
                        _evaluate_reaction_unpack, jobchunk
                    ):
                        #                    for result in _evaluate_reaction_unpack(jobchunk):
                        for reaction, products in results:
                            num_mols += len(products)
                            num_rxns += 1
                            if max_mols is not None and num_mols > max_mols:
                                return
                            if max_rxns is not None and num_rxns > max_rxns:
                                return
                            for mol in products:
                                self._mol_lib.add(mol)
                            self._rxn_lib.add(reaction)
                            exhausted = False
            gen += 1
            self.refresh()


class PriorityQueueStrategy(ABC):
    __slots__ = ()

    @abstractmethod
    def __init__(
        self,
        network: ChemNetwork,
        num_procs: Optional[int] = None,
        blacklist_key: Optional[str] = None,
    ) -> None:
        ...

    @abstractmethod
    def expand(
        self,
        max_recipes: Optional[int] = None,
        batch_size: int = 1,
        mol_filter_local: Optional[MolFilter] = None,
        mol_filter: Optional[MolFilter] = None,
        recipe_filter: Optional[RecipeFilter] = None,
        recipe_ranker: Optional[RecipeRanker] = None,
        mc_local: Optional[MetaDataCalculatorLocal] = None,
        mc_update: Optional[MetaDataUpdate] = DefaultMetaDataUpdate(),
    ) -> None:
        ...

    @abstractmethod
    def blacklist_key(self) -> str:
        ...


class PriorityQueueStrategyBasic(PriorityQueueStrategy):
    __slots__ = ("_network", "_blacklist_key", "_blacklist_func")

    def __init__(
        self,
        network: ChemNetwork,
        num_procs: Optional[int] = None,
        blacklist_key: Optional[str] = None,
    ) -> None:
        if num_procs is not None:
            raise NotImplementedError(
                f"Parallel processing not yet supported on {type(self)}"
            )
        self._network = network
        if blacklist_key is not None:
            self._blacklist_key = blacklist_key
        else:
            self._blacklist_key = f"_blacklist_{id(self)}"
        self._blacklist_func = ~MolFilterMetaVal(self._blacklist_key, True)

    def _blacklist_mol(self, index: _MolIndex) -> None:
        self._network.mol_meta(index, self._blacklist_key, True)

    def _unblacklist_mol(self, index: _MolIndex) -> None:
        self._network.mol_meta(index, self._blacklist_key, False)

    def _is_blacklisted(self, index: _MolIndex) -> bool:
        return self._blacklist_func(self._network.mols[index])

    def blacklist_key(self) -> str:
        return self._blacklist_key

    def expand(
        self,
        max_recipes: Optional[int] = None,
        heap_size: int = 1,
        mol_filter_local: Optional[MolFilter] = None,
        mol_filter: Optional[MolFilter] = None,
        recipe_filter: Optional[RecipeFilter] = None,
        recipe_ranker: Optional[RecipeRanker] = None,
        mc_local: Optional[MetaDataCalculatorLocal] = None,
        mc_update: Optional[MetaDataUpdate] = DefaultMetaDataUpdate(),
    ) -> None:

        if mol_filter_local is None:
            mol_filter_local = self._blacklist_func
        else:
            mol_filter_local = mol_filter_local & self._blacklist_func

        important_key_set: MetaKeyPacket = MetaKeyPacket()
        if mol_filter is not None:
            important_key_set = important_key_set + mol_filter.meta_required
        important_key_set = important_key_set + mol_filter_local.meta_required
        if recipe_filter is not None:
            important_key_set = important_key_set + recipe_filter.meta_required
        if recipe_ranker is not None:
            important_key_set = important_key_set + recipe_ranker.meta_required
        if mc_local is not None:
            important_key_set = important_key_set + mc_local.meta_required

        exhausted: bool = False
