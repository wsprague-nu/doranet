"""
Contains classes which define and implement network expansion strategies.

Classes:

    ExpansionStrategy
      CartesianStrategy*
      CartesianStrategyParallel*
"""

import heapq
from abc import ABC, abstractmethod
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass
from itertools import chain, islice
from itertools import product as iterproduct
from itertools import repeat
from math import prod
from operator import floordiv
from typing import (
    Any,
    Callable,
    Collection,
    Generator,
    Hashable,
    Iterable,
    Mapping,
    Optional,
    Protocol,
    Sequence,
    TypeVar,
    Union,
    final,
)

from pickaxe_generic.containers import ObjectLibrary
from pickaxe_generic.datatypes import (
    DataPacket,
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
    MolFilterMetaExist,
    MolFilterMetaVal,
    RankValue,
    ReactionFilter,
    ReactionFilterBase,
    RecipeFilter,
    RecipeRanker,
    ReplaceBlacklist,
)
from pickaxe_generic.network import ChemNetwork, Recipe, _MolIndex, _OpIndex


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
        heap_size: int = 1,
        batch_size: Optional[int] = None,
        # mol_filter_local: Optional[MolFilter] = None,
        # mol_filter: Optional[MolFilter] = None,
        # recipe_filter: Optional[RecipeFilter] = None,
        recipe_ranker: Optional[RecipeRanker] = None,
        # mc_local: Optional[MetaDataCalculatorLocal] = None,
        # mc_update: Optional[MetaDataUpdate] = DefaultMetaDataUpdate(),
    ) -> None:
        ...


@dataclass(frozen=True)
class RecipeRankingJob:
    __slots__ = (
        "operator",
        "op_args",
        "recipe_ranker",
        "min_rankvalue",
    )

    operator: DataPacket[OpDatBase]
    op_args: tuple[tuple[DataPacket[MolDatBase], ...]]
    recipe_ranker: Optional[RecipeRanker]
    min_rankvalue: Optional[RankValue]


def calc_batch_split(
    size_bundle: Sequence[int], batch_size: int
) -> tuple[int, ...]:
    num_split = list(repeat(1, len(size_bundle)))
    split_size = tuple(size_bundle)
    while prod(split_size) > batch_size:
        max_index = max(range(len(split_size)), key=split_size.__getitem__)
        num_split[max_index] += 1
        split_size = tuple(
            -(num_mols // -splitnum)
            for num_mols, splitnum in zip(size_bundle, num_split)
        )
    return tuple(num_split)


def _generate_recipe_batches(
    mol_table: Sequence[Sequence[_MolIndex]],
    table_indices: Sequence[int],
    batch_size: Optional[int] = None,
    updated_mols: set[_MolIndex] = set(),
) -> Generator[tuple[tuple[_MolIndex, ...], ...], None, None]:
    num_args = len(mol_table)

    # get number of old and new mols for each argument
    # num_old: list[int] = []
    # num_new: list[int] = []
    # for mols, i_counter in zip(mol_table, table_indices):
    #     n_new_from_old = len(updated_mols.intersection(mols[:i_counter]))
    #     n_new = len(updated_mols) - i_counter + n_new_from_old
    #     n_old = i_counter - n_new_from_old
    #     num_old.append(n_old)
    #     num_new.append(n_new)
    updated_vals: tuple[tuple[_MolIndex, ...], ...] = tuple(
        tuple(updated_mols.intersection(mol_list)) for mol_list in mol_table
    )
    num_old = tuple(table_indices)
    num_new = tuple(
        len(col) + len(mol_list) - i_counter
        for col, mol_list, i_counter in zip(
            updated_vals, mol_table, table_indices
        )
    )
    num_tot = tuple(len(col) for col in mol_table)

    # for each argument, generate a bundle
    for i_bundle in range(num_args):
        # value_gens = (
        #     tuple(
        #         (i for i in mol_list[:i_counter] if i not in updated_mols)
        #         for mol_list, i_counter in zip(
        #             mol_table[:i_bundle], table_indices[:i_bundle]
        #         )
        #     )
        #     + ((i for i in mol_table[i_bundle]),)
        #     + tuple(
        #         (
        #             i
        #             for i in chain(
        #                 updated_mols.intersection(mol_list[:i_counter]),
        #                 mol_list[i_counter:],
        #             )
        #         )
        #         for mol_list, i_counter in zip(
        #             mol_table[i_bundle + 1 :], table_indices[i_bundle + 1 :]
        #         )
        #     )
        # )

        # # if no batch size, yield entire bundle as batch
        # if batch_size is None:
        #     yield tuple(tuple(mol_gen) for mol_gen in value_gens)
        #     return

        # # get bundle size
        # size_bundle = tuple(
        #     chain(
        #         num_old[:i_bundle],
        #         (num_new[i_bundle],),
        #         (
        #             n_new + n_old
        #             for n_new, n_old in zip(
        #                 num_new[i_bundle + 1 :], num_old[i_bundle + 1 :]
        #             )
        #         ),
        #     )
        # )
        if batch_size is None:
            yield tuple(
                tuple(mol_list[:i_counter])
                for mol_list, i_counter in zip(
                    mol_table[:i_bundle], table_indices
                )
            ) + (
                updated_vals[i_bundle]
                + tuple(mol_table[i_bundle][table_indices[i_bundle] :]),
            ) + tuple(
                tuple(mol_list) for mol_list in mol_table[i_bundle + 1 :]
            )
            return

        size_bundle = (
            num_old[:i_bundle] + (num_new[i_bundle],) + num_tot[i_bundle + 1 :]
        )
        if not all(size_bundle):
            continue

        batch_split = calc_batch_split(size_bundle, batch_size)
        chunk_sizes = tuple(
            -(num_mols // -splitnum)
            for num_mols, splitnum in zip(size_bundle, batch_split)
        )
        split_indices = (range(num_splits) for num_splits in batch_split)
        for batch_subindices in iterproduct(*split_indices):
            prev_column_mols = tuple(
                tuple(
                    mol_list[:i_counter][
                        subindex * chunk_size : (subindex + 1) * chunk_size
                    ]
                )
                for mol_list, i_counter, chunk_size, subindex in zip(
                    mol_table[:i_bundle],
                    table_indices[:i_bundle],
                    chunk_sizes[:i_bundle],
                    batch_subindices[:i_bundle],
                )
            )
            cur_column_mols = tuple(
                mol_table[i_bundle][
                    batch_subindices[i_bundle]
                    * chunk_sizes[i_bundle] : (batch_subindices[i_bundle] + 1)
                    * chunk_sizes[i_bundle]
                ]
            )
            next_column_mols = tuple(
                tuple(
                    mol_list[
                        subindex * chunk_size : (subindex + 1) * chunk_size
                    ]
                )
                for mol_list, chunk_size, subindex in zip(
                    mol_table[i_bundle + 1 :],
                    chunk_sizes[i_bundle + 1 :],
                    batch_subindices[i_bundle + 1 :],
                )
            )
            yield prev_column_mols + (cur_column_mols,) + next_column_mols


def assemble_recipe_batch_job(
    op_index: _OpIndex,
    batch: tuple[tuple[_MolIndex]],
    network: ChemNetwork,
    keyset: MetaKeyPacket,
) -> RecipeRankingJob:

    return RecipeRankingJob()


class PriorityQueueStrategyBasic(PriorityQueueStrategy):
    __slots__ = "_network"

    def __init__(
        self,
        network: ChemNetwork,
        num_procs: Optional[int] = None,
    ) -> None:
        if num_procs is not None:
            raise NotImplementedError(
                f"Parallel processing not yet supported on {type(self)}"
            )
        self._network = network

    def expand(
        self,
        max_recipes: Optional[int] = None,
        heap_size: Optional[int] = None,
        batch_size: Optional[int] = None,
        # mol_filter_local: Optional[MolFilter] = None,
        # mol_filter: Optional[MolFilter] = None,
        # recipe_filter: Optional[RecipeFilter] = None,
        recipe_ranker: Optional[RecipeRanker] = None,
        # mc_local: Optional[MetaDataCalculatorLocal] = None,
        # mc_update: Optional[MetaDataUpdate] = DefaultMetaDataUpdate(),
    ) -> None:

        # set keysets so that updated reactions may occur and parameters may be
        # passed to parallel processes
        mol_filter_local_keyset: MetaKeyPacket = MetaKeyPacket()
        recipe_keyset: MetaKeyPacket = MetaKeyPacket()
        reaction_keyset: MetaKeyPacket = MetaKeyPacket()
        if recipe_ranker is not None:
            recipe_keyset = recipe_keyset + recipe_ranker.meta_required
        # if mol_filter_local is not None:
        #     mol_filter_local_keyset = (
        #         mol_filter_local_keyset + mol_filter_local.meta_required
        #     )
        # if mol_filter is not None:
        #     recipe_keyset = recipe_keyset + mol_filter.meta_required
        # if recipe_filter is not None:
        #     recipe_keyset = recipe_keyset + recipe_filter.meta_required
        # if recipe_ranker is not None:
        #     recipe_keyset = recipe_keyset + recipe_ranker.meta_required
        # if mc_local is not None:
        #     reaction_keyset = mc_local.meta_required
        #     recipe_keyset = recipe_keyset + reaction_keyset
        # total_keyset = mol_filter_local_keyset + recipe_keyset + reaction_keyset

        # initialize loop variables
        network = self._network
        compat_indices_table = [
            [0 for _ in network.compat_table(i)]
            for i in range(len(network.ops))
        ]
        updated_mols_set: set[_MolIndex] = set()
        updated_ops_set: set[_OpIndex] = set()
        recipe_heap: list[tuple[Optional[RankValue], Recipe]] = []

        while (
            any(
                any(
                    i_tested < len(compat_mols)
                    for i_tested, compat_mols in zip(
                        op_index_table, network.compat_table(op_index)
                    )
                )
                for op_index, op_index_table in enumerate(compat_indices_table)
            )
            or len(updated_mols_set) != 0
            or len(updated_ops_set) != 0
        ):
            # raise error if operator metadata has been updated
            if len(updated_ops_set) != 0:
                raise NotImplementedError(
                    "Updating necessary operator metadata is not yet supported"
                )

            # for each operator, create recipe batches
            for opIndex, op in enumerate(network.ops):
                # for each argument, accumulate a total of old_mols and new_mols
                compat_table = network.compat_table(opIndex)
                compat_indices = compat_indices_table[opIndex]
                if not all(compat_table):
                    continue
                # generate recipe batches
                for batch in _generate_recipe_batches(
                    compat_table, compat_indices, batch_size, updated_mols_set
                ):
                    # assemble recipe ranking job
                    print(batch)

            # update compat_indices_table
            compat_indices_table = [
                [len(mol_list) for mol_list in network.compat_table(i)]
                for i in range(len(network.ops))
            ]

    def _add_recipe_to_heap(
        self,
        recipe_heap: list[tuple[Optional[RankValue], Recipe]],
        recipe: Recipe,
        value: Optional[RankValue],
        heap_size: Optional[int] = None,
    ) -> None:
        if heap_size is None or len(recipe_heap) < heap_size:
            heapq.heappush(recipe_heap, (value, recipe))
        elif value is None:
            return
        elif recipe_heap[0][0] is None or recipe_heap[0][0] < value:
            heapq.heappushpop(recipe_heap, (value, recipe))
