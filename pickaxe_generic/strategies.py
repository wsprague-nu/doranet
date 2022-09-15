"""
Contains classes which define and implement network expansion strategies.

Classes:

    ExpansionStrategy
      CartesianStrategy*
      CartesianStrategyParallel*
"""

import heapq
import operator
from abc import ABC, abstractmethod
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass, field
from functools import reduce
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
    Iterator,
    Mapping,
    Optional,
    Protocol,
    Sequence,
    TypeVar,
    Union,
    final,
    overload,
)

from pickaxe_generic.containers import ObjectLibrary
from pickaxe_generic.datatypes import (
    DataPacket,
    DataPacketE,
    Identifier,
    MolDatBase,
    OpDatBase,
    RxnDatBase,
)
from pickaxe_generic.filters import (
    MetaDataCalculatorLocal,
    MetaDataUpdate,
    MetaKeyPacket,
    MolFilter,
    MolFilterMetaExist,
    MolFilterMetaVal,
    RankValue,
    ReactionFilter,
    RecipeFilter,
    RecipeRanker,
    ReplaceBlacklist,
)
from pickaxe_generic.metadata import (
    LocalPropertyCalc,
    MetaDataResolverFunc,
    MetaUpdateResolver,
    PropertyCompositor,
    ReactionFilterBase,
    RxnAnalysisStep,
    as_rxn_analysis_step,
)
from pickaxe_generic.network import (
    ChemNetwork,
    Reaction,
    ReactionExplicit,
    Recipe,
    RecipeExplicit,
    _MolIndex,
    _OpIndex,
    recipe_from_explicit,
)


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
        beam_size: Optional[int] = 1,
        # mol_filter_local: Optional[MolFilter] = None,
        # mol_filter: Optional[MolFilter] = None,
        recipe_filter: Optional[RecipeFilter] = None,
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
        "recipe_filter",
        "recipe_ranker",
        "heap_size",
    )

    operator: DataPacket[OpDatBase]
    op_args: tuple[tuple[DataPacket[MolDatBase], ...], ...]
    recipe_filter: Optional[RecipeFilter]
    recipe_ranker: Optional[RecipeRanker]
    heap_size: Optional[int]


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
    batch: tuple[tuple[_MolIndex, ...], ...],
    network: ChemNetwork,
    keyset: MetaKeyPacket,
    recipe_ranker: Optional[RecipeRanker] = None,
    heap_size: Optional[int] = None,
    recipe_filter: Optional[RecipeFilter] = None,
) -> RecipeRankingJob:
    mol_data: Union[
        Generator[Generator[None, None, None], None, None],
        Generator[Generator[MolDatBase, None, None], None, None],
    ]
    mol_meta: Generator[
        Union[Sequence[Mapping], Generator[None, None, None]], None, None
    ]
    if keyset.live_operator:
        op_data = network.ops[op_index]
    else:
        op_data = None
    if keyset.operator_keys:
        op_meta = network.op_metas((op_index,), keyset.operator_keys)[0]
    else:
        op_meta = None
    op = DataPacket(op_index, op_data, op_meta)
    if keyset.live_molecule:
        mol_data = ((network.mols[i] for i in mol_list) for mol_list in batch)
    else:
        mol_data = ((None for _ in mol_list) for mol_list in batch)
    if keyset.molecule_keys:
        mol_meta = (
            network.mol_metas(mol_list, keyset.molecule_keys)
            for mol_list in batch
        )
    else:
        mol_meta = ((None for _ in mol_list) for mol_list in batch)
    mol_batch = tuple(
        tuple(
            DataPacket(i, mol, m_meta)
            for i, mol, m_meta in zip(i_col, mol_col, meta_col)
        )
        for i_col, mol_col, meta_col in zip(batch, mol_data, mol_meta)
    )
    return RecipeRankingJob(
        op, mol_batch, recipe_filter, recipe_ranker, heap_size
    )


@dataclass(frozen=True)
class ReactionJob:
    __slots__ = (
        "operator",
        "op_args",
    )

    operator: DataPacketE[OpDatBase]
    op_args: tuple[DataPacketE[MolDatBase], ...]


def assemble_reaction_job(
    recipe: Recipe, network: ChemNetwork, keyset: MetaKeyPacket
) -> ReactionJob:
    op_index = recipe.operator
    op_data = network.ops[op_index]

    if keyset.operator_keys:
        op_meta = network.op_metas((op_index,), keyset.operator_keys)[0]
    else:
        op_meta = None

    op = DataPacketE(op_index, op_data, op_meta)

    mol_meta: Iterable[Optional[Mapping]]
    if keyset.molecule_keys:
        mol_meta = network.mol_metas(recipe.reactants, keyset.molecule_keys)
    else:
        mol_meta = (None for _ in recipe.reactants)

    reactants = tuple(
        DataPacketE(i, network.mols[i], meta)
        for i, meta in zip(recipe.reactants, mol_meta)
    )

    return ReactionJob(op, reactants)


@dataclass(frozen=True)
class RecipePriorityItem:
    rank: Optional[RankValue]
    recipe: Recipe

    def __lt__(self, other: "RecipePriorityItem") -> bool:
        if other.rank is None:
            if self.rank is None:
                return self.recipe < other.recipe
            return False
        if self.rank is None:
            return True
        return self.rank < other.rank

    def __eq__(self, other: object) -> bool:
        if (
            isinstance(other, RecipePriorityItem)
            and self.rank == other.rank
            and self.recipe == other.recipe
        ):
            return True
        return False


def execute_recipe_ranking(
    job: RecipeRankingJob,
    min_val: Optional[RecipePriorityItem],
    recipes_tested: Collection[Recipe],
) -> "RecipeHeap":
    recipe_generator = (
        (RecipeExplicit(job.operator, reactants_data), recipe)
        for reactants_data, recipe in (
            (
                reactants_data,
                Recipe(
                    _OpIndex(job.operator.i),
                    tuple(reactant.i for reactant in reactants_data),
                ),
            )
            for reactants_data in iterproduct(*job.op_args)
        )
        if recipe not in recipes_tested
    )

    if job.recipe_ranker is None:
        if job.recipe_filter is None:
            return RecipeHeap.from_iter(
                (
                    recipe_item
                    for recipe_item in (
                        RecipePriorityItem(
                            None,
                            recipe,
                        )
                        for recipe in (
                            Recipe(
                                _OpIndex(job.operator.i),
                                tuple(reactant.i for reactant in reactants),
                            )
                            for reactants in iterproduct(*job.op_args)
                        )
                        if recipe not in recipes_tested
                    )
                    if min_val is None or not (recipe_item < min_val)
                ),
                maxsize=job.heap_size,
            )
        return RecipeHeap.from_iter(
            (
                recipe_item
                for recipe_item in (
                    RecipePriorityItem(None, recipe)
                    for recipe_explicit, recipe in recipe_generator
                    if job.recipe_filter(recipe_explicit)
                )
                if min_val is None or not (recipe_item < min_val)
            ),
            job.heap_size,
        )

    # ideally, pass the min_rank to recipe_ranker, but only if the heap is full;
    # strategy utilizing heap and parallel reduction is likely best
    # but difficult to implement
    if job.recipe_filter is None:
        rank_generator = (
            RecipePriorityItem(job.recipe_ranker(recipe_explicit), recipe)
            for recipe_explicit, recipe in recipe_generator
        )
    else:
        rank_generator = (
            RecipePriorityItem(job.recipe_ranker(recipe_explicit), recipe)
            for recipe_explicit, recipe in recipe_generator
            if job.recipe_filter(recipe_explicit)
        )
    priority_item_generator = (
        recipe_item
        for recipe_item in rank_generator
        if (min_val is None) or not (recipe_item < min_val)
    )
    return RecipeHeap.from_iter(priority_item_generator, maxsize=job.heap_size)


class RecipeHeap:
    __slots__ = ("_heap", "_maxsize", "_ordered")

    _heap: Optional[list[RecipePriorityItem]]
    _ordered: Optional[list[RecipePriorityItem]]
    _maxsize: Optional[int]

    def __init__(
        self,
        maxsize: Optional[int] = None,
        heaps: Optional[Collection["RecipeHeap"]] = None,
    ) -> None:
        if heaps is None:
            self._heap = None
            self._ordered = None
        elif maxsize is None:
            self._heap = None
            self._ordered = list(heapq.merge(*heaps))
        else:
            self._heap = None
            self._ordered = list(
                reversed(
                    tuple(
                        islice(
                            heapq.merge(
                                *(reversed(heap) for heap in heaps),
                                reverse=True,
                            ),
                            0,
                            maxsize,
                        )
                    )
                )
            )
        self._maxsize = maxsize

    @classmethod
    def from_iter(
        cls, data: Iterable[RecipePriorityItem], maxsize: Optional[int] = None
    ) -> "RecipeHeap":
        heap = RecipeHeap(maxsize)
        for item in data:
            heap.add_recipe(item)
        return heap

    @property
    def min(self) -> Optional[RecipePriorityItem]:
        if self._heap is None:
            if self._ordered is None or len(self._ordered) == 0:
                return None
            return min(self._ordered)
        return self._heap[0]

    def popvals(self, n: Optional[int]) -> tuple[RecipePriorityItem, ...]:
        if self._ordered is None:
            if self._heap is None:
                return tuple()
            self._ordered = sorted(self._heap)
        self._heap = None
        if n is None:
            vals = tuple(reversed(self._ordered))
            self._ordered = None
            return vals
        return tuple(
            self._ordered.pop() for _ in range(min(len(self._ordered), n))
        )

    @overload
    def __getitem__(self, item: slice) -> Sequence[RecipePriorityItem]:
        ...

    @overload
    def __getitem__(self, item: int) -> RecipePriorityItem:
        ...

    def __getitem__(self, item: Union[int, slice]):
        if self._ordered is None:
            if self._heap is None:
                raise IndexError(f"Invalid index {item=}; heap is empty")
            self._ordered = sorted(self._heap)
        return self._ordered[item]

    def __len__(self) -> int:
        if self._heap is None:
            if self._ordered is None:
                return 0
            return len(self._ordered)
        return len(self._heap)

    def __iter__(self) -> Iterator[RecipePriorityItem]:
        if self._ordered is None:
            if self._heap is None:
                return iter([])
            self._ordered = sorted(self._heap)
        return iter(self._ordered)

    def __reversed__(self) -> Iterator[RecipePriorityItem]:
        if self._ordered is None:
            if self._heap is None:
                return iter([])
            self._ordered = sorted(self._heap)
        return reversed(self._ordered)

    def add_recipe(self, recipe: RecipePriorityItem) -> None:
        if self._heap is None:
            if self._ordered is None:
                self._heap = []
            else:
                self._heap = list(self._ordered)
                heapq.heapify(self._heap)
                self._ordered = None
        if self._maxsize is None or len(self._heap) < self._maxsize:
            heapq.heappush(self._heap, recipe)
        elif self._heap[0] < recipe:
            heapq.heapreplace(self._heap, recipe)

    def __add__(self, other: "RecipeHeap") -> "RecipeHeap":
        if self._maxsize != other._maxsize:
            raise ValueError(
                f"Heap sizes do not match ({self._maxsize} != {other._maxsize})"
            )

        return RecipeHeap(self._maxsize, (self, other))


def execute_reaction(
    rxn_job: ReactionJob,
) -> Generator[tuple[ReactionExplicit, bool], None, None]:
    reactants = tuple(mol.item for mol in rxn_job.op_args)
    if rxn_job.operator.item is None or any(mol is None for mol in reactants):
        raise ValueError("ReactionJob has non-None item components!")
    product_packets = rxn_job.operator.item(reactants)  # type: ignore
    reactant_datapackets = tuple(
        DataPacketE(p.i, p.item, p.meta) for p in rxn_job.op_args
    )
    operator_datapacket = DataPacketE(
        rxn_job.operator.i, rxn_job.operator.item, rxn_job.operator.meta
    )
    for rxn_products in product_packets:
        product_datapackets = tuple(
            DataPacketE(-1, product, None) for product in rxn_products
        )
        rxn = ReactionExplicit(
            operator_datapacket,
            reactant_datapackets,
            product_datapackets,
            None,
        )
        yield rxn, True


def execute_reactions(
    rxn_jobs: Collection[ReactionJob],
    rxn_analysis: Optional[RxnAnalysisStep] = None,
) -> Iterable[tuple[ReactionExplicit, bool]]:
    rxn_generator = reduce(
        chain, (execute_reaction(rxn_job) for rxn_job in rxn_jobs)  # type: ignore
    )
    if rxn_analysis is None:
        return rxn_generator
    return rxn_analysis.execute(rxn_generator)


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
        beam_size: Optional[int] = 1,
        # mol_filter_local: Optional[MolFilter] = None,
        # mol_filter: Optional[MolFilter] = None,
        recipe_filter: Optional[RecipeFilter] = None,
        recipe_ranker: Optional[RecipeRanker] = None,
        mc_local: Optional[
            Union[
                RxnAnalysisStep,
                PropertyCompositor,
                ReactionFilterBase,
                LocalPropertyCalc,
            ]
        ] = None,
        mc_update: Optional[MetaUpdateResolver] = None,
        # mc_update: Optional[MetaDataUpdate] = DefaultMetaDataUpdate(),
    ) -> None:

        rxn_analysis_task: Optional[RxnAnalysisStep] = None
        if mc_local is not None:
            rxn_analysis_task = as_rxn_analysis_step(mc_local)
        if mc_update is None:
            mc_update = MetaUpdateResolver({}, {}, {})

        if heap_size is not None and beam_size is not None:
            ValueError(
                f"Heap size ({heap_size}) must be greater than beam size ({beam_size})"
            )

        # set keysets so that updated reactions may occur and parameters may be
        # passed to parallel processes
        recipe_keyset: MetaKeyPacket = MetaKeyPacket()
        reaction_keyset: MetaKeyPacket = MetaKeyPacket()
        if recipe_filter is not None:
            recipe_keyset = recipe_keyset + recipe_filter.meta_required
        if recipe_ranker is not None:
            recipe_keyset = recipe_keyset + recipe_ranker.meta_required
        if rxn_analysis_task is not None:
            reaction_keyset = reaction_keyset + rxn_analysis_task.meta_required
        total_keyset = recipe_keyset + reaction_keyset
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
        #    reaction_keyset = mc_local.meta_required
        # total_keyset = mol_filter_local_keyset + recipe_keyset + reaction_keyset

        # initialize loop variables
        network = self._network
        compat_indices_table = [
            [0 for _ in network.compat_table(i)]
            for i in range(len(network.ops))
        ]
        updated_mols_set: set[_MolIndex] = set()
        updated_ops_set: set[_OpIndex] = set()
        recipe_heap: RecipeHeap = RecipeHeap(maxsize=heap_size)
        recipes_tested: set[Recipe] = set()

        while (max_recipes is None or len(recipes_tested) < max_recipes) and (
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
            or len(recipe_heap) > 0
        ):
            for op_index in updated_ops_set:
                compat_indices_table[op_index] = [
                    0 for _ in range(len(compat_indices_table[op_index]))
                ]

            # for each operator, create recipe batches
            for opIndex, _ in enumerate(network.ops):
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
                    recipejob = assemble_recipe_batch_job(
                        _OpIndex(opIndex),
                        batch,
                        network,
                        recipe_keyset,
                        recipe_ranker,
                        heap_size,
                        recipe_filter,
                    )
                    cur_min = recipe_heap.min
                    new_recipe_heap = execute_recipe_ranking(
                        recipejob, cur_min, recipes_tested
                    )
                    recipe_heap = recipe_heap + new_recipe_heap

                continue

            # update compat_indices_table
            compat_indices_table = [
                [len(mol_list) for mol_list in network.compat_table(i)]
                for i in range(len(network.ops))
            ]

            # perform beam expansion
            recipes_to_be_expanded = recipe_heap.popvals(beam_size)

            if len(recipes_to_be_expanded) == 0:
                break

            reaction_jobs = tuple(
                assemble_reaction_job(
                    reciperank.recipe, network, reaction_keyset
                )
                for reciperank in recipes_to_be_expanded
            )

            # execute reactions
            for rxn, pass_filter in execute_reactions(
                reaction_jobs, rxn_analysis_task
            ):
                # add product mols to network
                products_indices = tuple(
                    network.add_mol(mol.item, None, pass_filter)
                    for mol in rxn.products
                    if mol.item is not None
                )

                # build reaction
                reactants_indices = tuple(
                    _MolIndex(mol.i) for mol in rxn.reactants
                )
                rxn_implicit = Reaction(
                    _OpIndex(rxn.operator.i),
                    reactants_indices,
                    products_indices,
                )

                # add reaction to network
                rxn_index = network.add_rxn(rxn_implicit)

                updated_mols_set = set()
                updated_ops_set = set()

                # update reactant metadata
                for m_dat in zip(reactants_indices, rxn.reactants):
                    if m_dat[1].meta is not None:
                        for key, value in m_dat[1].meta.items():
                            if key in mc_update.mol_updates:
                                cur_val = network.mol_meta(m_dat[0], key)
                                if cur_val is not None:
                                    value = mc_update.mol_updates[key](
                                        value, cur_val
                                    )
                            if network.mol_meta(m_dat[0], key) != value:
                                network.mol_meta(m_dat[0], key, value)
                                if key in total_keyset.molecule_keys:
                                    updated_mols_set.add(m_dat[0])

                # update product metadata
                for m_dat in zip(products_indices, rxn.products):
                    if m_dat[1].meta is not None:
                        for key, value in m_dat[1].meta.items():
                            if key in mc_update.mol_updates:
                                cur_val = network.mol_meta(m_dat[0], key)
                                if cur_val is not None:
                                    value = mc_update.mol_updates[key](
                                        value, cur_val
                                    )
                            if network.mol_meta(m_dat[0], key) != value:
                                network.mol_meta(m_dat[0], key, value)
                                if key in total_keyset.molecule_keys:
                                    updated_mols_set.add(m_dat[0])

                # update operator metadata
                if rxn.operator.meta is not None:
                    for key, value in rxn.operator.meta.items():
                        if key in mc_update.op_updates:
                            cur_val = network.op_meta(
                                rxn_implicit.operator, key
                            )
                            if cur_val is not None:
                                value = mc_update.op_updates[key](
                                    value,
                                    network.op_meta(rxn_implicit.operator, key),
                                )
                        if network.op_meta(rxn_implicit.operator, key) != value:
                            network.op_meta(rxn_implicit.operator, key, value)
                            if key in total_keyset.operator_keys:
                                updated_ops_set.add(rxn_implicit.operator)

                # update reaction metadata
                if rxn.reaction_meta is not None:
                    for key, value in rxn.reaction_meta.items():
                        if key in mc_update.rxn_updates:
                            cur_val = network.rxn_meta(rxn_index, key)
                            if cur_val is not None:
                                value = mc_update.rxn_updates[key](
                                    value, cur_val
                                )
                        if network.rxn_meta(rxn_index, key) != value:
                            network.rxn_meta(rxn_index, key, value)

            recipes_tested.update(
                (reciperank.recipe for reciperank in recipes_to_be_expanded)
            )

            if len(recipe_heap) == 0:
                compat_indices_table = [
                    [0 for _ in network.compat_table(i)]
                    for i in range(len(network.ops))
                ]

            continue
