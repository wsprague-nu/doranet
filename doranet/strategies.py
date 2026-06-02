"""Contains classes which define and implement network expansion strategies."""

import collections
import collections.abc
import concurrent.futures
import dataclasses
import functools
import heapq
import itertools
import math
import operator
import typing
from collections.abc import Collection, Iterable, Sequence
from typing import TYPE_CHECKING, Generic, TypeVar, final

if TYPE_CHECKING:
    from _typeshed import SupportsDunderLT

from doranet import interfaces, metadata
from doranet import network as pgnetworks

_T_comp = TypeVar("_T_comp", bound="SupportsDunderLT")


def _generate_recipes_from_compat_table(
    compat: dict[interfaces.Identifier, list[list[interfaces.Identifier]]],
    known_cache: typing.Optional[
        set[tuple[interfaces.Identifier, tuple[interfaces.Identifier, ...]]]
    ] = None,
    uid_prefilter: typing.Optional[
        collections.abc.Callable[
            [
                interfaces.Identifier,
                collections.abc.Sequence[interfaces.Identifier],
            ],
            bool,
        ]
    ] = None,
    add_to_cache: bool = True,
) -> collections.abc.Generator[
    tuple[interfaces.Identifier, tuple[interfaces.Identifier, ...]], None, None
]:
    gen_cache = set()
    for op_uid, compat_list in compat.items():
        reactants: tuple[interfaces.Identifier, ...]
        for reactants in itertools.product(*compat_list):
            recipe = (op_uid, reactants)
            if known_cache is not None:
                if recipe in known_cache:
                    continue
                elif add_to_cache:
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


T = typing.TypeVar("T")


def _chunk_generator(
    n: int, iterable: collections.abc.Iterable[T]
) -> collections.abc.Generator[collections.abc.Iterable[T], None, None]:
    it = iter(iterable)
    while True:
        chunk_it = itertools.islice(it, n)
        try:
            first_el = next(chunk_it)
        except StopIteration:
            return
        yield itertools.chain((first_el,), chunk_it)


@dataclasses.dataclass(frozen=True, slots=True)
class RecipeRankingJob:
    operator: interfaces.DataPacket[interfaces.OpDatBase]
    op_args: tuple[
        tuple[interfaces.DataPacket[interfaces.MolDatBase], ...], ...
    ]
    mol_filter: typing.Optional[interfaces.MolFilter]
    bundle_filter: typing.Optional[interfaces.BundleFilter]
    recipe_filter: typing.Optional[interfaces.RecipeFilter]
    recipe_ranker: typing.Optional[interfaces.RecipeRanker]
    heap_size: typing.Optional[int]


def calc_batch_split(
    size_bundle: collections.abc.Sequence[int], batch_size: int
) -> tuple[int, ...]:
    num_split = list(itertools.repeat(1, len(size_bundle)))
    split_size = tuple(size_bundle)
    while math.prod(split_size) > batch_size:
        max_index = max(range(len(split_size)), key=split_size.__getitem__)
        num_split[max_index] += 1
        split_size = tuple(
            -(num_mols // -splitnum)
            for num_mols, splitnum in zip(size_bundle, num_split, strict=False)
        )
    return tuple(num_split)


def _generate_recipe_batches(
    mol_table: collections.abc.Sequence[
        collections.abc.Sequence[interfaces.MolIndex]
    ],
    table_indices: collections.abc.Sequence[int],
    batch_size: typing.Optional[int] = None,
    updated_mols: typing.Optional[set[interfaces.MolIndex]] = None,
) -> collections.abc.Generator[
    tuple[tuple[interfaces.MolIndex, ...], ...], None, None
]:
    if updated_mols is None:
        updated_mols = set()
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
    updated_vals: tuple[tuple[interfaces.MolIndex, ...], ...] = tuple(
        tuple(updated_mols.intersection(mol_list)) for mol_list in mol_table
    )
    num_old = tuple(table_indices)
    num_new = tuple(
        len(col) + len(mol_list) - i_counter
        for col, mol_list, i_counter in zip(
            updated_vals, mol_table, table_indices, strict=False
        )
    )
    num_tot = tuple(len(col) for col in mol_table)

    # for each argument, generate a bundle
    for i_bundle in range(num_args):
        if batch_size is None:
            yield (
                tuple(
                    tuple(mol_list[:i_counter])
                    for mol_list, i_counter in zip(
                        mol_table[:i_bundle], table_indices, strict=False
                    )
                )
                + (
                    updated_vals[i_bundle]
                    + tuple(mol_table[i_bundle][table_indices[i_bundle] :]),
                )
                + tuple(
                    tuple(mol_list) for mol_list in mol_table[i_bundle + 1 :]
                )
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
            for num_mols, splitnum in zip(
                size_bundle, batch_split, strict=False
            )
        )
        split_indices = (range(num_splits) for num_splits in batch_split)
        for batch_subindices in itertools.product(*split_indices):
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
                    strict=False,
                )
            )
            cur_column_mols = tuple(
                mol_table[i_bundle][
                    batch_subindices[i_bundle] * chunk_sizes[i_bundle] : (
                        batch_subindices[i_bundle] + 1
                    )
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
                    strict=False,
                )
            )
            yield prev_column_mols + (cur_column_mols,) + next_column_mols


def assemble_recipe_batch_job(
    op_index: interfaces.OpIndex,
    batch: tuple[tuple[interfaces.MolIndex, ...], ...],
    network: interfaces.ChemNetwork,
    keyset: interfaces.MetaKeyPacket,
    recipe_ranker: typing.Optional[interfaces.RecipeRanker] = None,
    heap_size: typing.Optional[int] = None,
    recipe_filter: typing.Optional[interfaces.RecipeFilter] = None,
    mol_filter: typing.Optional[interfaces.MolFilter] = None,
    bundle_filter: typing.Optional[interfaces.BundleFilter] = None,
) -> RecipeRankingJob:
    mol_data: typing.Union[
        collections.abc.Generator[
            collections.abc.Generator[None, None, None], None, None
        ],
        collections.abc.Generator[
            collections.abc.Generator[interfaces.MolDatBase, None, None],
            None,
            None,
        ],
    ]
    mol_meta: collections.abc.Generator[
        typing.Union[
            collections.abc.Sequence[collections.abc.Mapping],
            collections.abc.Generator[None, None, None],
        ],
        None,
        None,
    ]
    op_data = network.ops[op_index] if keyset.live_operator else None
    if keyset.operator_keys:
        op_meta = network.ops.meta(op_index, keyset.operator_keys)
    else:
        op_meta = None
    op = interfaces.DataPacket(op_index, op_data, op_meta)
    if any(len(b) == 0 for b in batch):
        return RecipeRankingJob(
            operator=op,
            op_args=tuple(() for _ in batch),
            mol_filter=mol_filter,
            bundle_filter=bundle_filter,
            recipe_filter=recipe_filter,
            recipe_ranker=recipe_ranker,
            heap_size=heap_size,
        )
    if keyset.live_molecule:
        mol_data = ((network.mols[i] for i in mol_list) for mol_list in batch)
    else:
        mol_data = ((None for _ in mol_list) for mol_list in batch)
    if keyset.molecule_keys:
        mol_meta = (
            tuple(network.mols.meta(mol_list, keyset.molecule_keys))
            for mol_list in batch
        )
    else:
        mol_meta = ((None for _ in mol_list) for mol_list in batch)
    mol_batch = tuple(
        tuple(
            interfaces.DataPacket(i, mol, m_meta)
            for i, mol, m_meta in zip(i_col, mol_col, meta_col, strict=False)
        )
        for i_col, mol_col, meta_col in zip(
            batch, mol_data, mol_meta, strict=False
        )
    )
    return RecipeRankingJob(
        op,
        mol_batch,
        mol_filter,
        bundle_filter,
        recipe_filter,
        recipe_ranker,
        heap_size,
    )


@dataclasses.dataclass(frozen=True, slots=True)
class ReactionJob:
    operator: interfaces.DataPacketE[interfaces.OpDatBase]
    op_args: tuple[interfaces.DataPacketE[interfaces.MolDatBase], ...]


def assemble_reaction_job(
    recipe: interfaces.Recipe,
    network: interfaces.ChemNetwork,
    keyset: interfaces.MetaKeyPacket,
) -> ReactionJob:
    op_index = recipe.operator
    op_data = network.ops[op_index]

    if keyset.operator_keys:
        op_meta = network.ops.meta(op_index, keyset.operator_keys)
    else:
        op_meta = None

    op = interfaces.DataPacketE(op_index, op_data, op_meta)

    mol_meta: collections.abc.Iterable[typing.Optional[collections.abc.Mapping]]
    if keyset.molecule_keys:
        mol_meta = network.mols.meta(recipe.reactants, keyset.molecule_keys)
    else:
        mol_meta = (None for _ in recipe.reactants)

    reactants = tuple(
        interfaces.DataPacketE(i, network.mols[i], meta)
        for i, meta in zip(recipe.reactants, mol_meta, strict=False)
    )

    return ReactionJob(op, reactants)


@dataclasses.dataclass(frozen=True, slots=True)
class RecipePriorityItem:
    rank: typing.Optional[interfaces.RankValue]
    recipe: interfaces.Recipe

    def __lt__(self, other: "RecipePriorityItem") -> bool:
        if other.rank is None:
            if self.rank is None:
                return self.recipe < other.recipe
            return False
        if self.rank is None:
            return True
        return self.rank < other.rank

    def __eq__(self, other: object) -> bool:
        return (
            isinstance(other, RecipePriorityItem)
            and self.rank == other.rank
            and self.recipe == other.recipe
        )

    def __hash__(self) -> int:
        return hash(dataclasses.astuple(self))


def execute_recipe_ranking(
    job: RecipeRankingJob,
    min_val: typing.Optional[RecipePriorityItem],
    recipes_tested: collections.abc.Collection[interfaces.Recipe],
) -> "RecipeHeap":
    # filter molecules
    args_edited: tuple[
        tuple[interfaces.DataPacket[interfaces.MolDatBase], ...], ...
    ]
    if job.mol_filter is not None:
        args_edited = tuple(
            tuple(
                mol for mol in arg_mols if job.mol_filter(mol, job.operator, i)
            )
            for i, arg_mols in enumerate(job.op_args)
        )
    else:
        args_edited = job.op_args

    # perform bundle filtering
    bundles: collections.abc.Iterable[interfaces.RecipeBundle]
    if job.bundle_filter is not None:
        bundles = job.bundle_filter(
            interfaces.RecipeBundle(job.operator, args_edited)
        )
    else:
        bundles = (interfaces.RecipeBundle(job.operator, args_edited),)

    recipe_generator = (
        (interfaces.RecipeExplicit(job.operator, reactants_data), recipe)
        for reactants_data, recipe in itertools.chain.from_iterable(
            (
                (
                    (
                        reactants_data,
                        interfaces.Recipe(
                            interfaces.OpIndex(job.operator.i),
                            tuple(
                                interfaces.MolIndex(reactant.i)
                                for reactant in reactants_data
                            ),
                        ),
                    )
                    for reactants_data in itertools.product(*bundle.args)
                )
                for bundle in bundles
            )
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
                        for _, recipe in recipe_generator
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


@dataclasses.dataclass(slots=True)
class RecipeHeap:
    _heap: typing.Optional[list[RecipePriorityItem]]
    _ordered: typing.Optional[list[RecipePriorityItem]]
    _maxsize: typing.Optional[int]

    def __init__(
        self,
        maxsize: typing.Optional[int] = None,
        heaps: typing.Optional[collections.abc.Collection["RecipeHeap"]] = None,
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
                        itertools.islice(
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
        cls,
        data: collections.abc.Iterable[RecipePriorityItem],
        maxsize: typing.Optional[int] = None,
    ) -> "RecipeHeap":
        heap = RecipeHeap(maxsize)
        for item in data:
            heap.add_recipe(item)
        return heap

    @property
    def min(self) -> typing.Optional[RecipePriorityItem]:
        if self._heap is None:
            if self._ordered is None or len(self._ordered) == 0:
                return None
            return min(self._ordered)
        return self._heap[0]

    def popvals(
        self, n: typing.Optional[int]
    ) -> tuple[RecipePriorityItem, ...]:
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

    @typing.overload
    def __getitem__(self, item: slice) -> list[RecipePriorityItem]: ...

    @typing.overload
    def __getitem__(self, item: int) -> RecipePriorityItem: ...

    def __getitem__(self, item: typing.Union[int, slice]):
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

    def __iter__(self) -> collections.abc.Iterator[RecipePriorityItem]:
        if self._ordered is None:
            if self._heap is None:
                return iter([])
            self._ordered = sorted(self._heap)
        return iter(self._ordered)

    def __reversed__(self) -> collections.abc.Iterator[RecipePriorityItem]:
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


@final
@dataclasses.dataclass(slots=True)
class RecipeStore(Generic[_T_comp]):
    """Container which returns the largest value on demand."""

    _store: dict[interfaces.Recipe, _T_comp]
    _maxsize: int | None
    _maxvalue: tuple[interfaces.Recipe, _T_comp] | None
    _minvalue: tuple[interfaces.Recipe, _T_comp] | None
    _minlimit: _T_comp | None

    @classmethod
    def new(cls, maxsize: int | None) -> "RecipeStore":
        if maxsize is not None and maxsize < 1:
            raise ValueError(
                f"`maxsize` must be positive or `None` (was {maxsize})"
            )
        return RecipeStore(
            _store={},
            _maxsize=maxsize,
            _maxvalue=None,
            _minvalue=None,
            _minlimit=None,
        )

    def __set_maxvalue(self) -> None:
        if len(self._store) == 0:
            self._maxvalue = None
            return
        max_val = max(self._store.items(), key=operator.itemgetter(1, 0))
        self._maxvalue = max_val

    def __set_minvalue(self) -> None:
        if len(self._store) == 0:
            self._minvalue = None
            return
        min_val = min(self._store.items(), key=operator.itemgetter(1, 0))
        self._minvalue = min_val

    def min_val(self) -> tuple[interfaces.Recipe, _T_comp] | None:
        if self._minvalue is None:
            self.__set_minvalue()
        return self._minvalue

    def insert(self, item: interfaces.Recipe, value: _T_comp) -> None:
        if self._minlimit is not None and value < self._minlimit:
            return
        # if store is empty, set relevant values directly
        if len(self._store) == 0:
            self._store[item] = value
            self._maxvalue = item, value
            self._minvalue = item, value
            return

        # minvalue and maxvalue should never be None if nonempty
        assert self._minvalue is not None
        assert self._maxvalue is not None

        cur_value = self._store.get(item)
        # if item already exists, update it (if necessary) and return
        if cur_value is not None and cur_value > value:
            self._store[item] = value
            # if item is current minimum, refresh current minimum marker
            if self._minvalue is not None and self._minvalue[0] == item:
                self.__set_minvalue()
            # if item is or exceeds current maximum,
            # set new maximum to new value
            if self._maxvalue[0] == item or (
                self._maxvalue[1],
                self._maxvalue[0],
            ) < (value, item):
                self._maxvalue = item, value

            return

        # item does not already exist in store

        # if stack is not yet full, add directly and return
        if self._maxsize is None or len(self._store) < self._maxsize:
            self._store[item] = value
            # if item has lower value than current min, replace it
            if (value, item) < (self._minvalue[1], self._minvalue[0]):
                self._minvalue = (item, value)
            # if item has greater value than current max, replace it
            elif (self._maxvalue[1], self._maxvalue[0]) < (value, item):
                self._maxvalue = (item, value)
            return

        # stack is full

        # if item is lower than current minimum, eject it and return
        if (value, item) < (self._minvalue[1], self._minvalue[0]):
            return

        # eject lowest item, insert new item, and refresh
        self._minlimit = self._minvalue[1]
        del self._store[self._minvalue[0]]
        self._store[item] = value
        if (self._maxvalue[1], self._maxvalue[0]) < (value, item):
            self._maxvalue = (item, value)
        self.__set_minvalue()

    def extend(self, data: Iterable[tuple[interfaces.Recipe, _T_comp]]) -> None:
        for item, value in data:
            self.insert(item, value)

    def pop_max(self) -> tuple[interfaces.Recipe, _T_comp] | None:
        if len(self) == 0:
            return None
        assert self._maxvalue is not None
        result = self._maxvalue
        del self._store[result[0]]
        self.__set_maxvalue()
        if len(self._store) == 0:
            self._minvalue = None
        return result

    def pop_maxn(
        self, n: int | None
    ) -> Sequence[tuple[interfaces.Recipe, _T_comp]]:
        if len(self) == 0:
            return ()
        if n is None:
            results = sorted(
                self._store.items(),
                key=operator.itemgetter(1, 0),
                reverse=True,
            )
            self._store.clear()
            self._maxvalue = None
            self._minvalue = None
            return results
        elif n == 1:
            result = self.pop_max()
            assert result is not None
            return (result,)

        results = heapq.nlargest(
            n,
            self._store.items(),
            key=operator.itemgetter(1, 0),
        )
        for item, _ in results:
            del self._store[item]
            if self._minvalue is not None and item == self._minvalue[0]:
                self._minvalue = None
        self.__set_maxvalue()
        return results

    def reset(self) -> None:
        if len(self._store) > 0:
            raise ValueError(
                f"Can only reset when empty (size: {len(self._store)})"
            )
        self._minlimit = None

    @classmethod
    def from_iter(
        cls,
        data: Iterable[tuple[interfaces.Recipe, _T_comp]],
        maxsize: int | None,
    ) -> "RecipeStore[_T_comp]":
        new_store: RecipeStore[_T_comp] = RecipeStore.new(maxsize=maxsize)
        for item, value in data:
            new_store.insert(item, value)
        return new_store

    def __len__(self) -> int:
        return len(self._store)


def execute_reaction(
    rxn_job: ReactionJob,
) -> collections.abc.Generator[
    tuple[interfaces.ReactionExplicit, bool], None, None
]:
    reactants = tuple(mol.item for mol in rxn_job.op_args)
    if rxn_job.operator.item is None or any(mol is None for mol in reactants):
        raise ValueError("ReactionJob has non-None item components!")
    product_packets = rxn_job.operator.item(*reactants)
    reactant_datapackets = tuple(
        interfaces.DataPacketE(p.i, p.item, p.meta) for p in rxn_job.op_args
    )
    operator_datapacket = interfaces.DataPacketE(
        rxn_job.operator.i, rxn_job.operator.item, rxn_job.operator.meta
    )
    for rxn_products in product_packets:
        product_datapackets = tuple(
            interfaces.DataPacketE(-1, product, None)
            for product in rxn_products
        )
        rxn = interfaces.ReactionExplicit(
            operator_datapacket,
            reactant_datapackets,
            product_datapackets,
            None,
        )
        yield rxn, True


def execute_reactions(
    rxn_jobs: collections.abc.Collection[ReactionJob],
    rxn_analysis: typing.Optional[metadata.RxnAnalysisStep] = None,
) -> collections.abc.Iterable[tuple[interfaces.ReactionExplicit, bool]]:
    rxn_generator = itertools.chain.from_iterable(
        execute_reaction(rxn_job) for rxn_job in rxn_jobs
    )
    if rxn_analysis is None:
        return rxn_generator
    return rxn_analysis.execute(rxn_generator)


class PriorityQueueStrategyBasic(interfaces.PriorityQueueStrategy):
    __slots__ = ("_network",)

    def __init__(
        self,
        network: interfaces.ChemNetwork,
        num_procs: typing.Optional[int] = None,
    ) -> None:
        if num_procs is not None:
            raise NotImplementedError(
                f"Parallel processing not yet supported on {type(self)}"
            )
        self._network = network

    def expand(
        self,
        max_recipes: typing.Optional[int] = None,
        mol_filter: typing.Optional[interfaces.MolFilter] = None,
        bundle_filter: typing.Optional[interfaces.BundleFilter] = None,
        recipe_filter: typing.Optional[interfaces.RecipeFilter] = None,
        recipe_ranker: typing.Optional[interfaces.RecipeRanker] = None,
        reaction_plan: typing.Optional[
            typing.Union[
                metadata.RxnAnalysisStep,
                metadata.PropertyCompositor,
                metadata.ReactionFilterBase,
                metadata.LocalPropertyCalc,
            ]
        ] = None,
        global_hooks: typing.Optional[
            collections.abc.Sequence[interfaces.GlobalUpdateHook]
        ] = None,
        heap_size: typing.Optional[int] = None,
        beam_size: typing.Optional[int] = 1,
        batch_size: typing.Optional[int] = None,
        save_unreactive: bool = True,
    ) -> None:
        rxn_analysis_task: typing.Optional[metadata.RxnAnalysisStep] = None
        if reaction_plan is not None:
            rxn_analysis_task = metadata.as_rxn_analysis_step(reaction_plan)
            mc_update = rxn_analysis_task.resolver
        else:
            mc_update = metadata.MetaUpdateResolver({}, {}, {})

        if heap_size is not None and beam_size is not None:
            raise ValueError(
                f"""Heap size ({heap_size}) must be greater than beam size
                    ({beam_size})"""
            )

        # set keysets so that updated reactions may occur and parameters may be
        # passed to parallel processes
        recipe_keyset: interfaces.MetaKeyPacket = interfaces.MetaKeyPacket()
        reaction_keyset: interfaces.MetaKeyPacket = interfaces.MetaKeyPacket()
        if mol_filter is not None:
            recipe_keyset = recipe_keyset + mol_filter.meta_required
        if bundle_filter is not None:
            recipe_keyset = recipe_keyset + bundle_filter.meta_required
        if recipe_filter is not None:
            recipe_keyset = recipe_keyset + recipe_filter.meta_required
        if recipe_ranker is not None:
            recipe_keyset = recipe_keyset + recipe_ranker.meta_required
        if rxn_analysis_task is not None:
            reaction_keyset = reaction_keyset + rxn_analysis_task.meta_required
        total_keyset = recipe_keyset + reaction_keyset

        # initialize loop variables
        network = self._network
        compat_indices_table = [
            [0 for _ in network.compat_table(interfaces.OpIndex(i))]
            for i in range(len(network.ops))
        ]
        updated_mols_set: set[interfaces.MolIndex] = set()
        updated_ops_set: set[interfaces.OpIndex] = set()
        recipe_heap: RecipeHeap = RecipeHeap(maxsize=heap_size)
        recipes_tested: set[interfaces.Recipe] = set()

        while (max_recipes is None or max_recipes > 0) and (
            any(
                any(
                    i_tested < len(compat_mols)
                    for i_tested, compat_mols in zip(
                        op_index_table,
                        network.compat_table(interfaces.OpIndex(op_index)),
                        strict=False,
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
            cur_min = recipe_heap.min
            for opIndex, _ in enumerate(network.ops):
                # for each argument, accumulate a total of old_mols and new_mols
                compat_table = network.compat_table(interfaces.OpIndex(opIndex))
                compat_indices = compat_indices_table[opIndex]
                if not all(compat_table):
                    continue
                # generate recipe batches
                for batch in _generate_recipe_batches(
                    compat_table, compat_indices, batch_size, updated_mols_set
                ):
                    # assemble recipe ranking job
                    recipejob = assemble_recipe_batch_job(
                        interfaces.OpIndex(opIndex),
                        batch,
                        network,
                        recipe_keyset,
                        recipe_ranker,
                        heap_size,
                        recipe_filter,
                        mol_filter,
                        bundle_filter,
                    )
                    new_recipe_heap = execute_recipe_ranking(
                        recipejob, cur_min, recipes_tested
                    )
                    recipe_heap = recipe_heap + new_recipe_heap

                continue

            # update compat_indices_table
            compat_indices_table = [
                [
                    len(mol_list)
                    for mol_list in network.compat_table(interfaces.OpIndex(i))
                ]
                for i in range(len(network.ops))
            ]

            # perform beam expansion
            recipes_to_be_expanded = recipe_heap.popvals(beam_size)

            if len(recipes_to_be_expanded) == 0:
                break

            if max_recipes is not None:
                if len(recipes_to_be_expanded) > max_recipes:
                    recipes_to_be_expanded = recipes_to_be_expanded[
                        :max_recipes
                    ]
                max_recipes = max_recipes - len(recipes_to_be_expanded)

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
                if not save_unreactive and not pass_filter:
                    continue
                # add product mols to network
                products_indices = tuple(
                    network.add_mol(mol.item, None, pass_filter)
                    for mol in rxn.products
                    if mol.item is not None
                )

                # build reaction
                reactants_indices = tuple(
                    interfaces.MolIndex(mol.i) for mol in rxn.reactants
                )
                rxn_implicit = interfaces.Reaction(
                    interfaces.OpIndex(rxn.operator.i),
                    reactants_indices,
                    products_indices,
                )

                # add reaction to network
                rxn_index = network.add_rxn(rxn=rxn_implicit)

                updated_mols_set = set()
                updated_ops_set = set()

                # update reactant metadata
                key: collections.abc.Hashable
                for m_dat in zip(
                    reactants_indices, rxn.reactants, strict=False
                ):
                    if m_dat[1].meta is not None:
                        cur_vals = network.mols.meta(m_dat[0], m_dat[1].meta)
                        for key, v in m_dat[1].meta.items():
                            value = v
                            if key in mc_update.mol_updates and key in cur_vals:
                                value = mc_update.mol_updates[key](
                                    value, cur_vals[key]
                                )
                            if key not in cur_vals or cur_vals[key] != value:
                                network.mols.set_meta(m_dat[0], {key: value})
                                if key in total_keyset.molecule_keys:
                                    updated_mols_set.add(m_dat[0])

                # update product metadata
                for m_dat in zip(products_indices, rxn.products, strict=False):
                    if m_dat[1].meta is not None:
                        cur_vals = network.mols.meta(m_dat[0], m_dat[1].meta)
                        for key, v in m_dat[1].meta.items():
                            value = v
                            if key in mc_update.mol_updates and key in cur_vals:
                                value = mc_update.mol_updates[key](
                                    value, cur_vals[key]
                                )
                            if key not in cur_vals or cur_vals[key] != value:
                                network.mols.set_meta(m_dat[0], {key: value})
                                if key in total_keyset.molecule_keys:
                                    updated_mols_set.add(m_dat[0])

                # update operator metadata
                if rxn.operator.meta is not None:
                    cur_vals = network.ops.meta(
                        rxn_implicit.operator, rxn.operator.meta
                    )
                    for key, v in rxn.operator.meta.items():
                        value = v
                        if key in mc_update.op_updates and key in cur_vals:
                            value = mc_update.op_updates[key](
                                value,
                                cur_vals[key],
                            )
                        if key not in cur_vals or cur_vals[key] != value:
                            network.ops.set_meta(
                                rxn_implicit.operator, {key: value}
                            )
                            if key in total_keyset.operator_keys:
                                updated_ops_set.add(rxn_implicit.operator)

                # update reaction metadata
                if rxn.reaction_meta is not None:
                    cur_vals = network.rxns.meta(rxn_index, rxn.reaction_meta)
                    for key, v in rxn.reaction_meta.items():
                        value = v
                        if key in mc_update.rxn_updates and key in cur_vals:
                            value = mc_update.rxn_updates[key](
                                value, cur_vals[key]
                            )
                        if key not in cur_vals or cur_vals[key] != value:
                            network.rxns.set_meta(rxn_index, {key: value})

            recipes_tested.update(
                (reciperank.recipe for reciperank in recipes_to_be_expanded)
            )

            if len(recipe_heap) == 0:
                compat_indices_table = [
                    [0 for _ in network.compat_table(interfaces.OpIndex(i))]
                    for i in range(len(network.ops))
                ]

            # run global hooks
            if global_hooks is not None:
                fake_network = pgnetworks.ChemNetworkFacadeMetaTrigger(
                    self._network, total_keyset
                )
                end_on_completion: bool = False
                for hook_func in global_hooks:
                    rval = hook_func(fake_network)
                    if (
                        rval
                        is interfaces.GlobalHookReturnValue.STOP_SHORTCIRCUIT
                    ):
                        return
                    if rval is interfaces.GlobalHookReturnValue.STOP:
                        end_on_completion = True
                if end_on_completion:
                    return

            continue


class CartesianStrategyUpdated:
    def __init__(
        self,
        network: interfaces.ChemNetwork,
        engine: interfaces.NetworkEngine,
        #        gen_key: collections.abc.Hashable = "generation",
    ) -> None:
        self._network = network
        self._engine = engine

    #        self._gen_key = gen_key

    def expand(
        self,
        num_iter: typing.Optional[int] = None,
        max_recipes: typing.Optional[int] = None,
        mol_filter: typing.Optional[interfaces.MolFilter] = None,
        bundle_filter: typing.Optional[interfaces.BundleFilter] = None,
        recipe_filter: typing.Optional[interfaces.RecipeFilter] = None,
        reaction_plan: typing.Optional[
            typing.Union[
                metadata.RxnAnalysisStep,
                metadata.PropertyCompositor,
                metadata.ReactionFilterBase,
                metadata.LocalPropertyCalc,
            ]
        ] = None,
        global_hooks: typing.Optional[
            collections.abc.Sequence[interfaces.GlobalUpdateHook]
        ] = None,
        save_unreactive: bool = True,
        # max_gen: typing.Optional[int] = None,
    ):
        engine = self._engine
        p_strat = engine.strat.pq(self._network)
        # mol_filter: typing.Optional[interfaces.MolFilter] = None
        # calc: typing.Optional[metadata.MolPropertyFromRxnCalc] = None
        # resolver: typing.Optional[metadata.MetaUpdateResolver] = None
        global_hook: list[interfaces.GlobalUpdateHook]
        global_hook = [] if global_hooks is None else list(global_hooks)
        # if max_gen is not None:
        #     mol_filter = engine.filter.mol.meta_func(
        #         self._gen_key, self.gen_test(max_gen)
        #     )
        #     calc = engine.meta.generation(self._gen_key)
        #     resolver = metadata.MetaUpdateResolver({self._gen_key: min},
        #                                            {}, {})
        if num_iter is not None:
            global_hook.append(engine.hook.max_iter(num_iter))
        p_strat.expand(
            max_recipes,
            # mol_filter,
            # reaction_plan=calc,
            # mc_update=resolver,
            beam_size=None,
            global_hooks=global_hook,
            mol_filter=mol_filter,
            bundle_filter=bundle_filter,
            recipe_filter=recipe_filter,
            reaction_plan=reaction_plan,
            save_unreactive=save_unreactive,
        )


def execute_reactions_single(
    rxn_job: ReactionJob, rxn_analysis: metadata.RxnAnalysisStep | None = None
) -> tuple[tuple[interfaces.ReactionExplicit, bool], ...]:
    results = tuple(execute_reactions((rxn_job,), rxn_analysis))
    return results


def execute_reactions_parallel(
    rxn_jobs: Collection[ReactionJob],
    executor: concurrent.futures.ProcessPoolExecutor | None,
    rxn_analysis: metadata.RxnAnalysisStep | None,
) -> Iterable[tuple[interfaces.ReactionExplicit, bool]]:
    if executor is None:
        return execute_reactions(rxn_jobs, rxn_analysis)
    exec_result = executor.map(
        functools.partial(
            execute_reactions_single,
            rxn_analysis=rxn_analysis,
        ),
        rxn_jobs,
    )
    result = itertools.chain.from_iterable(exec_result)
    return result


def _execute_recipe_job(
    job: RecipeRankingJob, recipes_tested: Collection[interfaces.Recipe]
) -> Iterable[interfaces.Recipe]:
    """Execute a single recipe processing job."""
    args_edited: tuple[
        tuple[interfaces.DataPacket[interfaces.MolDatBase], ...], ...
    ]
    if job.mol_filter is not None:
        args_edited = tuple(
            tuple(
                mol for mol in arg_mols if job.mol_filter(mol, job.operator, i)
            )
            for i, arg_mols in enumerate(job.op_args)
        )
    else:
        args_edited = job.op_args

    # perform bundle filtering
    bundles: collections.abc.Iterable[interfaces.RecipeBundle]
    if job.bundle_filter is not None:
        bundles = job.bundle_filter(
            interfaces.RecipeBundle(job.operator, args_edited)
        )
    else:
        bundles = (interfaces.RecipeBundle(job.operator, args_edited),)

    recipe_generator = (
        (interfaces.RecipeExplicit(job.operator, reactants_data), recipe)
        for reactants_data, recipe in itertools.chain.from_iterable(
            (
                (
                    (
                        reactants_data,
                        interfaces.Recipe(
                            interfaces.OpIndex(job.operator.i),
                            tuple(
                                interfaces.MolIndex(reactant.i)
                                for reactant in reactants_data
                            ),
                        ),
                    )
                    for reactants_data in itertools.product(*bundle.args)
                )
                for bundle in bundles
            )
        )
        if recipe not in recipes_tested
    )

    if job.recipe_filter is None:
        return (
            recipe_item.recipe
            for recipe_item in (
                RecipePriorityItem(
                    None,
                    recipe,
                )
                for _, recipe in recipe_generator
            )
        )
    return (
        recipe_item.recipe
        for recipe_item in (
            RecipePriorityItem(None, recipe)
            for recipe_explicit, recipe in recipe_generator
            if job.recipe_filter(recipe_explicit)
        )
    )


def _execute_recipe_ranking_linear(
    job: RecipeRankingJob,
    min_val: typing.Optional[RecipePriorityItem],
    recipes_tested: collections.abc.Collection[interfaces.Recipe],
) -> Iterable[RecipePriorityItem]:
    # filter molecules
    args_edited: tuple[
        tuple[interfaces.DataPacket[interfaces.MolDatBase], ...], ...
    ]
    if job.mol_filter is not None:
        args_edited = tuple(
            tuple(
                mol for mol in arg_mols if job.mol_filter(mol, job.operator, i)
            )
            for i, arg_mols in enumerate(job.op_args)
        )
    else:
        args_edited = job.op_args

    # perform bundle filtering
    bundles: collections.abc.Iterable[interfaces.RecipeBundle]
    if job.bundle_filter is not None:
        bundles = job.bundle_filter(
            interfaces.RecipeBundle(job.operator, args_edited)
        )
    else:
        bundles = (interfaces.RecipeBundle(job.operator, args_edited),)

    recipe_generator = (
        (interfaces.RecipeExplicit(job.operator, reactants_data), recipe)
        for reactants_data, recipe in itertools.chain.from_iterable(
            (
                (
                    (
                        reactants_data,
                        interfaces.Recipe(
                            interfaces.OpIndex(job.operator.i),
                            tuple(
                                interfaces.MolIndex(reactant.i)
                                for reactant in reactants_data
                            ),
                        ),
                    )
                    for reactants_data in itertools.product(*bundle.args)
                )
                for bundle in bundles
            )
        )
        if recipe not in recipes_tested
    )

    if job.recipe_ranker is None:
        if job.recipe_filter is None:
            return (
                recipe_item
                for recipe_item in (
                    RecipePriorityItem(
                        None,
                        recipe,
                    )
                    for _, recipe in recipe_generator
                )
                if min_val is None or not (recipe_item < min_val)
            )
        return (
            recipe_item
            for recipe_item in (
                RecipePriorityItem(None, recipe)
                for recipe_explicit, recipe in recipe_generator
                if job.recipe_filter(recipe_explicit)
            )
            if min_val is None or not (recipe_item < min_val)
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
    return priority_item_generator


def _execute_recipe_ranking_parallel(
    jobs: Collection[RecipeRankingJob],
    min_val: RecipePriorityItem | None,
    recipes_tested: Collection[interfaces.Recipe],
    executor: concurrent.futures.Executor | None,
) -> Iterable[RecipePriorityItem]:
    if executor is None:
        result_gen = itertools.chain.from_iterable(
            _execute_recipe_ranking_linear(
                job=job, min_val=min_val, recipes_tested=recipes_tested
            )
            for job in jobs
        )
    else:
        result_gen = itertools.chain.from_iterable(
            executor.map(
                functools.partial(
                    _execute_recipe_ranking_linear,
                    min_val=min_val,
                    recipes_tested=recipes_tested,
                ),
                jobs,
            )
        )
    return result_gen


@final
@dataclasses.dataclass(frozen=True, slots=True)
class CartesianStrategyV3:
    _network: interfaces.ChemNetwork

    @classmethod
    def new(cls, network: interfaces.ChemNetwork) -> "CartesianStrategyV3":
        """Attach new Cartesian Strategy based on V3 methods."""
        return CartesianStrategyV3(_network=network)

    def _expand(
        self,
        num_iter: int | None,
        max_recipes: int | None,
        mol_filter: interfaces.MolFilter | None,
        bundle_filter: interfaces.BundleFilter | None,
        recipe_filter: interfaces.RecipeFilter | None,
        reaction_plan: metadata.RxnAnalysisStep
        | metadata.PropertyCompositor
        | metadata.ReactionFilterBase
        | metadata.LocalPropertyCalc
        | None,
        global_hooks: Sequence[interfaces.GlobalUpdateHook] | None,
        save_unreactive: bool,
        executor: concurrent.futures.ProcessPoolExecutor | None,
    ) -> None:
        rxn_analysis_task: metadata.RxnAnalysisStep | None = None
        if reaction_plan is not None:
            rxn_analysis_task = metadata.as_rxn_analysis_step(reaction_plan)
            mc_update = rxn_analysis_task.resolver
        else:
            mc_update = metadata.MetaUpdateResolver({}, {}, {})

        # set keysets so that updated reactions may occur and parameters may
        # be passed to parallel processes
        recipe_keyset = interfaces.MetaKeyPacket()
        reaction_keyset = interfaces.MetaKeyPacket()
        if mol_filter is not None:
            recipe_keyset = recipe_keyset + mol_filter.meta_required
        if bundle_filter is not None:
            recipe_keyset = recipe_keyset + bundle_filter.meta_required
        if recipe_filter is not None:
            recipe_keyset = recipe_keyset + recipe_filter.meta_required
        if rxn_analysis_task is not None:
            reaction_keyset = reaction_keyset + rxn_analysis_task.meta_required
        total_keyset = recipe_keyset + reaction_keyset

        # initialize loop variables
        network = self._network
        compat_indices_tracker = [
            [0 for _ in network.compat_table(interfaces.OpIndex(i))]
            for i in range(len(network.ops))
        ]
        updated_mols_set: set[interfaces.MolIndex] = set()
        updated_ops_set: set[interfaces.OpIndex] = set()
        recipe_queue: collections.deque[interfaces.Recipe] = collections.deque()
        recipes_tested: set[interfaces.Recipe] = set()

        # set this to ensure iter count is correct
        _beam_size: int | None = None

        # set ending conditions
        remain_recipes = max_recipes
        cur_iter = 0
        while (
            (num_iter is None or cur_iter < num_iter)
            and (remain_recipes is None or remain_recipes > 0)
            and (
                any(
                    any(
                        i_tested < len(compat_mols)
                        for i_tested, compat_mols in zip(
                            op_index_table,
                            network.compat_table(interfaces.OpIndex(op_index)),
                            strict=False,
                        )
                    )
                    for op_index, op_index_table in enumerate(
                        compat_indices_tracker
                    )
                )
                or len(updated_mols_set) != 0
                or len(updated_ops_set) != 0
                or len(recipe_queue) > 0
            )
        ):
            cur_iter += 1
            for op_index in updated_ops_set:
                compat_indices_tracker[op_index] = [
                    0 for _ in range(len(compat_indices_tracker[op_index]))
                ]

            # for each operator, create recipe batches
            for opIndex in range(len(network.ops)):
                # for each argument, accumulate a total of
                # old mols and new mols
                compat_table = network.compat_table(interfaces.OpIndex(opIndex))
                compat_indices = compat_indices_tracker[opIndex]

                # if any of the entries is zero-length, no recipes possible
                if not all(compat_table):
                    continue

                # generate recipe batches
                for batch in _generate_recipe_batches(
                    compat_table, compat_indices, None, updated_mols_set
                ):
                    # assemble recipes and add to deque
                    recipejob = assemble_recipe_batch_job(
                        interfaces.OpIndex(opIndex),
                        batch,
                        network,
                        recipe_keyset,
                        None,
                        None,
                        recipe_filter,
                        mol_filter,
                        bundle_filter,
                    )
                    new_recipes = _execute_recipe_job(recipejob, recipes_tested)
                    recipe_queue.extend(new_recipes)

            # update compat_indices_tracker
            compat_indices_tracker = [
                [
                    len(mol_list)
                    for mol_list in network.compat_table(interfaces.OpIndex(i))
                ]
                for i in range(len(network.ops))
            ]

            # perform beam expansion
            recipes_to_be_expanded = (
                tuple(recipe_queue.popleft() for _ in range(_beam_size))
                if _beam_size is not None
                else tuple(
                    recipe_queue.popleft() for _ in range(len(recipe_queue))
                )
            )

            # if no available, recipes, break main loop (exit)
            if len(recipes_to_be_expanded) == 0:
                break

            if remain_recipes is not None:
                if len(recipes_to_be_expanded) > remain_recipes:
                    recipes_to_be_expanded = recipes_to_be_expanded[
                        :remain_recipes
                    ]
                remain_recipes -= len(recipes_to_be_expanded)

            reaction_jobs = tuple(
                assemble_reaction_job(recipe, network, reaction_keyset)
                for recipe in recipes_to_be_expanded
            )

            for rxn, pass_filter in execute_reactions_parallel(
                reaction_jobs, executor, rxn_analysis_task
            ):
                if not save_unreactive and not pass_filter:
                    continue
                # add product mols to network
                products_indices = tuple(
                    network.add_mol(mol.item, None, pass_filter)
                    for mol in rxn.products
                    if mol.item is not None
                )

                # build reaction
                reactants_indices = tuple(
                    interfaces.MolIndex(mol.i) for mol in rxn.reactants
                )
                rxn_implicit = interfaces.Reaction(
                    interfaces.OpIndex(rxn.operator.i),
                    reactants_indices,
                    products_indices,
                )

                # add reaction to network
                rxn_index = network.add_rxn(rxn=rxn_implicit)

                updated_mols_set = set()
                updated_ops_set = set()

                # update reactant metadata
                key: collections.abc.Hashable
                for m_dat in zip(
                    reactants_indices, rxn.reactants, strict=False
                ):
                    if m_dat[1].meta is not None:
                        cur_vals = network.mols.meta(m_dat[0], m_dat[1].meta)
                        for key, v in m_dat[1].meta.items():
                            value = v
                            if key in mc_update.mol_updates and key in cur_vals:
                                value = mc_update.mol_updates[key](
                                    value, cur_vals[key]
                                )
                            if key not in cur_vals or cur_vals[key] != value:
                                network.mols.set_meta(m_dat[0], {key: value})
                                if key in total_keyset.molecule_keys:
                                    updated_mols_set.add(m_dat[0])

                # update product metadata
                for m_dat in zip(products_indices, rxn.products, strict=False):
                    if m_dat[1].meta is not None:
                        cur_vals = network.mols.meta(m_dat[0], m_dat[1].meta)
                        for key, v in m_dat[1].meta.items():
                            value = v
                            if key in mc_update.mol_updates and key in cur_vals:
                                value = mc_update.mol_updates[key](
                                    value, cur_vals[key]
                                )
                            if key not in cur_vals or cur_vals[key] != value:
                                network.mols.set_meta(m_dat[0], {key: value})
                                if key in total_keyset.molecule_keys:
                                    updated_mols_set.add(m_dat[0])

                # update operator metadata
                if rxn.operator.meta is not None:
                    cur_vals = network.ops.meta(
                        rxn_implicit.operator, rxn.operator.meta
                    )
                    for key, v in rxn.operator.meta.items():
                        value = v
                        if key in mc_update.op_updates and key in cur_vals:
                            value = mc_update.op_updates[key](
                                value,
                                cur_vals[key],
                            )
                        if key not in cur_vals or cur_vals[key] != value:
                            network.ops.set_meta(
                                rxn_implicit.operator, {key: value}
                            )
                            if key in total_keyset.operator_keys:
                                updated_ops_set.add(rxn_implicit.operator)

                # update reaction metadata
                if rxn.reaction_meta is not None:
                    cur_vals = network.rxns.meta(rxn_index, rxn.reaction_meta)
                    for key, v in rxn.reaction_meta.items():
                        value = v
                        if key in mc_update.rxn_updates and key in cur_vals:
                            value = mc_update.rxn_updates[key](
                                value, cur_vals[key]
                            )
                        if key not in cur_vals or cur_vals[key] != value:
                            network.rxns.set_meta(rxn_index, {key: value})

            recipes_tested.update(recipes_to_be_expanded)

            if len(recipe_queue) == 0:
                compat_indices_tracker = [
                    [0 for _ in network.compat_table(interfaces.OpIndex(i))]
                    for i in range(len(network.ops))
                ]

            # run global hooks
            if global_hooks is not None:
                fake_network = pgnetworks.ChemNetworkFacadeMetaTrigger(
                    self._network, total_keyset
                )
                end_on_completion: bool = False
                for hook_func in global_hooks:
                    rval = hook_func(fake_network)
                    if (
                        rval
                        is interfaces.GlobalHookReturnValue.STOP_SHORTCIRCUIT
                    ):
                        return
                    if rval is interfaces.GlobalHookReturnValue.STOP:
                        end_on_completion = True
                if end_on_completion:
                    return

            continue

    def expand(
        self,
        num_iter: int | None = None,
        max_recipes: int | None = None,
        mol_filter: interfaces.MolFilter | None = None,
        bundle_filter: interfaces.BundleFilter | None = None,
        recipe_filter: interfaces.RecipeFilter | None = None,
        reaction_plan: metadata.RxnAnalysisStep
        | metadata.PropertyCompositor
        | metadata.ReactionFilterBase
        | metadata.LocalPropertyCalc
        | None = None,
        global_hooks: Sequence[interfaces.GlobalUpdateHook] | None = None,
        save_unreactive: bool = True,
        num_workers: int = 1,
    ) -> None:
        """
        Expand the network according to certain parameters.

        Parameters
        ----------
        max_recipes : int | None (default: None)
            Maximum number of recipes to evaluate.  Value of `None` indicates
            no limit to recipes.
        mol_filter : interfaces.MolFilter | None (default: None)
            Filters to use to stop specific molecules from being evaluated.
        bundle_filter : interfaces.BundleFilter | None (default: None)
            Filters to use to stop specific combinations of mols from being
            evaluated.
        recipe_filter : interfaces.RecipeFilter | None (default: None)
            Filters to use to stop recipes from being ranked or evaluated.  May
            be a composition of filters using basic logic operators.
        reaction_plan : metadata.RxnAnalysisStep | metadata.PropertyCompositor |
                        metadata.ReactionFilterBase |
                        metadata.LocalPropertyCalc | None (default: None)
            Metadata calculator and reaction filter.  Should be a "Reaction
            Analysis Flow" as detailed in examples.
        global_hooks : Sequence[interfaces.GlobalUpdateHook] | None
                       (default: None)
            Hook functions to run after a loop has completed.  If any return a
            boolean, expansion will terminate before running any more reactions.
            Return value of True terminates immediately, whereas return value
            of False terminates after remaining hooks also run.  `None`
            indicates expansion will continue until all recipes have been
            exhausted.

        Other Parameters
        ----------------
        save_unreactive : bool (default: True)
            Store reactions rejected by reaction_plan in the network.  If
            False, reactions which are rejected will simply be deleted, along
            with their products, instead of stored.  Depending on the
            proportion of rejected reactions, this may save a lot of memory.
        num_workers : int (default: 1)
            Number of parallel workers to deploy.  May result in
            non-deterministic output.  0 indicates using all cores.
        """
        if num_workers < 0:
            raise ValueError(
                f"Number of workers must be >= 0 (was {num_workers})"
            )
        match num_workers:
            # use default cores
            case 0:
                with concurrent.futures.ProcessPoolExecutor(
                    max_workers=None
                ) as executor:
                    self._expand(
                        num_iter=num_iter,
                        max_recipes=max_recipes,
                        mol_filter=mol_filter,
                        bundle_filter=bundle_filter,
                        recipe_filter=recipe_filter,
                        reaction_plan=reaction_plan,
                        global_hooks=global_hooks,
                        save_unreactive=save_unreactive,
                        executor=executor,
                    )
            # use no multiprocessing
            case 1:
                self._expand(
                    num_iter=num_iter,
                    max_recipes=max_recipes,
                    mol_filter=mol_filter,
                    bundle_filter=bundle_filter,
                    recipe_filter=recipe_filter,
                    reaction_plan=reaction_plan,
                    global_hooks=global_hooks,
                    save_unreactive=save_unreactive,
                    executor=None,
                )
            # use fixed number of max cores
            case _:
                with concurrent.futures.ProcessPoolExecutor(
                    max_workers=num_workers
                ) as executor:
                    self._expand(
                        num_iter=num_iter,
                        max_recipes=max_recipes,
                        mol_filter=mol_filter,
                        bundle_filter=bundle_filter,
                        recipe_filter=recipe_filter,
                        reaction_plan=reaction_plan,
                        global_hooks=global_hooks,
                        save_unreactive=save_unreactive,
                        executor=executor,
                    )


@final
@dataclasses.dataclass(frozen=True, slots=True)
class PriorityQueueStrategyV2(interfaces.PriorityQueueStrategy):
    _network: interfaces.ChemNetwork

    @classmethod
    def new(cls, network: interfaces.ChemNetwork) -> "PriorityQueueStrategyV2":
        """Attach new Priority Queue Strategy based on V2 methods."""
        return PriorityQueueStrategyV2(_network=network)

    def _expand(
        self,
        max_recipes: int | None,
        mol_filter: interfaces.MolFilter | None,
        bundle_filter: interfaces.BundleFilter | None,
        recipe_filter: interfaces.RecipeFilter | None,
        recipe_ranker: interfaces.RecipeRanker | None,
        reaction_plan: metadata.RxnAnalysisStep
        | metadata.PropertyCompositor
        | metadata.ReactionFilterBase
        | metadata.LocalPropertyCalc
        | None,
        global_hooks: Sequence[interfaces.GlobalUpdateHook] | None,
        heap_size: int | None,
        beam_size: int | None,
        batch_size: int | None,
        save_unreactive: bool,
        executor: concurrent.futures.ProcessPoolExecutor | None,
    ) -> None:
        remain_recipes = max_recipes
        del max_recipes
        rxn_analysis_task: metadata.RxnAnalysisStep | None = None
        if reaction_plan is not None:
            rxn_analysis_task = metadata.as_rxn_analysis_step(reaction_plan)
            mc_update = rxn_analysis_task.resolver
        else:
            mc_update = metadata.MetaUpdateResolver({}, {}, {})

        # set keysets so that updated reactions may occur and parameters may be
        # passed to parallel processes
        recipe_keyset: interfaces.MetaKeyPacket = interfaces.MetaKeyPacket()
        reaction_keyset: interfaces.MetaKeyPacket = interfaces.MetaKeyPacket()
        if mol_filter is not None:
            recipe_keyset = recipe_keyset + mol_filter.meta_required
        if bundle_filter is not None:
            recipe_keyset = recipe_keyset + bundle_filter.meta_required
        if recipe_filter is not None:
            recipe_keyset = recipe_keyset + recipe_filter.meta_required
        if recipe_ranker is not None:
            recipe_keyset = recipe_keyset + recipe_ranker.meta_required
        if rxn_analysis_task is not None:
            reaction_keyset = reaction_keyset + rxn_analysis_task.meta_required
        total_keyset = recipe_keyset + reaction_keyset

        # initialize loop variables
        network = self._network
        compat_indices_tracker = [
            [0 for _ in network.compat_table(interfaces.OpIndex(i))]
            for i in range(len(network.ops))
        ]
        updated_mols_set: set[interfaces.MolIndex] = set()
        updated_ops_set: set[interfaces.OpIndex] = set()
        recipe_store: RecipeStore = RecipeStore.new(maxsize=heap_size)
        recipes_tested: set[interfaces.Recipe] = set()

        while (remain_recipes is None or remain_recipes > 0) and (
            any(
                any(
                    i_tested < len(compat_mols)
                    for i_tested, compat_mols in zip(
                        op_index_table,
                        network.compat_table(interfaces.OpIndex(op_index)),
                        strict=False,
                    )
                )
                for op_index, op_index_table in enumerate(
                    compat_indices_tracker
                )
            )
            or len(updated_mols_set) != 0
            or len(updated_ops_set) != 0
            or len(recipe_store) > 0
        ):
            # if recipe store is empty, reset it and the tracker indices
            if len(recipe_store) == 0:
                recipe_store.reset()
                updated_ops_set = set(
                    interfaces.OpIndex(i) for i in range(len(network.ops))
                )

            for op_index in updated_ops_set:
                compat_indices_tracker[op_index] = [
                    0 for _ in range(len(compat_indices_tracker[op_index]))
                ]

            # for each operator, create recipe batches

            # current minimum only applies when store is full
            cur_minv = recipe_store.min_val()
            cur_min = (
                RecipePriorityItem(rank=cur_minv[1], recipe=cur_minv[0])
                if cur_minv is not None
                and heap_size is not None
                and len(recipe_store) >= heap_size
                else None
            )
            for opIndex in range(len(network.ops)):
                # for each argument, accumulate a total of old_mols and new_mols
                compat_table = network.compat_table(interfaces.OpIndex(opIndex))
                compat_indices = compat_indices_tracker[opIndex]
                if not all(compat_table):
                    continue
                # generate recipe batches
                recipe_batch_jobs: list[RecipeRankingJob] = []
                for batch in _generate_recipe_batches(
                    compat_table, compat_indices, batch_size, updated_mols_set
                ):
                    # assemble recipe ranking job
                    recipejob = assemble_recipe_batch_job(
                        interfaces.OpIndex(opIndex),
                        batch,
                        network,
                        recipe_keyset,
                        recipe_ranker,
                        heap_size,
                        recipe_filter,
                        mol_filter,
                        bundle_filter,
                    )
                    recipe_batch_jobs.append(recipejob)
                new_recipes = _execute_recipe_ranking_parallel(
                    recipe_batch_jobs,
                    cur_min,
                    recipes_tested,
                    executor=executor,
                )
                for recipe in new_recipes:
                    recipe_store.insert(recipe.recipe, recipe.rank)

                continue

            # update compat_indices_table
            compat_indices_tracker = [
                [
                    len(mol_list)
                    for mol_list in network.compat_table(interfaces.OpIndex(i))
                ]
                for i in range(len(network.ops))
            ]

            if len(recipe_store) == 0:
                break

            # perform beam expansion
            recipes_to_be_expanded = [
                RecipePriorityItem(r[1], r[0])
                for r in recipe_store.pop_maxn(beam_size)
            ]

            if len(recipes_to_be_expanded) == 0:
                break

            if remain_recipes is not None:
                if len(recipes_to_be_expanded) > remain_recipes:
                    recipes_to_be_expanded = recipes_to_be_expanded[
                        :remain_recipes
                    ]
                remain_recipes = remain_recipes - len(recipes_to_be_expanded)

            reaction_jobs = tuple(
                assemble_reaction_job(
                    reciperank.recipe, network, reaction_keyset
                )
                for reciperank in recipes_to_be_expanded
            )

            # execute reactions
            for rxn, pass_filter in execute_reactions_parallel(
                reaction_jobs, executor, rxn_analysis_task
            ):
                if not save_unreactive and not pass_filter:
                    continue
                # add product mols to network
                products_indices = tuple(
                    network.add_mol(mol.item, None, pass_filter)
                    for mol in rxn.products
                    if mol.item is not None
                )

                # build reaction
                reactants_indices = tuple(
                    interfaces.MolIndex(mol.i) for mol in rxn.reactants
                )
                rxn_implicit = interfaces.Reaction(
                    interfaces.OpIndex(rxn.operator.i),
                    reactants_indices,
                    products_indices,
                )

                # add reaction to network
                rxn_index = network.add_rxn(rxn=rxn_implicit)

                updated_mols_set = set()
                updated_ops_set = set()

                # update reactant metadata
                key: collections.abc.Hashable
                for m_dat in zip(
                    reactants_indices, rxn.reactants, strict=False
                ):
                    if m_dat[1].meta is not None:
                        cur_vals = network.mols.meta(m_dat[0], m_dat[1].meta)
                        for key, v in m_dat[1].meta.items():
                            value = v
                            if key in mc_update.mol_updates and key in cur_vals:
                                value = mc_update.mol_updates[key](
                                    value, cur_vals[key]
                                )
                            if key not in cur_vals or cur_vals[key] != value:
                                network.mols.set_meta(m_dat[0], {key: value})
                                if key in total_keyset.molecule_keys:
                                    updated_mols_set.add(m_dat[0])

                # update product metadata
                for m_dat in zip(products_indices, rxn.products, strict=False):
                    if m_dat[1].meta is not None:
                        cur_vals = network.mols.meta(m_dat[0], m_dat[1].meta)
                        for key, v in m_dat[1].meta.items():
                            value = v
                            if key in mc_update.mol_updates and key in cur_vals:
                                value = mc_update.mol_updates[key](
                                    value, cur_vals[key]
                                )
                            if key not in cur_vals or cur_vals[key] != value:
                                network.mols.set_meta(m_dat[0], {key: value})
                                if key in total_keyset.molecule_keys:
                                    updated_mols_set.add(m_dat[0])

                # update operator metadata
                if rxn.operator.meta is not None:
                    cur_vals = network.ops.meta(
                        rxn_implicit.operator, rxn.operator.meta
                    )
                    for key, v in rxn.operator.meta.items():
                        value = v
                        if key in mc_update.op_updates and key in cur_vals:
                            value = mc_update.op_updates[key](
                                value,
                                cur_vals[key],
                            )
                        if key not in cur_vals or cur_vals[key] != value:
                            network.ops.set_meta(
                                rxn_implicit.operator, {key: value}
                            )
                            if key in total_keyset.operator_keys:
                                updated_ops_set.add(rxn_implicit.operator)

                # update reaction metadata
                if rxn.reaction_meta is not None:
                    cur_vals = network.rxns.meta(rxn_index, rxn.reaction_meta)
                    for key, v in rxn.reaction_meta.items():
                        value = v
                        if key in mc_update.rxn_updates and key in cur_vals:
                            value = mc_update.rxn_updates[key](
                                value, cur_vals[key]
                            )
                        if key not in cur_vals or cur_vals[key] != value:
                            network.rxns.set_meta(rxn_index, {key: value})

            recipes_tested.update(
                (reciperank.recipe for reciperank in recipes_to_be_expanded)
            )

            if len(recipe_store) == 0:
                compat_indices_tracker = [
                    [0 for _ in network.compat_table(interfaces.OpIndex(i))]
                    for i in range(len(network.ops))
                ]

            # run global hooks
            if global_hooks is not None:
                fake_network = pgnetworks.ChemNetworkFacadeMetaTrigger(
                    self._network, total_keyset
                )
                end_on_completion: bool = False
                for hook_func in global_hooks:
                    rval = hook_func(fake_network)
                    if (
                        rval
                        is interfaces.GlobalHookReturnValue.STOP_SHORTCIRCUIT
                    ):
                        return
                    if rval is interfaces.GlobalHookReturnValue.STOP:
                        end_on_completion = True
                if end_on_completion:
                    return

            continue

    def expand(
        self,
        max_recipes: int | None = None,
        mol_filter: interfaces.MolFilter | None = None,
        bundle_filter: interfaces.BundleFilter | None = None,
        recipe_filter: interfaces.RecipeFilter | None = None,
        recipe_ranker: interfaces.RecipeRanker | None = None,
        reaction_plan: metadata.RxnAnalysisStep
        | metadata.PropertyCompositor
        | metadata.ReactionFilterBase
        | metadata.LocalPropertyCalc
        | None = None,
        global_hooks: Sequence[interfaces.GlobalUpdateHook] | None = None,
        heap_size: int | None = None,
        beam_size: int | None = 1,
        batch_size: int | None = None,
        save_unreactive: bool = True,
        num_workers: int = 1,
    ) -> None:
        """
        Expand the network according to certain parameters.

        Parameters
        ----------
        max_recipes : int | None (default: None)
            Maximum number of recipes to evaluate.  Value of `None` indicates
            no limit to recipes.
        mol_filter : interfaces.MolFilter | None (default: None)
            Filters to use to stop specific molecules from being evaluated.
        bundle_filter : interfaces.BundleFilter | None (default: None)
            Filters to use to stop specific combinations of mols from being
            evaluated.
        recipe_filter : interfaces.RecipeFilter | None (default: None)
            Filters to use to stop recipes from being ranked or evaluated.  May
            be a composition of filters using basic logic operators.
        recipe_ranker : typing.Optional[RecipeRanker] (default: None)
            Function which evaluates recipes and assigns them a sortable rank.
        reaction_plan : metadata.RxnAnalysisStep | metadata.PropertyCompositor |
                        metadata.ReactionFilterBase |
                        metadata.LocalPropertyCalc | None (default: None)
            Metadata calculator and reaction filter.  Should be a "Reaction
            Analysis Flow" as detailed in examples.
        global_hooks : Sequence[interfaces.GlobalUpdateHook] | None
                       (default: None)
            Hook functions to run after a loop has completed.  If any return a
            boolean, expansion will terminate before running any more reactions.
            Return value of True terminates immediately, whereas return value
            of False terminates after remaining hooks also run.  `None`
            indicates expansion will continue until all recipes have been
            exhausted.

        Other Parameters
        ----------------
        heap_size : int | None (default: None)
            Maximum size of internal queue.  More values will use more memory,
            but be less likely to exhaust.  If queue is exhausted, all recipes
            and their ranks will be regenerated.  Value of `None` indicates no
            limit to queue size.
        beam_size : int | None (default: 1)
            Number of recipes to evaluate concurrently, as a set.  Products
            generated by recipes in this set are used to generate new recipes
            only after all recipes in this set have been evaluated.  Also
            affects parallelism (maximum # of cores used for reactions is
            limited by this number, though beam_size can be set to greater than
            the number of cores without any computational performance penalty).
            Value of `None` indicates that all recipes in the queue will be
            evaluated every cycle.  If not None, must be less than or equal to
            heap_size.
        batch_size : int | None (default: None)
            Maximum number of recipes to permit each parallel recipe ranking
            job to evaluate at once.  Affects how recipes are bundled.  Value
            of `None` indicates no limit on number of recipes.  Tune when using
            parallel processes.
        save_unreactive : bool (default: True)
            Store reactions rejected by reaction_plan in the network.  If
            False, reactions which are rejected will simply be deleted, along
            with their products, instead of stored.  Depending on the
            proportion of rejected reactions, this may save a lot of memory.
        num_workers : int (default: 1)
            Number of parallel workers to deploy.  May result in
            non-deterministic output.  0 indicates using all cores.
        """
        if num_workers < 0:
            raise ValueError(
                f"Number of workers must be >= 0 (was {num_workers})"
            )
        match num_workers:
            # use default cores
            case 0:
                with concurrent.futures.ProcessPoolExecutor(
                    max_workers=None
                ) as executor:
                    self._expand(
                        max_recipes=max_recipes,
                        mol_filter=mol_filter,
                        bundle_filter=bundle_filter,
                        recipe_filter=recipe_filter,
                        recipe_ranker=recipe_ranker,
                        reaction_plan=reaction_plan,
                        global_hooks=global_hooks,
                        heap_size=heap_size,
                        beam_size=beam_size,
                        batch_size=batch_size,
                        save_unreactive=save_unreactive,
                        executor=executor,
                    )
            # use no multiprocessing
            case 1:
                self._expand(
                    max_recipes=max_recipes,
                    mol_filter=mol_filter,
                    bundle_filter=bundle_filter,
                    recipe_filter=recipe_filter,
                    recipe_ranker=recipe_ranker,
                    reaction_plan=reaction_plan,
                    global_hooks=global_hooks,
                    heap_size=heap_size,
                    beam_size=beam_size,
                    batch_size=batch_size,
                    save_unreactive=save_unreactive,
                    executor=None,
                )
            # use fixed number of max cores
            case _:
                with concurrent.futures.ProcessPoolExecutor(
                    max_workers=num_workers
                ) as executor:
                    self._expand(
                        max_recipes=max_recipes,
                        mol_filter=mol_filter,
                        bundle_filter=bundle_filter,
                        recipe_filter=recipe_filter,
                        recipe_ranker=recipe_ranker,
                        reaction_plan=reaction_plan,
                        global_hooks=global_hooks,
                        heap_size=heap_size,
                        beam_size=beam_size,
                        batch_size=batch_size,
                        save_unreactive=save_unreactive,
                        executor=executor,
                    )
