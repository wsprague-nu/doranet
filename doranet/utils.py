"""Contains classes which define and implement utility functions."""

import collections.abc
import dataclasses
import itertools
import typing

from doranet import interfaces


@dataclasses.dataclass(frozen=True, slots=True)
class RxnTrackerDepthFirstNetwork(interfaces.RxnTrackerNetwork):
    """Implements RxnTrackerNetwork interface."""

    network: interfaces.ChemNetwork

    def _getchains(
        self,
        cur_gen_mols: collections.abc.Collection[interfaces.MolIndex],
        prev_gens_mols: typing.Optional[set[interfaces.MolIndex]] = None,
        prev_gens_rxns: typing.Optional[set[interfaces.RxnIndex]] = None,
        reagent_table: typing.Optional[
            collections.abc.Container[interfaces.MolIndex]
        ] = None,
        fail_on_unknown_reagent: bool = False,
        depth: typing.Optional[int] = None,
    ) -> collections.abc.Generator[
        list[frozenset[interfaces.RxnIndex]], None, None
    ]:
        network = self.network
        n_mols = len(network.mols)
        # n_rxns = len(network.rxns)
        if depth is not None:
            if depth <= 0:
                return
            new_depth = depth - 1
        else:
            new_depth = None
        if len(cur_gen_mols) == 0:
            yield []
            return
        if prev_gens_mols is None:
            prev_gens_mols = set()
        if prev_gens_rxns is None:
            prev_gens_rxns = set()
        if reagent_table is None:
            reagent_table = set()
        rxnsets: list[list[interfaces.RxnIndex]] = []
        for mol in cur_gen_mols:
            if mol in reagent_table:
                continue
            elif mol >= n_mols:
                rxnsets.append([])
                continue
            newrxnset: list[interfaces.RxnIndex] = [
                rxn
                for rxn in network.producers(mol)
                if rxn not in prev_gens_rxns
                and all(
                    mol not in prev_gens_mols
                    for mol in network.rxns[rxn].reactants
                )
            ]
            rxnsets.append(newrxnset)
        if not fail_on_unknown_reagent:
            rxnsets = [rxnset for rxnset in rxnsets if len(rxnset) > 0]
            if len(rxnsets) == 0:  # these two lines
                return  # close the loophole
        if len(rxnsets) == 0:
            yield []
            return
        tested_combos: set[frozenset[interfaces.RxnIndex]] = set()
        for rxntuple in itertools.product(*rxnsets):
            rxncombo: frozenset[interfaces.RxnIndex] = frozenset(rxntuple)
            if rxncombo in tested_combos:
                continue
            else:
                tested_combos.add(frozenset(rxncombo))
            required_reagents = set(
                mol
                for mol in itertools.chain(
                    *(network.rxns[rxn].reactants for rxn in rxncombo)
                )
                if mol not in reagent_table
            )
            if len(required_reagents) == 0:
                yield [rxncombo]
                continue
            for path in self._getchains(
                required_reagents,
                prev_gens_mols.union(cur_gen_mols),
                prev_gens_rxns.union(rxncombo),
                reagent_table,
                fail_on_unknown_reagent,
                new_depth,
            ):
                path.append(rxncombo)
                yield path

    def getParentChains(
        self,
        target: interfaces.MolIndex,
        reagent_table: collections.abc.Container[interfaces.MolIndex] = tuple(),
        fail_on_unknown_reagent: bool = False,
        max_depth: typing.Optional[int] = None,
    ) -> collections.abc.Generator[
        list[frozenset[interfaces.RxnIndex]], None, None
    ]:
        if fail_on_unknown_reagent and not reagent_table:
            raise ValueError(
                "reagent table must be specified if fail_on_unknown_reagent is "
                "True"
            )
        return (
            path
            for path in self._getchains(
                [target],
                reagent_table=reagent_table,
                fail_on_unknown_reagent=fail_on_unknown_reagent,
                depth=max_depth,
            )
        )


T = typing.TypeVar("T")


def logreduce(
    function: collections.abc.Callable[[T, T], T],
    iterable: collections.abc.Iterable[T],
) -> T:
    """
    logreduce(function, iterable[, initial]) -> value.

    Apply a function of two arguments (satisfying the associative property)
    cumulatively to the items of a sequence or iterable, from left to right, so
    as to reduce the iterable to a single value.  For example, reduce(lambda x,
    y: x+y, [1, 2, 3, 4, 5]) calculates (((1+2)+(3+4))+5).  If initial is
    present, it is placed before the items of the iterable in the calculation,
    and serves as a default when the iterable is empty.

    Memory: maximum of log2(n) objects of iterable stored vs 2 for reduce()
    Speed: log2(n) calls to function() vs n for reduce()
    """
    i = iter(iterable)
    try:
        r_val, stop = _logreduce(function, i, 0)
    except StopIteration:
        raise TypeError(
            "logreduce() of empty iterable with no initial value"
        ) from None
    n = 0
    try:
        while not stop:
            n += 1
            new_val, stop = _logreduce(function, i, n - 1)
            r_val = function(r_val, new_val)
    except StopIteration:
        ...
    return r_val


def _logreduce(
    f: collections.abc.Callable[[T, T], T],
    i: collections.abc.Iterator[T],
    n: int,
) -> tuple[T, bool]:
    if n == 0:
        x = next(i)
        try:
            y = next(i)
        except StopIteration:
            return x, True
        return f(x, y), False
    x, stop = _logreduce(f, i, n - 1)
    if stop:
        return x, True
    try:
        y, stop = _logreduce(f, i, n - 1)
    except StopIteration:
        return x, True
    return f(x, y), False


# def getFigures(
#     target_smiles: str,
#     mol_smiles: tuple[str, ...],
#     num_gens: int,
#     engine: interfaces.NetworkEngine,
#     rxn_lib: interfaces.ObjectLibrary[interfaces.RxnDatBase],
#     job_name: str,
# ):
#     """
#     Create and display figures of reaction pathways.

#     If the target is found at the end of the expansion, create figures of the
#     reaction routes. Each figure contains all generations of one route. Also
#     prints reaction operators, reactants, and products. Figures are displayed,
#     also saved in the working folder. Based on the example code for using
#     RxnTrackerDepthFirst(rxn_lib).  Currently need to have the font file
#     'OpenSans-Regular.ttf' in the working folder.  May not work due to legacy
#     code source.

#     Parameters
#     ----------
#     target_smiles : str
#         Target molecule SMILES.
#     mol_smiles : tuple[str,...]
#         SMILES of starting chemicals.
#     num_gens : int
#         Number of expansion generations performed.
#     engine : NetworkEngine
#         Network Engine used to create data objects.
#     rxn_lib : ObjectLibrary[RxnDatBase]
#         Reaction library where reactions are contained.
#     job_name: str
#         Name of the job.
#     """
#     doubletracker = RxnTrackerDepthFirst(rxn_lib)
#     route = 1
#     for chain in doubletracker.getParentChains(
#         engine.Mol(target_smiles).uid,
#         reagent_table=[engine.Mol(smiles).uid for smiles in mol_smiles],
#         max_depth=num_gens,
#     ):
#         g = 0
#         rxn_image_list = []
#         for gen in chain:
#             print(f"Generation {g}:")
#             for rxnid in gen:
#                 print(rxn_lib[rxnid])
#                 rxn_smiles_str = ""
#                 for molecule in rxn_lib[rxnid].reactants:
#                     if not isinstance(molecule, str):
#                         raise TypeError(
#                             f"""Identifier {molecule} of type {type(molecule)}
#                                 must be SMILES string."""
#                         )
#                     rxn_smiles_str += molecule
#                     rxn_smiles_str += "."
#                 rxn_smiles_str += ">>"
#                 for molecule in rxn_lib[rxnid].products:
#                     if not isinstance(molecule, str):
#                         raise TypeError(
#                             f"""Identifier {molecule} of type {type(molecule)}
#                                 must be SMILES string."""
#                         )
#                     rxn_smiles_str += molecule
#                     rxn_smiles_str += "."
#                 rxn1 = rdkit.Chem.rdchemreactions.ReactionFromSmarts(
#                     rxn_smiles_str, useSmiles=True
#                 )
#                 img1 = rdkit.Chem.Draw.ReactionToImage(rxn1)
#                 msg = f"Route {route}, Generation {g}"
#                 img_w, img_h = img1.size
#                 I1 = PIL.ImageDraw.Draw(img1)
#                 myFont = PIL.ImageFont.truetype("OpenSans-Regular.ttf", 25)
#                 new_box = I1.textbbox((0, 0), msg, font=myFont)
#                 I1.text(
#                     ((img_w - new_box[2]) / 2, 0),
#                     msg,
#                     font=myFont,
#                     fill="black",
#                 )
#                 rxn_image_list.append(img1)
#             g += 1
#         max_width = 0
#         total_hight = 0
#         for img in rxn_image_list:
#             total_hight += img.size[1]
#             max_width = max(max_width, img.size[0])
#         combined_image = PIL.Image.new("RGB",
#                                        (max_width, total_hight),
#                                        "white")
#         current_hight = 0
#         for img in rxn_image_list:
#             combined_image.paste(
#                 img, ((max_width - img.size[0]) // 2, current_hight)
#             )
#             current_hight += img.size[1]
#         IPython.display(combined_image)  # type: ignore
#         combined_image.save(f"{job_name} route {route}.png")
#         route += 1
