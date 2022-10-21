"""
Contains classes which define and implement utility functions.

Classes:

    RxnTracker
      RxnTrackerSingle*
      RxnTrackerDepthFirst*
"""

from abc import ABC, abstractmethod
from itertools import chain
from itertools import product as iterproduct
from typing import (
    Collection,
    Dict,
    Generator,
    Iterable,
    List,
    Optional,
    Sequence,
    Set,
)

import IPython.display
import PIL
import rdkit.Chem

from pickaxe_generic.containers import ObjectLibrary
from pickaxe_generic.datatypes import Identifier, RxnDatBase


class RxnTracker(ABC):
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

    @abstractmethod
    def getParentChains(
        self,
        target: Identifier,
        reagent_table: Sequence[Identifier] = tuple(),
        fail_on_unknown_reagent: bool = False,
    ) -> Iterable[Iterable[Iterable[Identifier]]]:
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


class RxnTrackerSingle(RxnTracker):
    """Implements RxnTracker interface; only compatible with reactions
    involving a single reactant and product.  DEVELOPMENT ONLY"""

    _mol_lookup: Dict[Identifier, list[Identifier]]

    def __init__(self, rxn_lib: ObjectLibrary[RxnDatBase]) -> None:
        self._mol_lookup = {}
        self._rxn_lib = rxn_lib
        for rxnid in rxn_lib.ids():
            product_mol = sorted(rxn_lib[rxnid].products)[0]
            if product_mol not in self._mol_lookup:
                self._mol_lookup[product_mol] = []
            self._mol_lookup[product_mol].append(rxnid)

    def _getchains(
        self, cur_mol: Identifier, cur_mols: Optional[Set[Identifier]] = None
    ):
        if cur_mols is None:
            cur_mols = {cur_mol}
        noReactions = True
        if cur_mol not in self._mol_lookup:
            yield list()
        else:
            for rxnid in self._mol_lookup[cur_mol]:
                reactant = sorted(self._rxn_lib[rxnid].reactants)[0]
                if reactant in cur_mols:
                    continue
                noReactions = False
                for rxnpath in self._getchains(
                    reactant, cur_mols.union({reactant})
                ):
                    rxnpath.append(rxnid)
                    yield rxnpath
            if noReactions:
                yield list()

    def getParentChains(
        self,
        target: Identifier,
        reagent_table: Sequence[Identifier] = None,
        fail_on_unknown_reagent: bool = None,
        max_depth: Optional[int] = None,
    ) -> Iterable[Iterable[Iterable[RxnDatBase]]]:
        if reagent_table is not None or fail_on_unknown_reagent is not None:
            raise NotImplementedError("Arguments besides target not supported.")
        return ([path] for path in self._getchains(target))


class RxnTrackerDepthFirst(RxnTracker):
    """Implements RxnTracker interface; stores lookups as a hash table within
    the object.  Will eventually deprecate this functionality when the
    ObjectLibrary interface is updated to include native search functionality.
    """

    _mol_lookup: Dict[Identifier, list[Identifier]]
    _rxn_lib: ObjectLibrary[RxnDatBase]

    def __init__(self, rxn_lib: ObjectLibrary[RxnDatBase]) -> None:
        self._mol_lookup = {}
        self._rxn_lib = rxn_lib
        for rxnid in rxn_lib.ids():
            for product_mol in rxn_lib[rxnid].products:
                if product_mol not in self._mol_lookup:
                    self._mol_lookup[product_mol] = []
                self._mol_lookup[product_mol].append(rxnid)

    def _getchains(
        self,
        cur_gen_mols: Collection[Identifier],
        prev_gens_mols: Set[Identifier] = None,
        prev_gens_rxns: Set[Identifier] = None,
        reagent_table: Optional[Iterable[Identifier]] = None,
        fail_on_unknown_reagent: bool = False,
        depth: Optional[int] = None,
    ) -> Generator[list[frozenset[Identifier]], None, None]:
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
        rxnsets: List[List[Identifier]] = []
        for mol in cur_gen_mols:
            if mol in reagent_table:
                continue
            elif mol not in self._mol_lookup:
                rxnsets.append([])
                continue
            newrxnset: List[Identifier] = [
                rxn
                for rxn in self._mol_lookup[mol]
                if rxn not in prev_gens_rxns
                and all(
                    mol not in prev_gens_mols
                    for mol in self._rxn_lib[rxn].reactants
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
        tested_combos = set()
        for rxntuple in iterproduct(*rxnsets):
            rxncombo: frozenset[Identifier] = frozenset(rxntuple)
            if rxncombo in tested_combos:
                continue
            else:
                tested_combos.add(frozenset(rxncombo))
            required_reagents = set(
                mol
                for mol in chain(
                    *(self._rxn_lib[rxn].reactants for rxn in rxncombo)
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
        target: Identifier,
        reagent_table: Sequence[Identifier] = tuple(),
        fail_on_unknown_reagent: bool = False,
        max_depth: Optional[int] = None,
    ) -> Generator[list[frozenset[Identifier]], None, None]:
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


def getFigures(target_smiles, mol_smiles, num_gens, engine, rxn_lib, job_name):
    """If the target is found at the end of the expansion, create figures of the
       reaction routes. Each figure contains all generations of one route. Also
       prints reaction operators, reactants, and products. Figures are displayed,
       also saved in the working folder. Based on the example code for using
       RxnTrackerDepthFirst(rxn_lib).
       Currently need to have the font file 'OpenSans-Regular.ttf' in the working
       folder.

    Parameters
        ----------
        target_smiles : str, target molecule smiles
        mol_smiles : tuple, smiles of starting chemicals
        num_gens : int, number of expansion generations
        engine : example: engine, created with engine = create_engine()
        rxn_lib : example: rxn_lib, created with mol_lib, op_lib, rxn_lib = engine.Libs()
        job_name: str, name of the job
    """
    doubletracker = RxnTrackerDepthFirst(rxn_lib)
    route = 1
    for chain in doubletracker.getParentChains(
        engine.Mol(target_smiles).uid,
        reagent_table=[engine.Mol(smiles).uid for smiles in mol_smiles],
        max_depth=num_gens,
    ):
        g = 0
        rxn_image_list = []
        for gen in chain:
            print(f"Generation {g}:")
            for rxnid in gen:
                print(rxn_lib[rxnid])
                rxn_smiles_str = ""
                for molecule in rxn_lib[rxnid].reactants:
                    rxn_smiles_str += molecule
                    rxn_smiles_str += "."
                rxn_smiles_str += ">>"
                for molecule in rxn_lib[rxnid].products:
                    rxn_smiles_str += molecule
                    rxn_smiles_str += "."
                rxn1 = rdkit.Chem.rdchemreactions.ReactionFromSmarts(
                    rxn_smiles_str, useSmiles=True
                )
                img1 = rdkit.Chem.Draw.ReactionToImage(rxn1)
                msg = f"Route {route}, Generation {g}"
                img_w, img_h = img1.size
                I1 = PIL.ImageDraw.Draw(img1)
                myFont = PIL.ImageFont.truetype("OpenSans-Regular.ttf", 25)
                new_box = I1.textbbox((0, 0), msg, font=myFont)
                I1.text(
                    ((img_w - new_box[2]) / 2, 0),
                    msg,
                    font=myFont,
                    fill="black",
                )
                rxn_image_list.append(img1)
            g += 1
        max_width = 0
        total_hight = 0
        for img in rxn_image_list:
            total_hight += img.size[1]
            max_width = max(max_width, img.size[0])
        combined_image = PIL.Image.new("RGB", (max_width, total_hight), "white")
        current_hight = 0
        for img in rxn_image_list:
            combined_image.paste(
                img, ((max_width - img.size[0]) // 2, current_hight)
            )
            current_hight += img.size[1]
        IPython.display(combined_image)
        combined_image.save(f"{job_name} route {route}.png")
        route += 1
