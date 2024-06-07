"""Post-processing code for DORAnet."""

import collections.abc
import copy
import dataclasses
import io
import json
import math
import textwrap
import time
import typing
from collections import deque
from collections.abc import Iterable

import matplotlib.pyplot as plt
import networkx as nx
import networkx.exception as nxe
import numpy as np
import rdkit.Chem.rdmolfiles
import rdkit.Chem.rdmolops
from PIL import Image, ImageChops, ImageDraw, ImageFont
from rdkit import Chem
from rdkit.Chem import Descriptors, Draw

import doranet as dn
from doranet import interfaces, metadata

try:
    from pgthermo.properties import Hf
except ImportError:
    print("pgthermo not installed. All Hf will be 0")

    def Hf(SMILES):
        return 0


import re
from multiprocessing import Pool
from pathlib import Path

from rdkit.Chem.rdMolDescriptors import CalcMolFormula

from doranet.modules.synthetic.Reaction_Smarts_Forward import op_smarts

# post-processing for doranet network files

reaxys_not_supported: Iterable[str] = (
    "O",
    "O=O",
    "[H][H]",
    "Br",
    "[Br][Br]",
    "Cl",
    "[Cl][Cl]",
    "O=S(O)O",
    "N",
    "O=S(=O)(O)O",
    "O=NO",
    "N#N",
    "O=[N+]([O-])O",
    "NO",
    "C#N",
    "S",
    "O=S=O",
    "[C-]#[O+]",
    "O=C=O",
)
temp_ms = []
for i in reaxys_not_supported:  # clean up
    temp_ms.append(
        rdkit.Chem.rdmolfiles.MolToSmiles(
            rdkit.Chem.rdmolfiles.MolFromSmiles(i)
        )
    )
reaxys_not_supported = set(temp_ms)


def clean_SMILES(smiles):
    mol = rdkit.Chem.rdmolfiles.MolFromSmiles(smiles)
    rdkit.Chem.rdmolops.RemoveStereochemistry(mol)
    cpd_smiles = rdkit.Chem.rdmolfiles.MolToSmiles(mol)
    return cpd_smiles


def pretreat_networks(
    total_generations,
    starters,
    helpers,
    network_forward_name=False,
    network_forward_obj=False,
    network_retro_name=False,
    network_retro_obj=False,
    job_name="default_job_name",
    remove_pure_helpers_rxns=True,
    sanitize=True,
    calculate_dH=True,
    transform_enols_flag=True,
):
    """
    Load pickaxe network files and unpack them into a json file.

    Will remove rxn if it's between helpers
    Will remove duplicates, remove unconnected rxns

    Parameters
    ----------
                    network_forward_name: A string
                    network_retro_name:  A string

    Returns
    -------
                    Create a json file with reactions as strings
    """
    engine = dn.create_engine()
    print("Loading network files, it may take a while for large networks")

    if network_forward_name:
        print("loading network_forward file")
        network_forward = engine.network_from_file(network_forward_name)
        print(
            "number of reations in network_forward", len(network_forward.rxns)
        )
    if network_forward_obj:
        print("use network_forward from memory")
        network_forward = network_forward_obj
        print(
            "number of reations in network_forward", len(network_forward.rxns)
        )

    if network_retro_name:
        print("loading network_retro file")
        network_retro = engine.network_from_file(network_retro_name)
        print("number of reactions in network_retro:", len(network_retro.rxns))
    if network_retro_obj:
        print("use network_retro from memory")
        network_retro = network_retro_obj
        print("number of reactions in network_retro:", len(network_retro.rxns))

    print("generating reaction strings")
    whole_rxns_list = list()
    H_dict = dict()  # save calculated Hf

    def ngHCal(
        _mol_num, network_name
    ):  # calculate Hf of all moleclues in network
        _smiles = network_name.mols[_mol_num].uid
        if _smiles in H_dict:
            return H_dict[_smiles]
        else:
            _enthalpy_f = Hf(_smiles)
            if _enthalpy_f is None:
                print("Hf error", _smiles, _enthalpy_f)
            H_dict[_smiles] = _enthalpy_f
            return _enthalpy_f

    def RxnHCal(
        operatorNum, reas, pros, network_name, dirction, return_zero, rxn_type
    ):  # calculate dH of all reactions in network and generate a string
        if dirction == "f":
            reactants_stoi = network_name.ops.meta(
                operatorNum, keys=["reactants_stoi"]
            )["reactants_stoi"]
            products_stoi = network_name.ops.meta(
                operatorNum, keys=["products_stoi"]
            )["products_stoi"]
        if dirction == "r":
            reactants_stoi = network_name.ops.meta(
                operatorNum, keys=["products_stoi"]
            )["products_stoi"]
            products_stoi = network_name.ops.meta(
                operatorNum, keys=["reactants_stoi"]
            )["reactants_stoi"]
        Hcorr = network_name.ops.meta(
            operatorNum, keys=["enthalpy_correction"]
        )["enthalpy_correction"]
        dH = 0.0
        if rxn_type == "Enzymatic":
            dH = -9999
        elif rxn_type == "Catalytic" and not return_zero:
            for idx, mol in enumerate(pros):
                dH = dH + ngHCal(mol, network_name) * products_stoi[idx]
            for idx, mol in enumerate(reas):
                dH = dH - ngHCal(mol, network_name) * reactants_stoi[idx]
            if Hcorr is not None:
                dH = dH + Hcorr
        dH = round(dH, 4)
        return str(dH) + "$" + str(reactants_stoi) + "$" + str(products_stoi)

    if network_forward_name or network_forward_obj:
        for rxn in network_forward.rxns:
            reactants = rxn.reactants
            products = rxn.products
            operator = rxn.operator
            if (
                network_forward.ops.meta(operator, keys=["Reaction_type"])[
                    "Reaction_type"
                ]
                == "Enzymatic"
            ):
                dH = RxnHCal(
                    operator,
                    reactants,
                    products,
                    network_forward,
                    "f",
                    False,
                    "Enzymatic",
                )
            elif (
                network_forward.ops.meta(operator, keys=["Reaction_type"])[
                    "Reaction_type"
                ]
                == "Catalytic"
            ):
                if calculate_dH:
                    dH = RxnHCal(
                        operator,
                        reactants,
                        products,
                        network_forward,
                        "f",
                        False,
                        "Catalytic",
                    )
                else:
                    dH = RxnHCal(
                        operator,
                        reactants,
                        products,
                        network_forward,
                        "f",
                        True,
                        "Catalytic",
                    )
            rxn_string = str()
            for i in reactants:
                rxn_string = rxn_string + network_forward.mols[i].uid + "."
            rxn_string = rxn_string[:-1]
            rxn_string = (
                rxn_string
                + ">"
                + network_forward.ops.meta(operator, keys=["name"])["name"]
                + ">"
                + str(dH)
                + ">"
            )
            for i in products:
                rxn_string = rxn_string + network_forward.mols[i].uid + "."
            rxn_string = rxn_string[:-1]
            whole_rxns_list.append(rxn_string)

    if network_retro_name or network_retro_obj:
        for rxn in network_retro.rxns:
            reactants = rxn.products
            products = rxn.reactants
            operator = rxn.operator
            if (
                network_retro.ops.meta(operator, keys=["Reaction_type"])[
                    "Reaction_type"
                ]
                == "Enzymatic"
            ):
                dH = RxnHCal(
                    operator,
                    reactants,
                    products,
                    network_retro,
                    "r",
                    False,
                    "Enzymatic",
                )
            elif (
                network_retro.ops.meta(operator, keys=["Reaction_type"])[
                    "Reaction_type"
                ]
                == "Catalytic"
            ):
                if calculate_dH:
                    dH = RxnHCal(
                        operator,
                        reactants,
                        products,
                        network_retro,
                        "r",
                        False,
                        "Catalytic",
                    )
                else:
                    dH = RxnHCal(
                        operator,
                        reactants,
                        products,
                        network_retro,
                        "r",
                        True,
                        "Catalytic",
                    )
            rxn_string = str()
            for i in reactants:
                rxn_string = rxn_string + network_retro.mols[i].uid + "."
            rxn_string = rxn_string[:-1]
            rxn_string = (
                rxn_string
                + ">"
                + network_retro.ops.meta(operator, keys=["name"])["name"]
                + ">"
                + str(dH)
                + ">"
            )
            for i in products:
                rxn_string = rxn_string + network_retro.mols[i].uid + "."
            rxn_string = rxn_string[:-1]
            whole_rxns_list.append(rxn_string)

    def transform_enols(input_rxn_string):
        def smilesHfCal(smiles, return_Hf):
            if return_Hf:
                return Hf(smiles)
            else:
                return 0

        if input_rxn_string.split(">")[1] == "Keto-enol Tautomerization":
            return input_rxn_string
        to_return = ">".join(input_rxn_string.split(">")[:2]) + ">"
        old_dH = float(input_rxn_string.split(">")[2].split("$")[0])
        stoi = "$" + "$".join(input_rxn_string.split(">")[2].split("$")[1:])
        pros = input_rxn_string.split(">")[3].split(".")
        patt_enol = rdkit.Chem.rdmolfiles.MolFromSmarts("[C]=[C]-[OH]")
        rxn = rdkit.Chem.rdChemReactions.ReactionFromSmarts(
            "[C+0:1]=[C+0:2][O+0H:3]>>[*:1][*:2]=[*:3]"
        )
        pros_string = str()
        dH = 0
        for pro in pros:
            pro_mol = rdkit.Chem.rdmolfiles.MolFromSmiles(pro)
            if (
                len(pro_mol.GetSubstructMatches(patt_enol)) == 1
            ):  # ignore rare cases with multiple enol -OH
                comb = (pro_mol,)
                products_sets = rxn.RunReactants(comb)
                new_smiles = rdkit.Chem.rdmolfiles.MolToSmiles(
                    products_sets[0][0]
                )
                if (
                    smilesHfCal(new_smiles, calculate_dH) is not None
                ):  # new keto may not be supported by the enthalpy calculator
                    pros_string = pros_string + new_smiles + "."
                    dH = (
                        dH
                        + smilesHfCal(new_smiles, calculate_dH)
                        - smilesHfCal(pro, calculate_dH)
                    )
                else:
                    pros_string = pros_string + pro + "."
            else:
                pros_string = pros_string + pro + "."
        pros_string = pros_string[:-1]
        new_dH = old_dH + dH
        new_dH = round(new_dH, 4)
        return to_return + str(new_dH) + stoi + ">" + pros_string

    if transform_enols_flag:
        new_whole_rxns_list = list()
        for i in whole_rxns_list:
            new_whole_rxns_list.append(transform_enols(i))
        whole_rxns_list = new_whole_rxns_list

    data = set(whole_rxns_list)
    print("reaction strings generation finished")

    #### sanitization
    starters_set = set()
    for i in starters:
        starters_set.add(
            rdkit.Chem.rdmolfiles.MolToSmiles(
                rdkit.Chem.rdmolfiles.MolFromSmiles(i)
            )
        )

    helpers_set = set()
    for i in helpers:
        if (
            rdkit.Chem.rdmolfiles.MolToSmiles(
                rdkit.Chem.rdmolfiles.MolFromSmiles(i)
            )
            not in starters_set
        ):
            helpers_set.add(
                rdkit.Chem.rdmolfiles.MolToSmiles(
                    rdkit.Chem.rdmolfiles.MolFromSmiles(i)
                )
            )

    mol_smiles = starters_set.union(helpers_set)

    consumers_dict = dict()
    for rxn in data:
        reas = rxn.split(">")[0].split(".")
        for rea in reas:
            if rea not in consumers_dict:
                consumers_dict[rea] = {rxn}
            else:
                consumers_dict[rea].add(rxn)

    def find_consumers(_smiles):
        if _smiles in consumers_dict:
            return list(consumers_dict[_smiles])
        else:
            return []

    my_list = list()
    for i in data:
        reas = i.split(">")[0].split(".")
        name = i.split(">")[1]
        dH = i.split(">")[2].split("$")[0]
        rea_stoi = i.split(">")[2].split("$")[1]
        pro_stoi = i.split(">")[2].split("$")[2]
        pros = i.split(">")[3].split(".")

        helpers_reacting_flag = True  # check if the rxn is just between helpers
        for rea in reas:
            if rea not in helpers_set:
                helpers_reacting_flag = False

        # sort to remove dup later
        rea_stoi = eval(rea_stoi)
        pro_stoi = eval(pro_stoi)

        reas_sorted, rea_stoi = zip(
            *sorted(zip(reas, rea_stoi, strict=False)), strict=False
        )
        pros_sorted, pro_stoi = zip(
            *sorted(zip(pros, pro_stoi, strict=False)), strict=False
        )

        if len(reas) != len(reas_sorted) or len(pros) != len(pros_sorted):
            print("sorting error", i)

        midd = dH + "$" + str(rea_stoi) + "$" + str(pro_stoi)

        my_string = str()
        left = str()
        right = str()

        for j in reas_sorted:
            left = left + j + "."
        left = left[:-1]
        for j in pros_sorted:
            right = right + j + "."
        right = right[:-1]

        my_string = left + ">" + name + ">" + midd + ">" + right

        if left != right and (
            helpers_reacting_flag is False or remove_pure_helpers_rxns is False
        ):  # and float(dH) < 0:    # set dH threshold if needed
            my_list.append(my_string)

    my_list = list(set(my_list))

    # remove dead end reactions
    if sanitize:
        consumers_dict2 = dict()
        for rxn in my_list:
            reas = rxn.split(">")[0].split(".")
            for rea in reas:
                if rea not in consumers_dict2:
                    consumers_dict2[rea] = {rxn}
                else:
                    consumers_dict2[rea].add(rxn)

        def find_consumers2(_smiles):
            if _smiles in consumers_dict2:
                return list(consumers_dict2[_smiles])
            else:
                return []

        gen = 0
        good_mol_set = set()
        for mol in mol_smiles:
            good_mol_set.add(mol)

        while gen < total_generations:
            temp_set = set()
            for mol in good_mol_set:
                rxn_list = find_consumers2(mol)
                for rxn in rxn_list:
                    reas = rxn.split(">")[0].split(".")
                    good_rxn_flag = True
                    for rea in reas:
                        if rea not in good_mol_set:
                            good_rxn_flag = False
                    if good_rxn_flag is True:
                        pros = rxn.split(">")[3].split(".")
                        for pro in pros:
                            temp_set.add(pro)
            good_mol_set.update(temp_set)
            gen += 1

        reduced_rxns = set()
        for rxn in my_list:
            reas = rxn.split(">")[0].split(".")
            good_rxn_flag = True
            for rea in reas:
                if rea not in good_mol_set:
                    good_rxn_flag = False
            if good_rxn_flag is True:
                reduced_rxns.add(rxn)

        my_list = list(reduced_rxns)
        print("unconnected reactions removed")

    print("total number of reactions after pretreatment:", len(my_list))

    with open(
        f"{job_name}_network_pretreated.json", "w", encoding="utf-8"
    ) as convert_file:
        convert_file.write(json.dumps(my_list))


def pathway_finder(
    starters,
    helpers,
    target,
    search_depth,
    max_num_rxns,
    min_rxn_atom_economy,
    job_name="default_job_name",
    consider_name_difference=True,
):
    """
    Search for pathways within pretreated reaction network.

    Parameters
    ----------
                    starters: a list of SMILES strings
                    helpers: A dict
                    target: A SMILES string
                    search_depth: A int
                    max_num_rxns: A int
                    min_rxn_atom_economy: A float
                    job_name: A string

    Returns
    -------
                    Create a txt file with found pathways
                    Create a txt file with all reactions in pathways, can be
                      used as input file for reaxys batch query

    """
    start_time = time.time()
    target_smiles = rdkit.Chem.rdmolfiles.MolToSmiles(
        rdkit.Chem.rdmolfiles.MolFromSmiles(target)
    )
    generation = search_depth  # number of generations to backtrack
    pathway_max_length = max_num_rxns  # max number of rxns allowed in a pathway
    min_atom_economy = min_rxn_atom_economy  # 0~1, apply to each rxn in path
    # min_atom_economy is not suitable if a rxn produces 2 products and both are
    # useful

    mol_smiles = tuple(helpers.keys()) + tuple(starters)
    temp_ms = list()
    for i in mol_smiles:  # clean up
        temp_ms.append(
            rdkit.Chem.rdmolfiles.MolToSmiles(
                rdkit.Chem.rdmolfiles.MolFromSmiles(i)
            )
        )
    mol_smiles = temp_ms

    with open(
        f"{job_name}_network_pretreated.json", encoding="utf-8"
    ) as f:  # open json file
        data = json.load(f)
    starters_set = set(mol_smiles)

    print("total number of reactions", len(data))

    rxn_dict = (
        dict()
    )  # { rxn_idx: {name: xx   reas: []  pros: []  mid: dH$()$()   } }
    # rxn_idx = 0
    for rxn_idx, rxn in enumerate(data):
        pros = rxn.split(">")[3].split(".")
        reas = rxn.split(">")[0].split(".")
        name = rxn.split(">")[1]
        mid = rxn.split(">")[2]
        # mid = str(round(float(rxn.split(">")[2].split("$")[0]), 4)) +
        rxn_dict[rxn_idx] = {
            "name": name,
            "reas": reas,
            "pros": pros,
            "mid": mid,
        }
        # rxn_idx += 1

    producers_dict = dict()
    producers_dict_len = dict()
    consumers_dict = dict()

    # economy cal overhead: 61.66 seconds for 532 k rxns
    for rxn_idx in rxn_dict:
        # build consumers_dict, producers_dict for
        # find_consumers(), find_producers()
        products = rxn_dict[rxn_idx]["pros"]
        # also build atom econmy in producer_dict:
        # {smiles: {rxn_idx: economy, rxn_idx2: economy ...}}
        reas = rxn_dict[rxn_idx]["reas"]
        pro_stoi = eval(rxn_dict[rxn_idx]["mid"].split("$")[2])
        products_weight = dict()  # {pro: weight}
        total_weight = 0
        for _idx, rea in enumerate(reas):
            if rea not in consumers_dict:
                consumers_dict[rea] = {rxn_idx: 0}
            else:
                consumers_dict[rea][rxn_idx] = 0
        for idx, pro in enumerate(products):
            to_add = (
                Descriptors.MolWt(rdkit.Chem.rdmolfiles.MolFromSmiles(pro))
                * pro_stoi[idx]
            )
            if pro not in products_weight:
                products_weight[pro] = to_add
                total_weight += to_add
            else:  # if a rxn produces multiple same products
                products_weight[pro] += to_add
                total_weight += to_add
            if pro not in producers_dict:
                producers_dict[pro] = {rxn_idx: 0}
            else:
                producers_dict[pro][rxn_idx] = 0
        for pro in products:
            producers_dict[pro][rxn_idx] = products_weight[pro] / total_weight

    for smiles in producers_dict:  # number of producers, used to reduce forks
        producers_dict_len[smiles] = len(producers_dict[smiles])

    distance_to_startrs_dict = dict()  # min distance to starters
    for smiles in starters_set:
        distance_to_startrs_dict[smiles] = 0

    def find_consumers(_smiles):
        if _smiles in consumers_dict:
            return list(consumers_dict[_smiles].keys())
        else:
            return []

    def find_producers_num(_smiles):
        if _smiles in producers_dict_len:
            return producers_dict_len[_smiles]
        else:
            return 0

    def find_producers(_smiles):
        if _smiles in producers_dict:
            return list(producers_dict[_smiles].keys())  # [int1, int2, ...]
        else:
            return []

    # find min distance to starters
    expanded_nodes = copy.deepcopy(starters_set)
    # set to the generation number during pickaxe expansion,
    # could be more but not less
    expand_gen = generation
    i = 0
    while i < expand_gen:
        temp_set = set()
        for smiles in expanded_nodes:
            consumers = find_consumers(smiles)
            for consumer in consumers:
                valid_consumer = True
                pros = rxn_dict[consumer]["pros"]
                reas = rxn_dict[consumer]["reas"]
                for rea in reas:
                    if rea not in expanded_nodes:
                        valid_consumer = False
                if valid_consumer is True:
                    for pro in pros:
                        temp_set.add(pro)
                        if pro not in distance_to_startrs_dict:
                            distance_to_startrs_dict[pro] = i + 1
        i += 1
        expanded_nodes.update(temp_set)

    def find_min_distance(_smiles):
        if _smiles in distance_to_startrs_dict:
            return distance_to_startrs_dict[_smiles]  # return int
        else:
            return 999  # not produced by starters (retro network)

    work_list = deque()
    # [finished:1,dead end:2,to expand:0,
    # [(mol,int),(mol,int)],rxn_idx1,rxn_idx2,...]
    work_list.append([0, [(target_smiles, 0)]])
    gen_limit = generation

    keep_expanding_flag = True
    print("searching for pathways.")
    print("if it is taking too long, try adjusting pruning parameters")
    print()
    no_pathway = False
    while keep_expanding_flag is True:
        try:
            first_path = work_list.popleft()
            if first_path[0] == 0:
                nodes = first_path[1]
                node_idx = 0
                node = nodes[node_idx][0]  # smiles of expanding node
                node_gen = nodes[node_idx][1]
                producers = find_producers(node)  # list of int
                path_len = len(
                    first_path[2:]
                )  # number of rxns in current pathway
                for producer in producers:  # int, index of rxn
                    if producers_dict[node][producer] > min_atom_economy:
                        to_append = list(first_path)  # copy.deepcopy or list?
                        to_append[0] = 2  # default dead end, change later
                        to_append.append(
                            producer
                        )  # add producer rxn to this pathway
                        to_append[1] = to_append[1][1:]  # remove expended node
                        deadend_flag = False
                        for rea in rxn_dict[producer][
                            "reas"
                        ]:  # check if a lower gen mol already included
                            flag1 = True
                            if (
                                node_gen + 1 == gen_limit
                                and rea not in starters_set
                            ) or find_min_distance(
                                rea
                            ) + node_gen + 1 > gen_limit:
                                deadend_flag = True
                                break
                            for node_idx2, t in enumerate(to_append[1]):
                                if t[0] == rea:
                                    to_append[1][node_idx2] = (
                                        rea,
                                        max(node_gen + 1, t[1]),
                                    )  # keep the higer gen number, not lower
                                    flag1 = False
                            if flag1 is True and rea not in starters_set:
                                to_append[1].append((rea, node_gen + 1))
                        if deadend_flag is False:
                            Y_num_of_producers = list()
                            for _nodes in to_append[1]:
                                Y_num_of_producers.append(
                                    find_producers_num(_nodes[0])
                                )  # sort nodes by its producers, faster
                            to_append[1] = [
                                x
                                for _, x in sorted(
                                    zip(
                                        Y_num_of_producers,
                                        to_append[1],
                                        strict=False,
                                    )
                                )
                            ]
                            if not to_append[1]:
                                to_append[0] = 1  # finished
                                work_list.append(to_append)
                            if (
                                to_append[1]
                                and path_len + 1 < pathway_max_length
                            ):
                                to_append[0] = 0  # to expand next time
                                work_list.appendleft(to_append)
            else:
                work_list.appendleft(first_path)
                keep_expanding_flag = False
        except IndexError:
            no_pathway = True
            print("no pathway found! try adjusting pruning parameters.")
            break

    def index_to_path(_rxn_idx):
        result_string = str()
        for rea in rxn_dict[_rxn_idx]["reas"]:
            result_string = result_string + rea + "."
        result_string = result_string[:-1]
        result_string = (
            result_string
            + ">"
            + rxn_dict[_rxn_idx]["name"]
            + ">"
            + rxn_dict[_rxn_idx]["mid"]
            + ">"
        )
        for pro in rxn_dict[_rxn_idx]["pros"]:
            result_string = result_string + pro + "."
        result_string = result_string[:-1]
        return result_string

    to_save = []
    all_rxns = set()

    for i in list(work_list):
        path = i[2:]
        path = list(set(path))
        # path.sort()
        path_str = []
        for j in path:
            path_str.append(index_to_path(j))
        path_str.sort()
        if path_str not in to_save:
            to_save.append(path_str)
            # to_save_index.append(path)

    def get_rxn_smiles(_rxn_idx):
        result_string = str()
        for rea in rxn_dict[_rxn_idx]["reas"]:
            result_string = result_string + rea + "."
        result_string = result_string[:-1]
        result_string = result_string + ">>"
        for pro in rxn_dict[_rxn_idx]["pros"]:
            result_string = result_string + pro + "."
        result_string = result_string[:-1]
        return result_string

    def find_cyc(
        _pathway, _starters, _helpers
    ):  # return bool, if there're cycles in the pathway
        G = nx.DiGraph()
        # rxn_node = 0
        for rxn_node, rxn in enumerate(
            _pathway["SMILES"]
        ):  # add nodes and edges to graph
            reas = rxn.split(">>")[0].split(".")
            pros = rxn.split(">>")[1].split(".")
            G.add_node(rxn_node)
            for rea in reas:  # add nodes and edges for non-reagents
                if rea not in _helpers:
                    G.add_node(rea)
                    G.add_edge(rea, rxn_node)
                    for pro in pros:
                        if pro not in _helpers and pro not in _starters:
                            G.add_node(pro)
                            G.add_edge(rxn_node, pro)
            # rxn_node += 1
        try:
            nx.find_cycle(G, orientation="original")
            return True
        except nxe.NetworkXNoCycle:
            return False

    print("pathway search finished, removing loops if there's any.")
    end_time = time.time()
    elapsed_time = end_time - start_time
    print()
    print("time used:", "{:.2f}".format(elapsed_time), "seconds")

    # remove duplicates due to rxn name differences in forward
    # and retro reaction rules, remove pathways with cycles
    if no_pathway is False:
        pathways_list = list()
        # [{final_score:,eco:,pathy_by:,inter_by:{},SMILES:[],Nmaes:[],dH:[]}]
        pathway_num = 1
        for path_list in to_save:
            temp_dict = dict()
            temp_dict["number"] = str(pathway_num)
            pathway_num += 1
            temp_dict["final_score"] = "empty"
            temp_dict["atomic_economy"] = "empty"
            temp_dict["pathway_by-product"] = "empty"
            temp_dict["intermediate_by-product"] = "empty"

            clean_smiles_list = list()
            clean_name_list = list()
            clean_dH_list = list()
            stoi_list = list()
            for each_rxn_string in path_list:
                clean_smiles_list.append(
                    each_rxn_string.split(">")[0]
                    + ">>"
                    + each_rxn_string.split(">")[3]
                )
                clean_name_list.append(each_rxn_string.split(">")[1])
                clean_dH_list.append(
                    each_rxn_string.split(">")[2].split("$")[0]
                )

                stoi_list.append(
                    each_rxn_string.split(">")[2].split("$")[1]
                    + "$"
                    + each_rxn_string.split(">")[2].split("$")[2]
                )

            temp_dict["SMILES"] = clean_smiles_list
            temp_dict["name"] = clean_name_list
            temp_dict["enthalpy"] = clean_dH_list
            temp_dict["stoi"] = stoi_list

            pathways_list.append(temp_dict)

        path_rxnlist_set = list()
        path_dH_set = list()

        with open(
            f"{job_name}_pathways.txt", "w", encoding="utf-8"
        ) as f_result:
            new_ranking = 1

            clean_starters = set()
            for i in starters:
                clean_starters.add(
                    rdkit.Chem.rdmolfiles.MolToSmiles(
                        rdkit.Chem.rdmolfiles.MolFromSmiles(i)
                    )
                )
            clean_helpers = set()
            for i in helpers:
                clean_helpers.add(
                    rdkit.Chem.rdmolfiles.MolToSmiles(
                        rdkit.Chem.rdmolfiles.MolFromSmiles(i)
                    )
                )

            for path in pathways_list:
                if (
                    consider_name_difference
                    or path["SMILES"] not in path_rxnlist_set
                    or path["enthalpy"] not in path_dH_set
                ) and find_cyc(path, clean_starters, clean_helpers) is False:
                    f_result.write("pathway number " + str(new_ranking) + "\n")
                    f_result.write("place holder " + path["final_score"] + "\n")
                    f_result.write(
                        "place holder " + path["atomic_economy"] + "\n"
                    )
                    f_result.write(
                        "place holder " + path["pathway_by-product"] + "\n"
                    )
                    f_result.write(
                        "reaction SMILES stoichiometry "
                        + str(path["stoi"])
                        + "\n"
                    )
                    f_result.write(
                        "reaction SMILES, name, and enthalpy:" + "\n"
                    )

                    for i in path["SMILES"]:
                        f_result.write(i + "\n")
                    for i in path["name"]:
                        f_result.write(i + "\n")
                    for i in path["enthalpy"]:
                        f_result.write(i + "\n")

                    f_result.write("\n")
                    new_ranking += 1
                path_rxnlist_set.append(path["SMILES"])
                path_dH_set.append(path["enthalpy"])

        min_num = 999
        max_num = 0
        for i in pathways_list:
            min_num = min(min_num, len(i["SMILES"]))
            max_num = max(max_num, len(i["SMILES"]))
            for j in i["SMILES"]:
                all_rxns.add(j)
        print()
        print(new_ranking - 1, "pathways found!")
        print(len(all_rxns), "reactions found in all pathways")
        print("min number of reactions in a pathway:", min_num)
        print("max number of reactions in a pathway:", max_num)
        print()
        print("pruning parameters for this search:")
        print("search_depth:", search_depth)
        print("max_num_rxns:", max_num_rxns)
        print("min_rxn_atom_economy:", min_rxn_atom_economy)

        # generate reaxys batch query file   https://service.elsevier.com/app/answers/detail/a_id/26151/supporthub/reaxys/p/10958/#:~:text=Reaxys%20searches%20batch%20queries%20sequentially,individual%20query%20in%20the%20batch.

        rxn_set = set()

        def SMILES_rm_keku(_SMILES):  # remove reaxys unsuported mols,
            # kekulize mols, return a reaction SMILEs
            reas = _SMILES.split(">>")[0].split(".")
            pros = _SMILES.split(">>")[1].split(".")
            new_rxn = str()
            reas_set = set()
            pros_set = set()
            for rea in reas:
                if rea not in reaxys_not_supported and rea not in reas_set:
                    new_rxn = (
                        new_rxn
                        + rdkit.Chem.rdmolfiles.MolToSmiles(
                            rdkit.Chem.rdmolfiles.MolFromSmiles(rea),
                            kekuleSmiles=True,
                        )
                        + "."
                    )
                    reas_set.add(rea)
            new_rxn = new_rxn[:-1] + ">>"
            for pro in pros:
                if pro not in reaxys_not_supported and pro not in pros_set:
                    new_rxn = (
                        new_rxn
                        + rdkit.Chem.rdmolfiles.MolToSmiles(
                            rdkit.Chem.rdmolfiles.MolFromSmiles(pro),
                            kekuleSmiles=True,
                        )
                        + "."
                    )
                    pros_set.add(pro)
            if new_rxn[-1] == ".":
                new_rxn = new_rxn[:-1]
            return new_rxn

        for path in pathways_list:
            SMILES = path["SMILES"]
            for rxn in SMILES:
                rxn_set.add(SMILES_rm_keku(rxn))

        max_rxns_per_batch = 1000  # reaxys limit 1000
        rxn_set = list(rxn_set)
        rxn_set.sort()

        # if len(rxn_set) < max_rxns_per_batch:
        if True:
            with open(
                f"{job_name}_reaxys_batch_query.txt", "w", encoding="utf-8"
            ) as f_result:
                for rxn in rxn_set:
                    print_flag = True
                    if print_flag is True:
                        f_result.write("SMILES='" + rxn + "'")
                        f_result.write("\n")

        # else:
        if len(rxn_set) > max_rxns_per_batch:
            num_files = math.ceil(len(rxn_set) / max_rxns_per_batch)
            # print(num_files, "files generated")

            filenames = [
                f"{job_name}_reaxys_batch_query_part_" + str(i)
                for i in range(1, num_files + 1)
            ]
            output_dict = dict()

            clean_list = list(rxn_set)
            for idx, filename in enumerate(filenames):
                if idx < num_files - 1:
                    output_dict[filename] = clean_list[
                        idx * max_rxns_per_batch : (idx + 1)
                        * max_rxns_per_batch
                    ]
                else:
                    output_dict[filename] = clean_list[
                        idx * max_rxns_per_batch :
                    ]

            for key in output_dict:
                with open(key + ".txt", "w", encoding="utf-8") as f_result:
                    for query in output_dict[key]:
                        f_result.write("SMILES='" + query + "'")
                        f_result.write("\n")


@typing.final
@dataclasses.dataclass(frozen=True, slots=True)
class EnthalpyCalculator(
    metadata.MolPropertyCalc[float]
):  # Calculate Hf for molecules
    Enthalpy_key: collections.abc.Hashable

    @property
    def key(self) -> collections.abc.Hashable:
        return self.Enthalpy_key

    @property
    def meta_required(self) -> interfaces.MetaKeyPacket:
        return interfaces.MetaKeyPacket(molecule_keys={self.Enthalpy_key})

    @property
    def resolver(self) -> metadata.MetaDataResolverFunc[float]:
        return metadata.TrivialMetaDataResolverFunc

    def __call__(
        self,
        data: interfaces.DataPacketE[interfaces.MolDatBase],
        prev_value: typing.Optional[float] = None,
    ) -> typing.Optional[float]:
        if prev_value is not None:
            return prev_value
        item = data.item
        if data.meta is not None and self.Enthalpy_key in data.meta:
            return None
        if not isinstance(item, interfaces.MolDatRDKit):
            raise NotImplementedError(
                f"Not been implemented for molecule type {type(item)}"
            )
        _enthalpy_f = Hf(item.uid)
        if _enthalpy_f is None:
            # print("None Enthalpy returned by molecule:", item.uid)
            # tell user which molecule is causing problems
            return float("nan")
        return _enthalpy_f


@typing.final
@dataclasses.dataclass(frozen=True)
class MaxEnthalpyFilter(metadata.ReactionFilterBase):
    __slots__ = ("max_H", "H_key")
    max_H: float
    H_key: collections.abc.Hashable

    def __call__(self, recipe: interfaces.ReactionExplicit) -> bool:
        dH = 0.0
        if recipe.operator.meta is None:
            return False
        for idx, mol in enumerate(recipe.products):
            if mol.meta is None or mol.meta[self.H_key] == float("nan"):
                return False
            dH = (
                dH
                + mol.meta[self.H_key]
                * recipe.operator.meta["products_stoi"][idx]
            )
        for idx, mol in enumerate(recipe.reactants):
            if mol.meta is None or mol.meta[self.H_key] == float("nan"):
                return False
            dH = (
                dH
                - mol.meta[self.H_key]
                * recipe.operator.meta["reactants_stoi"][idx]
            )
        if recipe.operator.meta["enthalpy_correction"] is not None:
            dH = dH + recipe.operator.meta["enthalpy_correction"]
        if dH / recipe.operator.meta["number_of_steps"] < self.max_H:
            return True
        return False

    @property
    def meta_required(self) -> interfaces.MetaKeyPacket:
        return interfaces.MetaKeyPacket(
            molecule_keys={self.H_key},
            operator_keys={
                "reactants_stoi",
                "products_stoi",
                "enthalpy_correction",
                "number_of_steps",
            },
        )


@typing.final
@dataclasses.dataclass(frozen=True)
class Ring_Issues_Filter(metadata.ReactionFilterBase):
    def __call__(self, recipe: interfaces.ReactionExplicit) -> bool:
        if (
            recipe.operator.meta is not None
            and recipe.operator.meta["ring_issue"] is True
            and recipe.operator.meta["enthalpy_correction"] is None
        ):
            reactants_dict: dict[str, int | float] = dict()
            products_dict: dict[str, int | float] = dict()
            pattern = r"([A-Z][a-z]*)(\d*)"
            for idx, mol in enumerate(recipe.reactants):
                if not isinstance(mol.item, dn.interfaces.MolDatRDKit):
                    raise NotImplementedError(
                        f"""Ring_Issues_Filter not implemented for molecule type
                            {type(mol.item)}"""
                    )
                smiles = CalcMolFormula(mol.item.rdkitmol)
                matches: list[tuple[str, str]] = re.findall(pattern, smiles)
                for match in matches:
                    element, count_s = match
                    count = int(count_s) if count_s else 1
                    reactants_dict[element] = (
                        reactants_dict.get(element, 0)
                        + count * recipe.operator.meta["reactants_stoi"][idx]
                    )
            for idx, mol in enumerate(recipe.products):
                if not isinstance(mol.item, dn.interfaces.MolDatRDKit):
                    raise NotImplementedError(
                        f"""Ring_Issues_Filter not implemented for molecule type
                            {type(mol.item)}"""
                    )
                if (
                    "." in mol.item.smiles
                ):  # if there're fragments in a mol, indicates invalid rxn
                    return False
                smiles = CalcMolFormula(mol.item.rdkitmol)
                matches = re.findall(pattern, smiles)
                for match in matches:
                    element, count_s = match
                    count = int(count_s) if count_s else 1
                    products_dict[element] = (
                        products_dict.get(element, 0)
                        + count * recipe.operator.meta["products_stoi"][idx]
                    )
            if reactants_dict != products_dict:
                return False
        return True

    @property
    def meta_required(self) -> interfaces.MetaKeyPacket:
        return interfaces.MetaKeyPacket(
            operator_keys={
                "reactants_stoi",
                "products_stoi",
                "ring_issue",
                "enthalpy_correction",
            }
        )


@typing.final
@dataclasses.dataclass(frozen=True, slots=True)
class No_helpers_reaction_Filter(metadata.ReactionFilterBase):  # mol.item.uid
    intermediate_smiles: str

    def __call__(self, recipe: interfaces.ReactionExplicit) -> bool:
        for mol in recipe.reactants:
            if mol.item.uid == self.intermediate_smiles:
                return True
        return False


def Byproduct_index(path_idx, _input_list):
    my_num_gens = 1
    mol_smiles = _input_list
    engine = dn.create_engine()
    network = engine.new_network()

    for smiles in mol_smiles:
        network.add_mol(engine.mol.rdkit(smiles))

    for smarts in op_smarts:
        if smarts.kekulize_flag is False:
            network.add_op(
                engine.op.rdkit(smarts.smarts),
                meta={
                    "name": smarts.name,
                    "reactants_stoi": smarts.reactants_stoi,
                    "products_stoi": smarts.products_stoi,
                    "enthalpy_correction": smarts.enthalpy_correction,
                    "ring_issue": smarts.ring_issue,
                    "kekulize_flag": smarts.kekulize_flag,
                    "Retro_Not_Aromatic": smarts.Retro_Not_Aromatic,
                    "number_of_steps": smarts.number_of_steps,
                },
            )
        if smarts.kekulize_flag is True:
            network.add_op(
                engine.op.rdkit(
                    smarts.smarts,
                    kekulize=True,
                ),
                meta={
                    "name": smarts.name,
                    "reactants_stoi": smarts.reactants_stoi,
                    "products_stoi": smarts.products_stoi,
                    "enthalpy_correction": smarts.enthalpy_correction,
                    "ring_issue": smarts.ring_issue,
                    "kekulize_flag": smarts.kekulize_flag,
                    "Retro_Not_Aromatic": smarts.Retro_Not_Aromatic,
                    "number_of_steps": smarts.number_of_steps,
                },
            )

    strat = engine.strat.cartesian(network)

    H_calc = EnthalpyCalculator("Enthalpy_f")
    thermo_filter = MaxEnthalpyFilter(15, "Enthalpy_f")
    reaction_plan = Ring_Issues_Filter() >> H_calc >> thermo_filter

    strat.expand(
        num_iter=my_num_gens, reaction_plan=reaction_plan, save_unreactive=False
    )
    end_number = 0
    for mol in network.mols:
        if network.reactivity[network.mols.i(mol.uid)] is True:
            end_number += 1
    return (path_idx, end_number - len(_input_list))


def Inter_byproduct(
    path_idx, _intermediate, _smiles_list
):  # _intermediate included in _smiles_list
    my_num_gens = 1
    mol_smiles = _smiles_list
    engine = dn.create_engine()
    network = engine.new_network()

    for smiles in mol_smiles:
        network.add_mol(engine.mol.rdkit(smiles))

    for smarts in op_smarts:
        if smarts.kekulize_flag is False:
            network.add_op(
                engine.op.rdkit(smarts.smarts),
                meta={
                    "name": smarts.name,
                    "reactants_stoi": smarts.reactants_stoi,
                    "products_stoi": smarts.products_stoi,
                    "enthalpy_correction": smarts.enthalpy_correction,
                    "ring_issue": smarts.ring_issue,
                    "kekulize_flag": smarts.kekulize_flag,
                    "Retro_Not_Aromatic": smarts.Retro_Not_Aromatic,
                    "number_of_steps": smarts.number_of_steps,
                },
            )
        if smarts.kekulize_flag is True:
            network.add_op(
                engine.op.rdkit(
                    smarts.smarts,
                    kekulize=True,
                ),
                meta={
                    "name": smarts.name,
                    "reactants_stoi": smarts.reactants_stoi,
                    "products_stoi": smarts.products_stoi,
                    "enthalpy_correction": smarts.enthalpy_correction,
                    "ring_issue": smarts.ring_issue,
                    "kekulize_flag": smarts.kekulize_flag,
                    "Retro_Not_Aromatic": smarts.Retro_Not_Aromatic,
                    "number_of_steps": smarts.number_of_steps,
                },
            )

    strat = engine.strat.cartesian(network)

    H_calc = EnthalpyCalculator("Enthalpy_f")
    thermo_filter = MaxEnthalpyFilter(15, "Enthalpy_f")
    reaction_plan = (
        No_helpers_reaction_Filter(_intermediate)
        >> Ring_Issues_Filter()
        >> H_calc
        >> thermo_filter
    )

    strat.expand(
        num_iter=my_num_gens, reaction_plan=reaction_plan, save_unreactive=False
    )
    end_number = 0

    for mol in network.mols:
        if network.reactivity[network.mols.i(mol.uid)] is True:
            end_number += 1

    return (path_idx, _intermediate, end_number - len(_smiles_list))


# pathway ranking
def pathway_ranking(
    starters,
    helpers,
    target,
    weights,
    num_process,
    reaxys_result_name=False,
    job_name="default_job_name",
    cool_reactions=False,
):
    start_time = time.time()
    target_smiles = rdkit.Chem.rdmolfiles.MolToSmiles(
        rdkit.Chem.rdmolfiles.MolFromSmiles(target)
    )

    just_starters = set()
    for i in starters:  # clean up
        just_starters.add(
            rdkit.Chem.rdmolfiles.MolToSmiles(
                rdkit.Chem.rdmolfiles.MolFromSmiles(i)
            )
        )

    mol_smiles = tuple(helpers.keys()) + tuple(starters)
    temp_ms = []
    for i in mol_smiles:  # clean up
        temp_ms.append(
            rdkit.Chem.rdmolfiles.MolToSmiles(
                rdkit.Chem.rdmolfiles.MolFromSmiles(i)
            )
        )
    mol_smiles = temp_ms

    helpers = tuple(helpers.keys())
    temp_ms2 = []
    for i in helpers:  # clean up
        temp_ms2.append(
            rdkit.Chem.rdmolfiles.MolToSmiles(
                rdkit.Chem.rdmolfiles.MolFromSmiles(i)
            )
        )
    helpers = temp_ms2

    ## read pathway txt file
    with open(f"{job_name}_pathways.txt", encoding="utf-8") as f:
        # with open(f'{job_name}_pathways.txt', encoding = 'cp1252') as f:
        lines = f.readlines()

    clean_list = list()
    for i in lines:
        if i != "\n":
            clean_list.append(i.strip())

    pathways_list = list()
    # [{final_score:,eco:,pathy_by:,inter_by:{},SMILES:[], Nmaes:[],dH:[]}]
    pathway_marker = list()
    # pathway_num = 1
    for idx, i in enumerate(clean_list):
        if "pathway number" in i:
            pathway_marker.append(idx)

    for idx, marker in enumerate(pathway_marker):
        temp_dict = dict()
        temp_dict["stoi"] = eval(clean_list[marker + 4][30:])  # list of str

        if idx + 1 < len(pathway_marker):  # 1, 2, 3
            next_marker = pathway_marker[idx + 1]
        else:
            next_marker = len(clean_list)

        temp_dict["SMILES"] = clean_list[
            marker + 6 : marker + 6 + int((next_marker - (marker + 6)) / 3)
        ]
        temp_dict["name"] = clean_list[
            marker + 6 + int((next_marker - (marker + 6)) / 3) : marker
            + 6
            + int((next_marker - (marker + 6)) / 3) * 2
        ]
        temp_dict["enthalpy"] = clean_list[
            marker + 6 + int((next_marker - (marker + 6)) / 3) * 2 : marker
            + 6
            + int((next_marker - (marker + 6)) / 3) * 3
        ]
        pathways_list.append(temp_dict)

    # recreate rxn string format
    data = list()
    for record in pathways_list:
        temp_path_list = list()
        for idx, i in enumerate(record["SMILES"]):
            rxn_string = (
                i.split(">>")[0]
                + ">"
                + record["name"][idx]
                + ">"
                + record["enthalpy"][idx]
                + "$"
                + record["stoi"][idx]
                + ">"
                + i.split(">>")[1]
            )
            temp_path_list.append(rxn_string)
        data.append(temp_path_list)

    starters = set(mol_smiles)

    # remove pathway if the target somehow is also a reactant in a reaction
    temp_data = list()
    for pathway in data:
        error_flag = False
        for rxn in pathway:
            reas = rxn.split(">")[0].split(".")
            if target_smiles in reas:
                error_flag = True
        if error_flag is False:
            temp_data.append(pathway)
    data = temp_data
    print(
        "Start ranking pathways. Number of pathways after sanitization:",
        len(data),
    )

    # Enthalpy ###############################
    if weights["reaction_dH"] != 0:
        print("reaction_dH ranking working")
        H_list = []  # store max H for each pathway
        for path in data:
            path_max_H = -9999
            for rxn in path:
                name = rxn.split(">")[1]
                dH = rxn.split(">")[2].split("$")[0]
                if "2-step" in name or "&" in name:
                    dH = float(dH) / 2
                path_max_H = max(path_max_H, float(dH))
            H_list.append(path_max_H)
        min_H = min(H_list)
        max_H = max(H_list)
        diff_H = -(max_H - min_H)
        H_score_list = []
        for i in H_list:
            if diff_H != 0:
                H_score_list.append((i - min_H) / diff_H + 1)
            else:
                H_score_list.append((i - min_H) / 0.001 + 1)
        print("reaction_dH ranking finished")
    else:
        H_score_list = [0] * len(data)

    # Atom economy #############################
    eco_list = list()
    eco_score_list = list()

    def path_eco(
        _path,
    ):  # assuming no circular loops; if a middle rxn produces a starter
        # or intermediat consumed in upstream, it's considered recycled;
        _path = list(_path)
        left_dict = dict()  # smiles: number of mol
        right_dict = dict()
        target_weight = Descriptors.MolWt(
            rdkit.Chem.rdmolfiles.MolFromSmiles(target_smiles)
        )
        #  target_mol = float()
        idx_remove = 0
        for idx, rxn in enumerate(
            _path
        ):  # find rxn for final target, add reas in left,
            # by-product in right, remove rxn
            pros = rxn.split(">")[3].split(".")
            reas = rxn.split(">")[0].split(".")
            rea_stoi = eval(rxn.split(">")[2].split("$")[1])
            pro_stoi = eval(rxn.split(">")[2].split("$")[2])

            if target_smiles in pros:
                idx_remove = idx
                # target_idx = pros.index(target_smiles)
                target_stoi = 0
                for idx2, pro in enumerate(pros):
                    if pro == target_smiles:
                        target_stoi += pro_stoi[idx2]
                for idx2, pro in enumerate(pros):
                    if pro != target_smiles:
                        if pro not in right_dict:
                            right_dict[pro] = pro_stoi[idx2] / target_stoi
                        else:
                            right_dict[pro] = (
                                right_dict[pro] + pro_stoi[idx2] / target_stoi
                            )
                for idx2, rea in enumerate(reas):
                    if rea not in left_dict:
                        left_dict[rea] = rea_stoi[idx2] / target_stoi
                    else:
                        left_dict[rea] = (
                            left_dict[rea] + rea_stoi[idx2] / target_stoi
                        )
                break
        del _path[idx_remove]

        repeat_flag = True
        # explored_rxns = set()
        timeout = 100  # stop if encounter loops
        timecount = 0
        while (
            repeat_flag is True and len(_path) > 0
        ):  # repeat until left only contain starters;
            # check if all parents are starters
            timecount += 1
            popkey = 99
            for (
                i
            ) in left_dict:  # find idx for a intermediate molecule to produce
                if i not in starters:
                    popkey = i
                    break
            # num_mol = left_dict.pop(popkey)
            num_mol = left_dict[popkey]
            for idx, rxn in enumerate(_path):
                pros = rxn.split(">")[3].split(".")
                reas = rxn.split(">")[0].split(".")
                rea_stoi = eval(rxn.split(">")[2].split("$")[1])
                pro_stoi = eval(rxn.split(">")[2].split("$")[2])
                # if popkey in pros and rxn not in explored_rxns:
                # find a rxn that product the molecule
                if popkey in pros:
                    # explored_rxns.add(rxn)
                    # popmol_idx = pros.index(popkey)
                    popmol_stoi = 0
                    for idx2, pro in enumerate(pros):
                        if pro == popkey:
                            popmol_stoi += pro_stoi[idx2]
                    for idx2, pro in enumerate(pros):
                        # if pro != popkey: # add all reactants to left
                        if pro not in right_dict:
                            right_dict[pro] = (
                                num_mol / popmol_stoi * pro_stoi[idx2]
                            )
                        else:
                            right_dict[pro] = (
                                right_dict[pro]
                                + num_mol / popmol_stoi * pro_stoi[idx2]
                            )
                    for idx2, rea in enumerate(reas):
                        if rea not in left_dict:
                            left_dict[rea] = (
                                num_mol / popmol_stoi * rea_stoi[idx2]
                            )
                        else:
                            left_dict[rea] = (
                                left_dict[rea]
                                + num_mol / popmol_stoi * rea_stoi[idx2]
                            )

                    # move used rxn to end of list to avoid
                    # code loop if a loop in pathway
                    del _path[idx]
                    _path.append(rxn)
                    break

            key_list = list()
            for key in left_dict:  # balance left and right dict
                if key not in starters and key in right_dict:
                    key_list.append(key)
            for key in key_list:
                if left_dict[key] == right_dict[key]:
                    # print(left_dict)
                    left_dict.pop(key)
                    # print(left_dict)
                    right_dict.pop(key)
                elif left_dict[key] < right_dict[key]:
                    right_dict[key] = right_dict[key] - left_dict[key]
                    left_dict.pop(key)
                else:
                    left_dict[key] = left_dict[key] - right_dict[key]
                    right_dict.pop(key)

            done_flag = True
            for key in left_dict:
                if key not in starters:  # check if left only contain starters
                    done_flag = False
            if done_flag is True or timecount == timeout:
                repeat_flag = False

        if timecount == timeout:
            return None
        by_product_weight = 0
        for i in right_dict:
            by_product_weight += (
                Descriptors.MolWt(rdkit.Chem.rdmolfiles.MolFromSmiles(i))
                * right_dict[i]
            )
        to_return = target_weight / (target_weight + by_product_weight)
        return to_return

    if weights["atom_economy"] != 0:
        print("atom_economy ranking working")
        none_eco = 0
        min_eco = 1
        max_eco = 0
        for idx, path in enumerate(data):
            try:
                atom_eco = path_eco(path)
            except KeyError:
                atom_eco = 0
                print("Atom economy calculation error for pathway", idx + 1)
            # print(idx)
            if atom_eco is None:
                none_eco += 1
                eco_list.append(None)
            else:
                min_eco = min(min_eco, atom_eco)
                max_eco = max(max_eco, atom_eco)
                eco_list.append(atom_eco)

        for idx, i in enumerate(eco_list):
            if i is None:
                eco_list[idx] = min_eco

        print("min_eco", min_eco)
        print("max_eco", max_eco)
        # print("none_eco", none_eco)
        diff_eco = max_eco - min_eco
        for i in eco_list:
            if diff_eco != 0:
                eco_score_list.append((i - min_eco) / diff_eco)
            else:
                eco_score_list.append((i - min_eco) / 0.001)
        print("atom_economy ranking finished")

    else:
        eco_score_list = [0] * len(data)
        eco_list = [0] * len(data)

    # Number of steps ################################
    if weights["number_of_steps"] != 0:
        print("number_of_steps ranking working")
        step_list = []
        for path in data:
            step_list.append(len(path))
        min_step = min(step_list)
        max_step = max(step_list)
        diff_steps = -(max_step - min_step)
        if diff_steps == 0:  # if all pathways have same number of rxns
            diff_steps = 1  # value doesn't matter

        step_score_list = []
        for i in step_list:
            step_score_list.append((i - min_step) / diff_steps + 1)
        print("number_of_steps ranking finished")
    else:
        step_score_list = [0] * len(data)

    # By-product ##############################
    if weights["by_product_number"] != 0:
        print("by_product_number ranking working")
        by_pro_list = []
        intermediate_by_product_dict_list = []

        soup_list = list()  # [(idx, {SMILES,}), ]
        Inter_soup_list = list()  # [(idx, interSMILES, {SMILES,}), ]

        for idx, path in enumerate(data):
            mol_soup = set()  # for pathway_byproduct
            # mol_bypro_dict = dict()   # for intermediate by_product
            for rxn in path:
                products = rxn.split(">")[3].split(".")
                reas = rxn.split(">")[0].split(".")
                mol_soup.update(set(products), set(reas))

            soup_list.append((idx, mol_soup))
            # num_by_pro = Byproduct_index( list(mol_soup) )
            for inter_mol in mol_soup:
                if inter_mol not in helpers or inter_mol in just_starters:
                    Inter_soup_list.append((idx, inter_mol, mol_soup))

        with Pool(processes=num_process) as pool:
            results = [
                pool.apply_async(Byproduct_index, args=(work[0], work[1]))
                for work in soup_list
            ]  # workers pool
            by_pro_list_tuples = [
                r.get() for r in results
            ]  # [(path_idx, byproduct num),]

        with Pool(processes=num_process) as pool:
            results = [
                pool.apply_async(
                    Inter_byproduct, args=(work[0], work[1], work[2])
                )
                for work in Inter_soup_list
            ]  # workers pool
            inter_by_pro_list_tuples = [
                r.get() for r in results
            ]  # [(path_idx, SMILES, inter-byproduct num),]

        by_pro_list_tuples.sort()
        for i in by_pro_list_tuples:
            by_pro_list.append(i[1])

        inter_by_pro_list_tuples.sort()
        inter_by_pro_list_tuples.append((-1, "END", -1))
        path_index = 0
        mol_bypro_dict = dict()
        for i in inter_by_pro_list_tuples:
            if i[0] == path_index:
                mol_bypro_dict[i[1]] = i[2]
            else:
                intermediate_by_product_dict_list.append(mol_bypro_dict)
                mol_bypro_dict = dict()
                mol_bypro_dict[i[1]] = i[2]
                path_index += 1

        min_by_pro = min(by_pro_list)
        max_by_pro = max(by_pro_list)
        diff_by_pro = -(max_by_pro - min_by_pro)
        # diff_by_pro = 1

        by_pro_score_list = []

        for i in by_pro_list:
            if diff_by_pro != 0:
                by_pro_score_list.append((i - min_by_pro) / diff_by_pro + 1)
            else:
                by_pro_score_list.append((i - min_by_pro) / 0.001 + 1)
        print("by_product_number ranking finished")

    else:
        by_pro_score_list = [0] * len(data)
        by_pro_list = [0] * len(data)
        intermediate_by_product_dict_list = []
        for path in data:
            mol_soup = set()
            mol_bypro_dict = dict()
            for rxn in path:
                products = rxn.split(">")[3].split(".")
                reas = rxn.split(">")[0].split(".")
                mol_soup.update(set(products), set(reas))
            for inter_mol in mol_soup:
                if inter_mol not in helpers or inter_mol in just_starters:
                    mol_bypro_dict[inter_mol] = 0

            intermediate_by_product_dict_list.append(mol_bypro_dict)

    # if produce salt #################################
    salt_name_list = {
        "Hydrosulfide Anion Substitution with Alkyl Halides",
        "Hydrosulfide Anion Substitution with Alkyl Halides, Intramolecular",
        "Preparation of Sulfides",
        "Preparation of Sulfides, Intramolecular",
        "Thiols from Thiourea",
        "Sulfurous Acid Salts Addition to Oxiranes",
        "Organosulfates Reduction",
        "Alkenes Addition by Bisulfites",
        "Alkylation of Sulfinic Acids with Halides",
        "Reduction of Sulfonyl Chlorides to Sulfinic Acids",
        "Strecker Sulfite Alkylation",
    }

    if weights["salt_score"] != 0:
        print("salt_score ranking working")
        salt_rxn_num_list = []
        for path in data:
            salt_rxn_num = 0
            for rxn in path:
                name = rxn.split(">")[1]
                if name in salt_name_list:
                    salt_rxn_num += 1
            salt_rxn_num_list.append(salt_rxn_num)

        min_salt = min(salt_rxn_num_list)
        max_salt = max(salt_rxn_num_list)
        diff_salt = -(max_salt - min_salt)
        if diff_salt == 0:  # if all pathways have same number
            diff_salt = 1  # value doesn't matter

        salt_score_list = []
        for i in salt_rxn_num_list:
            salt_score_list.append((i - min_salt) / diff_salt + 1)
        print("salt_score ranking finished")
    else:
        salt_score_list = [0] * len(data)

    # coolness ############################
    if cool_reactions and weights["coolness"] != 0:
        print("coolness ranking working")
        cool_rxn_num_list = []
        for path in data:
            cool_rxn_num = 0
            for rxn in path:
                name = rxn.split(">")[1]
                if name in cool_reactions:
                    cool_rxn_num = cool_rxn_num + cool_reactions[name]
            cool_rxn_num_list.append(cool_rxn_num)

        min_cool = min(cool_rxn_num_list)
        max_cool = max(cool_rxn_num_list)
        diff_cool = max_cool - min_cool
        if diff_cool == 0:  # if all pathways have same number
            diff_cool = 1  # value doesn't matter

        cool_score_list = []
        for i in cool_rxn_num_list:
            cool_score_list.append((i - min_cool) / diff_cool)
        print("coolness ranking finished")
    else:
        cool_score_list = [0] * len(data)

    # reaxys score #############################
    if weights["in_reaxys"] != 0:
        print("in_reaxys ranking working")
        # load 2 files
        reaxys_batch = np.genfromtxt(
            f"{job_name}_reaxys_batch_query.txt",
            comments="?",
            dtype=str,
            delimiter=",",
            skip_header=0,
        )
        reaxys_result = np.genfromtxt(
            reaxys_result_name,
            comments="?",
            dtype=int,
            delimiter=",",
            skip_header=0,
            usecols=(2),
        )
        rxn_list = list()  # rxns, used as keys. value is number of hits
        for i in reaxys_batch:
            rxn_list.append(i[8:-1])

        in_reaxys = list()  # rxns in the set are reported in reaxys

        for idx, i in enumerate(reaxys_result):
            if i != 0:
                in_reaxys.append(rxn_list[idx])

        print("Number of reactions reported in Reaxys", len(in_reaxys))

        reaxys_set = set(in_reaxys)

        def SMILES_rm_keku(
            _SMILES,
        ):  # remove reaxys unsuported mols, kekulize mols,
            # return a reaction SMILES
            reas = _SMILES.split(">")[0].split(".")
            pros = _SMILES.split(">")[3].split(".")
            new_rxn = str()
            reas_set = set()
            pros_set = set()
            for rea in reas:
                if rea not in reaxys_not_supported and rea not in reas_set:
                    new_rxn = (
                        new_rxn
                        + rdkit.Chem.rdmolfiles.MolToSmiles(
                            rdkit.Chem.rdmolfiles.MolFromSmiles(rea),
                            kekuleSmiles=True,
                        )
                        + "."
                    )
                    reas_set.add(rea)
            new_rxn = new_rxn[:-1] + ">>"
            for pro in pros:
                if pro not in reaxys_not_supported and pro not in pros_set:
                    new_rxn = (
                        new_rxn
                        + rdkit.Chem.rdmolfiles.MolToSmiles(
                            rdkit.Chem.rdmolfiles.MolFromSmiles(pro),
                            kekuleSmiles=True,
                        )
                        + "."
                    )
                    pros_set.add(pro)
            if new_rxn[-1] == ".":
                new_rxn = new_rxn[:-1]
            return new_rxn

        reaxys_percent_list = []

        for path in data:
            total_num_rxns = len(path)
            in_reaxys_num = 0
            for rxn in path:
                if SMILES_rm_keku(rxn) in reaxys_set:
                    in_reaxys_num += 1
            reaxys_percent_list.append(in_reaxys_num / total_num_rxns)

        min_reaxys = min(reaxys_percent_list)
        max_reaxys = max(reaxys_percent_list)
        diff_reaxys = max_reaxys - min_reaxys
        if diff_reaxys == 0:  # if all pathways have same number
            diff_reaxys = 1  # value doesn't matter

        reaxys_score_list = []
        for i in reaxys_percent_list:
            reaxys_score_list.append((i - min_reaxys) / diff_reaxys)
        print("in_reaxys ranking finished")

    else:
        reaxys_score_list = [0] * len(data)

    # final score
    weight = list()

    #                dH  steps  by-pro   eco       salt      reaxys
    # weight =    [     0.0,  0.4  ,   0.2    ,0.1  ,   0.2 ,    0.3  ]
    weight.append(weights["reaction_dH"])
    weight.append(weights["number_of_steps"])
    weight.append(weights["by_product_number"])
    weight.append(weights["atom_economy"])
    weight.append(weights["salt_score"])
    weight.append(weights["in_reaxys"])
    weight.append(weights["coolness"])
    #
    final_score = list()

    for i in range(len(data)):
        final_score.append(
            H_score_list[i] * weight[0]
            + step_score_list[i] * weight[1]
            + by_pro_score_list[i] * weight[2]
            + eco_score_list[i] * weight[3]
            + salt_score_list[i] * weight[4]
            + reaxys_score_list[i] * weight[5]
            + cool_score_list[i] * weight[6]
        )

    # print("len(final_score)", len(final_score))
    print("min score", min(final_score))
    print("max score", max(final_score))
    print()

    (
        final_score,
        data,
        eco_list,
        by_pro_list,
        intermediate_by_product_dict_list,
        H_score_list,
        step_score_list,
        by_pro_score_list,
        eco_score_list,
        salt_score_list,
        reaxys_score_list,
        cool_score_list,
    ) = zip(
        *sorted(
            zip(
                final_score,
                data,
                eco_list,
                by_pro_list,
                intermediate_by_product_dict_list,
                H_score_list,
                step_score_list,
                by_pro_score_list,
                eco_score_list,
                salt_score_list,
                reaxys_score_list,
                cool_score_list,
                strict=False,
            )
        ),
        strict=False,
    )

    def clean_rxn_strings(
        _path,
    ):  # for a path, return a list of clean rxn smiles, names, dH of each rxn
        clean_path = list()
        names = list()
        dH_list = list()
        for rxn in _path:
            dH = rxn.split(">")[2].split("$")[0]
            name = rxn.split(">")[1]
            if "2-step" in name or "&" in name:
                dH = float(dH) / 2

            rxn_str = rxn.split(">")[0] + ">>" + rxn.split(">")[3]
            clean_path.append(rxn_str)
            names.append(rxn.split(">")[1])
            dH_list.append(str(dH))
        clean_path = clean_path + names + dH_list
        return clean_path

    k = 0
    with open(
        f"{job_name}_ranked_pathways.txt", "w", encoding="utf-8"
    ) as f_result:
        for idx, path in reversed(list(enumerate(data))):  # print result
            print_flag = True
            if print_flag is True:
                k += 1
                print("ranking", k)
                f_result.write("ranking " + str(k))
                f_result.write("\n")
                print("final score", final_score[idx])
                print(
                    "Max reaction enthalpy score",
                    H_score_list[idx],
                    " x ",
                    weight[0],
                    " = ",
                    H_score_list[idx] * weight[0],
                )
                print(
                    "Number of reactions score",
                    step_score_list[idx],
                    " x ",
                    weight[1],
                    " = ",
                    step_score_list[idx] * weight[1],
                )
                print(
                    "By-product score",
                    by_pro_score_list[idx],
                    " x ",
                    weight[2],
                    " = ",
                    by_pro_score_list[idx] * weight[2],
                )
                print(
                    "Pathway atom economy score",
                    eco_score_list[idx],
                    " x ",
                    weight[3],
                    " = ",
                    eco_score_list[idx] * weight[3],
                )
                print(
                    "Salt score",
                    salt_score_list[idx],
                    " x ",
                    weight[4],
                    " = ",
                    salt_score_list[idx] * weight[4],
                )
                print(
                    "Reaxys score",
                    reaxys_score_list[idx],
                    " x ",
                    weight[5],
                    " = ",
                    reaxys_score_list[idx] * weight[5],
                )
                print(
                    "Cool score",
                    cool_score_list[idx],
                    " x ",
                    weight[6],
                    " = ",
                    cool_score_list[idx] * weight[6],
                )

                f_result.write("final score " + str(final_score[idx]))
                f_result.write("\n")
                print("atom economy", eco_list[idx])
                f_result.write("atomic economy " + str(eco_list[idx]))
                f_result.write("\n")
                print("pathway by-product", by_pro_list[idx])
                f_result.write("pathway by-product " + str(by_pro_list[idx]))
                f_result.write("\n")
                print(
                    "intermediate by-product",
                    intermediate_by_product_dict_list[idx],
                )
                f_result.write(
                    "intermediate by-product "
                    + str(intermediate_by_product_dict_list[idx])
                )
                f_result.write("\n")
                f_result.write("reaction SMILES, name, and enthalpy:")
                f_result.write("\n")
                for i in clean_rxn_strings(path):
                    print(i)
                    f_result.write(i)
                    f_result.write("\n")
                print()
                f_result.write("\n")

    end_time = time.time()
    elapsed_time = end_time - start_time  # /60
    print()
    print("time used:", "{:.2f}".format(elapsed_time), "seconds")


# pathway visualization
def create_page(
    pathway_dict,
    page_number,
    reaxys_result_set,
    _my_start,
    _helpers,
    _reaxys_rxn_color,
    _normal_rxn_color,
):
    font_path = Path(__file__).parent / "OpenSans-Regular.ttf"

    def trim(im):  # trim white space around image
        bg = Image.new(im.mode, im.size, im.getpixel((0, 0)))
        diff = ImageChops.difference(im, bg)
        diff = ImageChops.add(diff, diff, 2.0, -100)
        bbox = diff.getbbox()
        if bbox:
            return im.crop(bbox)

    def add_margin(pil_img):  # add white margin on top and bottom of image
        width, height = pil_img.size
        new_width = width
        new_height = height + 60
        result = Image.new(
            pil_img.mode, (new_width, new_height), (255, 255, 255)
        )
        result.paste(pil_img, (0, 30))
        return result

    def add_text(img, msg):  # add a text on top-right of image
        img_w, img_h = img.size
        I1 = ImageDraw.Draw(img)
        myFont = ImageFont.truetype(str(font_path), 25)
        new_box = I1.textbbox((0, 0), msg, font=myFont)
        I1.text((img_w - new_box[2], 0), msg, font=myFont, fill=(0, 0, 255))
        return img

    def SMILES_rm_keku(
        _SMILES,
    ):  # remove reaxys unsuported mols, kekulize mols, return a reaction SMILES
        reas = _SMILES.split(">>")[0].split(".")
        pros = _SMILES.split(">>")[1].split(".")
        new_rxn = str()
        reas_set = set()
        pros_set = set()
        for rea in reas:
            if rea not in reaxys_not_supported and rea not in reas_set:
                new_rxn = (
                    new_rxn
                    + rdkit.Chem.rdmolfiles.MolToSmiles(
                        rdkit.Chem.rdmolfiles.MolFromSmiles(rea),
                        kekuleSmiles=True,
                    )
                    + "."
                )
                reas_set.add(rea)
        new_rxn = new_rxn[:-1] + ">>"
        for pro in pros:
            if pro not in reaxys_not_supported and pro not in pros_set:
                new_rxn = (
                    new_rxn
                    + rdkit.Chem.rdmolfiles.MolToSmiles(
                        rdkit.Chem.rdmolfiles.MolFromSmiles(pro),
                        kekuleSmiles=True,
                    )
                    + "."
                )
                pros_set.add(pro)
        if new_rxn[-1] == ".":
            new_rxn = new_rxn[:-1]
        return new_rxn

    help_dict = _helpers

    salt_name_set = {
        "Hydrosulfide Anion Substitution with Alkyl Halides",
        "Hydrosulfide Anion Substitution with Alkyl Halides, Intramolecular",
        "Preparation of Sulfides",
        "Preparation of Sulfides, Intramolecular",
        "Thiols from Thiourea",
        "Sulfurous Acid Salts Addition to Oxiranes",
        "Organosulfates Reduction",
        "Alkenes Addition by Bisulfites",
        "Alkylation of Sulfinic Acids with Halides",
        "Reduction of Sulfonyl Chlorides to Sulfinic Acids",
        "Strecker Sulfite Alkylation",
    }

    known_rxn_set = {  # name of rxns we want to be colored even if
        # they're not in returned reaxys result,
        # should only apply to specific rxns
        "Nitrogen Dioxide Disproportionation",
        "Andrussow Process",
        "Hydrogenation of Nitric Oxide",
    }

    my_start = _my_start

    rxn_list = pathway_dict["SMILES"]  # list of rxns in pathway
    name_list = pathway_dict["name"]

    for idx, _i in enumerate(name_list):
        name_list[idx] = textwrap.fill(
            name_list[idx], 25
        )  # wrape text so not too long

    title = (
        "Rank "
        + pathway_dict["ranking"]
        + "\n"
        + "Atomic Economy "
        + pathway_dict["atomic_economy"]
        + "\n"
        + "Pathway by-product "
        + pathway_dict["pathway_by-product"]
    )  # to display on each page

    split_list = pathway_dict["enthalpy"]
    dH_list = list()
    for _idx, i in enumerate(split_list):
        dH_list.append(str("%.1f" % float(i)) + " kcal/mol")

    for idx, _i in enumerate(dH_list):
        if "2-step" in name_list[idx] or "&" in name_list[idx]:
            dH_list[idx] = dH_list[idx] + " ×2"

    for idx, _i in enumerate(name_list):  # add dH to name
        if name_list[idx].replace("\n", " ") in salt_name_set:
            name_list[idx] = name_list[idx] + "(salt)"
        name_list[idx] = name_list[idx] + "\n" + dH_list[idx]

    G = nx.DiGraph(
        rankdir="TB",  # top to bottom layout
        # ranksep="0.8",
    )

    inter_py_pro = pathway_dict[
        "intermediate_by-product"
    ]  # intermediate by product number, dict, key smiles, value int
    inter_py_pro["O=C=O"] = ""  # for reagent node (CO2)
    # rxn_node = 0
    node_labels_dict = dict()
    node_labels_reagent = dict()

    reaxys_rxn_nodes = (
        set()
    )  # set of int, rxn nodes in reaxys_result_set, used for coloring edges

    def addjust_mol_image(mol):
        if rdkit.Chem.rdmolfiles.MolToSmiles(mol) == "C=O":
            mol = Chem.AddHs(mol)
        return Draw.MolToImage(mol)

    for rxn_node, rxn in enumerate(rxn_list):  # add nodes and edges to graph
        reas = rxn.split(">>")[0].split(".")
        pros = rxn.split(">>")[1].split(".")

        reagent_label = str()
        for rea in reas:  # to add reagents to label
            if rea in help_dict:
                reagent_label = reagent_label + "+" + help_dict[rea] + " "
        for pro in pros:
            if pro in help_dict:
                reagent_label = reagent_label + "-" + help_dict[pro] + " "

        G.add_node(rxn_node)
        if (
            SMILES_rm_keku(rxn) in reaxys_result_set
            or " ".join(name_list[rxn_node].replace("\n", " ").split()[:-2])
            in known_rxn_set
        ):
            reaxys_rxn_nodes.add(rxn_node)

        node_labels_dict[rxn_node] = name_list[rxn_node]
        node_labels_reagent[rxn_node] = reagent_label + "\n"
        for rea in reas:  # add nodes and edges for non-reagents
            if rea not in help_dict:
                G.add_node(
                    rea,
                    image1=add_text(
                        add_margin(
                            trim(
                                addjust_mol_image(
                                    rdkit.Chem.rdmolfiles.MolFromSmiles(rea)
                                )
                            )
                        ),
                        str(inter_py_pro[rea]),
                    ),
                )
                G.add_edge(rea, rxn_node, arror_margin=10)
                for pro in pros:
                    if pro not in help_dict:
                        G.add_node(
                            pro,
                            image1=add_text(
                                add_margin(
                                    trim(
                                        addjust_mol_image(
                                            rdkit.Chem.rdmolfiles.MolFromSmiles(
                                                pro
                                            )
                                        )
                                    )
                                ),
                                str(inter_py_pro[pro]),
                            ),
                        )
                        G.add_edge(rxn_node, pro, arror_margin=120)
        # rxn_node += 1

    fig, ax = plt.subplots(figsize=(25, 50))  # 25, 40

    G2 = copy.deepcopy(G)
    # my_start_being_produced = False

    for i in my_start:
        if i in G.nodes:
            my_start_being_produced = False
            starter_index = (
                0  # if starter usded multiple times, make them unique nodes
            )
            for node in G.nodes:
                if i in list(G.successors(node)):
                    my_start_being_produced = True

                if i in list(G.predecessors(node)):
                    G2.remove_edge(i, node)
                    G2.add_node(
                        i + str(starter_index),
                        image1=add_text(
                            add_margin(
                                trim(
                                    addjust_mol_image(
                                        rdkit.Chem.rdmolfiles.MolFromSmiles(i)
                                    )
                                )
                            ),
                            str(inter_py_pro[i]),
                        ),
                    )
                    G2.add_edge(
                        i + str(starter_index),
                        node,
                        arror_margin=10,
                        # minlen = 2
                    )
                    starter_index += 1

            if my_start_being_produced is False:
                G2.remove_node(i)

    pos = nx.nx_agraph.graphviz_layout(G2, prog="dot")

    edge_color_list = []  # for coloring edges,
    # black if not in reaxys, green if in reaxys
    for edge in G2.edges():
        if edge[0] in reaxys_rxn_nodes or edge[1] in reaxys_rxn_nodes:
            edge_color_list.append(_reaxys_rxn_color)  # "g"
        else:
            edge_color_list.append(_normal_rxn_color)  # "k"

    nx.draw_networkx_edges(  # draw edges
        G2,
        pos=pos,
        ax=ax,
        arrows=True,
        arrowsize=30,
        width=5,
        #  arrowstyle="-",
        min_source_margin=80,
        min_target_margin=80,  # 120,
        edge_color=edge_color_list,
    )

    label_props = dict(
        boxstyle="round,pad=0.3",
        fc="white",
        ec="white",
        #  lw=1,
    )

    nx.draw_networkx_labels(  # draw rxn nodes labels for reagents
        G2,
        pos=pos,
        ax=ax,
        font_size=20,
        font_color="r",
        labels=node_labels_reagent,
        bbox=label_props,
        verticalalignment="bottom",
    )

    nx.draw_networkx_labels(  # draw rxn nodes labels for rxn names
        G2,
        pos=pos,
        ax=ax,
        font_size=20,
        font_color="k",
        labels=node_labels_dict,
        bbox=label_props,
        verticalalignment="top",
    )

    # from networkx custom icon example
    # Transform from data coordinates (scaled between
    # xlim and ylim) to display coordinates
    tr_figure = ax.transData.transform

    # Transform from display to figure coordinates
    tr_axes = fig.transFigure.inverted().transform

    # Select the size of the image (relative to the X axis)
    # icon_size = (ax.get_xlim()[1] - ax.get_xlim()[0]) * 0.035
    icon_size = 0.1  # 0.12705   # icon size
    icon_center = icon_size / 2.0

    # Add the respective image to each node
    for n in G2.nodes:
        xf, yf = tr_figure(pos[n])
        xa, ya = tr_axes((xf, yf))
        # get overlapped axes and plot icon
        a = plt.axes([xa - icon_center, ya - icon_center, icon_size, icon_size])
        if "image1" in G2.nodes[n]:
            a.imshow(G2.nodes[n]["image1"])
        a.axis("off")

    ax.set_frame_on(False)  # don't show frame
    ax.set_title(title, fontsize=20, y=0.95)

    buf = io.BytesIO()
    plt.savefig(buf, format="pdf", bbox_inches="tight")
    buf.seek(0)
    plt.close()
    print("page done:", page_number, flush=True)
    return (page_number, buf)


def pathway_visualization(
    starters,
    helpers,
    num_process,
    reaxys_result_name=False,
    job_name="default_job_name",
    exclude_smiles=None,
    reaxys_rxn_color="blue",
    normal_rxn_color="black",
):
    from PyPDF2 import PdfFileMerger, PdfFileReader

    if exclude_smiles is None:
        exclude_smiles = []
    start_time = time.time()

    my_start = set()
    for i in starters:
        my_start.add(
            rdkit.Chem.rdmolfiles.MolToSmiles(
                rdkit.Chem.rdmolfiles.MolFromSmiles(i)
            )
        )

    if reaxys_result_name is not False:
        reaxys_batch = np.genfromtxt(
            f"{job_name}_reaxys_batch_query.txt",
            comments="?",
            dtype=str,
            delimiter=",",
            skip_header=0,
        )
        reaxys_result = np.genfromtxt(
            reaxys_result_name,
            comments="?",
            dtype=int,
            delimiter=",",
            skip_header=0,
            usecols=(2),
        )
        rxn_list = list()  # rxns, used as keys. value is number of hits
        for i in reaxys_batch:
            rxn_list.append(i[8:-1])
        in_reaxys = list()  # rxns in the set are reported in reaxys
        for idx, i in enumerate(reaxys_result):
            if i != 0:
                in_reaxys.append(rxn_list[idx])
    else:
        in_reaxys = list()

    reaxys_set = set(in_reaxys)

    with open(f"{job_name}_ranked_pathways.txt", encoding="utf-8") as f:
        lines = f.readlines()

    clean_list = list()
    for i in lines:
        if i != "\n":
            clean_list.append(i.strip())

    pathways_list = list()
    # [{final_score:,eco:,pathy_by:,inter_by:{},SMILES:[],Nmaes:[],dH:[]}]

    pathway_marker = list()
    pathway_num = 1

    for idx, i in enumerate(clean_list):
        if "ranking" in i:
            pathway_marker.append(idx)

    for idx, marker in enumerate(pathway_marker):
        temp_dict = dict()
        temp_dict["ranking"] = str(pathway_num)
        pathway_num += 1
        temp_dict["final_score"] = clean_list[marker + 1][12:]
        temp_dict["atomic_economy"] = (
            str("%.1f" % float(float(clean_list[marker + 2][15:]) * 100)) + "%"
        )
        temp_dict["pathway_by-product"] = clean_list[marker + 3][19:]
        temp_dict["intermediate_by-product"] = eval(clean_list[marker + 4][24:])

        if idx + 1 < len(pathway_marker):  # 1, 2, 3
            next_marker = pathway_marker[idx + 1]
        else:
            next_marker = len(clean_list)

        temp_dict["SMILES"] = clean_list[
            marker + 6 : marker + 6 + int((next_marker - (marker + 6)) / 3)
        ]

        all_reas = list()
        for j in temp_dict["SMILES"]:
            all_reas = all_reas + j.split(">>")[0].split(".")
        temp_dict["all_reactants"] = set(all_reas)

        temp_dict["name"] = clean_list[
            marker + 6 + int((next_marker - (marker + 6)) / 3) : marker
            + 6
            + int((next_marker - (marker + 6)) / 3) * 2
        ]
        temp_dict["enthalpy"] = clean_list[
            marker + 6 + int((next_marker - (marker + 6)) / 3) * 2 : marker
            + 6
            + int((next_marker - (marker + 6)) / 3) * 3
        ]

        pathways_list.append(temp_dict)

    print("Number of pathways: ", len(pathways_list))
    print("Number of reactions in reaxys: ", len(reaxys_set))

    work_list = list()  # list of (pathway_dict, page_number)

    exclude_set = set()
    for i in exclude_smiles:
        exclude_set.add(
            rdkit.Chem.rdmolfiles.MolToSmiles(
                rdkit.Chem.rdmolfiles.MolFromSmiles(i)
            )
        )
    for idx, i in enumerate(pathways_list):
        if not exclude_set & i["all_reactants"]:
            work_list.append((i, idx + 1))

    help_dict = dict()
    for key in helpers:
        if (
            rdkit.Chem.rdmolfiles.MolToSmiles(
                rdkit.Chem.rdmolfiles.MolFromSmiles(key)
            )
            not in my_start
        ):
            help_dict[
                rdkit.Chem.rdmolfiles.MolToSmiles(
                    rdkit.Chem.rdmolfiles.MolFromSmiles(key)
                )
            ] = helpers[key]

    print("Working on creating pages...")
    with Pool(processes=num_process) as pool:
        results = [
            pool.apply_async(
                create_page,
                args=(
                    work[0],
                    work[1],
                    reaxys_set,
                    my_start,
                    help_dict,
                    reaxys_rxn_color,
                    normal_rxn_color,
                ),
            )
            for work in work_list
        ]  # workers pool
        pages = [r.get() for r in results]

    print("Finished pages, writing to pdf")
    pages.sort()  # [(page_number, object)]
    merger = PdfFileMerger()
    for page in pages:
        merger.append(PdfFileReader(page[1]))

    filename_part2 = ""
    if exclude_smiles:
        filename_part2 = "excluded_some_molecules_"
    merger.write(f"{job_name}_{filename_part2}pathways_visualized.pdf")
    merger.close()

    end_time = time.time()
    elapsed_time = (end_time - start_time) / 60
    print()
    print("time used:", "{:.2f}".format(elapsed_time), "minutes")
