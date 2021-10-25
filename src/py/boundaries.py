#!/usr/bin/env python

from cwalks_analysis import load_cwalk_graph, load_files
import pandas as pd
import networkx as nx
from intervaltree import Interval, IntervalTree
from filtering import typical_chromosomes


def tad_boundaries(bedfile: str) -> dict:
    boundaries = pd.read_csv(bedfile, delimiter="\t", header=None)
    boundaries_chr = {x: y for x, y in boundaries.groupby(0)}  # dictionary with chr as keys
    keys, values = [], []
    for key in boundaries_chr.keys():
        left = boundaries_chr[key][1].tolist()
        right = boundaries_chr[key][2].tolist()
        sets_list = [(a, b) for a, b in zip(left, right)]
        keys.append(key)
        values.append(sets_list)
    return dict(zip(keys, values))


graphs, _ = load_files("./cwalks", load_cwalk_graph)  # load .txt cwalk graph
boundaries_dict = tad_boundaries("K562_Lieberman-raw_TADs.bed")

# prove that cwalks are in TAD boundaries
tree_dict = dict()  # ex. tree_dict["chr1"] will be an object of type IntervalTree
for key in boundaries_dict.keys():
    intervals = boundaries_dict[key]
    # Interval tree construction, separate for each chromosome
    tree_dict[key] = IntervalTree.from_tuples(intervals)

for graph in graphs:
    for cwalk in list(nx.connected_components(graph)):
        cwalk = list(cwalk)
        cwalk_boundaries = [Interval(node[0], node[1]) for node in cwalk]
        print(cwalk_boundaries)

