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


def main():
    graphs, _ = load_files("./cwalks", load_cwalk_graph)  # load .txt cwalk graph
    boundaries_dict = tad_boundaries("K562_Lieberman-raw_TADs.bed")  # list of tuples

    tree_dict = dict()  # ex. tree_dict["chr1"] will be an object of type IntervalTree
    for key in boundaries_dict.keys():
        intervals = boundaries_dict[key]
        # Interval tree construction, separate for each chromosome
        tree_dict[key] = IntervalTree.from_tuples(intervals)

    out_tad = 0
    in_tad = 0
    graphs_number = 0
    cwalks_number = []

    for graph in graphs:
        graphs_number += 1
        cwalks_number.append(nx.number_connected_components(graph))
        for cwalk in list(nx.connected_components(graph)):
            cwalk = list(cwalk)
            if cwalk[0][2] == "chrY":
                out_tad += 1
            else:
                cwalk_boundaries = [((node[0] + node[1])/2) for node in cwalk]
                itv_tad = [tree_dict[cwalk[0][2]][bound] for bound in cwalk_boundaries]
                identical = all(element == itv_tad[0] for element in itv_tad)
                if identical:
                    in_tad += 1
                else:
                    out_tad += 1

    print(f"In summary from each {graphs_number} graphs, we have {sum(cwalks_number)} cwalks from which {out_tad}"
          f" are outside TAD boundaries")


if __name__ == '__main__':
    main()

# In summary from each 37 graphs, we have 319577 cwalks from which 160346 are outside TAD boundaries
