#!/usr/bin/env python

from cwalks_analysis import load_cwalk_graph, load_files
import pandas as pd
import networkx as nx
from intervaltree import Interval, IntervalTree
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np


def barplot(hops: list):
    sns.set_style("whitegrid")
    fig, ax = plt.subplots(figsize=(12, 10))
    bar_4 = ax.bar(hops, many_tads_lst, bottom=np.add(np.add(one_tad_lst, two_tads_lst), three_tads_lst),
                   color="slategrey", edgecolor="black", width=1)
    bar_3 = ax.bar(hops, three_tads_lst, bottom=np.add(one_tad_lst, two_tads_lst), color="lightsteelblue",
                   edgecolor="black", width=1)
    bar_2 = ax.bar(hops, two_tads_lst, bottom=one_tad_lst, color="cornflowerblue", edgecolor="black", width=1)
    bar_1 = ax.bar(hops, one_tad_lst, color="royalblue", edgecolor="black", width=1)
    plt.xlabel("Number of hops")
    plt.ylabel("Percentage")
    ax.set_title("Inter-TAD C-walks", fontsize=18)
    plt.legend([bar_1, bar_2, bar_3, bar_4], ["1 TAD", "2 TAD", "3 TAD", "many TADs"], loc="upper right",
               prop={"size": 16})
    return plt.savefig("barplot.png")


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


def identical(itv):
    in_tad = 1
    while not all(itv[0] == element for element in itv):
        if all(element == itv[0] for element in itv):
            in_tad += 1
        else:
            a = itv[0]
            itv = [value for value in itv if value != a]
            in_tad += 1
    return in_tad


graphs, _ = load_files("./input_cwalks/txt", load_cwalk_graph)  # load .txt cwalk graph
boundaries_dict = tad_boundaries("K562_Lieberman-raw_TADs.bed")  # list of tuples

tree_dict = dict()  # ex. tree_dict["chr1"] will be an object of type IntervalTree
for key in boundaries_dict.keys():
    intervals = boundaries_dict[key]
    # Interval tree construction, separate for each chromosome
    tree_dict[key] = IntervalTree.from_tuples(intervals)


def cwalks_in_tad():
    count = 0
    cwalks_number = []
    one_tad, two_tads, three_tads, many_tads = 0, 0, 0, 0
    for graph in graphs:
        cwalks_number.append(nx.number_connected_components(graph))
        for cwalk in list(nx.connected_components(graph)):
            cwalk = list(cwalk)
            if all(cwalk[i][2] == cwalk[0][2] for i in range(0, len(cwalk))):  # inter or intra chrs cwalks
                if len(cwalk) == hop:
                    count += 1
                    if cwalk[0][2] == "chrY":
                        many_tads += 1
                    else:
                        cwalk_boundaries = [((node[0] + node[1]) / 2) for node in cwalk]
                        itv_tad = [tree_dict[cwalk[0][2]][bound] for bound in cwalk_boundaries]
                        in_tad = identical(itv_tad)
                        if in_tad == 1:
                            one_tad += 1
                        elif in_tad == 2:
                            two_tads += 1
                        elif in_tad == 3:
                            three_tads += 1
                        else:
                            many_tads += 1

    print(hop, count)  # number of cwalks of particular length
    print(one_tad, two_tads, three_tads, many_tads)
    return round((one_tad / count) * 100, 2), round((two_tads / count) * 100, 2), \
           round((three_tads / count) * 100, 2), round((many_tads / count) * 100, 2)


def main():
    one_tad_lst, two_tads_lst, three_tads_lst, many_tads_lst = [], [], [], []
    for hop in range(3, 12):
        one_tad, two_tads, three_tads, many_tads = cwalks_in_tad()
    
        one_tad_lst.append(one_tad)
        two_tads_lst.append(two_tads)
        three_tads_lst.append(three_tads)
        many_tads_lst.append(many_tads)
    
    print(one_tad_lst)  # list with percentage intra chrs cwalks, from 3 to 18
    print(two_tads_lst)
    print(three_tads_lst)
    print(many_tads_lst)
    barplot([hop for hop in range(3, 12)])


if __name__ == '__main__':
    main()
