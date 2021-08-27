#!/usr/bin/env python

import networkx as nx
import pandas as pd
from filtering import typical_chromosomes, collect_data
from intervaltree import IntervalTree
import pickle
import sys


def parse_positions(tsvfile: str, abs_threshold: int) -> zip:
    """ Returns lists of positions of aligns that are apart selected absolute threshold """
    df = pd.read_csv(tsvfile, sep='\t')
    df = df.where(df.abs_pos >= abs_threshold).dropna().reset_index(drop=True)
    return zip(df.chr.tolist(), df.pos_R1.tolist(), df.pos_R2.tolist())


def read_bedfile(bedfile: str, chromosome: str) -> list:
    """ Return lists of restrictions sites positions and chromosomes where they were found """
    df = pd.read_csv(bedfile, sep="\t", header=None)
    df = df.loc[df[0].isin(typical_chromosomes())].reset_index(drop=True)
    df = df.where(df[0] == chromosome).dropna().reset_index(drop=True)
    df[1] = df[1].astype(int)
    chr = df[0].tolist()
    pos = df[1].tolist()
    return pos, chr


def add_edge(u, v):
    """ Creating weighted edges """
    if G.has_edge(u, v):
        G[u][v]["weight"] += 1
    else:
        G.add_edge(u, v, weight=1)


def matching_edges(interval_tree_dict: dict, positions: zip):
    """ Graph construction:
    interval_tree: a dictionary storing the restriction intervals (as IntervalTree) for each chromosome
    positions: list of C-walk positions
    """
    for chr, position_R1, position_R2 in positions:

        left_edge = interval_tree_dict[chr][position_R1]
        right_edge = interval_tree_dict[chr][position_R2]

        right_edge = list(list(right_edge)[0])
        left_edge = list(list(left_edge)[0])
        right_edge[2], left_edge[2] = chr, chr

        if (left_edge and right_edge) and (left_edge != right_edge):  # prevention of empty sets and self-loops
            add_edge(tuple(left_edge), tuple(right_edge))  # ex. (77366342, 77367727, "chr1")


def cwalk(edges):
    """ Resolve the C-walk graph """
    P = nx.Graph()
    for u, v, a in edges:
        if (u not in P or P.degree[u] < 2) and (v not in P or P.degree[v] < 2):
            P.add_edge(u, v, weight=a["weight"])

    for cwalk in list(nx.connected_components(P)):
        if len(cwalk) < 3:
            for node in cwalk:
                P.remove_node(node)  # Remove cwalks that are include one hop
    return P


def save_as_bed(graph, experiment_name):
    numbers = list(range(1, nx.number_connected_components(graph) + 1))
    order = [2, 0, 1]
    for cwalk, number in zip(list(nx.connected_components(graph)), numbers):
        cwalk_length = range(1, len(cwalk) + 1)

        for node, length in zip(cwalk, cwalk_length):
            node = [node[i] for i in order]
            node.append(f"{length}_{number}")
            collect_data(node, f"{experiment_name}", "a")


def main():
    tree_dict = dict()  # ex. tree_dict["chr1"] will be an object of type IntervalTree
    for chr in typical_chromosomes():
        """ Interval tree construction, separate for each chromosome """
        restrictions, chromosomes = read_bedfile(sys.argv[2], chr)  # .bed file
        intervals = [(i, j) for i, j in zip(restrictions[:-1], restrictions[1:])]
        tree_dict[chr] = IntervalTree.from_tuples(intervals)

    """ Parse C-walk positions """
    positions = parse_positions(sys.argv[1], 500)  # .tsv file
    matching_edges(tree_dict, positions)

    """  C-walks construction """
    sorted_edges = sorted(G.edges(data=True), key=lambda x: x[2]["weight"], reverse=True)  # Sort edges by read-coverage
    P = cwalk(sorted_edges)
    save_as_bed(P, sys.argv[3])  # save as .bed outfile with cwalks
    pickle.dump(P, open(sys.argv[4], "wb"))  # save cwalks as .txt outfile in binary mode


if __name__ == '__main__':
    G = nx.Graph()
    main()
