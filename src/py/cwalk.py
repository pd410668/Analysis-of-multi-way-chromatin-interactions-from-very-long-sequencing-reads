#!/usr/bin/env python

import networkx as nx
import pandas as pd
from filtering import typical_chromosomes, collect_data
from intervaltree import IntervalTree
import pickle
import sys


def parse_positions(tsvfile: str, abs_threshold: int) -> zip:
    """ Returns lists of positions of aligns """
    df = pd.read_csv(tsvfile, sep="\t")
    df = df.where(df.abs_pos >= abs_threshold).dropna().reset_index(drop=True)
    return zip(df.chr_R1.tolist(), df.chr_R2.tolist(), df.pos_R1.astype(int).tolist(), df.pos_R2.astype(int).tolist())


def parse_bedfile(bedfile: str, organism: str) -> dict:
    """ Return lists of restrictions sites positions and chromosomes where they were found """
    df = pd.read_csv(bedfile, sep="\t", header=None)
    df = df.iloc[:, :2]
    df = df.loc[df[0].isin(typical_chromosomes(organism))].reset_index(drop=True)
    return {x: y for x, y in df.groupby(0)}


def add_weighted_edges(u, v, node1, node2):
    """ Creating weighted edges """
    if G.has_edge(u, v):
        G[u][v]["weight"] += 1
    else:
        G.add_edge(u, v, weight=1)
        nx.set_edge_attributes(G, {(u, v): {"node1": node1, "node2": node2}})


def matching_edges(interval_tree_dict: dict, positions: zip):
    """ 
    interval_tree: a dictionary storing the restriction intervals (as IntervalTree) for each chromosome
    positions: list of C-walk positions
    """
    for chr1, chr2, position_R1, position_R2 in positions:
        left_edge = interval_tree_dict[chr1][position_R1]
        right_edge = interval_tree_dict[chr2][position_R2]

        if left_edge and right_edge:  # prevention of empty sets
            left_edge = tuple(list(left_edge)[0])[0:2]
            right_edge = tuple(list(right_edge)[0])[0:2]
            add_weighted_edges(left_edge, right_edge, chr1, chr2)


def cwalk_construction(G):
    """ Resolve the C-walk graph """
    sorted_edges = sorted(G.edges(data=True), key=lambda x: x[2]["weight"], reverse=True)
    P = nx.Graph()
    for u, v, a in sorted_edges:
        node1 = a["node1"]
        node2 = a["node2"]
        if (u not in P or P.degree[u] < 2) and (v not in P or P.degree[v] < 2):
            P.add_edge(u, v, weight=a["weight"])
            nx.set_edge_attributes(P, {(u, v): {"node1": node1, "node2": node2}})

    for cwalk in list(nx.connected_components(P)):
        if len(cwalk) < 3:
            for node in cwalk:
                P.remove_node(node)  # Remove cwalks that are include one hop
    return P


def main():
    organism = sys.argv[1]
    restrictions_bedfile = sys.argv[2]

    restrictions_dict = parse_bedfile(restrictions_bedfile, organism)  # .bed file
    tree_dict = dict()  # ex. tree_dict["chr1"] will be an object of type IntervalTree
    for chr in typical_chromosomes(organism):
        """ Interval tree construction, separate for each chromosome """
        restrictions = restrictions_dict[chr][1].tolist()
        intervals = [(i, j) for i, j in zip(restrictions[:-1], restrictions[1:])]
        tree_dict[chr] = IntervalTree.from_tuples(intervals)

    positions = parse_positions(sys.argv[3], 500)  # .tsv file with selected absolute threshold
    matching_edges(tree_dict, positions)

    """  cwalks construction """
    P = cwalk_construction(G)
    P.remove_edges_from(nx.selfloop_edges(P))  # remove self-loops

    """ saving cwalks """
    pickle.dump(P, open(sys.argv[5], "wb"))  # save cwalks as .txt outfile in binary mode


if __name__ == '__main__':
    G = nx.Graph()
    main()
