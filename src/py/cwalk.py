#!/usr/bin/env python

import networkx as nx
import pandas as pd
from filtering import typical_chromosomes
import matplotlib.pyplot as plt
from intervaltree import IntervalTree
import pickle
import sys


def parse_positions(tsvfile: str) -> zip:
    """ Returns lists of positions of aligns that are apart => 500 """
    df = pd.read_csv(tsvfile, sep='\t')
    df = df.where(df.abs_pos >= 500).dropna().reset_index(drop=True)
    return zip(df.pos_R1.tolist(), df.pos_R2.tolist())


def read_bedfile(bedfile: str) -> list:
    """ Return lists of restrictions sites positions and chromosomes where they were found """
    df = pd.read_csv(bedfile, sep="\t", header=None)
    df = df.loc[df[0].isin(typical_chromosomes())].reset_index(drop=True)
    chr = df[0].tolist()
    pos = df[1].tolist()
    return pos, chr


def add_edge(u, v):
    """ Graph construction """
    if G.has_edge(u, v):
        G[u][v]["weight"] += 1
    else:
        G.add_edge(u, v, weight=1)


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


if __name__ == '__main__':
    
    """ Initiation basic dependencies """
    positions = parse_positions(sys.argv[1])  # .tsv file
    restrictions, chromosomes = read_bedfile(sys.argv[2])  # .bed file
    G = nx.Graph()

    """ Interval tree construction """
    intervals = [(i, j, chr) for i, j, chr in zip(restrictions[:-1], restrictions[1:], chromosomes[:-1]) if i <= j]
    tree = IntervalTree.from_tuples(intervals)

    for position_R1, position_R2 in positions:

        left_edge = tree[position_R1]
        right_edge = tree[position_R2]

        for i in list(left_edge):
            for j in list(right_edge):
                if i != j:  # prevention of self-loops
                    add_edge(i, j)

    """  C-walks construction """
    sorted_edges = sorted(G.edges(data=True), key=lambda x: x[2]["weight"], reverse=True)  # Sort edges by read-coverage
    P = cwalk(sorted_edges)
    pickle.dump(P, open(sys.argv[3], "wb"))  # save as .txt outfile in binary mode
