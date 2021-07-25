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
    df = pd.read_csv(bedfile, sep="\t", header=None)x
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


""" Initiation basic dependencies """

positions = parse_positions(sys.argv[1])  # .tsv file
restrictions, chromosomes = read_bedfile(sys.argv[2])  # .bed file
G = nx.Graph()

""" Interval tree construction """

intervals = [(i, j, chr) for i, j, chr in zip(restrictions[:-1], restrictions[1:], chromosomes[:-1]) if i <= j]
tree = IntervalTree.from_tuples(intervals)


def graph_construction(tsvpositions: zip):
    for position_R1, position_R2 in tsvpositions:

        right_edge = tree[position_R1]
        left_edge = tree[position_R2]

        if right_edge and left_edge:
            add_edge(repr(right_edge)[9:-1], repr(right_edge)[9:-1])  # ex. (50256114, 50257366, 'chr1')


graph_construction(positions)

"""
Find the node-pair for which its read coverage is maximal,
and define this coverage as R
"""

R = (max(dict(G.edges).items(), key=lambda x: x[1]["weight"]))[1].get("weight")

"""
Eliminate all other appearances of the nodes-pairs
if their coverage is less than 0.1 × R
"""

edge_weights = nx.get_edge_attributes(G, "weight")
G.remove_edges_from((edge for edge, weight in edge_weights.items() if weight < 0.1 * R))

""" Remove loops """

G.remove_edges_from(nx.selfloop_edges(G))

""" Write to .txt file graph in binary mode """

pickle.dump(G, open(sys.argv[3], "wb"))  # .txt outfile
