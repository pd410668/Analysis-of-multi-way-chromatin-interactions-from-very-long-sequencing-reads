import networkx as nx
<<<<<<< HEAD
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import sys

def parse_positions(tsvfile) -> zip:
    """ Returns lists of positions of aligns that are apart => 500 """
    df = pd.read_csv(tsvfile, sep='\t')
    clean_df = df.where(df.abs_pos >= 500).dropna().reset_index(drop=True)
    return zip(clean_df.pos_R1.tolist(), clean_df.pos_R2.tolist())

def read_headers(filename: str) -> list:
    """ return list of headers from fasta file """
    from Bio import SeqIO
    from filtering import typical_chromosomes
    return [record.id.split("-")[0] for record in SeqIO.parse(filename, "fasta") if record.id.split(":")[0] in typical_chromosomes()]

def intervals(headers):
    """ Get restriction intervals from fasta headers """
    restriction_intervals = []
    for i, j in zip(headers[:-1], headers[1:]):
        res_pos_beg = int(i.split(":")[1])
        res_pos_end = int(j.split(":")[1])
        if res_pos_beg <= res_pos_end:
            interval = pd.Interval(res_pos_beg, res_pos_end - 1, closed="both")
            restriction_intervals.append(interval)
        else:
            interval = pd.Interval(res_pos_beg, res_pos_beg + 4, closed="both")
            restriction_intervals.append(interval)
        return restriction_intervals

headers = read_headers(sys.argv[1])  # fasta 
restriction_intervals = intervals(headers)

nodes = headers[:-1]          # list
edges = [[i] for i in nodes]  # list of lists

for pos_R1, pos_R2 in parse_positions(sys.argv[2]):  # .tsv
    for i in range(len(restriction_intervals)):
        if pos_R1 and pos_R2 in restriction_intervals[i]:
            edges[i].append(nodes[i])

print(edges) 

graph_df = pd.DataFrame({"from": nodes, "to": edges})
print(graph_df.explode("to"))

G = nx.from_pandas_edgelist(graph_df.explode("to"), source="from", target="to")
nx.draw(G)
plt.savefig("graph.png")

""" Histogram of how many edges each node has (degree distribution) """
# plt.hist([v for k, v in nx.degree(G)])
# plt.show()







=======
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

pickle.dump(G, open("graph.txt", "wb"))
>>>>>>> 8c88a83f444b5e0113c0165457dfeb7f6e195201
