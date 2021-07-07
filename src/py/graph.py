import networkx as nx
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

to_remove = []
for edge in edges:
    if len(edge) == 1:
        to_remove.append(edge)
edges = [edges.remove(i) for i in to_remove]

for rm in to_remove:
    if rm[0] in nodes:
        index = nodes.index(rm[0])
        nodes[index] = None

for pos_R1, pos_R2 in parse_positions(sys.argv[2]):  # .tsv
    for i in range(len(restriction_intervals)):
        if pos_R1 and pos_R2 in restriction_intervals[i]:
            edges[i].append(nodes[i])

graph_df = pd.DataFrame({"from": nodes, "to": edges}).dropna()

G = nx.from_pandas_edgelist(graph_df.explode("to"), source="from", target="to")
nx.draw(G)
plt.savefig("graph.png")

""" Histogram of how many edges each node has (degree distribution) """
# plt.hist([v for k, v in nx.degree(G)])
# plt.show()

