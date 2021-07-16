import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from filtering import typical_chromosomes
import sys


def parse_positions(tsvfile: str, chr: str) -> zip:
    """ Returns lists of positions of aligns of selected chromosome that are apart => 500 """
    df = pd.read_csv(tsvfile, sep='\t')
    clean_df = df.where(df.abs_pos >= 500).dropna().reset_index(drop=True)  # Do I need to rm aligns with abs_pos<500?
    df_chrs = {x: y for x, y in clean_df[["chr", "pos_R1", "pos_R2"]].groupby("chr")}
    chrs = df_chrs[chr].reset_index(drop=True)
    return zip(chrs.pos_R1.tolist(), chrs.pos_R1.tolist())


def read_headers(filename: str, chr: str) -> list:
    """ return list of headers of selected chromosome from fasta file """
    from Bio import SeqIO
    from filtering import typical_chromosomes
    return [record.id.split("-")[0] for record in SeqIO.parse(filename, "fasta") if record.id.split(":")[0] == chr]


def add_edge(u, v):
    """ Graph construction """
    if G.has_edge(u, v):
        G[u][v]["weight"] += 1
    else:
        G.add_edge(u, v, weight=1)


def graph_construction(positions: zip):
    right_edge = ""
    left_edge = ""
    for position_R1, position_R2 in positions:
        for i, j in zip(headers[:-1], headers[1:]):  # Iterating through all intervals
            res_pos_beg = int(i.split(":")[1])  # extracting left boundary value of interval
            res_pos_end = int(j.split(":")[1])  # extracting right boundary value of interval

            if res_pos_beg <= res_pos_end:
                interval = pd.Interval(res_pos_beg, res_pos_end, closed="both")
            
                if position_R1 in interval:
                    right_edge = f"{i}-{j}"
                if position_R2 in interval:
                    left_edge = f"{i}-{j}"
                    
            """  I can consider those intervals as irrelevant
            else:   
                interval = pd.Interval(res_pos_beg, res_pos_beg + 4, closed="both")
            
                if position_R1 in interval:
                    right_edge = f"{i}-{j}"
                if position_R2 in interval:
                    left_edge = f"{i}-{j}"
            """ 

        if right_edge and left_edge:
            add_edge(right_edge, left_edge)


G = nx.Graph()
for chr in typical_chromosomes():
    headers = read_headers("DpnII_hg19_seq.fa", chr)
    positions = parse_positions("hs_k562_I_1.tsv", chr)
    graph_construction(positions)

    
"""
Find the node-pair for which its read coverage is maximal,
and define this coverage as R
"""

R = (max(dict(G.edges).items(), key=lambda x: x[1]["weight"]))[1].get("weight")
print(R)

""" 
Eliminate all other appearances of the nodes-pairs
if their coverage is less than 0.1 × R 
"""

edge_weights = nx.get_edge_attributes(G, "weight")
G.remove_edges_from((edge for edge, weight in edge_weights.items() if weight < 0.1 * R))

print(G.edges(data=True))
