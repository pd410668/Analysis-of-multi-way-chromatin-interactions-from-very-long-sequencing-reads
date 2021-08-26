#!/usr/bin/env python

import networkx as nx
import pandas as pd
from filtering import typical_chromosomes, collect_data
from intervaltree import IntervalTree
import sys


def parse_positions(tsvfile: str, abs_threshold: int) -> zip:
    """ Returns lists of positions of aligns that are apart selected absolute threshold """
    df = pd.read_csv(tsvfile, sep='\t')
    df = df.where(df.abs_pos >= abs_threshold).dropna().reset_index(drop=True)
    return zip(df.pos_R1.tolist(), df.pos_R2.tolist())


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

        
def matching_edges(interval_tree):
    """ Graph construction """
    for position_R1, position_R2 in positions:

        left_edge = interval_tree[position_R1]
        right_edge = interval_tree[position_R2]

        if (left_edge and right_edge) and (left_edge != right_edge):  # prevention of empty sets and self-loops
            add_edge(tuple(list(left_edge)[0]), tuple(list(right_edge)[0]))  # ex. (77366342.0, 77367727.0, 'chr1')
            # print(tuple(list(left_edge)[0]), tuple(list(right_edge)[0]))
            

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


def save_as_bed(graph):
    order = [2, 0, 1]
    for cwalk in list(nx.connected_components(graph)):
        cwalk_length = range(1, len(cwalk) + 1)
        for node, length in zip(cwalk, cwalk_length):
            node = [node[i] for i in order]
            node.append(length)
            collect_data(node, "cwalk.bed", "a")
            
            
if __name__ == '__main__':

    G = nx.Graph()

    for chr in typical_chromosomes():
        positions = parse_positions("hs_k562_I_1.tsv", 500)  # .tsv file
        restrictions, chromosomes = read_bedfile("DpnII_hg19.bed", chr)  # .bed file

        """ Interval tree construction """  # separate for each chromosome
        intervals = [(i, j, chr) for i, j, chr in zip(restrictions[:-1], restrictions[1:], chromosomes[:-1]) if i <= j]
        tree = IntervalTree.from_tuples(intervals)
        matching_edges(tree)

    """  C-walks construction """
    sorted_edges = sorted(G.edges(data=True), key=lambda x: x[2]["weight"], reverse=True)  # Sort edges by read-coverage
    P = cwalk(sorted_edges)
    save_as_bed(P)  # save as .bed outfile with cwalks
    pickle.dump(P, open("cwalk.txt", "wb"))  # save as .txt outfile in binary mode
