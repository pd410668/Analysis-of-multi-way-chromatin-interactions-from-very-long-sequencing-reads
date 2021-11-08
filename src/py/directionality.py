#!/usr/bin/env python

import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
import seaborn as sns
from filtering import typical_chromosomes, collect_data
import sys
import scipy.stats
import statistics
import matplotlib.ticker as ticker
from cwalks_analysis import load_cwalk_graph, load_files


def permutation(path: list) -> list:
    import itertools
    import random
    return list(random.choice(list(itertools.permutations(path))))


def directionality(path: list) -> bool:
    path = [node[0] for node in path]  # the distance is between the firsts points of the restriction intervals
    dist = [b - a for a, b in zip(path[:-1], path[1:])]
    same_sign_dist = [(a * b >= 0) for a, b in zip(dist[:-1], dist[1:])]
    # 0 if nodes are in the same restriction interval
    return same_sign_dist


def avg_directionality(signs: list) -> float:
    return round(statistics.mean(signs), 2)


def hist(data, name, avg):
    sns.set_style("whitegrid")
    fig, ax = plt.subplots(figsize=(10, 8))
    plt.hist(data, color="cornflowerblue", edgecolor="black")
    if len(data) >= 1000:
        ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, pos: "{0:g}".format(x / 1000) + "k"))
    ax.set_title(f"Distribution of directivity of cwalks with {name} length. \n "
                 f"Average value of directivity is {avg}.", fontsize=18)
    return plt.savefig(f"directivity_{name}.png"), plt.close()


graphs, _ = load_files("./cwalks", load_cwalk_graph)  # load .txt cwalk graph

duplicate_cwalk_length = []
for graph in graphs:
    for cwalk in list(nx.connected_components(graph)):
        duplicate_cwalk_length.append(len(cwalk))
cwalk_length = list(set(duplicate_cwalk_length))  # list with all possible length of cwalks in each graph

for each in cwalk_length:
    print(each, duplicate_cwalk_length.count(each))

for length in range(min(cwalk_length), max(cwalk_length) + 1):
    print(length)

    avg_signs_shuffle = []
    avg_signs = []

    for graph in graphs:
        for cwalk in list(nx.connected_components(graph)):
            if len(cwalk) == length:
                cwalk = list(cwalk)
                
                sign = directionality(cwalk)
                avg_sign = avg_directionality(sign)
                
                random.shuffle(cwalk)
                sign_shuffle = directionality(cwalk)
                avg_sign_shuffle = avg_directionality(sign_shuffle)
                
                avg_signs.append(avg_sign)
                avg_signs_shuffle.append(avg_sign_shuffle)

    if avg_signs:

        print(round(statistics.mean(avg_signs), 2))  # 0.33
        print(round(statistics.mean(avg_signs_shuffle), 2))  # 0.33
        hist(avg_signs_shuffle, length, round(statistics.mean(avg_signs_shuffle), 2))
        
"""
3
0.34
0.33
4
0.33
0.33
5
0.33
0.33
6
0.33
0.33
7
0.33
0.33
8
0.33
0.32
9
0.33
0.33
10
0.33
0.34
11
0.34
0.33
12
0.33
0.35
13
0.34
0.34
14
0.32
0.33
15
0.34
0.32
16
0.35
0.31
17
0.32
0.36
18
0.34
0.34
19
0.33
0.32
20
0.33
0.35
21
0.35
0.36
22
0.31
0.37
23
0.35
0.34
24
0.41
0.42
25
0.43
0.3
26
27
28
0.38
0.27
29
0.35
0.39
30
0.39
0.39
"""
