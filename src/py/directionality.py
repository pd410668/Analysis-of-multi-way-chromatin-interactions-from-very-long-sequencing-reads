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
    directions_avg = []
    directions = []
    numbers = []
    for graph in graphs:
        signs = [directionality(list(cwalk)) for cwalk in list(nx.connected_components(graph)) if len(cwalk) == length]
        avg_signs = [avg_directionality(sign) for sign in signs]
        directions_avg.extend(avg_signs)
        flat_signs = [sign for sublist in signs for sign in sublist]
        directions.extend(flat_signs)

    avg_sign = avg_directionality(directions)  # list of average value of directivity
    hist(directions_avg, length, avg_sign)
