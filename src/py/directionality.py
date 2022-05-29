#!/usr/bin/env python

import networkx as nx
import matplotlib.pyplot as plt
import seaborn as sns
import sys
import statistics
import matplotlib.ticker as ticker
from analysis import load_cwalk_graph, load_files
import random


def directionality(path: list) -> bool:
    path = [node[0] for node in path]  # the distance is between the firsts points of the restriction intervals
    dist = [b - a for a, b in zip(path[:-1], path[1:])]
    same_sign_dist = [(a * b >= 0) for a, b in zip(dist[:-1], dist[1:])]
    # 0 if nodes are in the same restriction interval
    return same_sign_dist


def avg_directionality(signs: list) -> float:
    return round(statistics.mean(signs), 2)


def hist(data, avg, len, name):
    sns.set_style("whitegrid")
    fig, ax = plt.subplots(figsize=(10, 8))
    plt.hist(data, color="royalblue", edgecolor="black")
    ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, pos: "{0:g}".format(x / 1000) + "k"))
    ax.set_title(f"Distribution of directionality of cwalks of {len} length.\n "
                 f"Average value of directivity is {avg}.", fontsize=18)
    return plt.savefig(f"directionality_{len}_{name}.png"), plt.close()


def main():
    graphs, _ = load_files(sys.argv[1], load_cwalk_graph)  # load .txt cwalk graph # "./cwalks"

    duplicate_cwalk_length = []
    for graph in graphs:
        for cwalk in list(nx.connected_components(graph)):
            cwalk = list(cwalk)
            if all(cwalk[i][2] == cwalk[0][2] for i in range(0, len(cwalk))):  # intra chrs cwalks
                duplicate_cwalk_length.append(len(cwalk))
    cwalk_length = list(set(duplicate_cwalk_length))  # list with all possible length of cwalks in each graph

    for length in range(min(cwalk_length), max(cwalk_length) + 1):

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

            hist(avg_signs_shuffle, round(statistics.mean(avg_signs_shuffle), 2), length, "shuffle")
            hist(avg_signs, round(statistics.mean(avg_signs), 2), length, "data")


if __name__ == '__main__':
    main()
