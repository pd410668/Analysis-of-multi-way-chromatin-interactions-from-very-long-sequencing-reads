#!/usr/bin/env python

import networkx as nx
import matplotlib.pyplot as plt
import seaborn as sns
import sys
import statistics
import matplotlib.ticker as ticker
from cwalk_analysis import load_cwalk_graph, load_files
import random
import pandas as pd
import random


def directionality(path: list) -> bool:
    path = [(node[0]+node[1])/2 for node in path]  # the distance is between the middle point of the restriction interval
    dist = [b - a for a, b in zip(path[:-1], path[1:])]
    same_sign_dist = [(a * b >= 0) for a, b in zip(dist[:-1], dist[1:])]
    # 0 if nodes are in the same restriction interval
    # if the same node occurs more than one time it is irrelevant in the directionality concept
    # as well as the circle
    return same_sign_dist


def avg_directionality(signs: list) -> float:
    return round(statistics.mean(signs), 2)


def hist(data, data_random, avg, avg_random, length):
    sns.set_style("whitegrid")
    fig, axs = plt.subplots(1, 2, sharex=True, sharey=True, tight_layout=True, figsize=(15, 6))
    fig.suptitle(f"Distribution of directionality of c-walks of {length} length from {sys.argv[1]} cells", fontsize=14)
    axs[0].hist(data, color="tab:blue", edgecolor="black")
    axs[0].set_title(f"c-walks from data \n average value of directionality is {avg}", fontsize=14)
    axs[1].hist(data_random, color="tab:red", edgecolor="black")
    axs[1].set_title(f"randomize c-walks \n average value of directionality is {avg_random}", fontsize=14)
    axs[0].set_ylabel("number of c-walks", fontsize=14)
    fig.supxlabel("value of directionality", fontsize=14)
    axs[1].yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, pos: "{0:g}".format(x / 1000) + "k"))
    axs[0].yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, pos: "{0:g}".format(x / 1000) + "k"))
    return plt.savefig(f"{sys.argv[3]}"), plt.close()


def main():
    graphs, _ = load_files(sys.argv[2], load_cwalk_graph)  # load .txt cwalk graph
    cwalks = []

    for graph in graphs:
        S = [graph.subgraph(c).copy() for c in nx.connected_components(graph)]

        for subgraph in S:
            sorted_cwalk = sorted(subgraph.edges(data=True), key=lambda x: x[2]["index"])

            edges = []
            nodes = []
            for next_edge in sorted_cwalk[0:]:
                v1, v2, _ = next_edge
                edges.append([v1, v2])
                for i in range(1, len(edges)):
                    for edge in edges:
                        v1, v2 = edge
                        if v1 == (edges[i][0] or edges[i][1]) and v1 not in nodes:
                            nodes.append(v1)
                        elif v2 not in nodes:
                            nodes.append(v2)
            if nodes[0] != edges[0][0]:
                nodes.insert(0, edges[0][0])
            elif nodes[0] != edges[0][1]:
                nodes.insert(0, edges[0][1])

            if nodes[-1] == nodes[0]:  # remove repetitive nodes (in case of circle)
                nodes = nodes[1:]  # remove first node in circle

            if all(nodes[i][2] == nodes[0][2] for i in range(0, len(nodes))):  # only intra-chrs cwalks
                cwalks.append(nodes)


    avg_signs = []
    avg_shuffle_signs = []

    for cwalk in cwalks:

        avg_signs.append(avg_directionality(directionality(cwalk)))

        random.shuffle(cwalk)
        avg_shuffle_signs.append(avg_directionality(directionality(cwalk)))

    hist(avg_signs, avg_shuffle_signs, round(statistics.mean(avg_signs), 2),
         round(statistics.mean(avg_shuffle_signs), 2), "all")

    from scipy.stats import wilcoxon
    res = wilcoxon(x=avg_signs, y=avg_shuffle_signs, zero_method="wilcox", alternative="less")
    print(f"The value of the test statistic: {res.statistic} and p-value: {res.pvalue}")


if __name__ =='__main__':
    main()
