#!/usr/bin/env python

import networkx as nx
import pandas as pd
import pickle
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.ticker as ticker
from matplotlib.ticker import MaxNLocator
import numpy as np
import statistics
import sys


import matplotlib.pylab as pylab
params = {'legend.fontsize': 'x-large',
         'axes.labelsize': 'x-large',
         'axes.titlesize':'x-large',
         'xtick.labelsize':'xx-large',
         'ytick.labelsize':'xx-large'}
pylab.rcParams.update(params)


def load_cwalk_graph(cwalk_graph):
    import pickle
    return pickle.load(open(cwalk_graph, "rb"))


def load_files(dictionary, parser):
    """
    load all available files from a folder
    and return lists of objects
    """
    import os
    files, labels = [], []
    for file in os.listdir(dictionary):
        labels.append(file)
        files.append(parser(f"{dictionary}/{file}"))
    return files, labels


def fractions(one_chr, two_chr, many_chr, name):
    sns.set_style("whitegrid")
    plt.figure(figsize=(12, 10))
    bar1 = plt.bar([i for i in range(3, 16)], one_chr, width=1, edgecolor="black", color="tab:blue")
    bar2 = plt.bar([i for i in range(3, 16)], two_chr, bottom=one_chr, width=1, edgecolor="black", color="tab:cyan")
    bar3 = plt.bar([i for i in range(3, 16)], many_chr,
                   bottom=np.add(one_chr, two_chr), width=1, edgecolor="black", color="tab:gray")
    plt.xlabel("Number of hops", fontsize=16)
    plt.ylabel("Percentage [%]", fontsize=16)
    plt.title(f"Fraction of c-walks in each class in {sys.argv[1]} cells", fontsize=18)
    plt.legend([bar1, bar2, bar3], ["one chromosome (class I)",
                                    "two chromosomes (class II)",
                                    "many chromosomes (class III)"])
    return plt.savefig(f"{name}"), plt.close() 


def barh(graphs, labels, name):
    data = []
    for graph in graphs:
        data.append(nx.number_connected_components(graph))
    sns.set_style("whitegrid")
    fig, ax = plt.subplots(figsize=(12, 16))
    ax.barh([label[:-11] for label in labels], data, color="tab:blue", edgecolor="black", height=1)
    ax.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, pos: "{0:g}".format(x / 1000) + "k"))
    ax.set_title(f"Number of c-walks in each graph in {sys.argv[1]} cells", fontsize=18)
    return plt.savefig(f"{name}"), plt.close()


def stats(graphs, name):
    i, j, k = 0, 0, 0
    lengths = []
    for graph in graphs:
        for cwalk in list(nx.connected_components(graph)):
            cwalk = list(cwalk)
            k += 1
            if all(cwalk[i][2] == cwalk[0][2] for i in range(0, len(cwalk))):
                lengths.append(len(cwalk))
                i += 1
            else:
                j += 1
                

    print(f"Number of intra-chromosomal cwalks: {i}")
    print(f"Number of inter-chromosomal cwalks: {j}")
    print(f"Number of all cwalks: {k}")
    print(f"Percentage of intra-chromosomal cwalks is {round((i / k) * 100, 3)}%")
    print(f"Average: {round(statistics.mean(lengths), 2)}")
    print(f"Median: {round(statistics.median(lengths), 2)}")
    print(f"Mode: {statistics.mode(lengths)}")
    print(f"Standard deviation: {round(statistics.stdev(lengths), 2)}")
    sns.set_style("whitegrid")
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))
    plt.suptitle(f"Distribution of c-walks length in {sys.argv[1]} cells", fontsize=20)
    ax1.hist(lengths, color="tab:blue", edgecolor="black")
    ax2.hist([x for x in lengths if x <= 6], color="tab:blue", edgecolor="black", rwidth=1)
    ax1.yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, pos: "{0:g}".format(x / 1000) + "k"))
    ax2.yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, pos: "{0:g}".format(x / 1000) + "k"))
    ax2.xaxis.set_major_locator(MaxNLocator(integer=True))
    ax1.set_xlabel("c-walks length", fontsize=16)
    ax1.set_ylabel("Number of c-walks", fontsize=16)
    ax2.set_xlabel("c-walks length", fontsize=16)
    ax2.set_ylabel("Number of c-walks", fontsize=16)
    return plt.savefig(f"{name}"), plt.close()


def identical(cwalk):
    span = 1
    chrs = [node[2] for node in cwalk]
    while not all(chrs[0] == element for element in chrs):
        if not all(chrs[0] == element for element in chrs):
            span += 1
            repeat = chrs[0]
            chrs = [chr for chr in chrs if chr != repeat]
    return span


def counting(graphs, length):
    i = 0
    one, two, many = 0, 0, 0
    for graph in graphs:
        for cwalk in list(nx.connected_components(graph)):
            cwalk = list(cwalk)
            if len(cwalk) == length:
                i += 1
                span = identical(cwalk)
                if span == 1:
                    one += 1
                elif span == 2:
                    two += 1
                elif span >= 3:
                    many += 1

    return [round((one / i) * 100, 2), round((two / i) * 100, 2), round((many / i) * 100, 2)]


def main():
    graphs, labels = load_files(sys.argv[2], load_cwalk_graph)  # load cwalk graphs
    one_chr, two_chr, many_chr = [], [], []
    for hop in range(3, 16):
        one, two, many = counting(graphs, hop)
        one_chr.append(one)  # fraction of inter-chrs cwalks of particular length
        two_chr.append(two)
        many_chr.append(many)

    """ Plotting """
    stats(graphs, sys.argv[3])
    fractions(one_chr, two_chr, many_chr, sys.argv[4])
    barh(graphs, labels, sys.argv[5])


if __name__ == '__main__':
    main()
