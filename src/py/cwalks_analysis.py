#!/usr/bin/env python

import matplotlib.ticker as ticker
from matplotlib.ticker import MaxNLocator
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
import seaborn as sns
from filtering import typical_chromosomes, collect_data
import statistics
import numpy as np
import sys


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


def histogram(data):
    sns.set_style("whitegrid")
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(18, 8))
    plt.suptitle("Distribution of c-walks length", fontsize=20)
    ax1.hist(data, color="firebrick", edgecolor="black")
    ax2.hist([x for x in data if x <= 6], color="firebrick", edgecolor="black")
    ax1.yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, pos: "{0:g}".format(x / 1000000) + "M"))
    ax2.yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, pos: "{0:g}".format(x / 1000000) + "M"))
    ax2.xaxis.set_major_locator(MaxNLocator(integer=True))
    ax1.set_xlabel("C-walk length", fontsize=16)
    ax1.set_ylabel("Frequency of occurrence", fontsize=16)
    ax2.set_xlabel("C-walk length", fontsize=16)
    ax2.set_ylabel("Frequency of occurrence", fontsize=16)
    return plt.savefig(f"{sys.argv[2]}"), plt.close()


def barh(data, labels):
    sns.set_style("whitegrid")
    fig, ax = plt.subplots(figsize=(12, 10))
    ax.barh([label[:-11] for label in labels], data, color="tab:blue", edgecolor="black")
    ax.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, pos: "{0:g}".format(x / 1000) + "k"))
    ax.set_title("Number of c-walks in each graph", fontsize=18)
    return plt.savefig(f"{sys.argv[3]}"), plt.close()


def main():
    graphs, labels = load_files(sys.argv[1], load_cwalk_graph)  # "./cwalks"

    cwalks_length, cwalks_number = [], []
    for graph in graphs:
        cwalks_number.append(nx.number_connected_components(graph))
        for cwalk in list(nx.connected_components(graph)):
            cwalks_length.append(len(list(cwalk)))

    # Save as .tsv basic statistics
    collect_data(["average", "median", "mode", "standard deviation"], "cwalks.tsv", "w")
    collect_data([round(statistics.mean(cwalks_length), 2),
                  round(statistics.median(cwalks_length), 2),
                  statistics.mode(cwalks_length),
                  round(statistics.stdev(cwalks_length), 2)], sys.argv[4], "a")  # "cwalks.tsv"

    histogram(cwalks_length)
    barh(cwalks_number, labels)


if __name__ == '__main__':
    main()
