#!/usr/bin/env python

import networkx as nx
import matplotlib.pyplot as plt
import seaborn as sns
import sys
import statistics
import matplotlib.ticker as ticker
from cwalk_analysis import load_cwalk_graph, load_files
import random
from filtering import typical_chromosomes
import pandas as pd
import random

import matplotlib.pylab as pylab

params = {'legend.fontsize': 'x-large',
          'axes.labelsize': 'x-large',
          'axes.titlesize': 'x-large',
          'xtick.labelsize': 'xx-large',
          'ytick.labelsize': 'xx-large'}
pylab.rcParams.update(params)


def parse_tf(tf: str) -> list:
    """
    load transcription factor binding sites
    return: list of dicts (dict for each tf) where chomosomes are keys and values are peaks within them
    """
    tf = pd.read_csv(tf, sep='\t', header=None)
    tf = tf.iloc[:, 0:3]
    tf[3] = ((tf[1] + tf[2]) / 2).astype(float)
    tf.columns = tf.columns.map(str)
    tf_dict = tf.groupby("0")["3"].agg(list).to_dict()
    return tf_dict


def mirror_peaks(chr_sizes: str, peaks_dict: dict) -> dict:
    """
    Mirror image of peaks
    load .tsv file with chromosomes sizes
    return: dict where chromosomes are keys and values are "mirror peaks" within them
    """
    df_sizes = pd.read_csv(chr_sizes, sep='\t', header=None)
    df_sizes = df_sizes.loc[df_sizes[0].isin(typical_chromosomes("human"))].reset_index(drop=True)  # sys instead of str
    df_sizes.columns = df_sizes.columns.map(str)
    sizes_dict = df_sizes.groupby("0")["1"].agg(list).to_dict()

    keys, mirrors = [], []
    for key in peaks_dict.keys():
        mirror = [sizes_dict[key][0] - peak for peak in peaks_dict[key]]
        keys.append(key)
        mirrors.append(mirror)
    return dict(zip(keys, mirrors))


def count_boundaries(subgraph) -> list:
    # edge boundaries for one cwalk
    # subgraph is one cwalk
    edge_boundaries = []
    sorted_cwalk = sorted(subgraph.edges(data=True), key=lambda x: x[2]["index"])
    if all(i[2] == j[2] for i, j, _ in sorted_cwalk):  # only inter-chrs cwalks
        for next_edge in sorted_cwalk[0:]:
            v1, v2, _ = next_edge

            edge_bound = (pd.Interval(min([v1[0], v1[0], v2[0], v2[0]]),
                                      max([v1[0], v1[0], v2[0], v2[0]]), closed="both"))
            edge_boundaries.append((edge_bound, v1[2]))
    return edge_boundaries


def peak_in_edges(bounds: list, chromosome: str, tf_dict: dict) -> int:
    # number of edges which has at least one peak in one cwalk
    isin = 0
    for edge in bounds:
        if chromosome != "chrY":
            for peak in tf_dict[chromosome]:
                if peak in edge:
                    isin += 1
                    break
    return isin


def histogram(labels, fractions, mirror_fractions):

    sns.set_style("whitegrid")
    fig, axs = plt.subplots(4, 3, sharex=True, sharey=True, tight_layout=True, figsize=(16, 18))

    if statistics.mean(fractions[0]) < statistics.mean(mirror_fractions[0]): 
        axs[0, 0].hist(fractions[0], color="tab:blue", edgecolor="black")
        axs[0, 0].hist(mirror_fractions[0], color="tab:red", edgecolor="black")
    else:
        axs[0, 0].hist(mirror_fractions[0], color="tab:red", edgecolor="black")
        axs[0, 0].hist(fractions[0], color="tab:blue", edgecolor="black")

    axs[0, 0].set_title(f"Overlap with {labels[0]} binding sites", fontsize=14)

    if statistics.mean(fractions[1]) < statistics.mean(mirror_fractions[1]):
        axs[0, 1].hist(fractions[1], color="tab:blue", edgecolor="black")
        axs[0, 1].hist(mirror_fractions[1], color="tab:red", edgecolor="black")
    else:
        axs[0, 1].hist(mirror_fractions[1], color="tab:red", edgecolor="black")
        axs[0, 1].hist(fractions[1], color="tab:blue", edgecolor="black")

    axs[0, 1].set_title(f"Overlap with {labels[1]} binding sites", fontsize=14)

    if statistics.mean(fractions[2]) < statistics.mean(mirror_fractions[2]):
        axs[0, 2].hist(fractions[2], color="tab:blue", edgecolor="black")
        axs[0, 2].hist(mirror_fractions[2], color="tab:red", edgecolor="black")
    else:
        axs[0, 2].hist(mirror_fractions[2], color="tab:red", edgecolor="black")
        axs[0, 2].hist(fractions[2], color="tab:blue", edgecolor="black")

    axs[0, 2].set_title(f"Overlap with {labels[2]} binding sites", fontsize=14)

    if statistics.mean(fractions[3]) < statistics.mean(mirror_fractions[3]):
        axs[1, 0].hist(fractions[3], color="tab:blue", edgecolor="black")
        axs[1, 0].hist(mirror_fractions[3], color="tab:red", edgecolor="black")
    else:
        axs[1, 0].hist(mirror_fractions[3], color="tab:red", edgecolor="black")
        axs[1, 0].hist(fractions[3], color="tab:blue", edgecolor="black")

    axs[1, 0].set_title(f"Overlap with {labels[3]} binding sites", fontsize=14)

    if statistics.mean(fractions[4]) < statistics.mean(mirror_fractions[4]):
        axs[1, 1].hist(fractions[4], color="tab:blue", edgecolor="black")
        axs[1, 1].hist(mirror_fractions[4], color="tab:red", edgecolor="black")
    else:
        axs[1, 1].hist(mirror_fractions[4], color="tab:red", edgecolor="black")
        axs[1, 1].hist(fractions[4], color="tab:blue", edgecolor="black")

    axs[1, 1].set_title(f"Overlap with {labels[4]} binding sites", fontsize=14)

    if statistics.mean(fractions[5]) < statistics.mean(mirror_fractions[5]):
        axs[1, 2].hist(fractions[5], color="tab:blue", edgecolor="black")
        axs[1, 2].hist(mirror_fractions[5], color="tab:red", edgecolor="black")
    else:
        axs[1, 2].hist(mirror_fractions[5], color="tab:red", edgecolor="black")
        axs[1, 2].hist(fractions[5], color="tab:blue", edgecolor="black")

    axs[1, 2].set_title(f"Overlap with {labels[5]} binding sites", fontsize=14)

    if statistics.mean(fractions[6]) < statistics.mean(mirror_fractions[6]):
        axs[2, 0].hist(fractions[6], color="tab:blue", edgecolor="black")
        axs[2, 0].hist(mirror_fractions[6], color="tab:red", edgecolor="black")
    else:
        axs[2, 0].hist(mirror_fractions[6], color="tab:red", edgecolor="black")
        axs[2, 0].hist(fractions[6], color="tab:blue", edgecolor="black")

    axs[2, 0].set_title(f"Overlap with {labels[6]} binding sites", fontsize=14)

    if statistics.mean(fractions[7]) < statistics.mean(mirror_fractions[7]):
        axs[2, 1].hist(fractions[7], color="tab:blue", edgecolor="black")
        axs[2, 1].hist(mirror_fractions[7], color="tab:red", edgecolor="black")
    else:
        axs[2, 1].hist(mirror_fractions[7], color="tab:red", edgecolor="black")
        axs[2, 1].hist(fractions[7], color="tab:blue", edgecolor="black")

    axs[2, 1].set_title(f"Overlap with {labels[7]} binding sites", fontsize=14)

    if statistics.mean(fractions[8]) < statistics.mean(mirror_fractions[8]):
        axs[2, 2].hist(fractions[8], color="tab:blue", edgecolor="black")
        axs[2, 2].hist(mirror_fractions[8], color="tab:red", edgecolor="black")
    else:
        axs[2, 2].hist(mirror_fractions[8], color="tab:red", edgecolor="black")
        axs[2, 2].hist(fractions[8], color="tab:blue", edgecolor="black")

    axs[2, 2].set_title(f"Overlap with {labels[8]} binding sites", fontsize=14)
    
    if statistics.mean(fractions[9]) < statistics.mean(mirror_fractions[9]):
        axs[3, 0].hist(fractions[9], color="tab:blue", edgecolor="black")
        axs[3, 0].hist(mirror_fractions[9], color="tab:red", edgecolor="black")
    else:
        axs[3, 0].hist(mirror_fractions[9], color="tab:red", edgecolor="black")
        axs[3, 0].hist(fractions[9], color="tab:blue", edgecolor="black")
        
    axs[3, 0].set_title(f"Overlap with {labels[9]} binding sites", fontsize=14)
    
    if statistics.mean(fractions[10]) < statistics.mean(mirror_fractions[10]):
        axs[3, 1].hist(fractions[10], color="tab:blue", edgecolor="black")
        axs[3, 1].hist(mirror_fractions[10], color="tab:red", edgecolor="black")
    else:
        axs[3, 1].hist(mirror_fractions[10], color="tab:red", edgecolor="black")
        axs[3, 1].hist(fractions[10], color="tab:blue", edgecolor="black")
        
    axs[3, 1].set_title(f"Overlap with {labels[10]} binding sites", fontsize=14)

    fig.delaxes(axs[3, 2])
    fig.supxlabel("Fraction of edges in c-walk having a peak", fontsize=18)

    axs[0, 1].yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, pos: "{0:g}".format(x / 1000) + "k"))
    axs[1, 1].yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, pos: "{0:g}".format(x / 1000) + "k"))
    axs[2, 1].yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, pos: "{0:g}".format(x / 1000) + "k"))
    axs[3, 1].yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, pos: "{0:g}".format(x / 1000) + "k"))

    return plt.savefig(sys.argv[4]), plt.close()


def counting(boundaries, tf_dict, tf_mirror_dict):
    boundaries = filter(None, boundaries)  # remove empty lists
    fractions, mirror_fractions = [], []
    for cwalk_bounds in boundaries:
        # iter over each cwalk edge bounds (edge bounds are pandas itvs)
        bounds = [i[0] for i in cwalk_bounds]
        chromosome = [i[1] for i in cwalk_bounds][0]
        isin = peak_in_edges(bounds, chromosome, tf_dict)
        mirror_isin = peak_in_edges(bounds, chromosome, tf_mirror_dict)
        # fraction of edges in cwalks which has at least one peak, normalized by cwalk length (number of nodes)
        fraction = round(isin / (len(cwalk_bounds) + 1), 2)
        mirror_fraction = round(mirror_isin / (len(cwalk_bounds) + 1), 2)

        fractions.append(fraction)
        mirror_fractions.append(mirror_fraction)

    return fractions, mirror_fractions


def main():
    graphs, _ = load_files(sys.argv[1], load_cwalk_graph)  # load folder with cwalks

    lst_tf_peaks, labels = load_files((sys.argv[3], parse_tf)
    lst_tf_mirror = []
    for tf_dict in lst_tf_peaks:
        lst_tf_mirror.append(mirror_peaks((sys.argv[2], tf_dict))

    boundaries = []
    for graph in graphs:
        S = [graph.subgraph(c).copy() for c in nx.connected_components(graph)]
        for subgraph in S:  # subgraph is one cwalk
            edge_bounds = count_boundaries(subgraph)  # edge bounds for one cwalk
            boundaries.append(edge_bounds)

    tf_labels = []
    tfs_fractions = []
    tfs_mirror_fractions = []
    for label, tf_dict, tf_mirror_dict in zip(labels, lst_tf_peaks, lst_tf_mirror):
        print(label[:-14])
        tf_labels.append(label[:-14])
        fractions, mirror_fractions = counting(boundaries, tf_dict, tf_mirror_dict)

        print(f"Average value of {label[:-14]} distribution is {statistics.mean(fractions)}")
        print(f"Average value of {label[:-14]} with randomize peaks distribution is {statistics.mean(mirror_fractions)}")

        from scipy.stats import wilcoxon
        if fractions and mirror_fractions:

            print("two-sided")
            res = wilcoxon(x=fractions, y=mirror_fractions, zero_method="zsplit", alternative="two-sided")
            print(f"The value of the test statistic: {res.statistic} and p-value: {res.pvalue}")

            res = wilcoxon(x=fractions, y=mirror_fractions, zero_method="zsplit", alternative="greater")
            print("greater")
            print(f"The value of the test statistic: {res.statistic} and p-value: {res.pvalue}")

            res = wilcoxon(x=fractions, y=mirror_fractions, zero_method="zsplit", alternative="less")
            print("less")
            print(f"The value of the test statistic: {res.statistic} and p-value: {res.pvalue}")

            print("wilcox")
            res = wilcoxon(x=fractions, y=mirror_fractions, zero_method="wilcox", alternative="greater")
            print("greater")
            print(f"The value of the test statistic: {res.statistic} and p-value: {res.pvalue}")

            res = wilcoxon(x=fractions, y=mirror_fractions, zero_method="wilcox", alternative="less")
            print("less")
            print(f"The value of the test statistic: {res.statistic} and p-value: {res.pvalue}")

        zeros = [each for each in fractions if each == 0.0]
        mirror_zeros = [each for each in mirror_fractions if each == 0.0]

        if fractions and mirror_fractions:
            print(
                f"Percentage of c-walks from data which edges has no peaks within is {round(len(zeros) / (len(fractions)) * 100, 2)}%")
            print(
                f"Percentage of c-walks which edges has no random peaks within is {round(len(mirror_zeros) / (len(mirror_fractions)) * 100, 2)}%")

        tfs_fractions.append(fractions)
        tfs_mirror_fractions.append(mirror_fractions)

    histogram(tf_labels, tfs_fractions, tfs_mirror_fractions)


if __name__ == '__main__':
    main()
