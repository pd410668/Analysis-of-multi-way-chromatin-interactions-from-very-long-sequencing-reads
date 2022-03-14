#!/usr/bin/env python

import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
import seaborn as sns
from filtering import typical_chromosomes, collect_data
import sys
import scipy.stats
import statistics
from cwalks_analysis import load_cwalk_graph, load_files


def parse_tf(tf: str) -> dict:
    """
    load transcription factor binding sites
    return: dict where chomosomes are keys and values are peaks within them
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


def counting(path: list, peaks_dict: dict, chromosome: str) -> int:
    """ How many times peaks are in one cwalk from single chromosome """
    cuts = 0
    temp_cut = 0
    for node in path:
        itv = pd.Interval(node[0], node[1], closed="both")
        for peak in peaks_dict[chromosome]:  # element in list of peaks in current chromosome
            if peak in itv:
                temp_cut += 1
        if temp_cut >= 1:
            cuts += 1
    return cuts


def histogram(x, y, name):
    sns.set_style("whitegrid")
    fig, axs = plt.subplots(1, 2, sharex=True, sharey=True, tight_layout=True, figsize=(16, 9))
    fig.suptitle(f"Overlap with {name} binding sites", fontsize=20)
    # axs[0].hist(x, color="tab:blue", edgecolor="black")
    sns.histplot(x, ax=axs[0], color="tab:blue", edgecolor="black")
    axs[0].set_title("Peaks from CHIP-seq data", fontsize=14)
    axs[0].set_xlabel("Fraction of edges in cwalk having a peak", fontsize=14)
    axs[0].set_ylabel("Number of c-walks", fontsize=14)
    # axs[1].hist(y, color="tab:green", edgecolor="black")
    sns.histplot(y, ax=axs[1], color="tab:green", edgecolor="black")
    axs[1].set_title("Randomize peaks", fontsize=14)
    axs[1].set_xlabel("Fraction of edges in cwalk having a peak", fontsize=14)
    return plt.savefig(f"charts/{name}_histplot.png"), plt.close()


def main(label: str):
    graphs, _ = load_files("./cwalks", load_cwalk_graph)
    mirror_peaks_dict = mirror_peaks("hg19.chrom.sizes.tsv", tf_peaks_dict)  # load chromosomes sizes

    normalized_peaks = []
    normalized_random_peaks = []
    for chr in tf_peaks_dict.keys():
        for graph in graphs:
            for cwalk in list(nx.connected_components(graph)):
                cwalk = list(cwalk)
                if all(cwalk[i][2] == cwalk[0][2] for i in range(0, len(cwalk))):  # taking into consideration only intra-chrs cwlks
                    cut_tf = counting(list(cwalk), tf_peaks_dict, chr)
                    cut_mirror = counting(list(cwalk), mirror_peaks_dict, chr)

                    # normalization
                    norm_cut_tf = round(cut_tf / len(list(cwalk)), 2)
                    norm_cut_mirror = round(cut_mirror / len(list(cwalk)), 2)

                    normalized_peaks.append(norm_cut_tf)
                    normalized_random_peaks.append(norm_cut_mirror)

    # cwalks with no cuts are considered as irrelevant
    norm_peaks_limited = [peak for peak in normalized_peaks if peak != 0]
    norm_random_peaks_limited = [peak for peak in normalized_random_peaks if peak != 0]
    histogram(norm_peaks_limited, norm_random_peaks_limited, label)


if __name__ == '__main__':
    Tfs, labels = load_files("./chip-seq_hg19_K562", parse_tf)
    labels = [label[:-14] for label in labels]
    for tf_peaks_dict, label in zip(Tfs, labels):
        main(label)

