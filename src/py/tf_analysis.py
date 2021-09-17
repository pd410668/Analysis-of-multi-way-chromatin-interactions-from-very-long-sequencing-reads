#!/usr/bin/env python

import pickle
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
import seaborn as sns
import random
from filtering import typical_chromosomes
import numpy as np
import sys


def load_cwalk_graph(cwalk_graph):
    return pickle.load(open(cwalk_graph, "rb"))


def parse_tf(tf: str) -> dict:
    """
    load transcription factor binding sites
    return: dict where chomosomes are keys and values are peaks within them
    """
    tf = pd.read_csv(tf, sep='\t', header=None)
    tf = tf.iloc[:, 0:3]
    tf[3] = ((tf[1] + tf[2]) / 2).astype(float)  # changed to float
    tf.columns = tf.columns.map(str)
    tf_dict = tf.groupby("0")["3"].agg(list).to_dict()
    return tf_dict


def mirror_peaks(chr_sizes: str) -> dict:
    """ Mirror image of peaks """
    df_sizes = pd.read_csv(chr_sizes, sep='\t', header=None)
    df_sizes = df_sizes.loc[df_sizes[0].isin(typical_chromosomes())].reset_index(drop=True)
    df_sizes.columns = df_sizes.columns.map(str)
    sizes_dict = df_sizes.groupby("0")["1"].agg(list).to_dict()
    
    keys, mirrors = [], []
    for key in peaks_dict.keys():
        mirror = [sizes_dict[key][0] - peak for peak in peaks_dict[key]]
        keys.append(key)
        mirrors.append(mirror)
    return dict(zip(keys, mirrors))


def counting(path: set) -> int:
    """ how many times peaks are in one cwalk from single chromosome """
    cut = 0
    path = list(path)  # [(21404672, 21405189, 'chr14'), (21403394, 21404672, 'chr14')]
    for node in path:
        if node[2] == chr:
            itv = pd.Interval(node[0], node[1], closed="both")
            for peak in peaks_dict[chr]:  # element in list of peaks in current chromosome
                if peak in itv:
                    cut += 1
    return cut


def tf_histplot(x, y):
    sns.set_style("whitegrid")
    fig, axs = plt.subplots(1, 2, sharex=True, sharey=True, tight_layout=True, figsize=(16, 9))
    axs[0].hist(x, color="tab:blue", edgecolor="black")
    axs[0].set_title("Peaks from data", fontsize=15)
    axs[0].set_xlabel("Number of peaks within one cwalk", fontsize=12)
    axs[0].set_ylabel("Number of c-walks", fontsize=15)
    axs[1].hist(y, color="tab:green", edgecolor="black")
    axs[1].set_title("Randomize peaks", fontsize=15)
    axs[1].set_xlabel("Number of peaks within one cwalk", fontsize=15)
    return plt.savefig("tf_histplot.png")


if __name__ == '__main__':

    P = load_cwalk_graph(sys.argv[1])  # load .txt cwalk graph
    peaks_dict = parse_tf(sys.argv[2])  # load tf binding sites

    normalized_peaks = []
    for chr in typical_chromosomes():
        if chr in peaks_dict.keys():
            for cwalk in list(nx.connected_components(P)):
                cwalk_len = len(cwalk)
                cut = counting(cwalk)
                if cut != 0:  # excluded cwalks with no cuts
                    cut = cut / cwalk_len  # normalization
                    normalized_peaks.append(cut)

    random_peaks = np.random.uniform(low=min(normalized_peaks), high=max(normalized_peaks),
                                     size=(len(normalized_peaks),))  # randomly create peaks for histogram

    tf_histplot(normalized_peaks, random_peaks)  # histograms creation for comparison




