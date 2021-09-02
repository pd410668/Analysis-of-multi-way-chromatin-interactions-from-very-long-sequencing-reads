#!/usr/bin/env python

import pickle
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
import seaborn as sns
import random
from filtering import typical_chromosomes


def load_cwalk_graph(cwalk_graph):
    return pickle.load(open(cwalk_graph, "rb"))


def parse_tf(tf: str) -> list:
    """
    load transcription factor binding sites
    return: list of chromosomes and peaks peaks within them
    """
    tf = pd.read_csv(tf, sep='\t', header=None)
    tf = tf.iloc[:, 0:3]
    tf[3] = ((tf[1] + tf[2]) / 2).astype(int)
    return tf[0].tolist(), tf[3].tolist()


def count_peaks(path: set, peaks: list) -> int:
    """ how many times peaks are in one cwalk """
    cut = 0
    path = list(path)
    for node in path:
        itv = pd.Interval(node[0], node[1], closed="both")
        for peak in peaks:
            if peak in itv:
                cut += 1
    return cut


def normalizing(graph, peaks):
    cuts = []
    cwalk_length = []
    for cwalk in list(nx.connected_components(graph)):
        cut = count_peaks(cwalk, peaks)
        cuts.append(cut)
        cwalk_length.append(len(cwalk))
    return [i / j for i, j in zip(cuts, cwalk_length)]


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


def main():
    random_peaks = random.randint(min(tf_peaks), max(tf_peaks) + 1)  # random peaks
    normalize = normalizing(P, tf_peaks)
    random_normalize = normalizing(P, random_peaks)
    tf_histplot(normalize, random_normalize)


if __name__ == '__main__':
    import sys
    P = load_cwalk_graph(sys.argv[1])  # load .txt cwalk graph
    chrs, tf_peaks = parse_tf(sys.argv[2])  # load CTCF binding sites
    main()
