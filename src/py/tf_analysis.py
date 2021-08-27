#!/usr/bin/env python

import pandas as pd
import pickle
import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
import random


def parse_tf(tf: str) -> list :
    tf = pd.read_csv(tf, sep='\t', header=None)  # load transcription factor binding sites
    tf = tf.iloc[:, 0:3]
    tf[3] = ((tf[1] + tf[2])/2).astype(int)
    return tf[3].tolist()


def count_cuts(path: set, peaks: list) -> int:
    """ how many times peaks are in one cwalk """
    cut = 0
    for node in path:
        itv = pd.Interval(node[0], node[1], closed="both")
        for peak in peaks:
            if peak in itv:
                cut += 1
    return cut


def normalizing_cuts(graph, peaks):
    cuts = []
    cwalk_length = []
    for cwalk in list(nx.connected_components(graph)):
        cut = count_cuts(cwalk, peaks)
        cuts.append(cut)
        cwalk_length.append(len(cwalk))
    return [i/j for i, j in zip(cuts, cwalk_length)]


P = pickle.load(open("cwalk.txt", "rb"))  # Load cwalks.txt file

tf_peaks = parse_tf("wgEncodeAwgTfbsUtaK562CtcfUniPk.narrowPeak.gz")
random_peaks = random.randint(min(tf_peaks), max(tf_peaks) + 1)  # random peaks

normalize = normalizing_cuts(P, tf_peaks)
random_normalize = normalizing_cuts(P, random_peaks)


""" Plotting """
plt.hist(normalize)
plt.show()
plt.savefig("normalize.png")

plt.hist(random_normalize)
plt.show()
plt.savefig("random_normalize.png")
