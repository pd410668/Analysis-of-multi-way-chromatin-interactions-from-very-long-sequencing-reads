#!/usr/bin/env python

import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
import seaborn as sns
from filtering import typical_chromosomes, collect_data
import sys
import scipy.stats
import statistics


def load_cwalk_graph(cwalk_graph):
    import pickle
    return pickle.load(open(cwalk_graph, "rb"))


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


def load_tfs(dictionary):
    """
    load all available .narrowPeak.gz files
    and return list of dicts
    """
    import os
    tfs = []
    for file in os.listdir(dictionary):  # "./chip-seq_hg19_K562"
        tf = f"{dictionary}/{file}"
        entries = parse_tf(tf)
        tfs.append(entries)
    return tfs


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


def counting(path: list, peaks_dict: dict) -> int:
    """ How many times peaks are in one cwalk from single chromosome """
    cut_tf = 0
    for node in path:
        itv = pd.Interval(node[0], node[1], closed="both")
        for peak in peaks_dict[chr]:  # element in list of peaks in current chromosome
            if peak in itv:
                cut_tf += 1
    return cut_tf


def tf_histplot(x, y, name):
    sns.set_style("whitegrid")
    fig, axs = plt.subplots(1, 2, sharex=True, sharey=True, tight_layout=True, figsize=(16, 9))
    fig.suptitle("Max-TF", fontsize=20)
    # axs[0].hist(x, color="tab:blue", edgecolor="black")
    sns.histplot(x, ax=axs[0], kde=True, color="tab:blue", edgecolor="black")
    axs[0].set_title("Peaks from data", fontsize=14)
    axs[0].set_xlabel("Number of peaks within one cwalk", fontsize=14)
    axs[0].set_ylabel("Number of c-walks", fontsize=14)
    # axs[1].hist(y, color="tab:green", edgecolor="black")
    sns.histplot(y, ax=axs[1], kde=True, color="tab:green", edgecolor="black")
    axs[1].set_title("Randomize peaks", fontsize=14)
    axs[1].set_xlabel("Number of peaks within one cwalk", fontsize=14)
    return plt.savefig(f"{name}_tf_histplot.png")


if __name__ == '__main__':

    P = load_cwalk_graph("hs_k562_I_1_cwalks.txt")  # load .txt cwalk graph
    tf_peaks_dict = parse_tf("chip-seq_hg19_K562/wgEncodeAwgTfbsHaibK562MaxV0416102UniPk.narrowPeak.gz")  # load tf binding sites
    mirror_peaks_dict = mirror_peaks("hg19.chrom.sizes.tsv", tf_peaks_dict)  # load chromosomes sizes

    """ create supportive .tsv file with field names """
    collect_data(["cwalk", "cwalk_length", "cuts_number"], "cwalks_binding_sites.tsv", "w")

    normalized_peaks = []
    normalized_random_peaks = []
    for chr in tf_peaks_dict.keys():
        for cwalk in list(nx.connected_components(P)):
            cwalk = list(cwalk)
            cwalk_len = len(cwalk)

            cut_tf = counting(cwalk, tf_peaks_dict)
            cut_mirror = counting(cwalk, mirror_peaks_dict)

            # save supportive .tsv file
            collect_data([cwalk, cwalk_len, cut_tf], "cwalks_binding_sites/cwalks_binding_sites.tsv", "a")  # sys instead of str

            # normalization
            norm_cut_tf = round(cut_tf / cwalk_len, 2)
            norm_cut_mirror = round(cut_mirror / cwalk_len, 2)

            normalized_peaks.append(norm_cut_tf)
            normalized_random_peaks.append(norm_cut_mirror)

    tf_histplot(normalized_peaks, normalized_random_peaks, "charts/standard_Max")  # non readable

    """ Cwalks where cuts are found insignificant number of times (average is between (0, 1] """
    # cwalks with no cuts are considered as irrelevant
    norm_peaks_limited = [peak for peak in normalized_peaks if 0 < peak <= 1]
    norm_random_peaks_limited = [peak for peak in normalized_random_peaks if 0 < peak <= 1]

    tf_histplot(norm_peaks_limited, norm_random_peaks_limited, "charts/limited_Max")

    """ Cwalks where cuts are found significant more times (average is greater than 1) """
    norm_peaks_unlimited = [peak for peak in normalized_peaks if peak > 1]
    norm_random_peaks_unlimited = [peak for peak in normalized_random_peaks if peak > 1]

    tf_histplot(norm_peaks_unlimited, norm_random_peaks_unlimited, "charts/unlimited_Max")

    # Statistical tests

    # The Shapiro-Wilk test tests the null hypothesis that the data was drawn from a normal distribution
    sw_norm_peaks_limited = scipy.stats.shapiro(norm_peaks_limited)
    sw_norm_random_peaks_limited = scipy.stats.shapiro(norm_random_peaks_limited)
    print(sw_norm_peaks_limited, sw_norm_random_peaks_limited)

    sw_norm_peaks_unlimited = scipy.stats.shapiro(norm_peaks_unlimited)
    sw_norm_random_peaks_unlimited = scipy.stats.shapiro(norm_random_peaks_unlimited)
    print(sw_norm_peaks_unlimited, sw_norm_random_peaks_unlimited)

    """ Perform the Mann-Whitney U rank test on two independent samples """
    _, mw_pv_limited = scipy.stats.mannwhitneyu(norm_peaks_limited, norm_random_peaks_limited,
                                                alternative="greater")
    _, mw_pv_unlimited = scipy.stats.mannwhitneyu(norm_peaks_unlimited, norm_random_peaks_unlimited,
                                                  alternative="greater")  # can be performed t-test instead
    print(mw_pv_limited) 
    print(mw_pv_unlimited)  

    print(statistics.mean(norm_peaks_unlimited))
    print(statistics.mean(norm_random_peaks_unlimited))

    print(statistics.variance(norm_peaks_unlimited))
    print(statistics.variance(norm_random_peaks_unlimited))

    print(len(norm_peaks_unlimited))
    print(len(norm_random_peaks_unlimited))

    t_test = scipy.stats.ttest_ind(norm_peaks_unlimited, norm_random_peaks_unlimited, equal_var=False)
    print(t_test)

    _, all_pw = scipy.stats.mannwhitneyu(normalized_peaks, normalized_random_peaks,
                                                  alternative="two-sided")
    print(all_pw)