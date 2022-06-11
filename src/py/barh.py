#!/usr/bin/env python

from os import listdir
from os.path import isfile, join
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import matplotlib.ticker as ticker
import sys


def load_tsvfiles():
    """ load all available .tsv files """
    return [sample for sample in listdir(sys.argv[1]) if isfile(join(sys.argv[1], sample))]


def count_aligns(infiles):
    """
    Extract information about the number of all aligns,
    aligns with abs_pos >= 500,
    rejected aligns
    and its labels.
    """
    global rejected, labels
    raw_data, data = [], []
    for file in infiles:
        df = pd.read_csv(f"{sys.argv[1]}/{file}", sep='\t')  # load single .tsv file
        df_clean = df.where(df.abs_pos >= 500).dropna()
        df_clean.reset_index(drop=True, inplace=True)
        raw_data.append(len(df.index))
        data.append(len(df_clean.index))
        rejected = [a - b for (a, b) in zip(raw_data, data)]
        labels = [file[:-4] for file in infiles]
    return [data, rejected, labels]


def barh(data, rejected, labels, name):
    sns.set_style("whitegrid")
    fig, ax = plt.subplots(figsize=(12, 16))
    y_pos = np.arange(len(labels))
    bar_1 = ax.barh(y_pos, data, color="tab:blue", edgecolor="black", height=1)
    bar_2 = ax.barh(y_pos, rejected, left=data, color="tab:red", edgecolor="black", height=1)
    ax.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, pos: "{0:g}".format(x / 1000000) + "M"))
    plt.yticks(y_pos, labels)
    ax.set_title("Number of useful and rejected alignments in each sample", fontsize=18)
    plt.legend([bar_1, bar_2], ["Useful", "Rejected"], loc="upper right", prop={"size": 16})
    return plt.savefig(f"{name}"), plt.close()


if __name__ == '__main__':
    counted = count_aligns(load_tsvfiles())
    barh(counted[0], counted[1], counted[2], sys.argv[2])
