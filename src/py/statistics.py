#!/usr/bin/env python

from os import listdir
from os.path import isfile, join
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as st
import seaborn as sns
import matplotlib.ticker as ticker
import sys

""" load all available .tsv files """
infiles = [sample for sample in listdir(sys.argv[1]) if isfile(join(sys.argv[1], sample))]

def count_aligns(infiles):
    """
    Extract information about the number of all aligns,
    aligns with abs_pos >= 500,
    rejected aligns
    and its labels.
    """
    count = []
    count_clean = []
    for file in infiles:
        df = pd.read_csv(f"{sys.argv[1]}/{file}", sep='\t')
        df_clean = df.where(df.abs_pos >= 500).dropna()
        df_clean.reset_index(drop=True, inplace=True)
        count.append(len(df.index))
        count_clean.append(len(df_clean.index))
        rejected = [x1 - x2 for (x1, x2) in zip(count, count_clean)]
        labels = [file[:-4] for file in infiles]
    return count, count_clean, rejected, labels

def barh(data, labels, name):
    sns.set_style("whitegrid")
    fig, ax = plt.subplots(figsize=(12, 10))
    y_pos = np.arange(len(labels))
    ax.barh(y_pos, data, color="royalblue", edgecolor="black", height=1)
    ax.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, pos: "{0:g}".format(x / 1000000) + "M"))
    plt.yticks(y_pos, labels)
    ax.set_title("Number of useful alignments in each sample", fontsize=18)
    return plt.savefig(f"{name}"), plt.close()

def displot(data, name):
    sns.set_style("whitegrid")
    sns.displot(data, kind="kde", legend=False)
    plt.legend(["all alignments", "alignments with condition", "rejected alignments"])
    return plt.savefig(f"{name}"), plt.close()

counted = count_aligns(infiles)
COUNT_all = counted[0]
COUNT_clean = counted[1]
COUNT_rejected = counted[2]
LABELS = counted[3]

""" The quantity of useful sub-samples has a normal distribution """
print(st.shapiro(COUNT_all))       # ShapiroResult(statistic=0.7983113527297974, pvalue=9.660232535679825e-06)
print(st.shapiro(COUNT_clean))     # ShapiroResult(statistic=0.8580167293548584, pvalue=0.00020296304137445986)
print(st.shapiro(COUNT_rejected))  # ShapiroResult(statistic=0.6690419912338257, pvalue=5.432401195548664e-08)

""" Create 95% confidence interval (using normal distribution) for population """
print(st.norm.interval(alpha=0.95, loc=np.mean(COUNT_all), scale=st.sem(COUNT_all)))             #  (1390578.197220105, 2233960.06593779)
print(st.norm.interval(alpha=0.95, loc=np.mean(COUNT_clean), scale=st.sem(COUNT_clean)))         #  (766191.0162746658, 1168563.7205674392)
print(st.norm.interval(alpha=0.95, loc=np.mean(COUNT_rejected), scale=st.sem(COUNT_rejected)))   #  (572586.2517241447, 1117197.2745916448)

""" Plotting """
barh(COUNT_clean, LABELS, sys.argv[2])
displot(counted, sys.argv[3])

