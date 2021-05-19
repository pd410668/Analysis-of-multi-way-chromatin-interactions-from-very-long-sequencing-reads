#!/usr/bin/env python

"""
statistics.py taking as input tsv file as outfile from filtering.py
and return bar plot and hist plot
usage:
chmod 777 statistics.py
./statistics.py statistics_experiment.tsv 
"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import sys

def RvsR_barplot(df, name_experiment):
    sns.set(rc={'figure.figsize': (8, 6)})
    sns.set_style("whitegrid")
    ax = sns.barplot(x=df["RvsR"].value_counts().index, y=df["RvsR"].value_counts(), color="royalblue")
    ax.set(xlabel="Combinations", ylabel="number of matches", title=f"Sample from {name_experiment} human cells")
    return plt.savefig(f"bar_plot_{name_experiment}"), plt.close()

def distances_histplot(df, name_experiment):
    distances = df["position"]
    comparisions = df.reset_index(level=0, inplace=True)
    sns.set(rc={'figure.figsize': (8, 6)})
    sns.set_style("whitegrid")
    ax = sns.histplot(data=df, y=comparisions, x=distances, color="royalblue")
    ax.set(xlabel="Position", ylabel="Number of comparisons", title=f"Sample from {name_experiment} human cells")
    return plt.savefig(f"hist_plot_{name_experiment}"), plt.close()

def main():
    experiment = sys.argv[1]
    name_experiment = experiment[10:-4]
    df = pd.read_csv(f"statistics_{name_experiment}.tsv", sep='\t')
    df.columns = ["seqname", "position", "strand_1", "strand_2", "RvsR"]

    RvsR_barplot(df, name_experiment)
    distances_histplot(df, name_experiment)

if __name__ == '__main__':
    main()
