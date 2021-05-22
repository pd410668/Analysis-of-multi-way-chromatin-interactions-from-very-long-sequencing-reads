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

def RvsR_barplot(df, name):
    sns.set(rc={'figure.figsize': (8, 6)})
    sns.set_style("whitegrid")
    ax = sns.barplot(x=df["RvsR"].value_counts().index, y=df["RvsR"].value_counts(), color="royalblue")
    ax.set(xlabel="Combinations", ylabel="number of matches", title=f"Sample from {name} human cells")
    return plt.savefig(f"RvsR_{name}.png"), plt.close()

def distances_histplot(df, name):
    distances = df["position"]
    df.reset_index(level=0, inplace=True)
    comparisions = df["index"]
    sns.set(rc={'figure.figsize': (8, 6)})
    sns.set_style("whitegrid")
    ax = sns.histplot(data=df, y=comparisions, x=distances, color="darkblue")
    ax.set(xlabel="Position", ylabel="Number of comparisons", title=f"Sample from {name} human cells")
    return plt.savefig(f"distances_{name}.png"), plt.close()

def main():
    input = sys.argv[1]
    name = sys.argv[1][37:-4]

    df = pd.read_csv(sys.argv[1], sep='\t')
    df.columns = ["seqname", "position", "strand_1", "strand_2", "RvsR"]

    RvsR_barplot(df, name)
    distances_histplot(df, name)

if __name__ == '__main__':
    main()
