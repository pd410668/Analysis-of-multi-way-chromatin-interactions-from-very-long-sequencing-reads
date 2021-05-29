#!/usr/bin/env python

"""
statistics.py taking as input tsv file as outfile from filtering.py
and return bar plot and hist plot
usage:
chmod 777 statistics.py
./statistics.py statistics_experiment.tsv histplots_name barplot_name
"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.ticker as ticker
import sys


def main():
    def distancess_histplot(df, experiment_name):
        sns.set_style("whitegrid")
        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex=False, sharey=False, figsize=(16, 12))
        plt.subplots_adjust(bottom=0.1, top=0.9, hspace=0.25, wspace=0.3, right=0.93, left=0.15)

        if "R1 vs R1" in RvsR_keys:
            R1vsR1 = df["R1 vs R1"].reset_index(drop=True)
            sns.histplot(ax=ax1, data=df, x=R1vsR1["abs_pos"], color="red", edgecolor="black")
        else:
            fig.delaxes(ax1)

        if "R2 vs R2" in RvsR_keys:
            R2vsR2 = df["R2 vs R2"].reset_index(drop=True)
            sns.histplot(ax=ax2, data=df, x=R2vsR2["abs_pos"], color="red", edgecolor="black")
        else:
            fig.delaxes(ax2)

        if "R1 vs R2" in RvsR_keys:
            R1vsR2 = df["R1 vs R2"].reset_index(drop=True)
            sns.histplot(ax=ax3, data=df, x=R1vsR2["abs_pos"], color="red", edgecolor="black")
            ax3.yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, pos: "{:,.0f}".format(x/1000) + 'K'))
        else:
            fig.delaxes(ax3)

        if "R2 vs R1" in RvsR_keys:
            R2vsR1 = df["R2 vs R1"].reset_index(drop=True)
            sns.histplot(ax=ax4, data=df, x=R2vsR1["abs_pos"], color="red", edgecolor="black")
        else:
            fig.delaxes(ax4)

        ax1.set_title("R1 vs R1", fontsize=18)
        ax2.set_title("R2 vs R2", fontsize=18)
        ax3.set_title("R1 vs R2", fontsize=18)
        ax4.set_title("R2 vs R1", fontsize=18)

        fig.suptitle(f"Sample from {experiment_name} human cells", fontsize=20)

        axes = [ax1, ax2, ax3]
        for ax in axes:
            ax.set_xlabel("Absolute value comparing each two aligns", fontsize=15)
            ax.set_ylabel("Number of compared aligns", fontsize=15)
        return plt.savefig(f"{sys.argv[2]}"), plt.close() 

    def RvsR_barplot(df, experiment_name):
        sns.set(rc={'figure.figsize': (8, 6)})
        sns.set_style("whitegrid")
        ax = sns.barplot(x=df["RvsR"].value_counts().index, y=df["RvsR"].value_counts(), color="royalblue", edgecolor="black")
        ax.set(xlabel="Combinations", ylabel="number of matches", title=f"Sample from {experiment_name} human cells")
        ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, pos: "{:,.0f}".format(x / 1000) + 'K'))
        return plt.savefig(f"{sys.argv[3]}"), plt.close() 


    df = pd.read_csv(sys.argv[1], sep='\t')
    experiment_name = sys.argv[1][25:-4]

    # Initiation basic dependencies
    df_RvsR = {x: y for x, y in df.groupby("RvsR")}
    RvsR_keys = list(df_RvsR.keys())

    # Plotting
    RvsR_barplot(df, experiment_name)
    distancess_histplot(df_RvsR, experiment_name)

if __name__ == '__main__':
    main()
