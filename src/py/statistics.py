#!/usr/bin/env python

"""
statistics.py taking as input tsv file as outfile from filtering.py
and return bar plot and hist plot
usage:
chmod 777 statistics.py
./statistics.py statistics_experiment.tsv histplots_name barplot_name log10_name
"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.ticker as ticker
import sys

def main():
    def distances_histplot(df, df2, keys, experiment_name, name):
        sns.set_style("whitegrid")
        fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex=False, sharey=False, figsize=(16, 12))
        plt.subplots_adjust(bottom=0.1, top=0.9, hspace=0.25, wspace=0.3, right=0.93, left=0.15)
        if "R1 vs R1" in keys:
            R1vsR1 = df["R1 vs R1"].reset_index(drop=True)
            sns.histplot(ax=ax1, data=df, x=R1vsR1["abs_pos"], color="red", edgecolor="black")
            ax1.set_title("R1 vs R1", fontsize=18)
        elif "forward vs forward" in keys:
            S1vsS1 = df["forward vs forward"].reset_index(drop=True)
            sns.histplot(ax=ax1, data=df, x=S1vsS1["abs_pos"], color="forestgreen", edgecolor="black")
            ax1.yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, pos: "{0:g}".format(x / 1000) + 'K'))
            ax1.set_title("forward vs forward", fontsize=18)
        else:
            fig.delaxes(ax1)
        if "R2 vs R2" in keys:
            R2vsR2 = df["R2 vs R2"].reset_index(drop=True)
            sns.histplot(ax=ax2, data=df, x=R2vsR2["abs_pos"], color="red", edgecolor="black")
            ax2.set_title("R2 vs R2", fontsize=18)
        elif "reverse vs reverse" in keys:
            S1vsS1 = df["reverse vs reverse"].reset_index(drop=True)
            sns.histplot(ax=ax2, data=df, x=S1vsS1["abs_pos"], color="forestgreen", edgecolor="black")
            ax2.yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, pos: "{0:g}".format(x / 1000) + 'K'))
            ax2.set_title("reverse vs reverse", fontsize=18)
        else:
            fig.delaxes(ax2)
        if "R1 vs R2" in keys:
            R1vsR2 = df2["R1 vs R2"].reset_index(drop=True)
            sns.histplot(ax=ax3, data=df2, x=R1vsR2["abs_pos"], color="red", edgecolor="black")
            ax3.yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, pos: "{0:g}".format(x / 1000) + 'K'))
            ax3.set_title("R1 vs R2", fontsize=18)
        elif "forward vs reverse" in keys:
            S1vsS1 = df2["forward vs reverse"].reset_index(drop=True)
            sns.histplot(ax=ax3, data=df2, x=S1vsS1["abs_pos"], color="forestgreen", edgecolor="black")
            ax3.yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, pos: "{0:g}".format(x / 1000) + 'K'))
            ax3.set_title("forward vs reverse", fontsize=18)
        else:
            fig.delaxes(ax3)
        if "R2 vs R1" in keys:
            R2vsR1 = df2["R2 vs R1"].reset_index(drop=True)
            sns.histplot(ax=ax4, data=df2, x=R2vsR1["abs_pos"], color="red", edgecolor="black")
            ax4.set_title("R2 vs R1", fontsize=18)
        elif "reverse vs forward" in keys:
            S1vsS1 = df2["reverse vs forward"].reset_index(drop=True)
            sns.histplot(ax=ax4, data=df2, x=S1vsS1["abs_pos"], color="forestgreen", edgecolor="black")
            ax4.yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, pos: "{0:g}".format(x / 1000) + 'K'))
            ax4.set_title("reverse vs forward", fontsize=18)
        else:
            fig.delaxes(ax4)
        fig.suptitle(f"Sample from {experiment_name} human cells", fontsize=20)

        axes = [ax1, ax2, ax3, ax4]
        for ax in axes:
            ax.set_xlabel("Absolute value comparing each two aligns", fontsize=15)
            ax.set_ylabel("Number of compared aligns", fontsize=15)
        return plt.savefig(f"{sys.argv[2]}_{name}"), plt.close()

    def feature_barplot(df, experiment_name, name):
        sns.set_style("whitegrid")
        plt.title(f"Sample from {experiment_name} human cells")
        if name == "RvsR":
            fig = plt.figure(figsize=(10, 5))
            fig.suptitle(f"Sample from {experiment_name} human cells")
            plt.subplots_adjust(wspace=0.3, top=0.85)
            fig.add_subplot(121)
            ax1 = sns.barplot(x=df[name].value_counts().index, y=df[name].value_counts(), color="royalblue", edgecolor="black")
            ax1.set(xlabel="Combinations", ylabel="Number of compared aligns [mln]")
            ax1.yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, pos: "{0:g}".format(x / 1000000) + 'M'))
            fig.add_subplot(122)
            ax2 = sns.barplot(x=df[name].value_counts().iloc[1:].index, y=df[name].value_counts().iloc[1:], color="royalblue", edgecolor="black")
            ax2.set(xlabel="Combinations", ylabel="Number of compared aligns")
        elif name == "strand_1vs2":
            fig = plt.figure(figsize=(8, 5))
            fig.suptitle(f"Sample from {experiment_name} human cells")
            ax = sns.barplot(x=df[name].value_counts().index, y=df[name].value_counts(), color="forestgreen", edgecolor="black")
            ax.set(xlabel="Combinations", ylabel="Number of compared aligns [mln]")
            ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, pos: "{0:g}".format(x / 1000) + 'K'))
        return plt.savefig(f"{sys.argv[3]}_{name}"), plt.close()

    def log_scale_histplot(df, experiment_name, name):
        sns.set(rc={'figure.figsize': (8, 6)})
        sns.set_style("whitegrid")
        ax = sns.histplot(data=df, x=DIST_500_inf["abs_pos"], log_scale=True, color="red", edgecolor="black")
        ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, pos: "{0:g}".format(x / 1000) + 'K'))
        ax.set(xlabel="Absolute value of position comparing each two aligns [log10]",
               ylabel="Number of compared aligns",
               title=f"Sample from {experiment_name} human cells")
        return plt.savefig(f"{sys.argv[4]}_{name}"), plt.close()

    # load .tsv file
    df = pd.read_csv(sys.argv[1], sep='\t')
    experiment_name = sys.argv[1][26:-4]

    """ Initiation basic dependencies """

    DIST_no_limit = df
    DIST_0_1500 = df.where(df.abs_pos <= 1500)
    DIST_600_1600 = df.where(df.abs_pos >= 600).where(df.abs_pos <= 1600)
    DIST_500_inf = df.where(df.abs_pos >= 500)
    DIST_0_10000 = df.where(df.abs_pos <= 10000)

    df_RvsR_0_1500 = {x: y for x, y in DIST_0_1500.groupby("RvsR")}
    df_RvsR_600_1600 = {x: y for x, y in DIST_600_1600.groupby("RvsR")}
    df_RvsR_0_10000 = {x: y for x, y in DIST_0_10000.groupby("RvsR")}

    df_strand_0_1500 = {x: y for x, y in DIST_0_1500.groupby("strand_1vs2")}
    df_strand_600_1600 = {x: y for x, y in DIST_600_1600.groupby("strand_1vs2")}
    df_strand_0_10000 = {x: y for x, y in DIST_0_10000.groupby("strand_1vs2")}

    RvsR_keys = list(df_RvsR_0_1500.keys())
    strand_keys = list(df_strand_0_1500.keys())

    """ Plotting """

    # RvsR
    feature_barplot(DIST_no_limit, experiment_name, "RvsR")
    distances_histplot(df_RvsR_0_1500, df_RvsR_0_1500, RvsR_keys, experiment_name, "0_1500_R")
    distances_histplot(df_RvsR_600_1600, df_RvsR_600_1600, RvsR_keys, experiment_name, "600_1600_R")
    distances_histplot(df_RvsR_0_10000, df_RvsR_600_1600, RvsR_keys, experiment_name, "0_10000_R")

    # strand 1vs2
    feature_barplot(DIST_no_limit, experiment_name, "strand_1vs2")
    distances_histplot(df_strand_0_1500, df_strand_0_1500, strand_keys, experiment_name, "0_1500_S")
    distances_histplot(df_strand_600_1600, df_strand_600_1600, strand_keys, experiment_name, "600_1600_S")
    distances_histplot(df_strand_0_10000, df_strand_600_1600, strand_keys, experiment_name, "0_10000_S")

    # log10-scale
    log_scale_histplot(DIST_no_limit, experiment_name, "log10")

if __name__ == '__main__':
    main()
