#!/usr/bin/env python

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.ticker as ticker
import sys

""" load .tsv file """
df = pd.read_csv(sys.argv[1], sep='\t')
experiment_name = sys.argv[1][21:-4]

""" Initiation basic dependencies """
DIST_no_limit = df
DIST_500_inf = df.where(df.abs_pos >= 500)
DIST_0_5000 = df.where(df.abs_pos <= 5000)
DIST_500_10000 = df.where(df.abs_pos >= 500).where(df.abs_pos <= 10000)

df_RvsR_0_5000 = {x: y for x, y in DIST_0_5000.groupby("RvsR")}
df_RvsR_500_10000 = {x: y for x, y in DIST_500_10000.groupby("RvsR")}

df_strand_0_5000 = {x: y for x, y in DIST_0_5000.groupby("strand_1vs2")}
df_strand_500_10000 = {x: y for x, y in DIST_500_10000.groupby("strand_1vs2")}

RvsR_keys = list(df_RvsR_0_5000.keys())
strand_keys = list(df_strand_0_5000.keys())


def main():
    def distances_histplot(df, keys, experiment_name, name):
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
            ax1.yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, pos: "{0:g}".format(x / 1000) + "k"))
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
            ax2.yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, pos: "{0:g}".format(x / 1000) + "k"))
            ax2.set_title("reverse vs reverse", fontsize=18)
        else:
            fig.delaxes(ax2)
        if "R1 vs R2" in keys:
            R1vsR2 = df["R1 vs R2"].reset_index(drop=True)
            sns.histplot(ax=ax3, data=df, x=R1vsR2["abs_pos"], color="red", edgecolor="black")
            ax3.yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, pos: "{0:g}".format(x / 1000) + "k"))
            ax3.set_title("R1 vs R2", fontsize=18)
        elif "forward vs reverse" in keys:
            S1vsS1 = df["forward vs reverse"].reset_index(drop=True)
            sns.histplot(ax=ax3, data=df, x=S1vsS1["abs_pos"], color="forestgreen", edgecolor="black")
            ax3.yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, pos: "{0:g}".format(x / 1000) + "k"))
            ax3.set_title("forward vs reverse", fontsize=18)
        else:
            fig.delaxes(ax3)
        if "R2 vs R1" in keys:
            R2vsR1 = df["R2 vs R1"].reset_index(drop=True)
            sns.histplot(ax=ax4, data=df, x=R2vsR1["abs_pos"], color="red", edgecolor="black")
            ax4.set_title("R2 vs R1", fontsize=18)
        elif "reverse vs forward" in keys:
            S1vsS1 = df["reverse vs forward"].reset_index(drop=True)
            sns.histplot(ax=ax4, data=df, x=S1vsS1["abs_pos"], color="forestgreen", edgecolor="black")
            ax4.yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, pos: "{0:g}".format(x / 1000) + 'K'))
            ax4.set_title("reverse vs forward", fontsize=18)
        else:
            fig.delaxes(ax4)
        fig.suptitle(f"Sample from {experiment_name} human cells", fontsize=20)

        axes = [ax1, ax2, ax3, ax4]
        for ax in axes:
            ax.set_xlabel("Absolute value comparing each two aligns", fontsize=15)
            ax.set_ylabel("Number of compared aligns", fontsize=15)
        return plt.savefig(f"{name}"), plt.close()

    def feature_barplot(df, experiment_name, name, save_as):
        sns.set_style("whitegrid")
        plt.title(f"Sample from {experiment_name} human cells")
        if name == "RvsR":
            fig = plt.figure(figsize=(10, 5))
            fig.suptitle(f"Sample from {experiment_name} human cells")
            plt.subplots_adjust(wspace=0.3, top=0.85)
            fig.add_subplot(121)
            ax1 = sns.barplot(x=df[name].value_counts().index, y=df[name].value_counts(), color="royalblue", edgecolor="black")
            ax1.set(xlabel="Combinations", ylabel="Number of compared aligns")
            ax1.yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, pos: "{0:g}".format(x / 1000000) + "M"))
            fig.add_subplot(122)
            ax2 = sns.barplot(x=df[name].value_counts().iloc[1:].index, y=df[name].value_counts().iloc[1:], color="royalblue", edgecolor="black")
            ax2.set(xlabel="Combinations", ylabel="Number of compared aligns")
        elif name == "strand_1vs2":
            fig = plt.figure(figsize=(8, 5))
            fig.suptitle(f"Sample from {experiment_name} human cells")
            ax = sns.barplot(x=df[name].value_counts().index, y=df[name].value_counts(), color="forestgreen", edgecolor="black")
            ax.set(xlabel="Combinations", ylabel="Number of compared aligns")
            ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, pos: "{0:g}".format(x / 1000) + "k"))
        return plt.savefig(f"{save_as}"), plt.close()

    def log_scale_histplot(df, experiment_name, name):
        sns.set(rc={'figure.figsize': (12, 8)})
        sns.set_style("whitegrid")
        ax = sns.histplot(data=df, x=DIST_500_inf["abs_pos"], log_scale=True, color="red", edgecolor="black")
        ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, pos: "{0:g}".format(x / 1000) + "k"))
        ax.set(xlabel="Absolute value of position comparing each two aligns [log10]",
               ylabel="Number of compared aligns",
               title=f"Sample from {experiment_name} human cells")
        return plt.savefig(f"{name}"), plt.close()

    
    """ Plotting """
    # RvsR
    feature_barplot(DIST_no_limit, experiment_name, "RvsR", sys.argv[2])  # "RvsR"
    distances_histplot(df_RvsR_0_5000, RvsR_keys, experiment_name, sys.argv[3])  # "0_5000_R"
    distances_histplot(df_RvsR_500_10000, RvsR_keys, experiment_name, sys.argv[4])  # "500_10000_R"

    # strand 1vs2
    feature_barplot(DIST_no_limit, experiment_name, "strand_1vs2", sys.argv[5])  # "strand_1vs2"
    distances_histplot(df_strand_0_5000, strand_keys, experiment_name, sys.argv[6])  # "0_5000_S"
    distances_histplot(df_strand_500_10000, strand_keys, experiment_name, sys.argv[7])  # "500_10000_S"

    # log10-scale
    log_scale_histplot(DIST_500_10000, experiment_name, sys.argv[8])  # log10_500_1000"

if __name__ == '__main__':
    main()

