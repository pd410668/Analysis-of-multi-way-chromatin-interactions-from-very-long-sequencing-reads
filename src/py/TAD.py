#!/usr/bin/env python

from cwalk import parse_bedfile
from cwalk_analysis import load_cwalk_graph, load_files
import pandas as pd
import networkx as nx
from intervaltree import Interval, IntervalTree
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.ticker as ticker
import numpy as np
from cwalk import parse_bedfile
from filtering import typical_chromosomes
import sys

import matplotlib.pylab as pylab
params = {'legend.fontsize': 'x-large',
         'axes.labelsize': 'x-large',
         'axes.titlesize':'x-large',
         'xtick.labelsize':'xx-large',
         'ytick.labelsize':'xx-large'}
pylab.rcParams.update(params)


def boundaries(bedfile: str):
    # parse bedfile with TAD domains boundaries
    # return DataFrame with "chrom", "start", "end", "mode", "size" columns
    df = pd.read_csv(bedfile, delimiter="\t")
    return df[["chrom", "start", "end", "mode", "size"]]


def identical(itv: list) -> int:
    # return number of visited intervals from list of intervals
    span = 1
    while not all(itv[0] == element for element in itv):
        if not all(itv[0] == element for element in itv):
            span += 1
            a = itv[0]
            itv = [value for value in itv if value != a]
    return span


def chrs_sizes(chr_sizes: str) -> dict:
    # dict where chrs are keys and values are list with size of this chr
    df_sizes = pd.read_csv(chr_sizes, sep="\t", header=None)
    df_sizes = df_sizes.loc[df_sizes[0].isin(typical_chromosomes(sys.argv[1]))].reset_index(
        drop=True)
    df_sizes.columns = df_sizes.columns.map(str)
    return df_sizes.groupby("0")["1"].agg(list).to_dict()


def random(cwalk, chrs_dict, rest_dict):
    # return random cwalk (mirror image relative to the center of the chromosome)
    center = [(((i[0] + i[1]) / 2), i[2]) for i in cwalk]  # ex. (139285246.5, 'chr3')
    reflection = [(chrs_dict[i[1]][0] - i[0], i[1]) for i in center]
    itv = [(rest_dict[i[1]][i[0]], i[1]) for i in reflection]
    itv = [i for i in itv if len(i[0]) != 0]  # remove empty sets
    reflected = [(list(i[0])[0][0], list(i[0])[0][1], i[1]) for i in itv]
    return reflected


def tad_tree(boundaries):
    # interval tree where keys are chrs and values are zip object (start, end, mode)
    chrs_dict = {x: y for x, y in boundaries.groupby("chrom")}
    keys, values = [], []
    for key in chrs_dict.keys():
        left = chrs_dict[key]["start"].tolist()
        right = chrs_dict[key]["end"].tolist()
        mode = chrs_dict[key]["mode"].tolist()
        sets_list = [(a, b, c) for a, b, c in zip(left, right, mode)]
        keys.append(key)
        values.append(sets_list)
    boundaries_dict = dict(zip(keys, values))
    tree = dict()  # ex. tree_dict["chr1"] will be an object of type IntervalTree
    for key in boundaries_dict.keys():
        intervals = boundaries_dict[key]
        # Interval tree construction, separate for each chromosome
        tree[key] = IntervalTree.from_tuples(intervals)
    return tree


def counting(tad_dict, graphs, length):
    i = 0
    general_in, general_out = 0, 0
    active, passive, two, three, many = 0, 0, 0, 0, 0  # number of cwalks in eah TAD type
    for graph in graphs:
        for cwalk in list(nx.connected_components(graph)):
            cwalk = list(cwalk)
            if len(cwalk) == length and all(
                    cwalk[i][2] == cwalk[0][2] for i in range(0, len(cwalk))):  # only intra-chrs c-walks
                i += 1  # number of intra-chrs c-walks

                # list of the middle of restriction intervals of each node in cwalk
                cwalk_boundaries = [((node[0] + node[1]) / 2) for node in cwalk]

                # found list of TAD bounds (itv) for each node in cwalk
                itv_tad = [tad_dict[cwalk[0][2]][bound] for bound in cwalk_boundaries]

                in_tad = identical(itv_tad)  # number of TAD visited by cwalk
                modes = [list(i)[0][2] for i in list(itv_tad)]  # list of TAD mode for each node in cwalk

                if in_tad == 1:
                    general_in += 1  # number of cwalks in one TAD (active or passive)
                else:
                    general_out += 1  # number of cwalks in more than one TAD
                if in_tad == 1 and modes[0] == "active" and all(modes[0] == element for element in modes):
                    active += 1  # number of cwalks in one active TAD
                elif in_tad == 1 and modes[0] == "passive" and all(modes[0] == element for element in modes):
                    passive += 1  # number of cwalks in one passive TAD
                elif in_tad == 2:
                    two += 1  # number of cwalks in two TADs (active or passive)
                elif in_tad == 3:
                    three += 1  # number of cwalks in three TADs (active or passive)
                elif in_tad > 3:
                    many += 1  # number of cwalks in many TADs (active or passive)

    return active, passive, two, three, many, i, general_in, general_out


def random_counting(tad_dict, graphs, length, chrs_dict, rest_dict):
    i = 0
    general_in, general_out = 0, 0
    active, passive, two, three, many = 0, 0, 0, 0, 0  # number of random cwalks in eah TAD type
    for graph in graphs:
        for cwalk in list(nx.connected_components(graph)):
            cwalk = list(cwalk)

            if len(cwalk) == length and all(
                    cwalk[i][2] == cwalk[0][2] for i in range(0, len(cwalk))):  # only intra-chrs c-walks

                random_cwalk = random(cwalk, chrs_dict, rest_dict)

                if len(cwalk) == len(random_cwalk):
                    i += 1  # number of random intra-chrs c-walks

                    # list of the middle of restriction intervals of each node in cwalk
                    cwalk_boundaries = [((node[0] + node[1]) / 2) for node in random_cwalk]

                    # found list of TAD bounds (itv) for each node in cwalk
                    itv_tad = [tad_dict[random_cwalk[0][2]][bound] for bound in cwalk_boundaries]

                    in_tad = identical(itv_tad)  # number of TAD visited by cwalk
                    modes = [list(i)[0][2] for i in list(itv_tad)]  # list of TAD mode for each node in random cwalk

                    if in_tad == 1:
                        general_in += 1  # number of random cwalks in one TAD (active or passive)
                    else:
                        general_out += 1  # number of random cwalks in more than one TAD
                    if in_tad == 1 and modes[0] == "active" and all(modes[0] == element for element in modes):
                        active += 1  # number of cwalks in one active TAD
                    elif in_tad == 1 and modes[0] == "passive" and all(modes[0] == element for element in modes):
                        passive += 1  # number of cwalks in one passive TAD
                    elif in_tad == 2:
                        two += 1  # number of random cwalks in two TADs (active or passive)
                    elif in_tad == 3:
                        three += 1  # number of random cwalks in three TADs (active or passive)
                    elif in_tad > 3:
                        many += 1  # number of random cwalks in many TADs (active or passive)

    return active, passive, two, three, many, i, general_in, general_out


def tad_count(graphs, tree_dict, name):
    active, passive, two_active, two_passive, three_active, three_passive, other = 0, 0, 0, 0, 0, 0, 0
    labels = ["one active", "one passive", "two active", "two passive", "three active", "three passive", "other"]
    for graph in graphs:
        for cwalk in list(nx.connected_components(graph)):
            cwalk = list(cwalk)
            if all(cwalk[i][2] == cwalk[0][2] for i in range(0, len(cwalk))):  # only intra chrs cwalks
                cwalk_boundaries = [((node[0] + node[1]) / 2) for node in cwalk]
                itv_tad = [tree_dict[cwalk[0][2]][bound] for bound in cwalk_boundaries]
                in_tad = identical(itv_tad)
                modes = [list(i)[0][2] for i in list(itv_tad)]
                if in_tad == 1 and modes[0] == "active" and all(modes[0] == element for element in modes):
                    active += 1
                elif in_tad == 1 and modes[0] == "passive" and all(modes[0] == element for element in modes):
                    passive += 1
                elif in_tad == 2 and modes[0] == "active" and all(modes[0] == element for element in modes):
                    two_active += 1
                elif in_tad == 2 and modes[0] == "passive" and all(modes[0] == element for element in modes):
                    two_passive += 1
                elif in_tad == 3 and modes[0] == "active" and all(modes[0] == element for element in modes):
                    three_active += 1
                elif in_tad == 3 and modes[0] == "passive" and all(modes[0] == element for element in modes):
                    three_passive += 1
                else:
                    other += 1
    data = [active, passive, two_active, two_passive, three_active, three_passive, other]

    sns.set_style("whitegrid")
    fig, ax = plt.subplots(figsize=(15, 12))
    ax.bar(labels, data, color="tab:blue", edgecolor="black", width=1)
    plt.xticks(rotation=30)
    ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, pos: "{0:g}".format(x / 1000) + "k"))
    plt.ylabel("Number of c-walks", fontsize=15)
    ax.set_title(f"Number of c-walks based on location in TAD types in {sys.argv[1]} cells", fontsize=20)
    return plt.savefig(f"{name}"), plt.close()


def chart_comparison(general_in, general_out, name):
    labels = ["data", "random"]
    sns.set_style("whitegrid")
    fig, ax = plt.subplots(figsize=(10, 8))
    bar1 = plt.bar(labels, general_in, label="in a single TAD", color="tab:blue", edgecolor="black", width=1)
    bar2 = plt.bar(labels, general_out, bottom=general_in, label="in more than one TAD", color="tab:red", edgecolor="black", width=1)
    ax.set_title(f"Fraction of c-walks contained in one TAD and c-walks in more than one TAD \n"
                 f"Comparison c-walks from data and random in {sys.argv[1]} cells", fontsize=16)
    plt.ylabel("Percentage")
    plt.legend([bar1, bar2], ["in a single TAD", "in more than one TAD"], loc="upper right")
    return plt.savefig(f"{name}"), plt.close()


def tad_fractions(one_tad_active, one_tad_passive, two_tad, three_tad, many_tad, name):
    sns.set_style("whitegrid")
    plt.figure(figsize=(12, 10))
    bar1 = plt.bar([i for i in range(3, 16)], one_tad_active, width=1, edgecolor="black", color="royalblue")
    bar2 = plt.bar([i for i in range(3, 16)], one_tad_passive, bottom=one_tad_active, width=1, edgecolor="black",
                   color="cornflowerblue")
    bar3 = plt.bar([i for i in range(3, 16)], two_tad,
                   bottom=np.add(one_tad_active, one_tad_passive), width=1, edgecolor="black", color="tab:cyan")
    bar4 = plt.bar([i for i in range(3, 16)], three_tad,
                   bottom=np.add(np.add(one_tad_active, one_tad_passive), two_tad), width=1, edgecolor="black",
                   color="lightsteelblue")
    bar5 = plt.bar([i for i in range(3, 16)], many_tad,
                   bottom=np.add(np.add(np.add(one_tad_active, one_tad_passive), two_tad), three_tad), width=1,
                   edgecolor="black", color="slategrey")
    plt.xlabel("Number of hops")
    plt.ylabel("Percentage")
    plt.title(f"Intra- and Inter-TAD c-walks in {sys.argv[1]} cells", fontsize=18)
    plt.legend([bar1, bar2, bar3, bar4, bar5], ["1 TAD active", "1 TAD passive", "2 TAD", "3 TAD", "many TADs"],
               loc="upper right",
               prop={"size": 16})
    return plt.savefig(f"{name}"), plt.close()


def hop_count(active, passive, two, three, many, name):
    active = active[:6]
    passive = passive[:6]
    two = two[:6]
    three = three[:6]
    many = many[:6]
    sns.set_style("whitegrid")
    fig, ax = plt.subplots(figsize=(10, 8))
    bar1 = plt.bar([i for i in range(3, 9)], active, width=1, edgecolor="black", color="royalblue")
    bar2 = plt.bar([i for i in range(3, 9)], passive, bottom=active, width=1, edgecolor="black",
                   color="cornflowerblue")
    bar3 = plt.bar([i for i in range(3, 9)], two,
                   bottom=np.add(active, passive), width=1, edgecolor="black", color="tab:cyan")
    bar4 = plt.bar([i for i in range(3, 9)], three,
                   bottom=np.add(np.add(active, passive), two), width=1, edgecolor="black",
                   color="lightsteelblue")
    bar5 = plt.bar([i for i in range(3, 9)], many,
                   bottom=np.add(np.add(np.add(active, passive), two), three), width=1,
                   edgecolor="black", color="slategrey")
    ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, pos: "{0:g}".format(x / 1000) + "k"))
    plt.xlabel("Number of hops")
    plt.ylabel("Number of c-walks")
    plt.title(f"Intra- and inter-TAD c-walks in {sys.argv[1]} cells", fontsize=18)
    plt.legend([bar1, bar2, bar3, bar4, bar5], ["1 TAD active", "1 TAD passive", "2 TAD", "3 TAD", "many TADs"],
               loc="upper right",
               prop={"size": 12})
    return plt.savefig(f"{name}"), plt.close()


def domains(boundaries, name):
    active_boundaries = boundaries.loc[boundaries["mode"] == "active"].dropna().reset_index(drop=True)
    passive_boundaries = boundaries.loc[boundaries["mode"] == "passive"].dropna().reset_index(drop=True)
    sns.set_style("whitegrid")
    plt.figure(figsize=(10, 8))
    sns.distplot(active_boundaries["size"], hist=False, color="red")
    sns.distplot(passive_boundaries["size"], hist=False, color="black")
    plt.title(f"{sys.argv[1]} domains size distribution")
    plt.legend(labels=["Active", "Inactive"])
    plt.xlabel("Domain size [Mb]")
    return plt.savefig(f"{name}"), plt.close()


def main():
    graphs, _ = load_files(sys.argv[2], load_cwalk_graph)  # load .txt cwalk graph
    tad_dict = tad_tree(boundaries(sys.argv[3]))  # interval tree dict with TAD domains
    chrs_dict = chrs_sizes(sys.argv[4])

    # plot domains size distribution
    domains(boundaries(sys.argv[3]), sys.argv[11])

    # Return dict where chrs are keys and values are list of restriction itvs
    restrictions_dict = parse_bedfile(sys.argv[5], sys.argv[1])

    # Interval tree construction, separate for each chromosome
    rest_dict = dict()  # ex. tree_dict["chr1"] will be an object of type IntervalTree
    for chr in typical_chromosomes(sys.argv[1]):
        restrictions = restrictions_dict[chr][1].tolist()
        restrictions[0] = 0
        restrictions[-1] = chrs
_dict[chr][0]
        intervals = [(i, j) for i, j in zip(restrictions[:-1], restrictions[1:])]
        rest_dict[chr] = IntervalTree.from_tuples(intervals)

    lst_active, lst_passive, lst_two, lst_three, lst_many, lst_i, lst_general_in, lst_general_out = \
        [], [], [], [], [], [], [], []
    lst_rdm_active, lst_rdm_passive, lst_rdm_two, lst_rdm_three, lst_rdm_many, lst_rdm_i, lst_rdm_general_in, lst_rdm_general_out = \
        [], [], [], [], [], [], [], []

    count_active, count_passive, count_two, count_three, count_many, count_general_in, count_general_out = \
        [], [], [], [], [], [], []
    rdm_count_active, rdm_count_passive, rdm_count_two, rdm_count_three, rdm_count_many, rdm_count_general_in, rdm_count_general_out = \
        [], [], [], [], [], [], []

    for length in range(3, 30):
        active, passive, two, three, many, i, general_in, general_out = counting(tad_dict, graphs, length)
        # list of % of intra-chrs cwalks of each length
        if i != 0:
            count_general_in.append(general_in)
            count_general_out.append(general_out)

            lst_general_in.append(round((general_in / i) * 100, 2))  # % of intra-chrs cwalks in a single TAD
            lst_general_out.append(round((general_out / i) * 100, 2))  # % of intra-chrs cwalks in more than one TAD
            lst_i.append(i)  # number of cwalks of each length

            rdm_active, rdm_passive, rdm_two, rdm_three, rdm_many, rdm_i, rdm_general_in, rdm_general_out = \
                random_counting(tad_dict, graphs, length, chrs_dict, rest_dict)

            rdm_count_general_in.append(rdm_general_in)
            rdm_count_general_out.append(rdm_general_out)

            lst_rdm_general_in.append(round((rdm_general_in / i) * 100, 2))  # % of intra-chrs cwalks in a single TAD
            lst_rdm_general_out.append(round((rdm_general_out / i) * 100, 2))  # % of intra-chrs cwalks in more than one TAD
            lst_rdm_i.append(rdm_i)  # number of random cwalks of each length

    """ Plotting """
    chart_comparison([round((sum(count_general_in) / sum(lst_i)) * 100, 2),
                      round((sum(rdm_count_general_in) / sum(lst_rdm_i)) * 100, 2)],
                     [round((sum(count_general_out) / sum(lst_i)) * 100, 2),
                      round((sum(rdm_count_general_out) / sum(lst_rdm_i)) * 100, 2)], sys.argv[9])

    print(f"Total number of c-walks: {sum(lst_i)}")
    print(f"Total number of random c-walks: {sum(lst_rdm_i)}")

    print(f"Percentage of intra-chrs cwalks in a single TAD: {round((sum(count_general_in) / sum(lst_i)) * 100, 2)}%")
    print(f"Percentage of intra-chrs cwalks in more than one TAD: {round((sum(count_general_out) / sum(lst_i)) * 100, 2)}%")

    print(f"Percentage of random intra-chrs cwalks in a single TAD: {round((sum(rdm_count_general_in) / sum(lst_rdm_i)) * 100, 2)}%")
    print(f"Percentage of random intra-chrs cwalks in more than one TAD: {round((sum(rdm_count_general_out) / sum(lst_rdm_i)) * 100, 2)}%")

    print(
        f"Number of c-walks from data inside one TAD: {[sum(count_general_in)][0]} and in more than one TAD:"
        f" {[sum(count_general_out)][0]}")
    print(
        f"Number of random c-walks inside one TAD: {[sum(rdm_count_general_in)][0]} and in more than one TAD: "
        f"{[sum(rdm_count_general_out)][0]}")
    print(f"We are considering sample consists with {sum(lst_i) + sum(lst_rdm_i)} c-walks "
          f"for which we examine two features (c-walks from data and random)."
          f" Significance level: 0.05")

    table = np.array([[[sum(count_general_in)][0], [sum(count_general_out)][0]],
                      [[sum(rdm_count_general_in)][0], [sum(rdm_count_general_out)][0]]])

    print("Fisher test:")
    from scipy.stats import fisher_exact

    oddsr, p = fisher_exact(table, alternative="two-sided")
    print(f"Two-sided Fisher: {oddsr} and p-value: {p}")

    oddsr, p = fisher_exact(table, alternative="greater")
    print(f"One-sided Fisher: {oddsr} and p-value: {p}")

    for length in range(3, 16):
        active, passive, two, three, many, i, general_in, general_out = counting(tad_dict, graphs, length)
        # list of % of intra-chrs cwalks of each length
        if i != 0:
            lst_active.append(round((active / i) * 100, 2))  # % of intra-chrs cwalks in a one active TAD
            lst_passive.append(round((passive / i) * 100, 2))  # % of intra-chrs cwalks in a one passive TAD
            lst_two.append(round((two / i) * 100, 2))  # % of intra-chrs cwalks in two TADs
            lst_three.append(round((three / i) * 100, 2))  # % of intra-chrs cwalks three TADs
            lst_many.append(round((many / i) * 100, 2))  # % of intra-chrs cwalks in many TADs

            # list of intra-chrs cwalks counts of each length
            count_active.append(active)
            count_passive.append(passive)
            count_two.append(two)
            count_three.append(three)
            count_many.append(many)

            rdm_active, rdm_passive, rdm_two, rdm_three, rdm_many, rdm_i, rdm_general_in, rdm_general_out = \
                random_counting(tad_dict, graphs, length, chrs_dict, rest_dict)
            # list of % of intra-chrs cwalks of each length
            lst_rdm_active.append(round((rdm_active / rdm_i) * 100, 2))  # % of intra-chrs cwalks in a one active TAD
            lst_rdm_passive.append(round((rdm_passive / rdm_i) * 100, 2))  # % of intra-chrs cwalks in a one passive TAD
            lst_rdm_two.append(round((rdm_two / rdm_i) * 100, 2))  # % of intra-chrs cwalks in two TADs
            lst_rdm_three.append(round((rdm_three / rdm_i) * 100, 2))  # % of intra-chrs cwalks three TADs
            lst_rdm_many.append(round((rdm_many / rdm_i) * 100, 2))  # % of intra-chrs cwalks in many TADs

            # list of intra-chrs random cwalks counts of each length
            rdm_count_active.append(active)
            rdm_count_passive.append(passive)
            rdm_count_two.append(two)
            rdm_count_three.append(three)
            rdm_count_many.append(many)

    """ Plotting """
    hop_count(count_active, count_passive, count_two, count_three, count_many, sys.argv[10])
    tad_count(graphs, tad_dict, sys.argv[8])
    tad_fractions(lst_active, lst_passive, lst_two, lst_three, lst_many, sys.argv[6])
    tad_fractions(lst_rdm_active, lst_rdm_passive, lst_rdm_two, lst_rdm_three, lst_rdm_many,
                  sys.argv[7])


if __name__ == '__main__':
    main()
