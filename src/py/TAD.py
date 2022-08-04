#!/usr/bin/env python

from analysis import load_cwalk_graph, load_files
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


def boundaries(bedfile: str) -> dict:
    df = pd.read_csv(bedfile, delimiter="\t")
    return df[["chrom", "start", "end", "mode"]]


def identical(itv):
    span = 1
    while not all(itv[0] == element for element in itv):
        if not all(itv[0] == element for element in itv):
            span += 1
            a = itv[0]
            itv = [value for value in itv if value != a]
    return span


def tree(boundaries):
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
    plt.title(f"Intra- and Inter-TAD C-walks in {sys.argv[1]} cells", fontsize=18)
    plt.legend([bar1, bar2, bar3, bar4, bar5], ["1 TAD active", "1 TAD passive", "2 TAD", "3 TAD", "many TADs"],
               loc="upper right",
               prop={"size": 16})
    return plt.savefig(f"{name}"), plt.close()


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
    fig, ax = plt.subplots(figsize=(12, 10))
    ax.bar(labels, data, color="tab:blue", edgecolor="black", width=1)
    ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, pos: "{0:g}".format(x / 1000) + "k"))
    plt.ylabel("Number of c-walks")
    ax.set_title(f"Number of c-walks in each TAD type in {sys.argv[1]} cells", fontsize=18)
    return plt.savefig(f"{name}"), plt.close()


def mirror(chr_sizes: str) -> float:
    df_sizes = pd.read_csv(chr_sizes, sep="\t", header=None)
    df_sizes = df_sizes.loc[df_sizes[0].isin(typical_chromosomes(sys.argv[1]))].reset_index(
        drop=True)  # sys instead of str
    df_sizes.columns = df_sizes.columns.map(str)
    return df_sizes.groupby("0")["1"].agg(list).to_dict()


def comparison(general_in, general_out, name):
    labels = ["data", "random"]
    sns.set_style("whitegrid")
    fig, ax = plt.subplots(figsize=(8, 6))
    bar1 = plt.bar(labels, general_in, label="in TAD", color="tab:blue", edgecolor="black", width=1)
    bar2 = plt.bar(labels, general_out, bottom=general_in, label="out TAD", color="tab:red", edgecolor="black", width=1)
    ax.set_title(f"Fraction of c-walks contained in one TAD and the others outside \n"
                 f"Comparison c-walks from data and random in {sys.argv[1]} cells", fontsize=16)
    plt.xlabel("Number of c-walks")
    plt.ylabel("Percentage")
    plt.legend([bar1, bar2], ["in TAD", "out TAD"], loc="upper right")
    return plt.savefig(f"{name}"), plt.close()


def hop_count(active, passive, two, three, many, name):
    active = active[:6]
    passive = passive[:6]
    two = two[:6]
    three = three[:6]
    many = many[:6]
    sns.set_style("whitegrid")
    fig, ax = plt.subplots(figsize=(12, 10))
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
               prop={"size": 16})
    return plt.savefig(f"{name}"), plt.close()


def random(cwalk, mirror_dict, rest_dict):
    # return random cwalk (mirror image relative to the center of the chromosome)
    center = [(((i[0] + i[1]) / 2), i[2]) for i in cwalk]  # ex. (139285246.5, 'chr3')
    reflection = [(mirror_dict[i[1]][0] - i[0], i[1]) for i in center]
    itv = [(rest_dict[i[1]][i[0]], i[1]) for i in reflection]
    itv = [i for i in itv if len(i[0]) != 0]  # remove empty sets
    reflected = [(list(i[0])[0][0], list(i[0])[0][1], i[1]) for i in itv]
    return reflected


def random_counting(mirror_dict, rest_dict, tree_dict, graphs, length):
    i = 0
    rdm_general_in, rdm_general_out = 0, 0
    active, passive, two, three, many = 0, 0, 0, 0, 0
    for graph in graphs:
        for cwalk in list(nx.connected_components(graph)):
            cwalk = list(cwalk)
            if len(cwalk) == length:
                if all(cwalk[i][2] == cwalk[0][2] for i in range(0, len(cwalk))):  # only intra chrs cwalks
                    reflected = random(cwalk, mirror_dict, rest_dict)
                    if len(reflected) == len(cwalk):
                        i += 1
                        cwalk_boundaries = [((node[0] + node[1]) / 2) for node in reflected]
                        itv_tad = [tree_dict[reflected[0][2]][bound] for bound in cwalk_boundaries]
                        in_tad = identical(itv_tad)
                        modes = [list(i)[0][2] for i in list(itv_tad)]
                        if in_tad == 1:
                            rdm_general_in += 1
                        else:
                            rdm_general_out += 1
                        if in_tad == 1 and modes[0] == "active" and all(modes[0] == element for element in modes):
                            active += 1
                        elif in_tad == 1 and modes[0] == "passive" and all(modes[0] == element for element in modes):
                            passive += 1
                        elif in_tad == 2:
                            two += 1
                        elif in_tad == 3:
                            three += 1
                        elif in_tad > 3:
                            many += 1

    return [round((active / i) * 100, 2), round((passive / i) * 100, 2), round((two / i) * 100, 2),
            round((three / i) * 100, 2), round((many / i) * 100, 2), i, rdm_general_in, rdm_general_out]


def counting(tree_dict, graphs, length):
    i = 0
    general_in, general_out = 0, 0
    active, passive, two, three, many = 0, 0, 0, 0, 0
    for graph in graphs:
        for cwalk in list(nx.connected_components(graph)):
            cwalk = list(cwalk)
            if len(cwalk) == length and all(cwalk[i][2] == cwalk[0][2] for i in range(0, len(cwalk))):
                i += 1
                cwalk_boundaries = [((node[0] + node[1]) / 2) for node in cwalk]
                itv_tad = [tree_dict[cwalk[0][2]][bound] for bound in cwalk_boundaries]
                in_tad = identical(itv_tad)
                modes = [list(i)[0][2] for i in list(itv_tad)]
                if in_tad == 1:
                    general_in += 1
                else:
                    general_out += 1
                if in_tad == 1 and modes[0] == "active" and all(modes[0] == element for element in modes):
                    active += 1
                elif in_tad == 1 and modes[0] == "passive" and all(modes[0] == element for element in modes):
                    passive += 1
                elif in_tad == 2:
                    two += 1
                elif in_tad == 3:
                    three += 1
                elif in_tad > 3:
                    many += 1

    return active, passive, two, three, many, i, general_in, general_out


def main():
    graphs, _ = load_files(sys.argv[2], load_cwalk_graph)  # load .txt cwalk graph
    df_boundaries = boundaries(sys.argv[3])  # dict where keys are chrs and values are TADs active boundaries
    mirror_dict = mirror(sys.argv[4])
    restrictions_dict = parse_bedfile(sys.argv[5], sys.argv[1])  # .bed file
    rest_dict = dict()  # ex. tree_dict["chr1"] will be an object of type IntervalTree
    for chr in typical_chromosomes(sys.argv[1]):  # sys.argv[1]
        """ Interval tree construction, separate for each chromosome """
        restrictions = restrictions_dict[chr][1].tolist()
        intervals = [(i, j) for i, j in zip(restrictions[:-1], restrictions[1:])]
        rest_dict[chr] = IntervalTree.from_tuples(intervals)

    tree_dict = tree(df_boundaries)

    one_tad_active, one_tad_passive, two_tad, three_tad, many_tad = [], [], [], [], []
    count = []
    count_active, count_passive, count_two, count_three, count_many = [], [], [], [], []
    rdm_one_tad_active, rdm_one_tad_passive, rdm_two_tad, rdm_three_tad, rdm_many_tad = [], [], [], [], []
    rdm_count = []
    general_in_lst, general_out_lst = [], []
    rdm_general_in_lst, rdm_general_out_lst = [], []
    for hop in range(3, 16):
        active, passive, two, three, many, i, general_in, general_out = counting(tree_dict, graphs, hop)
        one_tad_active.append(round((active / i) * 100, 2))  # fraction of cwalks inside one TAD of particular length
        one_tad_passive.append(round((passive / i) * 100, 2))
        two_tad.append(round((two / i) * 100, 2))
        three_tad.append(round((three / i) * 100, 2))
        many_tad.append(round((many / i) * 100, 2))
        count.append(i)

        count_active.append(active)
        count_passive.append(passive)
        count_two.append(two)
        count_three.append(three)
        count_many.append(many)

        general_in_lst.append(general_in)
        general_out_lst.append(general_out)

        rdm_active, rdm_passive, rdm_two, rdm_three, rdm_many, j, rdm_general_in, rdm_general_out = \
            random_counting(mirror_dict, rest_dict, tree_dict, graphs, hop)
        rdm_one_tad_active.append(rdm_active)  # fraction of radomm cwalks inside one TAD of particular length
        rdm_one_tad_passive.append(rdm_passive)
        rdm_two_tad.append(rdm_two)
        rdm_three_tad.append(rdm_three)
        rdm_many_tad.append(rdm_many)
        rdm_count.append(j)

        rdm_general_in_lst.append(rdm_general_in)
        rdm_general_out_lst.append(rdm_general_out)

    """ Plotting """
    tad_fractions(one_tad_active, one_tad_passive, two_tad, three_tad, many_tad, f"{sys.argv[6]}")
    tad_fractions(rdm_one_tad_active, rdm_one_tad_passive, rdm_two_tad, rdm_three_tad, rdm_many_tad, f"{sys.argv[7]}")
    tad_count(graphs, tree_dict, sys.argv[8])
    comparison([round((sum(general_in_lst) / sum(count)) * 100, 2),
                round((sum(rdm_general_in_lst) / sum(rdm_count)) * 100, 2)],
               [round((sum(general_out_lst) / sum(count)) * 100, 2),
                round((sum(rdm_general_out_lst) / sum(rdm_count)) * 100, 2)], sys.argv[9])
    hop_count(count_active, count_passive, count_two, count_three, count_many, sys.argv[10])

    print(f"Total number of c-walks: {sum(count)}")
    print(f"Total number of random c-walks: {sum(rdm_count)}")

    print(
        f"Number of c-walks from data inside one TAD: {[sum(general_in_lst)][0]} and others: {[sum(general_out_lst)][0]}")
    print(
        f"Number of random c-walks inside one TAD: {[sum(rdm_general_in_lst)][0]} and others: {[sum(rdm_general_out_lst)][0]}")
    print(f"We are considering sample consists with {sum(count) + sum(rdm_count)} c-walks "
          f"for which we examine two features (c-walks from data and random)."
          f" Significance level: 0.05")  # only intra-chrs

    from statsmodels.stats.contingency_tables import cochrans_q
    print("Cochran condition:")
    table = np.array([[[sum(general_in_lst)][0], [sum(general_out_lst)][0]],
                      [[sum(rdm_general_in_lst)][0], [sum(rdm_general_out_lst)][0]]])
    cochran = cochrans_q(table, return_object=True)
    print(cochran)
    print("Cochran condition is satisfied")

    print("Chi-square test:")
    from scipy.stats import chi2_contingency
    stat, p, dof, expected = chi2_contingency(table)  # Pearson's Chi-square test
    print(f"Chi-square statistics: {stat}")
    print(f"Observed: {table}")
    print(f"Expected: {expected}")
    print(f"Degrees of freedom: {dof}")
    print(f"The relationship is present. P-value: {p}")

    print("Fisher test:")
    from scipy.stats import fisher_exact
    oddsr, p = fisher_exact(table, alternative="two-sided")
    print(f"Two-sided Fisher: {oddsr} and p-value: {p}")

    oddsr, p = fisher_exact(table, alternative="greater")
    print(f"One-sided Fisher: {oddsr} and p-value: {p}")

    
if __name__ == '__main__':
    main()
