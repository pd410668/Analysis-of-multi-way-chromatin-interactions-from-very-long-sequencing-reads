#!/usr/bin/env python

import csv
import matplotlib.pyplot as plt
import seaborn as sns
import sys

"""
statistics.py taking as input txt file with statistics from filtering.py
and return bar plot and line plots
usage:
chmod 777 statistics.py
./statistics.py statistics.txt experiment_name
"""

def main():
    experiment = sys.argv[1]
    name_experiment = sys.argv[2]
    R_list = ["R1 vs R1", "R1 vs R2", "R2 vs R1", "R2 vs R2"]

    def load_txt(txt_file):
        alignments = []
        with open(txt_file) as file:
            infile = csv.reader(file)
            for line in infile:
                if line != ['[]']:
                    align = [line[0][2:-1], int(line[1]), line[2][2:-1], line[3][2:-1], line[4][2:-1], int(line[5][:2])]
                    alignments.append(align)
        return alignments

    def RvsR_barplot(data, experiment):
        plt.figure(figsize=(12, 10))
        sns.set_style("whitegrid")
        plt.bar(R_list, data, width=0.8,
                color="cornflowerblue", edgecolor="navy", linewidth=1)
        plt.title(f"Sample from {name_experiment} human cells", size=20, y=1.04)
        plt.xlabel("Combinations", size=15)
        plt.ylabel("number of matches", size=15)
        return plt.savefig(f"bar_plot_{name_experiment}.png"), plt.close()

    def distances_lineplot(data, experiment, R, linewidth):
        sns.set_style("whitegrid")
        plt.figure(figsize=(16, 12))
        plt.plot(data, color="cornflowerblue", linewidth=linewidth)
        plt.title(f"Sample from {name_experiment} human cells {R}", size=20, y=1.04)
        plt.ylabel("pairwise distance", size=15)
        plt.xlabel("number of comparisons", size=15)
        return plt.savefig(f"line_plot_{name_experiment}_{R}.png"), plt.clf()

    # load alignments from txt file
    alignments = load_txt(experiment)
    alignments = load_txt(experiment)

    # counting matches and collecting distances
    RvsR = [0, 0, 0, 0]
    distances = [[], [], [], []]
    for align in alignments:
        if align[4] == 'R1 vs R1':
            RvsR[0] += align[5]
            distances[0].append(align[1])
        if align[4] == 'R1 vs R2':
            RvsR[1] += align[5]
            distances[1].append(align[1])
        if align[4] == 'R2 vs R1':
            RvsR[2] += align[5]
            distances[2].append(align[1])
        if align[4] == 'R2 vs R2':
            RvsR[3] += align[5]
            distances[3].append(align[1])

    # Line plots for distances for all R vs R
    if distances[0]:
        distances_lineplot(distances[0], name_experiment, R_list[0], 2)
    if distances[1]:
        distances_lineplot(distances[1], name_experiment, R_list[1], 0.6)
    if distances[2]:
        distances_lineplot(distances[2], name_experiment, R_list[2], 2)
    if distances[3]:
        distances_lineplot(distances[3], name_experiment, R_list[3], 2)

    # Bar plot for matches RvsR
    RvsR_barplot(RvsR, name_experiment)

if __name__ == '__main__':
    main()




