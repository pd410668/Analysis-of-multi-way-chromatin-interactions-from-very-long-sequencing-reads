#!/usr/bin/env python

import functions
import csv

experiment = "hs_k562_I_1_R1"

def read_RvsR_statistics(in_txt_file):
    statistics = []
    with open(f"RvsR_{in_txt_file}.txt") as file:
        read = csv.reader(file)
        for stat in read:
            count_list = [int(stat[0][1:]), int(stat[1]), int(stat[2]), int(stat[3][:2])]
            statistics.append(count_list)
    summary_statistics = [0, 0, 0, 0]
    for stat in statistics:
        summary_statistics[0] += stat[0]
        summary_statistics[1] += stat[1]
        summary_statistics[2] += stat[2]
        summary_statistics[3] += stat[3]
    return summary_statistics

def read_distances_statistics(in_txt_file):
    statistics = []
    with open(f"distances_{in_txt_file}.txt") as file:
        read = csv.reader(file)
        for stat in read:
            count_list = [int(stat[0][1:]), int(stat[1]), int(stat[2]), int(stat[3][:2])]
            statistics.append(count_list)
    R11, R12, R21, R22 = [], [], [], []
    for stat in statistics:
        R11.append(stat[0])
        R12.append(stat[1])
        R21.append(stat[2])
        R22.append(stat[3])
        summary_R11 = [i for i in R11 if i != 0]
        summary_R12 = [i for i in R12 if i != 0]
        summary_R21 = [i for i in R21 if i != 0]
        summary_R22 = [i for i in R22 if i != 0]
        distances = [summary_R11, summary_R12, summary_R21, summary_R22]
    return distances

# saving chart (bar plot)
functions.bar_plot(read_RvsR_statistics(experiment), experiment)

# saving charts (line plots)
distances = read_distances_statistics(experiment)
functions.line_plot(distances[0], experiment, "R1 vs R1")
functions.line_plot(distances[1], experiment, "R1 vs R2")
functions.line_plot(distances[2], experiment, "R2 vs R1")
functions.line_plot(distances[3], experiment, "R2 vs R2")
