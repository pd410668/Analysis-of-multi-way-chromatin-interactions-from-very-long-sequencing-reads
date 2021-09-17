#!/usr/bin/env python

import sys
from tf import load_cwalk_graph
import networkx as nx
import statistics


def directionality(path: list) -> bool:
    path = [node[0] for node in path] # the distance is between the firsts points of the restriction intervals
    dist = [b - a for a, b in zip(path[:-1], path[1:])]
    same_sign_dist = [(a * b >= 0) for a, b in zip(dist[:-1], dist[1:])]
    # 0 if nodes are in the same restriction interval
    return same_sign_dist


def avg_directionality(signs: list) -> float:
    return round(statistics.mean(signs), 2)


P = load_cwalk_graph("hs_k562_I_1_cwalks.txt")  # load .txt cwalk graph

signs = [directionality(list(cwalk)) for cwalk in list(nx.connected_components(P))]  # list of lists of bool values
avg_signs = [avg_directionality(sign) for sign in signs]  # list of average value of directivity

plus, minus = 0, 0
for avg_sign in avg_signs:
    if avg_sign >= 0.5:
        plus += 1
    else:
        minus += 1

print(plus, minus)  # 1571 2832

# Student’s t-Test for independent samples
chunks_plus = [x for x in avg_signs if x >= 0.5]
chunks_minus = [x for x in avg_signs if x < 0.5]

var_plus = round(statistics.variance(chunks_plus), 2)  # 0.05
var_minus = round(statistics.variance(chunks_minus), 2)  # 0,02

t_test = scipy.stats.ttest_ind(chunks_plus, chunks_minus, equal_var=False)
print(t_test)  # Ttest_indResult(statistic=105.92401155048114, pvalue=0.0), the averages are significantly different
