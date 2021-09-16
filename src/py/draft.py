#!/usr/bin/env python

import networkx as nx
from tf import load_cwalk_graph
import statistics


def directionality(path: list) -> bool:
    path = [node[0] for node in path] # the distance is between the firsts points of the restriction intervals
    dist = [b - a for a, b in zip(path[:-1], path[1:])]
    same_sign_dist = [(a * b >= 0) for a, b in zip(dist[:-1], dist[1:])]
    # 0 if nodes are in the same restriction interval
    return same_sign_dist


def avg_directionality(signs: list) -> float:
    return round(statistics.mean(signs), 2)



path = [1, 5, 8, 11, 1, 16]
print(path[:-1])
print(path[1:])
dist = [b - a for a, b in zip(path[:-1], path[1:])]
print(dist)
same_sign_dist = [(a * b >= 0) for a, b in zip(dist[:-1], dist[1:])]
print(same_sign_dist)
print(avg_directionality(same_sign_dist))

"""
[1, 5, 8, 11, 1]
[5, 8, 11, 1, 16]
[4, 3, 3, -10, 15] - dist
[True, True, False, False]
0.5
"""
