#!/usr/bin/env python

import sys
from tf import load_cwalk_graph
import networkx as nx


def directivity(path: list):
    path = [node[0] for node in path]  # the distance is between the firsts points of the restriction intervals
    distance = [(a - b) for a, b in zip(path, path[1:])]
    distances.extend(distance)

    # Directivity for each walks separately
    direct = [1 if dist > 0 else 0 for dist in distance]
    directs.append(direct)

    return distances, directs


P = load_cwalk_graph("hs_k562_I_1_cwalks.txt")  # load .txt cwalk graph

distances = []
directs = []
for cwalk in list(nx.connected_components(P)):
    cwalk = list(cwalk)  # [(21404672, 21405189, 'chr14'), (21403394, 21404672, 'chr14')]
    distances, directs = directivity(cwalk)


# Directivity for all cwalks
plus, minus = 0, 0
for dist in distances:
    if dist > 0:
        plus += 1
    else:
        minus += 1

print(distances)  # [57, -648, -1278, 3010850, -628562, -30661, -855, -8394 ...]
print(directs)  # [[1, 1, 0], [1, 0], [0, 1], [0, 1, 0], [0, 1, 0, 1, 1, 1, 1, 0, 0, 0, 1]]
print(plus, minus)  # 7993 8019
