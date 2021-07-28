#!/usr/bin/env python

import networkx as nx
import pickle
import sys


def main():
    G = pickle.load(open(sys.argv[1], "rb"))  # Load graph.txt file
    sorted_edges = sorted(G.edges(data=True), key=lambda x: x[2]["weight"], reverse=True)  # Sort edges by read-coverage

    """ Resolve the C-walk graph """

    P = nx.Graph()
    for u, v, a in sorted_edges:
        if (u not in P or P.degree[u] < 2) and (v not in P or P.degree[v] < 2):
            P.add_edge(u, v, weight=a["weight"])

    """ Write to .txt file c-walk graph in binary mode """

    pickle.dump(P, open(sys.argv[2], "wb"))  # .txt outfile


if __name__ == '__main__':
    main()
