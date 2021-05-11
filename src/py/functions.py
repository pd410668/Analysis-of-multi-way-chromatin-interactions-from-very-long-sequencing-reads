#!/usr/bin/env python

import pysam
from prettytable import PrettyTable
import matplotlib.pyplot as plt
from math import isclose

def parse_bam(inbamfile, read):
    alignments = pysam.AlignmentFile(inbamfile, "rb")
    alignments_list = []
    for line in alignments.fetch(until_eof=True):
        if line.is_reverse:
            strand = "reverse"
        else:
            strand = "forward"
        name = line.qname
        if line.qname[0].isdigit():
            name = line.qname[2:]
        if line.qname[0].isalpha():
            name = line.qname
        chromosome = f"chr{line.rname + 1}"
        align = (name, f"R{read}", chromosome, line.pos, line.mapq, strand)
        alignments_list.append(align)
    return alignments_list

def save_as_table(rows, name):
    table = PrettyTable()
    table.field_names = ["seqname", "read", "chromosome", "position", "mapq", "strand"]
    table.add_rows(rows)
    data = table.get_string()
    with open(f"table_{name}.txt", "w") as file:
        outfile = file.write(data)
    return outfile

def in_proximity(align_1, align_2):
    # if two alignments are found on the same chromosome
    # and their positions are up to 10,000 bp apart
    return align_1[2] == align_2[2] and isclose(align_1[3], align_1[3], abs_tol=10000)

def collect_statistics(statistics, type, name):
    writer = open(f"{type}_{name}.txt", "a")
    outfile = writer.write((str(statistics) + "\n"))
    return outfile

def cleaning(rows):
    # filtering not mapped alignment and with quality below 30
    filtered_rows = []
    for align in rows:
        if align[4] >= 30:
            filtered_rows.append(align)
    return filtered_rows

def bar_plot(data, experiment):
    plt.bar(["R1 vs R1", "R1 vs R2", "R2 vs R1", "R2 vs R2"], data, width=0.7,
            color="cornflowerblue", edgecolor="navy", linewidth=1.2)
    plt.title(f"Sample from {experiment} human cells", size=20)
    plt.xlabel("Combinations", size=12)
    plt.ylabel("number of matches", size=12)
    return plt.savefig(f"bar_plot_{experiment}.png"), plt.clf()

def line_plot(data, experiment, RvR):
    plt.plot(data, color="cornflowerblue")
    plt.title(f"Sample from {experiment} human cells {RvR}", size=20)
    plt.ylabel("pairwise distance", size=12)
    plt.xlabel("number of comparisons", size=12)
    return plt.savefig(f"line_plot_{experiment}_{RvR}.png"), plt.clf()
