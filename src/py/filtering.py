#!/usr/bin/env python

import pysam
from math import isclose
from prettytable import PrettyTable
import argparse
import matplotlib.pyplot as plt
import sys

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

def save_as_table(rows):
    table = PrettyTable()
    table.field_names = ["seqname", "read", "chromosome", "position", "mapq", "strand"]
    table.add_rows(rows)
    data = table.get_string()
    with open("table.txt", "w") as file:
        outfile = file.write(data)
    return outfile

def bar_plot(data):
    plt.bar(["R1 vs R1", "R1 vs R2", "R2 vs R1", "R2 vs R2"], data)
    return plt.savefig("bar_plot.png"), plt.clf()

def line_plot(data):
    plt.plot(data)
    return plt.savefig("line_plot.png"), plt.clf()

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-1", help="input R1 BAM file")
    parser.add_argument("-2", help="input R2 BAM file")
    args = parser.parse_args()

    args.alignments_R1 = parse_bam("R1.bam", 1)
    args.alignments_R2 = parse_bam("R2.bam", 2)

    iter_1 = iter(args.alignments_R1)
    iter_2 = iter(args.alignments_R2)

    next1 = next(iter_1)
    next2 = next(iter_2)

    all_alignments = []
    rows = []
    while next1 is not None:
        seqname = next1[0]
        assert(seqname == next2[0])

        align_1 = []
        try:
            while next1[0] == seqname:
                align_1.append(next1)
                next1 = next(iter_1)
        except StopIteration:
            next1 = None

        align_2 = []
        try:
            while next2[0] == seqname:
                align_2.append(next2)
                next2 = next(iter_2)
        except StopIteration:
            next2 = None

        alignments_list = align_1 + align_2
        all_alignments.extend(alignments_list)

        def in_proximity(align_1, align_2):
            # if two alignments are found on the same chromosome
            # and their positions are up to 10,000 bp apart
            return align_1[2] == align_2[2] and isclose(align_1[3], align_1[3], abs_tol=10000)                   

        # take all alignment pairs
        for i in range(len(lr)):
            for j in range(i + 1, len(lr)):
                if in_proximity(lr[i], lr[j])
                    if lr[i][1] == "R1" and lr[i][1] == "R1":
                        stats[0] += 1
                    if lr[i][1] == "R1" and lr[i][1] == "R2":
                        stats[1] += 1
                    if lr[i][1] == "R2" and lr[i][1] == "R1":
                        stats[2] += 1
                    if lr[i][1] == "R2" and lr[i][1] == "R2":
                        stats[3] += 1
       
        if len(align_1) == 1:
            rows.append(list((align_1[0])))
        else:
            align_1_zip = zip(align_1, range(0, len(align_1)))
            for i, j in align_1_zip:
                seqn = f"{i[0]}.{j}"
                i = list(i)
                i[0] = seqn
                rows.append(i)

        if len(align_2) == 1:
            rows.append((align_2[0]))
        else:
            align_2_zip = zip(align_2, range(0, len(align_2)))
            for i, j in align_2_zip:
                seqn = f"{i[0]}.{j}"
                i = list(i)
                i[0] = seqn
                rows.append(i)   
        
    # filtering not mapped alignment and with quality below 30
    filtered_rows = []
    for align in rows:
        if align[4] >= 30:
            filtered_rows.append(align)

    save_as_table(filtered_rows)
    bar_plot(stats)
    
    chrs = [[] for _ in range(1, 25)]
    for i in range(1, 25):
        for align in all_alignments:
            if align[2] == f"chr{i}":
                chrs[i].append(align[3])

    line_plot(chrs[1])

if __name__ == '__main__':
    main()
