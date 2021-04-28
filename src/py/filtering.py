#!/usr/bin/env python

import pysam
from math import isclose
from prettytable import PrettyTable
import argparse
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
    table.field_names = ["seqname", "read pair", "chromosome", "position", "mapq", "strand"]
    table.add_rows(rows)
    data = table.get_string()
    with open("table.txt", "w") as file:
        outfile = file.write(data)
    return outfile

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-1", help="input R1 BAM file")
    parser.add_argument("-2", help="input R2 BAM file")
    args = parser.parse_args()

    args.alignments_R1 = parse_bam("R1.bam", 1)
    args.alignments_R2 = parse_bam("R2.bam", 2)

    equal = lambda x, y: x == y

    iter_1 = iter(args.alignments_R1)
    iter_2 = iter(args.alignments_R2)

    next1 = next(iter_1)
    next2 = next(iter_2)

    rows = []
    while next1 is not None:
        seqname = next1[0]
        assert(seqname == next2[0])

        l1 = []
        try:
            while next1[0] == seqname:
                l1.append(next1)
                next1 = next(iter_1)
        except StopIteration:
            next1 = None

        l2 = []
        try:
            while next2[0] == seqname:
                l2.append(next2)
                next2 = next(iter_2)
        except StopIteration:
            next2 = None


        lr = l1 + l2

        def in_proximity(a1, a2):
            # if two alignments are found on the same chromosome
            # and their positions are up to 10,000 bp apart
            return a1[2] == a2[2] and isclose(a1[3], a2[3], abs_tol=10000)                    

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

        # # another idea for statistics to collect:
        # in_proximity(l1[1], l2[1])
        # in_proximity(l1[1], l2[-1])
        # in_proximity(l1[-1], l2[1])
        # in_proximity(l1[-1], l2[-1])


        # append all the alignments to the final table
        rows.extend(lr)

    filtered_rows = []
    for align in rows:
        if align[4] >= 30:
            filtered_rows.append(align)

    save_as_table(filtered_rows)

    print(stats)

if __name__ == '__main__':
    main()
