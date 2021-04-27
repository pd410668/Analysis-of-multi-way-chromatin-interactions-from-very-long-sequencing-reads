#!/usr/bin/env python

import pysam
from math import isclose
from prettytable import PrettyTable
import argparse
import sys

def parse_bam(inbamfile, read):
    alignments = pysam.AlignmentFile(inbamfile, "rb")
    alignments_list = []
    strands = ["forward", "reverse"]
    for line in alignments.fetch(until_eof=True):
        if line.flag == 16:
            strand = strands[0]
        else:
            strand = strands[1]
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

        if len(l1) == 1 and len(l2) == 1:
            rows.append(list(l1[0]))
            rows.append(list(l2[0]))
        else:
            for i in range(0, len(l1)):
                for j in range(0, len(l2)):
                    if equal(l1[i][1], l2[j][1]):
                        if equal(l1[i][2], l2[j][2]):
                            if isclose(l1[i][3], l2[j][3], abs_tol=10000):
                                l1.pop(i)

            if len(l1) == 1:
                rows.append((l1[0]))
            else:
                l1_zip = zip(l1, range(0, len(l1)))
                for i, j in l1_zip:
                    seqn = f"{i[0]}.{j}"
                    i = list(i)
                    i[0] = seqn
                    rows.append(i)

            if len(l2) == 1:
                rows.append((l2[0]))
            else:
                l2_zip = zip(l2, range(0, len(l2)))
                for i, j in l2_zip:
                    seqn = f"{i[0]}.{j}"
                    i = list(i)
                    i[0] = seqn
                    rows.append(i)

    filtered_rows = []
    for align in rows:
        if align[4] >= 30:
            filtered_rows.append(align)

    save_as_table(filtered_rows)

if __name__ == '__main__':
    main()
