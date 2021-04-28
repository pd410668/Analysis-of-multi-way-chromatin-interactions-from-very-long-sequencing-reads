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

            for i in range(0, len(l1)):
                for j in range(0, len(l2)):
                    if equal(l1[i][1], l2[j][1]):
                        if equal(l1[i][2], l2[j][2]):
                            if isclose(l1[i][3], l2[j][3], abs_tol=10000):
                                l1.pop(i)

        # append all the alignments to the final table
        rows.extend(lr)

    filtered_rows = []
    for align in rows:
        if align[4] >= 30:
            filtered_rows.append(align)

    save_as_table(filtered_rows)

if __name__ == '__main__':
    main()
