#!/usr/bin/env python

import pysam
from math import isclose
import sys

"""
filtering.py taking as input two bam files R1, R2 and name of 
return statistics_experiment.txt
usage:
chmod 777 filtering.py
./filtering.py experiment_name
"""

def parse_bam(inbamfile, read):
    """
    :param inbamfile: load two bam files R1 and R2 and name of experimeent
    :param read: number of read
    :return: list od tuples
    """
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
        align = (name, f"R{read}", chromosome, line.pos, line.mapq, strand) # order of tuples
        alignments_list.append(align)
    return alignments_list

def cleaning(alignments):
    # filtering not mapped alignment and with quality below 30
    filtered_alignments = []
    for align in alignments:
        if align[4] >= 30:
            filtered_alignments.append(align)
    return filtered_alignments

def in_proximity(align_1, align_2):
    # if two alignments are found on the same chromosome
    # and their positions are up to 1,000 bp apart
    return align_1[2] == align_2[2] and isclose(align_1[3], align_1[3], abs_tol=1000)

def collect_statistics(statistics, name):
    # saving supportive txt file to make statistics
    writer = open(f"statistics_{name}.txt", "a")
    outfile = writer.write((str(statistics) + "\n"))
    return outfile

def main():
    experiment = sys.argv[1]
#     experiment_R2 = sys.argv[2]
#     name_experiment = sys.argv[3]

    alignments_R1 = parse_bam(f"{experiment}_R1.bowtie2.bam", 1)
    alignments_R2 = parse_bam(f"{experiment}_R2.bowtie2.bam", 2)

    iter_1 = iter(alignments_R1)
    iter_2 = iter(alignments_R2)

    next1 = next(iter_1)
    next2 = next(iter_2)

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

        all_alignments = align_1 + align_2
        filtered_alignments = cleaning(all_alignments)

        statistics = []
        R_counts = [0, 0, 0, 0]
        for i in range(len(filtered_alignments)):
            for j in range(i + 1, len(filtered_alignments)):
                if in_proximity(filtered_alignments[i], filtered_alignments[j]):
                    if abs(filtered_alignments[i][3] - filtered_alignments[j][3]) != 0:
                        if all_alignments[i][1] == "R1" and all_alignments[j][1] == "R1":
                            R_counts[0] += 1
                        if all_alignments[i][1] == "R1" and all_alignments[j][1] == "R2":
                            R_counts[1] += 1
                        if all_alignments[i][1] == "R2" and all_alignments[j][1] == "R1":
                            R_counts[2] += 1
                        if all_alignments[i][1] == "R2" and all_alignments[j][1] == "R2":
                            R_counts[3] += 1
                        for match in R_counts:
                            if match != 0:
                                statistics = [
                                    filtered_alignments[i][0],
                                    abs(filtered_alignments[i][3] - filtered_alignments[j][3]),
                                    filtered_alignments[i][5],
                                    filtered_alignments[j][5],
                                    f"{filtered_alignments[i][1]} vs {filtered_alignments[j][1]}",
                                    match
                                ]

                        collect_statistics(statistics, experiment)

if __name__ == '__main__':
    main()

