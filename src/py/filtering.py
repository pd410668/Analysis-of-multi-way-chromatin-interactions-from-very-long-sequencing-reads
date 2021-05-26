#!/usr/bin/env python

import pysam
import sys
import csv

"""
filtering.py taking as input two bam files R1, R2
and return experimen_name.tsv file
usage: 
chmod 777 filtering.py
./filtering.py experiment.R1 experiment.R2
"""

def parse_bam(inbamfile, read):
    # load bam file and create list of tuples
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

def collect_statistics(data, experiment_name, WvsA):
    # saving supportive tsv file to make statistics
    with open(f"{experiment_name}.tsv", WvsA, newline='') as outfile:
        tsv_output = csv.writer(outfile, delimiter='\t')
        tsv_output.writerow(data)

# initiation of basic dependencies
experiment_R1 = sys.argv[1]
experiment_R2 = sys.argv[2]
experiment_name = sys.argv[3]

fieldnames = ["seqname", "chr_R1", "pos_R1", "strand_R1", "chr_R2", "pos_R2", "strand_R1", "RvsR"]
collect_statistics(fieldnames, experiment_name, "w")

def main():
    alignments_R1 = parse_bam(f"{experiment_R1}", 1)
    alignments_R2 = parse_bam(f"{experiment_R2}", 2)

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

        for i in range(len(filtered_alignments)):
            for j in range(i + 1, len(filtered_alignments)):
                if filtered_alignments[i][2] == filtered_alignments[j][2]:
                    statistics = [
                        filtered_alignments[i][0],
                        filtered_alignments[i][2],  # chromosome R1
                        filtered_alignments[i][3],  # position R1
                        filtered_alignments[i][5],  # strand R1
                        filtered_alignments[j][2],  # chromosome R2
                        filtered_alignments[j][3],  # position R2
                        filtered_alignments[j][5],  # strand R2
                        f"{filtered_alignments[i][1]} vs {filtered_alignments[j][1]}"
                    ]
                    collect_statistics(statistics, experiment_name, "a")

if __name__ == '__main__':
    main()
