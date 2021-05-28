#!/usr/bin/env python

import pysam
import sys
import csv

"""
filtering.py taking as input two bam files R1, R2
and return experimen_name.tsv file
usage: 
chmod 777 filtering.py
./filtering.py experiment.R1 experiment.R2 experiment_name
"""

def parse_bam(inbamfile, read):
    # load bam two bam files and create list of tuples
    alignments = pysam.AlignmentFile(inbamfile, "rb")
    alignments_list = []
    for line in alignments:
        if line.is_reverse:
            strand = "reverse"
        else:
            strand = "forward"
        name = line.qname
        if line.qname[0].isdigit():
            name = line.qname[2:]
        if line.qname[0].isalpha():
            name = line.qname
        align = (name, f"R{read}", line.reference_name, line.pos, line.mapq, strand)  # order of tuples
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
    # and the absolute value of their position does not exceed 1,000 bp
    return align_1[2] == align_2[2] and abs(align_1[3] - align_2[3]) <= 1000

def collect_statistics(data, experiment_name, WvsA):
    # saving supportive tsv file to make statistics
    with open(f"{experiment_name}", WvsA, newline='') as outfile:
        tsv_output = csv.writer(outfile, delimiter='\t')
        tsv_output.writerow(data)

experiment_R1 = sys.argv[1]
experiment_R2 = sys.argv[2]
experiment_name = sys.argv[3]

# create .tsv file with field names
fieldnames = ["seqname", "chr_R1", "pos_R1", "strand_R1", "chr_R2", "pos_R2", "strand_R1", "RvsR", "abs_pos"]
collect_statistics(fieldnames, experiment_name, "w")

# List of atypical chromosomes found in sample
atypical_chrs = ["chrM",
                "chrUn_gl000232",
                "chrUn_gl000220",
                "chrUn_gl000224",
                "chrUn_gl000234",
                "chrUn_gl000231",
                "chrUn_gl000240",
                "chr7_gl000195_random",
                "chr4_gl000194_random",
                "chr4_gl000193_random",
                "chr19_gl000208_random"
                ]

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
                if in_proximity(filtered_alignments[i], filtered_alignments[j]):
                    if not filtered_alignments[i][2] in atypical_chrs:
                        statistics = [
                            filtered_alignments[i][0],
                            filtered_alignments[i][2],  # chromosome R1
                            filtered_alignments[i][3],  # position R1
                            filtered_alignments[i][5],  # strand R1
                            filtered_alignments[j][2],  # chromosome R2
                            filtered_alignments[j][3],  # position R2
                            filtered_alignments[j][5],  # strand R2
                            f"{filtered_alignments[i][1]} vs {filtered_alignments[j][1]}",
                            abs(filtered_alignments[i][3] - filtered_alignments[j][3])
                        ]
                        # appends aligns to created .tsv file
                        collect_statistics(statistics, experiment_name, "a")

if __name__ == '__main__':
    main()
