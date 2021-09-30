#!/usr/bin/env python

import pysam
import sys


def parse_bam(inbamfile, read):
    """ load bam two bam files and create list of tuples """
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
    """ filtering not mapped alignment and with quality below 30 """
    return [align for align in alignments if align[4] >= 30]


def collect_data(data, experiment_name, WvsA):
    """ saving supportive file to further analysis """
    import csv
    with open(f"{experiment_name}", WvsA, newline='') as outfile:
        tsv_output = csv.writer(outfile, delimiter='\t')
        tsv_output.writerow(data)


def typical_chromosomes(organism) -> list:
    """ return list of typical chromosomes for human or mouse """
    chrs = [f"chr{i}" for i in range(1, 23)] if organism == "human" else [f"chr{i}" for i in range(1, 20)]
    chrs.extend(["chrX", "chrY"])
    return chrs


if __name__ == '__main__':

    experiment_R1 = sys.argv[1]
    experiment_R2 = sys.argv[2]
    experiment_name = sys.argv[3]

    """ create .tsv file with field names """
    fieldnames = ["seqname", "chr_R1", "chr_R2", "pos_R1", "pos_R2", "strand_1vs2", "RvsR", "abs_pos"]
    collect_data(fieldnames, experiment_name, "w")

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
                if filtered_alignments[i][2] == filtered_alignments[j][2]:  # if two aligns are found on the same chromosome
                    if filtered_alignments[i][2] in typical_chromosomes("human"):  # removed aligns with atypical chromosomes
                        statistics = [
                            filtered_alignments[i][0],  # seqname
                            filtered_alignments[i][2],  # chromosome R1
                            filtered_alignments[j][2],  # chromosome R2
                            filtered_alignments[i][3],  # position R1
                            filtered_alignments[j][3],  # position R2
                            f"{filtered_alignments[i][5]} vs {filtered_alignments[j][5]}",  # strand_1vs2
                            f"{filtered_alignments[i][1]} vs {filtered_alignments[j][1]}",  # RvsR
                            abs(filtered_alignments[i][3] - filtered_alignments[j][3])      # abs_pos
                        ]
                        """ appends aligns to created .tsv file """
                        collect_data(statistics, experiment_name, "a")
