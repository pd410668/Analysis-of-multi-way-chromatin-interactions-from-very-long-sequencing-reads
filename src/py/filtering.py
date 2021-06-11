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
    return [align for align in alignments if align[4] >= 30]

def collect_statistics(data, experiment_name, WvsA):
    # saving supportive tsv file to make statistics
    with open(f"{experiment_name}", WvsA, newline='') as outfile:
        tsv_output = csv.writer(outfile, delimiter='\t')
        tsv_output.writerow(data)
        
def main():
    experiment_R1 = sys.argv[1]
    experiment_R2 = sys.argv[2]
    experiment_name = sys.argv[3]

    # create .tsv file with field names
    fieldnames = ["seqname", "chr", "pos_R1", "pos_R2", "strand_1vs2", "RvsR", "abs_pos"]
    collect_statistics(fieldnames, experiment_name, "w")

    # List of all atypical chromosomes
    atypical_chrs = ["chrM",
                    "chr1_gl000191_random",
                    "chr1_gl000192_random",
                    "chr4_gl000193_random",
                    "chr4_gl000194_random",
                    "chr7_gl000195_random",
                    "chr8_gl000196_random",
                    "chr8_gl000197_random",
                    "chr9_gl000198_random",
                    "chr9_gl000199_random",
                    "chr9_gl000200_random",
                    "chr9_gl000201_random",
                    "chr11_gl000202_random",
                    "chr17_gl000203_random",
                    "chr17_gl000204_random",
                    "chr17_gl000205_random",
                    "chr17_gl000206_random",
                    "chr18_gl000207_random",
                    "chr19_gl000208_random",
                    "chr19_gl000209_random",
                    "chr21_gl000210_random",
                    "chrUn_gl000211",
                    "chrUn_gl000212",
                    "chrUn_gl000213",
                    "chrUn_gl000214",
                    "chrUn_gl000215",
                    "chrUn_gl000216",
                    "chrUn_gl000217",
                    "chrUn_gl000218",
                    "chrUn_gl000219",
                    "chrUn_gl000243",
                    "chrUn_gl000244",
                    "chrUn_gl000245",
                    "chrUn_gl000246",
                    "chrUn_gl000247",
                    "chrUn_gl000248",
                    "chrUn_gl000249",
                    "chrUn_gl000220",
                    "chrUn_gl000221",
                    "chrUn_gl000222",
                    "chrUn_gl000223",
                    "chrUn_gl000224",
                    "chrUn_gl000225",
                    "chrUn_gl000226",
                    "chrUn_gl000227",
                    "chrUn_gl000228",
                    "chrUn_gl000229",
                    "chrUn_gl000230",
                    "chrUn_gl000231",
                    "chrUn_gl000232",
                    "chrUn_gl000233",
                    "chrUn_gl000234",
                    "chrUn_gl000235",
                    "chrUn_gl000236",
                    "chrUn_gl000237",
                    "chrUn_gl000238",
                    "chrUn_gl000239",
                    "chrUn_gl000240",
                    "chrUn_gl000241",
                    "chrUn_gl000242",
                    "chrUn_gl000243",
                    "chrUn_gl000244",
                    "chrUn_gl000245",
                    "chrUn_gl000246",
                    "chrUn_gl000247",
                    "chrUn_gl000248",
                    "chrUn_gl000249"
                    ]

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
                    if not filtered_alignments[i][2] in atypical_chrs:  # removed aligns with atypical chromosomes
                        statistics = [
                            filtered_alignments[i][0],  # seqname
                            filtered_alignments[i][2],  # chromosome
                            filtered_alignments[i][3],  # position R1
                            filtered_alignments[j][3],  # position R2
                            f"{filtered_alignments[i][5]} vs {filtered_alignments[j][5]}",  # strand_1vs2
                            f"{filtered_alignments[i][1]} vs {filtered_alignments[j][1]}",  # RvsR
                            abs(filtered_alignments[i][3] - filtered_alignments[j][3])  # abs_pos
                        ]
                        # appends aligns to created .tsv file
                        collect_statistics(statistics, experiment_name, "a")
                        
if __name__ == '__main__':
    main()
