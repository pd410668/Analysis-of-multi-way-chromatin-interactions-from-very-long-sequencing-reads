#!/usr/bin/env python

import functions
import argparse

def main():

    experiment = "hs_k562_I_1_R1"

    parser = argparse.ArgumentParser()
    parser.add_argument("-1", help="input R1 BAM file")
    parser.add_argument("-2", help="input R2 BAM file")
    args = parser.parse_args()

    args.alignments_R1 = functions.parse_bam(f"{experiment}_R1.bowtie2.bam", 1)
    args.alignments_R2 = functions.parse_bam(f"{experiment}_R2.bowtie2.bam", 2)

    iter_1 = iter(args.alignments_R1)
    iter_2 = iter(args.alignments_R2)

    next1 = next(iter_1)
    next2 = next(iter_2)

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

        # adding rows to table
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

        all_alignments = align_1 + align_2

        for i in range(len(all_alignments)):
            for j in range(i + 1, len(all_alignments)):
                if all_alignments[i][2] == all_alignments[j][2]: # same chromosome
                    statistics = [
                        all_alignments[i][1],
                        all_alignments[j][1],
                        abs(all_alignments[i][3] - all_alignments[j][3]),
                        all_alignments[i][5],
                        all_alignments[j][5]
                    ]

                # saving supportive txt file to make statistics
                functions.collect_statistics(statistics, "RvsR", experiment)

    # filtering not mapped alignment and with quality below 30
    # save to the table
    functions.cleaning(rows)
    functions.save_as_table(functions.cleaning(rows), experiment)

if __name__ == '__main__':
    main()
