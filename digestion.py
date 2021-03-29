"""
The script is taking as input fastq file. Next the digest function cuts at the point
of occurrence under the GATC sequence and it returns the digested fastq file.
Usage:
chmod 777 digestion.py
./digestion.py input.fastq output.fastq
"""

#!/usr/bin/env python

from Bio.SeqIO.QualityIO import FastqGeneralIterator
import sys

def division(i, seq, positions) -> str:
    if i == len(positions): sub_seq = seq[positions[i-1]-1:]
    elif i == 0: sub_seq = seq[:positions[i]-1]
    else: sub_seq = seq[positions[i-1]-1:positions[i]-1]
    return sub_seq

def digest(input_reads) -> zip:
    sub_title, sub_seq, sub_qual = [], [], []
    for title, seq, qual in input_reads:
        positions = [i for i in range(len(seq)) if seq.startswith("GATC", i)]
        if positions[0] == 0: positions.pop(0)
        if len(positions) != 0:
            for i in range(0, len(positions)+1):
                sub_title.append(f"{i}.{title}")
                sub_seq.append(division(i, seq, positions))
                sub_qual.append(division(i, qual, positions))
        else:
            sub_title.append(title), sub_seq.append(seq), sub_qual.append(qual)
    return zip(sub_title, sub_seq, sub_qual)

def fastq_write(output_zip):
    with open(sys.argv[2], 'w') as output:
        for title, seq, qual in output_zip:
            output.write(f"@{title}\n{seq}\n+\n{qual}\n")

reads = FastqGeneralIterator(sys.argv[1])
fastq_write(digest(reads))