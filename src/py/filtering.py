import pysam
import itertools

alignments_R1 = pysam.AlignmentFile("R1.bam", "rb")
alignments_R2 = pysam.AlignmentFile("R2.bam", "rb")
input_zip = itertools.zip_longest(alignments_R1, alignments_R2)
outfile = pysam.AlignmentFile("outfile.bam", "wb", template=alignments_R1)

for alignment_R1, alignment_R2 in input_zip:
    try:
        if alignment_R1.qname == alignment_R2.qname or alignment_R1.qname[2:] == alignment_R2.qname[2:] \
            or alignment_R1.qname[2:] == alignment_R2.qname or alignment_R1.qname == alignment_R2.qname[2:]:
            if alignment_R1.qname == alignment_R2.qname and alignment_R1.qname[0].isalpha():
                outfile.write(alignment_R1)
                outfile.write(alignment_R2)
            # if alignment_R1.qname[0].isdigit() and alignment_R2.qname[0].isalpha():

            if alignment_R2.qname[0].isdigit() and alignment_R1.qname[0].isalpha():
                outfile.write(alignment_R1)


    except Exception:
        if alignment_R1 == None:
            outfile.write(alignment_R2)
        if alignment_R2 == None:
            outfile.write(alignment_R1)
