SAMPLES=["hs_k562_I_2_R1", "hs_k562_I_2_R2"]

rule all:
	input: 
		expand("data/bam/k562_I/{sample}.bowtie2.bam", sample=SAMPLES) # data/bam/k562_I && II

rule digestion:
	input:
		"data/fastq/k562_I/{sample}.fastq" # data/fastq/k562_I & II & III
	output:
		"data/fastq_digested/{sample}.digested.fastq"
	shell:
		"src/py/./digestion.py {input} {output}"

rule bowtie2:
	input:
		index="data/Bowtie2Index/hg19",
		fastq="data/fastq_digested/{sample}.digested.fastq"
	output:
		"data/sam/{sample}.sam"
	shell:
		"bowtie2 -x {input.index}/hg19 -U {input.fastq} -S {output}"

rule samtools:
	input:
		"data/sam/{sample}.sam"
	output:
		"data/bam/k562_I/{sample}.bowtie2.bam" # data/bam/k562_I & II & III
	shell:
		"samtools view -u {input} -o {output}"
