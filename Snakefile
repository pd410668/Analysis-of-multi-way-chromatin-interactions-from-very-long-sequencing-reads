SAMPLES=["hs_k562_I_1_R1", "hs_k562_I_1_R2"]

rule all:
	input: 
		expand("data/bam/{sample}.bowtie2.sorted.bam", sample=SAMPLES)

rule digestion:
	input:
		"data/fastq/{sample}.fastq"
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
		"data/bam/{sample}.bowtie2.sorted.bam"
	shell:
		"samtools view -q 31 -F 4 -u {input} -o {output}"