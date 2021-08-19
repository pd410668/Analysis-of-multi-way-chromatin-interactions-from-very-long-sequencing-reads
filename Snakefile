	# experiment on k562 human cells

SAMPLES=["hs_k562_I_2_R1", "hs_k562_I_2_R2"] 
RES = list(set([i.rsplit('_R')[0] for i in SAMPLES]))

rule all:
	input:
		expand("data/cwalks/{res}_cwalks.txt", res=RES)
		
rule digestion:
	input:
		"data/fastq/k562/hs_k562/{sample}.fastq"
	output:
		"data/fastq_digested/{sample}.digested.fastq"
	shell:
		"src/py/digestion.py {input} {output}"

rule bowtie2:
	input:
		index = "data/Bowtie2Index/hg19",
		fastq = "data/fastq_digested/{sample}.digested.fastq"
	output:
		"data/sam/{sample}.sam"
	shell:
		"bowtie2 -x {input.index}/hg19 -U {input.fastq} -S {output}"

rule samtools:
	input:
		"data/sam/{sample}.sam"
	output:
		"data/bam/{sample}.bowtie2.bam"
	shell:
		"samtools view -u {input} -o {output}"

rule filtering:
	input:
		expand("data/bam/{sample}.bowtie2.bam", sample=SAMPLES)
	output:
		"data/supportive_filtering/{res}.tsv"
	shell:
		"src/py/filtering.py {input} {output}"

rule cwalk:
	input:
		tsvfile = "data/supportive_filtering/{res}.tsv",
		bedfile = "data/restriction_positions/DpnII_hg19.bed"
	output:
		"data/cwalks/{res}_cwalks.txt"
	shell:
		"src/py/cwalk.py {input.tsvfile} {input.bedfile} {output}"