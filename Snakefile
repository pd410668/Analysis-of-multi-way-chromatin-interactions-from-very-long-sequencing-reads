# experiment on k562 human cells

SAMPLES=["hs_k562_I_1_R1", "hs_k562_I_1_R2"] 
RES = list(set([i.rsplit('_R')[0] for i in SAMPLES]))

rule all:
	input:
		expand("data/cwalk_graph/{res}_cwalk_graph.txt", res=RES)
		
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

rule graph:
	input:
		tsvfile = "data/supportive_filtering/{res}.tsv",
		bedfile = "data/restriction_positions/DpnII_hg19.bed"
	output:
		"data/supportive_graph/{res}_graph.txt"
	shell:
		"src/py/graph.py {input.tsvfile} {input.bedfile} {output}"

rule cwalk:
	input:
		"data/supportive_graph/{res}_graph.txt"
	output:
		"data/cwalk_graph/{res}_cwalk_graph.txt"
	shell:
		"src/py/cwalk.py {input} {output}"