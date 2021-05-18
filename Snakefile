SAMPLES=["hs_k562_I_2_R1", "hs_k562_I_2_R2"]

rule all:
	input: 
		expand("data/supportive_filtering/k562_I/{sample}", sample=SAMPLES) # data/bam/k562_I && II

rule digestion:
	input:
		"data/fastq/{sample}.fastq" # data/fastq/k562_I & II & III
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

rule filtering:
	input :
		R1 = "data/bam/k562_I/{sample}.bowtie2.bam",
		R2 = "data/bam/k562_I/{sample}.bowtie2.bam"
		name = "hs_k562_I_2"
	output:
		"data/supportive_filtering/{sample}"
	shell:
		"src/py/./filtering.py {input.R1} {input.R2} {input.name} "

# rule statistics:
# 	input :
# 		"data/supportive_filtering/statistics_{sample}.txt"
# 	output:
# 		"data/test_files/filtering"
# 	shell:
# 		"src/py/./statistics.py {input}	