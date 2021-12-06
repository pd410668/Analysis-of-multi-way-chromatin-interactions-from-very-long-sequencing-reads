SAMPLES=["test.digested.R1", "test.digested.R2"]
RES = list(set([i.rsplit('_R')[0] for i in SAMPLES]))
EXP = list(set([i.rsplit('_I')[0] for i in SAMPLES]))

rule all:
	input:
		expand("data/bam/k562/{sample}.bowtie2.bam", sample=SAMPLES)
		# expand("data/cwalks/k562/bed/{res}_cwalks.bed", res=RES),
		# expand("data/cwalks/k562/txt/{res}_cwalks.txt", res=RES),
		
		# # charts
		# expand("data/analysis/charts/k562/{res}_RvsR.png", res=RES),
		# expand("data/analysis/charts/k562/{res}_0_5000_R.png", res=RES),
		# expand("data/analysis/charts/k562/{res}_500_10000_R.png", res=RES),
		# expand("data/analysis/charts/k562/{res}_strand_1vs2.png", res=RES),
		# expand("data/analysis/charts/k562/{res}_0_5000_S.png", res=RES),
		# expand("data/analysis/charts/k562/{res}_500_10000_S.png", res=RES),
		# expand("data/analysis/charts/k562/{res}_log10_500_1000.png", res=RES),
		# expand("data/analysis/charts/k562/human_barh.png")


rule digestion:  # working
	input:
		"data/fastq/k562/{sample}.fastq"
	output:
		"data/digested/k562/{sample}.digested.fastq"
	shell:
		"src/py/digestion.py {input} {output}"

rule bowtie2:  #  bedtools bamtofastq [OPTIONS] -i <BAM> -fq <FASTQ> - command to check (working)
	input:
		index = "data/Bowtie2Index/hg19",
		fastq = "data/digested/k562/{sample}.fastq"
	output:
		"data/sam/k562/{sample}.sam"
	shell:
		"bowtie2 -x {input.index}/hg19 -U {input.fastq} -S {output}"

rule samtools:
	input:
		"data/sam/k562/{sample}.sam"
	output:
		"data/bam/k562/{sample}.bowtie2.bam"
	shell:
		"samtools view -u {input} -o {output}"  # "bowtie2 -x {input.index}/hg19.zip -U {input.fastq} | samtools view -b -o {output}" - can be done faster

# rule filtering:
# 	input:
# 		expand("data/bam/k562/{sample}.bowtie2.bam", sample=SAMPLES)
# 	output:
# 		"data/supportive/k562/{res}.tsv"
# 	shell:
# 		"src/py/filtering.py {input} {output}"

# rule barh:
# 	input:
# 		"data/supportive/k562"
# 	output:
# 		"data/analysis/charts/k562/human_barh.png"
# 	shell:
# 		"src/py/barh.py {input} {output}"

# rule charts:
# 	input:
# 		"data/supportive/k562/{res}.tsv"
# 	output:
# 		RvsR_barplot = "data/analysis/charts/k562/{res}_RvsR.png",
# 		dist_0_5000_R = "data/analysis/charts/k562/{res}_0_5000_R.png",
# 		dist_500_10000_R = "data/analysis/charts/k562/{res}_500_10000_R.png",
# 		strand_1vs2_barplot = "data/analysis/charts/k562/{res}_strand_1vs2.png",
# 		dist_0_5000_S = "data/analysis/charts/k562/{res}_0_5000_S.png",
# 		dist_500_10000_S = "data/analysis/charts/k562/{res}_500_10000_S.png",
# 		log10 = "data/analysis/charts/k562/{res}_log10_500_1000.png"
# 	shell:
# 		"src/py/charts.py {input} {output.RvsR_barplot} {output.dist_0_5000_R} {output.dist_500_10000_R} \
# 		{output.strand_1vs2_barplot} {output.dist_0_5000_S} {output.dist_500_10000_S} {output.log10}"

# rule cwalk:
# 	input:
# 		tsvfile = "data/supportive/k562/{res}.tsv",
# 		bedfile = "data/restrictions/DpnII_hg19.bed"
# 	output:
# 		cwalk_bed = "data/cwalks/k562/bed/{res}_cwalks.bed",
# 		cwalk_txt = "data/cwalks/k562/txt/{res}_cwalks.txt"
# 	shell:
# 		"src/py/cwalk.py {input.tsvfile} {input.bedfile} {output.cwalk_bed} {output.cwalk_txt}"

# rule cwalks_analysis:
# 	input:
# 		"data/cwalks/k562/txt"
# 	output:
# 		"data/analysis/cwalks/cwalks"
# 	shell:
# 		"src/py/cwalks_analysis.py {input} {output}"
