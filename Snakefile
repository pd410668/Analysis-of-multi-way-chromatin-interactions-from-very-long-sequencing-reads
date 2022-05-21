SAMPLES=["hs_k562_I_1_R1", "hs_k562_I_1_R2"]  # human
# SAMPLES=["hsp_mesc_I_5_R1", "hsp_mesc_I_5_R2"]  # mouse

RES = list(set([i.rsplit("_R")[0] for i in SAMPLES]))

rule all:
	input:
		expand("data/cwalks/k562/bed/{res}_cwalks.bed", res=RES),  # human
		expand("data/cwalks/k562/txt/{res}_cwalks.txt", res=RES),  # human

		# expand("data/cwalks/mesc/bed/{res}_cwalks.bed", res=RES),  # mouse
		# expand("data/cwalks/mesc/txt/{res}_cwalks.txt", res=RES),  # mouse
		
		# human charts
		expand("data/charts/k562/{res}_RvsR.png", res=RES),
		expand("data/charts/k562/{res}_0_5000_R.png", res=RES),
		expand("data/charts/k562/{res}_500_10000_R.png", res=RES),
		expand("data/charts/k562/{res}_strand_1vs2.png", res=RES),
		expand("data/charts/k562/{res}_0_5000_S.png", res=RES),
		expand("data/charts/k562/{res}_500_10000_S.png", res=RES),
		expand("data/charts/k562/{res}_log10_500_1000.png", res=RES),
		expand("data/charts/k562/human_barh.png")

		# mouse charts
		# expand("data/charts/mm9/{res}_RvsR.png", res=RES),
		# expand("data/charts/mm9/{res}_0_5000_R.png", res=RES),
		# expand("data/charts/mm9/{res}_500_10000_R.png", res=RES),
		# expand("data/charts/mm9/{res}_strand_1vs2.png", res=RES),
		# expand("data/charts/mm9/{res}_0_5000_S.png", res=RES),
		# expand("data/charts/mm9/{res}_500_10000_S.png", res=RES),
		# expand("data/charts/mm9/{res}_log10_500_1000.png", res=RES),
		# expand("data/charts/mm9/mouse_barh.png")


rule digestion:
	input:
		"data/fastq/k562/{sample}.fastq"  # human
		# "data/fastq/mesc/{sample}.fastq"  # mouse
	output:
		"data/digested/k562/{sample}.digested.fastq"  # human
		# "data/digested/mesc/{sample}.digested.fastq"  # mouse
	shell:
		"src/py/digestion.py {input} {output}"

rule fastq2bam:
	input:
		index = "data/Bowtie2Index/hg19",  # human
		# index = "data/Bowtie2Index/mm9",  # mouse

		fastq = "data/digested/k562/{sample}.digested.fastq"  # human
		# fastq = "data/digested/mesc/{sample}.digested.fastq"  # mouse 
	output:
		"data/bam/k562/{sample}.bowtie2.bam"  # human
		# "data/bam/mesc/{sample}.bowtie2.bam"  # mouse
	shell:
		"bowtie2 -x {input.index}/hg19 -U {input.fastq} | samtools view -bS -> {output}" # human
		# "bowtie2 -x {input.index}/mm9 -U {input.fastq} | samtools view -bS -> {output}" # mouse


rule filtering:
	input:
		expand("data/bam/k562/{sample}.bowtie2.bam", sample=SAMPLES) # human
		# expand("data/bam/mesc/{sample}.bowtie2.bam", sample=SAMPLES)  # mouse
	output:
		"data/supportive/k562/{res}.tsv"  # humman
		# "data/supportive/mesc/{res}.tsv"  # mouse
	shell:
		"src/py/filtering.py human {input} {output}"   # human
		# "src/py/filtering.py mouse {input} {output}"   # mouse

rule barh:
	input:
		"data/supportive/k562"  # human
		# "data/supportive/mesc"  # mouse
	output:
		"data/charts/k562/human_barh.png"  # human
		# "data/charts/mesc/mouse_barh.png"  # mouse
	shell:
		"src/py/barh.py {input} {output}"

rule charts:
	input:
		"data/supportive/k562/{res}.tsv"
	output:
		RvsR_barplot = "data/charts/k562/{res}_RvsR.png",
		dist_0_5000_R = "data/charts/k562/{res}_0_5000_R.png",
		dist_500_10000_R = "data/charts/k562/{res}_500_10000_R.png",
		strand_1vs2_barplot = "data/charts/k562/{res}_strand_1vs2.png",
		dist_0_5000_S = "data/charts/k562/{res}_0_5000_S.png",
		dist_500_10000_S = "data/charts/k562/{res}_500_10000_S.png",
		log10 = "data/charts/k562/{res}_log10_500_1000.png"
	shell:
		"src/py/charts.py {input} {output.RvsR_barplot} {output.dist_0_5000_R} {output.dist_500_10000_R} \
		{output.strand_1vs2_barplot} {output.dist_0_5000_S} {output.dist_500_10000_S} {output.log10}"

rule cwalk:
	input:
		tsvfile = "data/supportive/k562/{res}.tsv",  # human
		# tsvfile = "data/supportive/mesc/{res}.tsv",  # mouse

		bedfile = "data/restrictions/DpnII_hg19.bed"  # human
		# bedfile = "data/restrictions/DpnII_mm9.bed"  # mouse

	output:
		cwalk_bed = "data/cwalks/k562/bed/{res}_cwalks.bed",  # human
		# cwalk_bed = "data/cwalks/k562/bed/{res}_cwalks.bed",  # mouse

		cwalk_txt = "data/cwalks/k562/txt/{res}_cwalks.txt"  # human
		# cwalk_txt = "data/cwalks/mesc/txt/{res}_cwalks.txt"  # mouse
	shell:
		"src/py/cwalk.py human {input.bedfile} {input.tsvfile} {output.cwalk_bed} {output.cwalk_txt}"  # human
		# "src/py/cwalk.py mouse {input.bedfile} {input.tsvfile} {output.cwalk_bed} {output.cwalk_txt}"  # mouse
