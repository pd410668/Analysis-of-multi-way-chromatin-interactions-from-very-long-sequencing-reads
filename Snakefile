# SAMPLES=["nsp_k562_V_8_R1", "nsp_k562_V_8_R2"]  # human 
SAMPLES=["ms_mesc_III_32_R1", "ms_mesc_III_32_R2"]  # mouse
RES = list(set([i.rsplit("_R")[0] for i in SAMPLES]))

rule all:
	input:
		# human
		# expand("data/charts/k562/{res}_RvsR.png", res=RES),
		# expand("data/charts/k562/{res}_0_5000_R.png", res=RES),
		# expand("data/charts/k562/{res}_500_10000_R.png", res=RES),
		# expand("data/charts/k562/{res}_strand_1vs2.png", res=RES),
		# expand("data/charts/k562/{res}_0_5000_S.png", res=RES),
		# expand("data/charts/k562/{res}_500_10000_S.png", res=RES),
		# expand("data/charts/k562/{res}_log10_500_1000.png", res=RES),
		# expand("data/charts/TAD/k562/tad_fraction_data_human.png"),
		# expand("data/charts/TAD/k562/tad_fraction_random_human.png"),
		# expand("data/charts/TAD/k562/tad_types_human.png"),
		# expand("data/charts/TAD/k562/tad_comparision_human.png"),
		# expand("data/charts/TAD/k562/tad_counts_human.png"),
		# expand("data/charts/TAD/k562/tad_doms_human.png" )
		# expand("data/cwalks/k562/txt/{res}_cwalks.txt", res=RES),  
		# expand("data/cwalks/k562/bed/{res}_cwalks.bed", res=RES),
		# expand("data/charts/cwalks/k562/hist_intra_human.png"),
		# expand("data/charts/cwalks/k562/human_fractions.png"),
		# expand("data/charts/cwalks/k562/cwalk_human_barh.png")
		# expand("data/charts/k562/human_barh.png"),
		# expand("data/charts/directionality/k562")
		# mouse
		# expand("data/charts/mesc/{res}_RvsR.png", res=RES),
		# expand("data/charts/mesc/{res}_0_5000_R.png", res=RES),
		# expand("data/charts/mesc/{res}_500_10000_R.png", res=RES),
		# expand("data/charts/mesc/{res}_strand_1vs2.png", res=RES),
		# expand("data/charts/mesc/{res}_0_5000_S.png", res=RES),
		# expand("data/charts/mesc/{res}_500_10000_S.png", res=RES),
		# expand("data/charts/mesc/{res}_log10_500_1000.png", res=RES),
		# expand("data/charts/TAD/mesc/tad_fraction_data_mouse.png"),
		# expand("data/charts/TAD/mesc/tad_fraction_random_mouse.png"),
		# expand("data/charts/TAD/mesc/tad_types_mouse.png"),
		# expand("data/charts/TAD/mesc/tad_comparision_mouse.png"),
		# expand("data/charts/TAD/mesc/tad_counts_mouse.png"),
		# expand("data/charts/TAD/mesc/tad_doms_mouse.png" )
		# expand("data/cwalks/mesc/bed/{res}_cwalks.bed", res=RES), 
		# expand("data/cwalks/mesc/txt/{res}_cwalks.txt", res=RES)
		# expand("data/charts/cwalks/k562/hist_intra_mouse.png"),
		# expand("data/charts/cwalks/k562/human_fractions.png"),
		# expand("data/charts/cwalks/k562/cwalk_mouse_barh.png")
		# expand("data/charts/mesc/mouse_barh.png"),
		# expand("data/charts/directionality/mesc"),

rule digestion:
	input:
		# "data/fastq/k562/{sample}.fastq"  # human
		"data/fastq/mesc/{sample}.fastq"  # mouse
	output:
		# "data/digested/k562/{sample}.digested.fastq"  # human
		"data/digested/mesc/{sample}.digested.fastq"  # mouse
	shell:
		"src/py/digestion.py {input} {output}"

rule fastq2bam:
	input:
		# index = "data/Bowtie2Index/hg19",  # human
		index = "data/Bowtie2Index/mm9",  # mouse

		# fastq = "data/digested/k562/{sample}.digested.fastq"  # human
		fastq = "data/digested/mesc/{sample}.digested.fastq"  # mouse 
	output:
		# "data/bam/k562/{sample}.bowtie2.bam"  # human
		"data/bam/mesc/{sample}.bowtie2.bam"  # mouse
	shell:
		# "bowtie2 -x {input.index}/hg19 -U {input.fastq} | samtools view -bS -> {output}" # human
		"bowtie2 -x {input.index}/mm9 -U {input.fastq} | samtools view -bS -> {output}" # mouse



# rule bowtie2:  #  bedtools bamtofastq [OPTIONS] -i <BAM> -fq <FASTQ> - command to check (working)
# 	input:
# 		index = "data/Bowtie2Index/hg19",
# 		fastq = "data/digested/k562/{sample}.digested.fastq"
# 	output:
# 		"data/sam/k562/{sample}.sam"
# 	shell:
# 		"bowtie2 -x {input.index}/hg19 -U {input.fastq} -S {output}"

# rule samtools:
# 	input:
# 		"data/sam/k562/{sample}.sam"
# 	output:
# 		"data/bam/k562/{sample}.bowtie2.2.bam"
# 	shell:
# 		"samtools view -u {input} -o {output}"


rule filtering:
	input:
		# expand("data/bam/k562/{sample}.bowtie2.bam", sample=SAMPLES) # human
		expand("data/bam/mesc/{sample}.bowtie2.bam", sample=SAMPLES)  # mouse
	output:
		# "data/supportive/k562/{res}.tsv"  # humman
		"data/supportive/mesc/{res}.tsv"  # mouse
	shell:
		# "src/py/filtering.py human {input} {output}"   # human
		"src/py/filtering.py mouse {input} {output}"   # mouse

rule barh:
	input:
		# "data/supportive/k562"  # human
		"data/supportive/mesc"  # mouse
	output:
		# "data/charts/k562/human_barh.png"  # human
		"data/charts/mesc/mouse_barh.png"  # mouse
	shell:
	 	# "src/py/barh.py human {input} {output}"  # human
		"src/py/barh.py mouse {input} {output}"  # mouse

rule charts:
	input:
		"data/supportive/k562/{res}.tsv"  # human
		# "data/supportive/mesc/{res}.tsv"  # mouse
	output:
	# human
		RvsR_barplot = "data/charts/k562/{res}_RvsR.png",
		dist_0_5000_R = "data/charts/k562/{res}_0_5000_R.png",
		dist_500_10000_R = "data/charts/k562/{res}_500_10000_R.png",
		strand_1vs2_barplot = "data/charts/k562/{res}_strand_1vs2.png",
		dist_0_5000_S = "data/charts/k562/{res}_0_5000_S.png",
		dist_500_10000_S = "data/charts/k562/{res}_500_10000_S.png",
		log10 = "data/charts/k562/{res}_log10_500_1000.png"
	# mouse
		# RvsR_barplot = "data/charts/mesc/{res}_RvsR.png",
		# dist_0_5000_R = "data/charts/mesc/{res}_0_5000_R.png",
		# dist_500_10000_R = "data/charts/mesc/{res}_500_10000_R.png",
		# strand_1vs2_barplot = "data/charts/mesc/{res}_strand_1vs2.png",
		# dist_0_5000_S = "data/charts/mesc/{res}_0_5000_S.png",
		# dist_500_10000_S = "data/charts/mesc/{res}_500_10000_S.png",
		# log10 = "data/charts/mesc/{res}_log10_500_1000.png"

	shell:
		"src/py/charts.py human {input} {output.RvsR_barplot} {output.dist_0_5000_R} {output.dist_500_10000_R} \
		{output.strand_1vs2_barplot} {output.dist_0_5000_S} {output.dist_500_10000_S} {output.log10}"
		# "src/py/charts.py mouse {input} {output.RvsR_barplot} {output.dist_0_5000_R} {output.dist_500_10000_R} \
		# {output.strand_1vs2_barplot} {output.dist_0_5000_S} {output.dist_500_10000_S} {output.log10}"

rule cwalk:
	input:
		# tsvfile = "data/supportive/k562/{res}.tsv",  # human
		tsvfile = "data/supportive/mesc/{res}.tsv",  # mouse

		# bedfile = "data/restrictions/DpnII_hg19.bed"  # human
		bedfile = "data/restrictions/DpnII_mm9.bed"  # mouse

	output:
		# cwalk_txt = "data/cwalks/k562/txt/{res}_cwalks.txt",  # human
		cwalk_txt = "data/cwalks/mesc/txt/{res}_cwalks.txt",  # mouse
		# cwalk_bed = "data/cwalks/k562/bed/{res}_cwalks.bed"  # human 
		cwalk_bed = "data/cwalks/mesc/bed/{res}_cwalks.bed"  # mouse
		
	shell:
		# "src/py/cwalk.py human {input.bedfile} {input.tsvfile} {output.cwalk_txt} {output.cwalk_bed}"  # human 
		"src/py/cwalk.py mouse {input.bedfile} {input.tsvfile} {output.cwalk_txt} {output.cwalk_bed}"  # mouse

rule analysis:
	input:
		"data/cwalks/k562/txt"  # human
		# "data/cwalks/mesc/txt"  # mouse
	output:
		stats_intra = "data/charts/cwalks/k562/hist_all_human.png",  # human
		fractions = "data/charts/cwalks/k562/human_fractions.png",  # human
		barh_all = "data/charts/cwalks/k562/cwalk_human_barh.png"  # human
		# stats = "data/analysis/cwalks/mesc/hist_all_mouse.png",  # mouse
		# fractions = "data/analysis/cwalks/mesc/mouse_fractions.png",  # mouse
		# barh = "data/analysis/cwalks/mesc/cwalk_mouse_barh.png"  # mouse
	shell:
		"src/py/analysis.py human {input} {output.stats_intra} {output.fractions} {output.barh_all}"  # human
		# "src/py/analysis.py mouse {input} {output.stats} {output.fractions} {output.barh}"  # mouse

rule directionality:
	input:
		"data/cwalks/k562/txt"  # human
		# "data/cwalks/mesc/txt"  # mouse
	output:
		"data/charts/directionality/k562"  # human 
		# "data/charts/directionality"  # mouse
	shell:
		"src/py/directionality.py human data/cwalks/k562/txt {output}"  # human
		# "src/py/directionality.py mouse data/cwalks/mesc/txt"  # mouse

rule tad:
	input:
		cwalk = "data/cwalks/k562/txt",  # human
		doms = "data/domains/human_doms.txt",  # human
		genome = "data/genomes/hg19.chrom.sizes",  # human
		rest = "data/restrictions/DpnII_hg19.bed"  # human
		# cwalk = "data/cwalks/mesc/txt",  # mouse
		# doms = "data/domains/mouse_doms.txt", # mouse
		# genome = "data/genomes/mm9.chrom.sizes",  # mouse
		# rest = "data/restrictions/DpnII_mm9.bed"  # mouse

	output:
		frac_data = "data/charts/TAD/k562/tad_fraction_data_human.png",  # human
		frac_random = "data/charts/TAD/k562/tad_fraction_random_human.png", # human
		types = "data/charts/TAD/k562/tad_types_human.png",  # human
		comparision = "data/charts/TAD/k562/tad_comparision_human.png",  # human
		counts = "data/charts/TAD/k562/tad_counts_human.png",  # human
		doms = "data/charts/TAD/k562/tad_doms_human.png"  # human
		# frac_data = "data/charts/TAD/mesc/tad_fraction_data_mouse.png",  # mouse
		# frac_random = "data/charts/TAD/mesc/tad_fraction_random_mouse.png", # mouse
		# types = "data/charts/TAD/mesc/tad_types_mouse.png",  # mouse
		# comparision = "data/charts/TAD/mesc/tad_comparision_mouse.png",  # mouse
		# counts = "data/charts/TAD/mesc/tad_counts_mouse.png",  # mouse
		# doms = "data/charts/TAD/mesc/tad_doms_mouse.png"  # mouse

	shell:
		"src/py/TAD.py human {input.cwalk} {input.doms} {input.genome} {input.rest} \
		{output.frac_data} {output.frac_random} {output.types} {output.comparision} {output.counts} {output.doms}"
