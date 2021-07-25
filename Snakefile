# experiment on k562 human cells

SAMPLES=["hs_k562_I_1_R1", "hs_k562_I_1_R2"] 
RES = list(set([i.rsplit('_R')[0] for i in SAMPLES]))
EXP = list(set([i.rsplit('_I')[0] for i in SAMPLES]))

rule all:
	input:
<<<<<<< HEAD
		expand("data/analysis/plots/{res}_RvsR.png", res=RES),
		expand("data/analysis/plots/{res}_0_5000_R.png", res=RES),
		expand("data/analysis/plots/{res}_500_10000_R.png", res=RES),
		expand("data/analysis/plots/{res}_strand_1vs2.png", res=RES),
		expand("data/analysis/plots/{res}_0_5000_S.png", res=RES),
		expand("data/analysis/plots/{res}_500_10000_S.png", res=RES),
		expand("data/analysis/plots/{res}_log10_500_1000.png", res=RES)
		expand("data/analysis/plots/{exp}_barh.png", exp=EXP),
		expand("data/supportive_graph/{res}_graph.txt", res=RES)
=======
		# expand("data/analysis/plots/{res}_RvsR.png", res=RES),
		# expand("data/analysis/plots/{res}_0_5000_R.png", res=RES),
		# expand("data/analysis/plots/{res}_500_10000_R.png", res=RES),
		# expand("data/analysis/plots/{res}_strand_1vs2.png", res=RES),
		# expand("data/analysis/plots/{res}_0_5000_S.png", res=RES),
		# expand("data/analysis/plots/{res}_500_10000_S.png", res=RES),
		# expand("data/analysis/plots/{res}_log10_500_1000.png", res=RES)
		# expand("data/analysis/plots/{exp}_barh.png", exp=EXP),
		expand("data/analysis/plots/{exp}_graph.png")
>>>>>>> 8c88a83f444b5e0113c0165457dfeb7f6e195201
		
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

rule charts:
	input:
		"data/supportive_filtering/{res}.tsv"
	output:
		RvsR_barplot = "data/analysis/plots/{res}_RvsR.png",
		dist_0_5000_R = "data/analysis/plots/{res}_0_5000_R.png",
		dist_500_10000_R = "data/analysis/plots/{res}_500_10000_R.png",
		strand_1vs2_barplot = "data/analysis/plots/{res}_strand_1vs2.png",
		dist_0_5000_S = "data/analysis/plots/{res}_0_5000_S.png",
		dist_500_10000_S = "data/analysis/plots/{res}_500_10000_S.png",
		log10 = "data/analysis/plots/{res}_log10_500_1000.png"
	shell:
		"src/py/charts.py {input} {output.RvsR_barplot} {output.dist_0_5000_R} {output.dist_500_10000_R} {output.strand_1vs2_barplot} {output.dist_0_5000_S} {output.dist_500_10000_S} {output.log10}"

rule statistics:
	input:
		expand("data/supportive_filtering", exp=EXP)
	output:
		barh = "data/analysis/plots/{exp}_barh.png",
		displot = "data/analysis/plots/{exp}_displot.png"
	shell:
		"src/py/statistics.py {input} {output.barh} {output.displot}"

<<<<<<< HEAD
rule graph:
	input:
		tsvfile = "data/supportive_filtering/{res}.tsv"
		bedfile = "data/restriction_positions/DpnII_hg19.bed"
	output:
		"data/supportive_graph/{res}_graph.txt"
	shell:
		"src/py/graph.py {input.tsvfile} {input.bedfile} {output}"  
=======
rule bedtools:
	input:
		ref = "data/reference_fasta/hg19.fa.gz",
		DpnII_hg19 = "data/rest_site_positions/bed/DpnII_hg19.bed"
	output:
		DpnII_hg19_seq = "data/rest_site_positions/seq/DpnII_hg19_seq.fa"
	shell:
		"bedtools getfasta -fi {input.ref} -bed {input.DpnII_hg19} -fo {output.DpnII_hg19_seq}"

rule graph:
	input:
		"data/rest_site_positions/seq/DpnII_hg19_seq.fa"
	output:
		expand("data/analysis/plots/{exp}_graph.png")
	shell:
		"src/py/graph.py {input} {output}"  
>>>>>>> 8c88a83f444b5e0113c0165457dfeb7f6e195201
