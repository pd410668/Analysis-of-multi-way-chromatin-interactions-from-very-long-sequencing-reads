# experiment on k562 human cells

SAMPLES=["hs_k562_I_1_R1", "hs_k562_I_1_R2"] 
RES = list(set([i.rsplit('_R')[0] for i in SAMPLES]))
EXP = list(set([i.rsplit('_I')[0] for i in SAMPLES]))

rule all:
	input:
		expand("data/cwalks/{res}_cwalks.bed", res=RES),
		expand("data/cwalks/{res}_cwalks.txt", res=RES),

		# plots
		expand("data/analysis/plots/{res}_RvsR.png", res=RES),
		expand("data/analysis/plots/{res}_0_5000_R.png", res=RES),
		expand("data/analysis/plots/{res}_500_10000_R.png", res=RES),
		expand("data/analysis/plots/{res}_strand_1vs2.png", res=RES),
		expand("data/analysis/plots/{res}_0_5000_S.png", res=RES),
		expand("data/analysis/plots/{res}_500_10000_S.png", res=RES),
		expand("data/analysis/plots/{res}_log10_500_1000.png", res=RES),
		expand("data/analysis/plots/{exp}_barh.png", exp=EXP),
		expand("data/analysis/plots/human_barh.png"),
		expand("data/analysis/plots/tf_histplot.png")
		
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
		fastq = "data/digested/{sample}.digested.fastq"
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
		"data/supportive/{res}.tsv"
	shell:
		"src/py/filtering.py {input} {output}"

rule charts:
	input:
		"data/supportive/{res}.tsv"
	output:
		RvsR_barplot = "data/analysis/plots/{res}_RvsR.png",
		dist_0_5000_R = "data/analysis/plots/{res}_0_5000_R.png",
		dist_500_10000_R = "data/analysis/plots/{res}_500_10000_R.png",
		strand_1vs2_barplot = "data/analysis/plots/{res}_strand_1vs2.png",
		dist_0_5000_S = "data/analysis/plots/{res}_0_5000_S.png",
		dist_500_10000_S = "data/analysis/plots/{res}_500_10000_S.png",
		log10 = "data/analysis/plots/{res}_log10_500_1000.png"
	shell:
		"src/py/charts.py {input} {output.RvsR_barplot} {output.dist_0_5000_R} {output.dist_500_10000_R} \
		{output.strand_1vs2_barplot} {output.dist_0_5000_S} {output.dist_500_10000_S} {output.log10}"

rule barh:
	input:
		"data/supportive/"
	output:
		"data/analysis/plots/human_barh.png"
	shell:
		"src/py/barh.py {input} {output}"

rule cwalk:
	input:
		tsv = "data/supportive/{res}.tsv",
		bed = "data/restrictions/DpnII_hg19.bed"
	output:
		cbed = "data/cwalks/{res}_cwalks.bed",
		ctxt = "data/cwalks/{res}_cwalks.txt"
	shell:
		"src/py/cwalk.py {input.tsv} {input.bed} {output.cbed} {output.ctxt}"

rule transcription_factor:
	input:
		ctxt = "data/cwalks/{res}_cwalks.txt",
		tf = "data/factors/wgEncodeAwgTfbsUtaK562CtcfUniPk.narrowPeak.gz"
	output:
		"data/analysis/plots/tf_histplot.png"
	shell:
		"src/py/barh.py {input.ctxt} {input.tf}"
