SAMPLES=["hs_k562_I_4_R1", "hs_k562_I_4_R2"]
MOUSE_SAMPLES = [] 
RES = list(set([i.rsplit('_R')[0] for i in SAMPLES]))
EXP = list(set([i.rsplit('_I')[0] for i in SAMPLES]))
ORG = ["human, mouse"]

rule all:
	input:
		expand("data/cwalks/human/bed/{res}_cwalks.bed", res=RES),
		expand("data/cwalks/human/txt/{res}_cwalks.txt", res=RES),
		# charts
		expand("data/analysis/charts/{res}_RvsR.png", res=RES),
		expand("data/analysis/charts/{res}_0_5000_R.png", res=RES),
		expand("data/analysis/charts/{res}_500_10000_R.png", res=RES),
		expand("data/analysis/charts/{res}_strand_1vs2.png", res=RES),
		expand("data/analysis/charts/{res}_0_5000_S.png", res=RES),
		expand("data/analysis/charts/{res}_500_10000_S.png", res=RES),
		expand("data/analysis/charts/{res}_log10_500_1000.png", res=RES),
		expand("data/analysis/charts/human_barh.png") # remove the graph every time you run the snakefile
		
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
		"data/supportive/human/{res}.tsv"
	shell:
		"src/py/filtering.py {input} {output}"

rule barh:
	input:
		# expand("data/supportive/{org}", org=ORG)
		"data/supportive/human"
	output:
		# "data/analysis/charts/{org}_barh.png"
		"data/analysis/charts/human_barh.png"
	shell:
		"src/py/barh.py {input} {output}"

rule charts:
	input:
		"data/supportive/human/{res}.tsv"
	output:
		RvsR_barplot = "data/analysis/charts/{res}_RvsR.png",
		dist_0_5000_R = "data/analysis/charts/{res}_0_5000_R.png",
		dist_500_10000_R = "data/analysis/charts/{res}_500_10000_R.png",
		strand_1vs2_barplot = "data/analysis/charts/{res}_strand_1vs2.png",
		dist_0_5000_S = "data/analysis/charts/{res}_0_5000_S.png",
		dist_500_10000_S = "data/analysis/charts/{res}_500_10000_S.png",
		log10 = "data/analysis/charts/{res}_log10_500_1000.png"
	shell:
		"src/py/charts.py {input} {output.RvsR_barplot} {output.dist_0_5000_R} {output.dist_500_10000_R} \
		{output.strand_1vs2_barplot} {output.dist_0_5000_S} {output.dist_500_10000_S} {output.log10}"


rule cwalk:
	input:
		tsvfile = "data/supportive/human/{res}.tsv",
		bedfile = "data/restrictions/DpnII_hg19.bed"
	output:
		cwalkbed = "data/cwalks/human/bed/{res}_cwalks.bed",
		cwalktxt = "data/cwalks/human/txt/{res}_cwalks.txt"
	shell:
		"src/py/cwalk.py {input.tsvfile} {input.bedfile} {output.cwalkbed} {output.cwalktxt}"

# rule transcription_factor:
# 	input:
# 		ctxt = "data/cwalks/{res}_cwalks.txt",
# 		tf = "data/factors/wgEncodeAwgTfbsUtaK562CtcfUniPk.narrowPeak.gz"
# 	output:
# 		"data/analysis/charts/tf_histplot.png"
# 	shell:
# 		"src/py/barh.py {input.ctxt} {input.tf}"
