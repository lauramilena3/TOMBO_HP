import os
import re
import pandas as pd
#======================================================
# Config files
#======================================================
configfile: "config.yaml"

#======================================================
# Global variables
#======================================================

OUTPUT_DIR=config['work_dir'].rstrip("/")
RAW_DATA_DIR =config['input_dir']
FLOWCELL=config['flowcell']
KIT=config['kit']
MODEL_GUPPY=config['model']
HACKED=config['hacked']
ALTERNATIVE_MODELS =config['alternative_models']

GENOME_dir=config['genome_dir']
if len(GENOME_dir)==0:
	GENOME_dir=OUTPUT_DIR + "/GENOMES"
GENOME_name,=glob_wildcards(GENOME_dir + "/{GENOME}"+ ".fasta")

SAMPLE_SHEET=config['sample_sheet']
SAMPLES=config['samples'].split()
CONTROL=config['control'].split()
BARCODES=SAMPLES+CONTROL
BARCODES = list( dict.fromkeys(BARCODES) )
BARCODES = list(filter(None, BARCODES))
MAPPING_TYPES= ["default", "loose"]

dir_list = ["RULES_DIR","ENVS_DIR","DB", "TOOLS", "GENOMES", "BASECALLED", "MAPPING", "QC", "DEMULTIPLEXED", "SINGLE", "TOMBO", "ASSEMBLY_DIR", "PLOTS_DIR", "RAW_NOTEBOOKS","NOTEBOOKS_DIR"]
dir_names = ["rules", "../envs", OUTPUT_DIR + "/db", OUTPUT_DIR + "/tools" ,GENOME_dir , OUTPUT_DIR + "/01_BASECALLED", OUTPUT_DIR + "/02_MAPPING", OUTPUT_DIR + "/03_QC", OUTPUT_DIR + "/04_DEMULTIPLEXED", OUTPUT_DIR + "/05_FAST5_SINGLE" , OUTPUT_DIR + "/06_TOMBO", OUTPUT_DIR + "/07_ASSEMBLY", OUTPUT_DIR + "/FIGURES_AND_TABLES", "notebooks", OUTPUT_DIR + "/NOTEBOOKS"]
dirs_dict = dict(zip(dir_list, dir_names))
MODEL_MEGALODON=config['megalodon_model']
#SAMPLES,=glob_wildcards(RAW_DATA_DIR + "/{{input.sample}}_" +".fast5")
print("outdir", OUTPUT_DIR)
print("samples", SAMPLES)
print("control", CONTROL)
print("barcodes", BARCODES)
print("genome_dir", GENOME_dir)
print("genome", GENOME_name)
#======================================================
# Rules
#======================================================
wildcard_constraints:
	barcode="barcode..",

rule all:
	input:
		#expand(dirs_dict["GUPPY"] + "/{barcode}/fastq/{barcode}.fastq",barcode=BARCODES),
#		cp fastq_runid_*{params.barcode_number}_0.fastq {output.basecalled}
#		directory(expand(dirs_dict["BASECALLED"] + "/{barcode}"), barcode=BARCODES),
#		expand(dirs_dict["TOMBO"] + "/resquiggled_checkpoint_{barcode}.txt", barcode=BARCODES),
#		expand(dirs_dict["BASECALLED"] + "/annotated_checkpoint_{barcode}.txt", barcode=BARCODES),
#		expand(dirs_dict["SINGLE"] + "/{barcode}", barcode=BARCODES),
#		expand(dirs_dict["MEGALODON"] + "/{barcode}", barcode=BARCODES),
#		expand(dirs_dict["DEEPSIGNAL"] + "/{barcode}_deepsignal-prob.tsv", barcode=BARCODES),
		expand(dirs_dict["TOMBO"] + "/{genome}_{barcode}.tombo_denovo/{genome}_{barcode}_denovo.tombo.stats", barcode=BARCODES, genome=GENOME_name),
		expand(dirs_dict["TOMBO"] + "/{genome}_{sample}_{control}.tombo_sampleCompare/{genome}_{sample}_{control}.tombo.stats", sample=SAMPLES, control=CONTROL, genome=GENOME_name),
		expand(dirs_dict["QC"] + "/{genome}_{barcode}_nanoQC", barcode=BARCODES, genome=GENOME_name),
		expand(dirs_dict["PLOTS_DIR"] + "/sampleCompare_{genome}_{sample}_{control}_histogram_pentanucleotide.pdf", sample=SAMPLES, control=CONTROL, genome=GENOME_name),
		expand(dirs_dict["PLOTS_DIR"] + "/denovo_{genome}_{sample}_histogram_pentanucleotide.pdf", sample=SAMPLES, genome=GENOME_name),

def input_modifications_batch(wildcards):
# Read counts
	sample_sheet=pd.read_csv(SAMPLE_SHEET, sep="\t")
	inputs=[]
	print(sample_sheet)
	for index,row in sample_sheet.iterrows():
		row_sample=row["sample"]
		row_control=row["control"]
		row_genome=row["genome"]
		inputs.extend(expand(dirs_dict["QC"] + "/{genome}_{sample}_{mapping}_nanoQC", sample=[row_sample], genome=[row_genome], mapping=MAPPING_TYPES)),
		inputs.extend(expand(dirs_dict["QC"] + "/{genome}_{control}_{mapping}_nanoQC", control=[row_control], genome=[row_genome], mapping=MAPPING_TYPES)),
		inputs.extend(expand(dirs_dict["PLOTS_DIR"] + "/{genome}/sampleCompare_{mapping}/sampleCompare_{genome}_{sample}_{control}_{mapping}_histogram_dinucleotide.pdf", sample=row_sample, control=row_control, genome=row_genome, mapping=[MAPPING_TYPES])),
		inputs.extend(expand(dirs_dict["PLOTS_DIR"] + "/{genome}/denovo_{mapping}/denovo_{genome}_{sample}_{mapping}_histogram_dinucleotide.pdf", sample=row_sample, genome=row_genome, mapping=[MAPPING_TYPES])),
		# inputs.extend(expand(dirs_dict["TOMBO"] + "/resquiggled_checkpoint_{genome}_{sample}_loose.txt",sample=row_sample, genome=row_genome)),

	return inputs

rule run_modifications_batch:
	input:
		input_modifications_batch,

rule demultiplex_run:
	input:
		expand(dirs_dict["DEMULTIPLEXED"] + "/{barcode}_checkpoint.txt", barcode=BARCODES),

# rule megalodon_run:
# 	input:
# 		expand(dirs_dict["MEGALODON"] + "/{barcode}", barcode=BARCODES),

# rule deepsignal_run:
# 	input:
# 		expand(dirs_dict["DEEPSIGNAL"] + "/{barcode}_deepsignal-prob.tsv", barcode=BARCODES),

rule tombo_run_denovo:
	input:
		expand(dirs_dict["TOMBO"] + "/{genome}_{barcode}.tombo_denovo/{genome}_{barcode}_denovo.tombo.stats", barcode=BARCODES, genome=GENOME_name),
		expand(dirs_dict["QC"] + "/{genome}_{barcode}_nanoQC", barcode=BARCODES, genome=GENOME_name),

rule tombo_run_sampleCompare:
	input:
		expand(dirs_dict["TOMBO"] + "/{genome}_{sample}_{control}.tombo_sampleCompare/{genome}_{sample}_{control}.tombo.stats", sample=SAMPLES, control=CONTROL, genome=GENOME_name),
		expand(dirs_dict["QC"] + "/{genome}_{barcode}_nanoQC", barcode=BARCODES, genome=GENOME_name),

rule tombo_run_alternative:
	input:
		expand(dirs_dict["TOMBO"] + "/{genome}_{sample}.tombo_alternative_{model}/{genome}_{sample}.tombo_alternative_{model}.{model}.tombo.stats", sample=SAMPLES, model=ALTERNATIVE_MODELS, genome=GENOME_name),

rule get_rerio_model:
	output:
		rerio_dir=directory(dirs_dict["TOOLS"]+ "/rerio"),
	params:
		tools_dir=dirs_dict["TOOLS"],
	conda:
		"envs/env1.yaml"
	message:
		"Get Deepbinner"
	shell:
		"""
		git clone https://github.com/nanoporetech/rerio {params.tools_dir}/rerio
		cd {params.tools_dir}/rerio
		./download_model.py basecall_models/res_dna_r941_min_modbases_5mC_CpG_v001
		"""

# if HACKED:

# 	rule guppy_basecalling:
# 		input:
# 			raw_data=RAW_DATA_DIR,
# 		output:
# 	#		demultiplexed_dir=directory(expand((dirs_dict["DEMULTIPLEXED"] + "/{barcode}"), barcode=BARCODES)),
# 			basecalled_summary=dirs_dict["BASECALLED"] + "/sequencing_summary.txt",
# 			basecalled_dir=directory(expand(dirs_dict["BASECALLED"] + "/{barcode}",barcode=BARCODES)),
# 			workspace_dir=directory(dirs_dict["BASECALLED"] + "/workspace"),
# 		params:
# 			flowcell=FLOWCELL,
# 			kit=KIT,
# 			model=MODEL_GUPPY,
# 			basecalled_dir=directory(dirs_dict["BASECALLED"]),
# 		conda:
# 			"envs/env1.yaml"
# 		message:
# 			"Basecalling single fast5 files with guppy"
# 		threads: 32
# 		shell:
# 			"""
# 			/home/krakenosh/my_scripts/guppy_V5.0.16/bin/guppy_basecaller -i {input.raw_data} -s {params.basecalled_dir} -q 0 -r -x 'cuda:0' -c /home/krakenosh/my_scripts/ont-guppy/data/{params.model} --barcode_kits {params.kit} --fast5_out --chunks_per_runner 160
# 			guppy_barcoder  -i {params.basecalled_dir}  -t {threads} --kit SQK16S-GXO -s {params.basecalled_dir} -r			
# 			"""


rule guppy_basecalling:
	input:
		raw_data=RAW_DATA_DIR,
	output:
#		demultiplexed_dir=directory(expand((dirs_dict["DEMULTIPLEXED"] + "/{barcode}"), barcode=BARCODES)),
		basecalled_summary=dirs_dict["BASECALLED"] + "/sequencing_summary.txt",
		basecalled_dir=directory(expand(dirs_dict["BASECALLED"] + "/{barcode}",barcode=BARCODES)),
		workspace_dir=directory(dirs_dict["BASECALLED"] + "/workspace"),
	params:
		flowcell=FLOWCELL,
		kit=KIT,
		model=MODEL_GUPPY,
		basecalled_dir=directory(dirs_dict["BASECALLED"]),
	conda:
		"envs/env1.yaml"
	message:
		"Basecalling single fast5 files with guppy"
	threads: 32
	shell:
		"""
		guppy_basecaller -i {input.raw_data} -s {params.basecalled_dir} -q 0 -r -x 'cuda:0 cuda:1' -c /opt/ont/guppy/data/{params.model} --barcode_kits {params.kit} --post_out --disable_qscore_filtering --chunks_per_runner 32 --gpu_runners_per_device 2 --num_callers 8
		# --chunks_per_runner 1 --gpu_runners_per_device 1 --num_callers 1
		# --chunks_per_runner 128
		"""

rule merge_fastq:
	input:
		basecalled_dir=dirs_dict["BASECALLED"] + "/{barcode}",
	output:
		merged_fastq=temp(dirs_dict["QC"] + "/{genome}_{barcode}_merged.fastq"),
		merged_fastq_porechopped=(dirs_dict["QC"] + "/{genome}_{barcode}_merged_porechop.fastq"),
	conda:
		"envs/env3.yaml"
	message:
		"Merging fastq files"
	threads: 2
	shell:
		"""
		cat {input.basecalled_dir}/*fastq > {output.merged_fastq}
		porechop -i {output.merged_fastq} -o {output.merged_fastq_porechopped} --threads {threads}
		"""

rule map_to_genomes_default:
	input:
		merged_fastq_porechopped=(dirs_dict["QC"] + "/{genome}_{barcode}_merged_porechop.fastq"),
		genome=GENOME_dir + "/{genome}.fasta",
	output:
		sam=dirs_dict["MAPPING"] + "/{genome}_{barcode}_default.sam",
	conda:
		"envs/env1.yaml"
	message:
		"Mapping reads to genomes with minimap2"
	threads: 4
	shell:
		"""
		# DEFAULT MAPPING
		minimap2 -t {threads} -k 15 -w 10 -B 4 -O 4,24 -ax map-ont {input.genome} {input.merged_fastq_porechopped} > {output.sam}
		"""

rule map_to_genomes_loose:
	input:
		merged_fastq_porechopped=(dirs_dict["QC"] + "/{genome}_{barcode}_merged_porechop.fastq"),
		genome=GENOME_dir + "/{genome}.fasta",
	output:
		sam=dirs_dict["MAPPING"] + "/{genome}_{barcode}_loose.sam",
	conda:
		"envs/env1.yaml"
	message:
		"Mapping reads to genomes with minimap2"
	threads: 4
	shell:
		"""
		# LOOSE MAPPING
		minimap2 -t {threads} -k 10 -w 3 -B 3 -O 2,8 -ax map-ont {input.genome} {input.merged_fastq_porechopped} > {output.sam}
		"""

rule genome_stats:
	input:
		merged_fastq_porechopped=(dirs_dict["QC"] + "/{genome}_{barcode}_merged_porechop.fastq"),
		sam=dirs_dict["MAPPING"] + "/{genome}_{barcode}_{mapping}.sam",
	output:
		bam=dirs_dict["MAPPING"] + "/{genome}_{barcode}_{mapping}_sorted.bam",
		plus_cov=dirs_dict["MAPPING"] + "/{genome}_{barcode}_{mapping}_coverage_plus.bedgraph",
		minus_cov=dirs_dict["MAPPING"] + "/{genome}_{barcode}_{mapping}_coverage_minus.bedgraph",
		mapped_list_forward=dirs_dict["MAPPING"] + "/{genome}_{barcode}_{mapping}_list_mapped_forward.txt",
		mapped_list_reverse=dirs_dict["MAPPING"] + "/{genome}_{barcode}_{mapping}_list_mapped_reverse.txt",
		mapped_fastq_forward=(dirs_dict["MAPPING"] + "/{genome}_{barcode}_{mapping}_mapped_forward.fastq"),
		mapped_fastq_reverse=(dirs_dict["MAPPING"] + "/{genome}_{barcode}_{mapping}_mapped_reverse.fastq"),
		mapped_list=dirs_dict["MAPPING"] + "/{genome}_{barcode}_{mapping}_fast5_list_mapped.txt",
		mapped_fastq=(dirs_dict["MAPPING"] + "/{genome}_{barcode}_{mapping}_mapped.fastq"),
	conda:
		"envs/env1.yaml"
	message:
		"Mapping statistics"
	threads: 1
	shell:
		"""
		samtools view -bS {input.sam} | samtools sort -o {output.bam}
		samtools index {output.bam}
		
		bedtools genomecov -ibam {output.bam} -bg -strand + > {output.plus_cov}
		bedtools genomecov -ibam {output.bam} -bg -strand - > {output.minus_cov}
		
		samtools view -F 0x10 {output.bam} | cut -f1 > {output.mapped_list_forward}
		samtools view -f 0x10 {output.bam} | cut -f1 > {output.mapped_list_reverse}
		cat {output.mapped_list_forward} {output.mapped_list_reverse} > {output.mapped_list}

		seqtk subseq {input.merged_fastq_porechopped} {output.mapped_list_forward} > {output.mapped_fastq_forward}
		seqtk subseq {input.merged_fastq_porechopped} {output.mapped_list_reverse} > {output.mapped_fastq_reverse}
		seqtk subseq {input.merged_fastq_porechopped} {output.mapped_list} > {output.mapped_fastq}
		"""

rule qualityCheckNanopore:
	input:
		merged_fastq_porechopped=(dirs_dict["QC"] + "/{genome}_{barcode}_merged_porechop.fastq"),
	output:
		nanoqc_dir=directory(dirs_dict["QC"] + "/{genome}_{barcode}_{mapping}_nanoQC"),
	message:
		"Performing nanoQC statistics"
	conda:
		"envs/env3.yaml"
	shell:
		"""
		nanoQC -o {output.nanoqc_dir} {input.merged_fastq_porechopped}
		"""

# rule asemblyFlye:
# 	input:
# 		mapped_fastq_forward=(dirs_dict["BASECALLED"] + "/{genome}_{barcode}_{mapping}_mapped_forward.fastq"),
# 		mapped_fastq_reverse=(dirs_dict["BASECALLED"] + "/{genome}_{barcode}_{mapping}_mapped_reverse.fastq"),
# 		merged_fastq_porechopped=(dirs_dict["BASECALLED"] + "/{genome}_{barcode}_merged_porechop.fastq"),
# 	output:
# 		scaffolds=dirs_dict["ASSEMBLY_DIR"] + "/flye_{sample}_{sampling}/assembly.fasta",
# 		scaffolds_final=dirs_dict["ASSEMBLY_DIR"] + "/{sample}_contigs_flye.{sampling}.fasta"
# 	message:
# 		"Assembling Nanopore reads with Flye"
# 	params:
# 		assembly_dir=dirs_dict["ASSEMBLY_DIR"] + "/flye_{sample}_{sampling}",
# 		genome_size=config["genome_size"],
# 		metagenomic_flag=METAGENOME_FLAG,
# 	conda:
# 		dirs_dict["ENVS_DIR"] + "/env1.yaml"
# 	benchmark:
# 		dirs_dict["BENCHMARKS"] +"/asemblyFlye/{sample}_{sampling}.tsv"
# 	threads: 4
# 	shell:
# 		"""
# 		flye --nano-raw {input.mapped_fastq_forward} --out-dir {params.assembly_dir} --genome-size {params.genome_size} --threads {threads} 
# 		cp {output.scaffolds} {output.scaffolds_final}
# 		"""

rule demultiplexing:
	input:
		basecalled_summary=dirs_dict["BASECALLED"] + "/sequencing_summary.txt",
		workspace_dir=dirs_dict["BASECALLED"] + "/workspace",
		basecalled_dir=dirs_dict["BASECALLED"] + "/{barcode}",
		mapped_list=dirs_dict["MAPPING"] + "/{genome}_{barcode}_loose_fast5_list_mapped.txt",
	output:
		demultiplexed_dir=directory(dirs_dict["DEMULTIPLEXED"] + "/{genome}_{barcode}"),
		length_list=dirs_dict["DEMULTIPLEXED"] + "/{genome}_{barcode}_fast5_list_length.txt",
		demultiplexed_list=dirs_dict["DEMULTIPLEXED"] + "/{genome}_{barcode}_fast5_list_pass.txt",
	message:
		"Demultiplexing single fast5 files"
	params:
		min_read_length=config['min_read_length']
	conda:
		"envs/env1.yaml"
	threads: 32
	shell:
		"""
		cat {input.basecalled_dir}/*fastq | sed -n '1~4s/^@/>/p;2~4p' | sed 's/\s.*$//' |
			awk '$0 ~ ">" {{print c; c=0;printf substr($0,2,100) "\t"; }} $0 !~ ">" {{c+=length($0);}} END {{ print c; }}' |
			awk '$2>{params.min_read_length}' | cut -f1 > {output.length_list}
		comm -12  <(sort {input.mapped_list}) <(sort {output.length_list}) > {output.demultiplexed_list}

		split -d -n l/40 {output.demultiplexed_list} {output.demultiplexed_list}_split
		parallel fast5_subset -i {input.workspace_dir} -s {output.demultiplexed_dir}_parallel{{}} -l {output.demultiplexed_list}_split{{}} -f split{{}} -n 100000 -t 4 ::: {{00..39}}
		
		mkdir {output.demultiplexed_dir}
		mv {output.demultiplexed_dir}_parallel*/*fast5 {output.demultiplexed_dir}
		rm -rf {output.demultiplexed_dir}_parallel* {output.demultiplexed_list}_split*
		"""

rule multi_to_single_fast5:
	input:
		demultiplexed_dir=(dirs_dict["DEMULTIPLEXED"] + "/{genome}_{barcode}"),
	output:
		single_data=directory(dirs_dict["SINGLE"]+ "/{genome}_{barcode}"),
	conda:
		"envs/env1.yaml"
	message:
		"Converting multi fast5 to single fast5"
	threads: 32
	shell:
		"""
		multi_to_single_fast5 --input_path {input.demultiplexed_dir} --save_path {output.single_data} -t {threads}
		"""

# rule annotate_tombo:
# 	input:
# 		basecalled_summary=dirs_dict["BASECALLED"] + "/sequencing_summary.txt",
# 		demultiplexed_dir=directory(dirs_dict["DEMULTIPLEXED"] + "/{barcode}"),
# 		basecalled_dir=directory(dirs_dict["BASECALLED"] + "/pass/{barcode}"),
# 	output:
# #		demultiplexed_dir=directory(expand((dirs_dict["DEMULTIPLEXED"] + "/{barcode}"), barcode=BARCODES)),
# 		annotated=(dirs_dict["BASECALLED"] + "/annotated_checkpoint_{barcode}.txt"),
# 	params:
# 		flowcell=FLOWCELL,
# 		kit=KIT,
# 	conda:
# 		"envs/env2.yaml"
# 	message:
# 		"Annotating fast5 files with fastq basecalls"
# 	threads: 16
# 	shell:
# 		"""
# 		tombo preprocess annotate_raw_with_fastqs --fast5-basedir {input.demultiplexed_dir} --fastq-filenames {input.basecalled_dir}/*fastq --sequencing-summary-filenames {input.basecalled_summary} --overwrite --processes {threads}
# 		touch {output.annotated}
# 		"""

rule resquiggle_tombo:
	input:
		#demultiplexed_dir=dirs_dict["DEMULTIPLEXED"] + "/{barcode}",
		single_data=(dirs_dict["SINGLE"]+ "/{genome}_{barcode}"),
		genome=GENOME_dir + "/{genome}.fasta",
	output:
		resquiggled=(dirs_dict["TOMBO"] + "/resquiggled_checkpoint_{genome}_{barcode}_default.txt"),
	threads: 16
	conda:
		"envs/env1.yaml"
	shell:
		"""
		tombo resquiggle --dna {input.single_data} {input.genome} --processes {threads} --overwrite --ignore-read-locks --corrected-group "default"
		touch {output.resquiggled}
		"""

rule resquiggle_tombo_loose:
	input:
		#demultiplexed_dir=dirs_dict["DEMULTIPLEXED"] + "/{barcode}",
		resquiggled=(dirs_dict["TOMBO"] + "/resquiggled_checkpoint_{genome}_{barcode}.txt"),
		single_data=(dirs_dict["SINGLE"]+ "/{genome}_{barcode}"),
		genome=GENOME_dir + "/{genome}.fasta",
	output:
		resquiggled=(dirs_dict["TOMBO"] + "/resquiggled_checkpoint_{genome}_{barcode}_loose.txt"),
	threads: 16
	conda:
		"envs/env1_loose.yaml"
	shell:
		"""
		tombo resquiggle --dna {input.single_data} {input.genome} --processes {threads} --overwrite --ignore-read-locks --corrected-group "loose"
		touch {output.resquiggled}
		"""

rule tombo_sample_compare:
	input:
		sample=(dirs_dict["SINGLE"] + "/{genome}_{sample}"),
		control=(dirs_dict["SINGLE"] + "/{genome}_{control}"),
		resquiggled=(dirs_dict["TOMBO"] + "/resquiggled_checkpoint_{genome}_{sample}_{mapping}.txt"),
		resquiggled2=(dirs_dict["TOMBO"] + "/resquiggled_checkpoint_{genome}_{control}_{mapping}.txt"),
		genome=GENOME_dir + "/{genome}.fasta",
	output:
		stats=dirs_dict["TOMBO"] + "/{genome}_{sample}_{control}_{mapping}.tombo_sampleCompare/{genome}_{sample}_{control}_{mapping}.tombo.stats" ,
		significant=dirs_dict["TOMBO"] + "/{genome}_{sample}_{control}_{mapping}.tombo_sampleCompare/{genome}_{sample}_{control}_{mapping}_tombo_results.significant_regions.fasta",
	params:
		name="{genome}_{sample}_{control}_{mapping}",
		tombo_results_dir=directory(dirs_dict["TOMBO"] + "/{genome}_{sample}_{control}_{mapping}.tombo_sampleCompare"),
		meme=dirs_dict["TOMBO"] + "/{genome}_{sample}_{control}_{mapping}.tombo_sampleCompare/{genome}_{sample}_{control}_{mapping}_tombo_sampleCompare_results.motif_detection.meme",
	conda:
		"envs/env1.yaml"
	message:
		"Detecting modified bases with Tombo sample compare"
	threads: 8
	shell:
		"""
		rm -r {params.tombo_results_dir} || true
		mkdir {params.tombo_results_dir}
		cd {params.tombo_results_dir}
		# tombo filter clear_filters --fast5-basedirs {input.sample}
		# tombo filter clear_filters --fast5-basedirs {input.control}
		tombo detect_modifications model_sample_compare --fast5-basedirs {input.sample} --control-fast5-basedirs {input.control} --statistics-file-basename {params.name} --per-read-statistics-basename {params.name} --processes {threads} --corrected-group {wildcards.mapping}
		tombo text_output browser_files --fast5-basedirs {input.sample} --control-fast5-basedirs {input.control} --statistics-filename {output.stats} --genome-fasta {input.genome} --browser-file-basename {params.name} --file-types coverage valid_coverage fraction dampened_fraction signal signal_sd
		tombo text_output signif_sequence_context --statistics-filename {output.stats} --genome-fasta {input.genome} --num-regions 100 --num-bases 10 --sequences-filename {output.significant}
		"""

rule tombo_denovo:
	input:
		sample=(dirs_dict["SINGLE"] + "/{genome}_{barcode}"),
		resquiggled=(dirs_dict["TOMBO"] + "/resquiggled_checkpoint_{genome}_{barcode}_{mapping}_.txt"),
		genome=GENOME_dir + "/{genome}.fasta",
	output:
		stats=dirs_dict["TOMBO"] + "/{genome}_{barcode}_{mapping}.tombo_denovo/{genome}_{barcode}_{mapping}_denovo.tombo.stats" ,
		significant=dirs_dict["TOMBO"] + "/{genome}_{barcode}_{mapping}.tombo_denovo/{genome}_{barcode}_{mapping}_tombo_denovo_results.significant_regions.fasta",
	params:
		name="{genome}_{barcode}_denovo",
		tombo_results_dir=directory(dirs_dict["TOMBO"] + "/{genome}_{barcode}_{mapping}.tombo_denovo"),
		meme=dirs_dict["TOMBO"] + "/{genome}_{barcode}_{mapping}.tombo_denovo/{genome}_{barcode}_{mapping}_tombo_denovo_results.motif_detection.meme",
	wildcard_constraints:
		control="barcode..",
		sample="barcode..",
	conda:
		"envs/env1.yaml"
	message:
		"Detecting modified bases with Tombo de novo"
	threads: 8
	shell:
		"""
		rm -r {params.tombo_results_dir} || true
		mkdir {params.tombo_results_dir}
		cd {params.tombo_results_dir}
		# tombo filter clear_filters --fast5-basedirs {input.sample}
		tombo detect_modifications de_novo --fast5-basedirs {input.sample} --statistics-file-basename {params.name} --per-read-statistics-basename {params.name} --processes {threads}
		tombo text_output browser_files --fast5-basedirs {input.sample} --statistics-filename {output.stats} --genome-fasta {input.genome} --browser-file-basename {params.name} --file-types coverage valid_coverage fraction dampened_fraction signal signal_sd
		tombo text_output signif_sequence_context --statistics-filename {output.stats} --genome-fasta {input.genome} --num-regions 100 --num-bases 10 --sequences-filename {output.significant}
		"""

rule tombo_alternative:
	input:
		sample=(dirs_dict["SINGLE"] + "/{genome}_{sample}"),
		resquiggled=(dirs_dict["TOMBO"] + "/resquiggled_checkpoint_{genome}_{sample}_{mapping}.txt"),
		genome=GENOME_dir + "/{genome}.fasta",
	output:
		stats=dirs_dict["TOMBO"] + "/{genome}_{sample}_{mapping}.tombo_alternative_{model}/{genome}_{sample}_{mapping}.tombo_alternative_{model}.{model}.tombo.stats" ,
		significant=dirs_dict["TOMBO"] + "/{genome}_{sample}_{mapping}.tombo_alternative_{model}/{genome}_{sample}_{mapping}_tombo_alternative_{model}_results.significant_regions.fasta",
	params:
		name="{genome}_{sample}_{mapping}.tombo_alternative_{model}",
		tombo_results_dir=(dirs_dict["TOMBO"] + "/{genome}_{sample}_{mapping}.tombo_alternative_{model}"),
		readstats="{genome}_{sample}_{mapping}.tombo_alternative_{model}_per_read" ,
	wildcard_constraints:
		control="barcode..",
		sample="barcode..",
	conda:
		"envs/env1.yaml"
	message:
		"Detecting modified bases with Tombo de novo"
	threads: 8
	shell:
		"""
		rm -r {params.tombo_results_dir} || true
		mkdir {params.tombo_results_dir}
		cd {params.tombo_results_dir}
		tombo detect_modifications alternative_model --alternate-bases {wildcards.model} --fast5-basedirs {input.sample} --statistics-file-basename {params.name} --per-read-statistics-basename {params.readstats} --processes {threads}
		tombo text_output browser_files --fast5-basedirs {input.sample} --statistics-filename {output.stats} --genome-fasta {input.genome} --browser-file-basename {params.name} --file-types coverage valid_coverage fraction dampened_fraction signal signal_sd
		tombo text_output signif_sequence_context --statistics-filename {output.stats} --genome-fasta {input.genome} --num-regions 10000 --num-bases 10 --sequences-filename {output.significant}
		"""

rule parse_tombo_results_sampleCompare:
	input:
		stats_sampleCompare=dirs_dict["TOMBO"] + "/{genome}_{sample}_{control}_{mapping}.tombo_sampleCompare/{genome}_{sample}_{control}_{mapping}.tombo.stats" ,
	output:
		modfrac_png= dirs_dict["PLOTS_DIR"] + "/{genome}/sampleCompare_{mapping}/sampleCompare_{genome}_{sample}_{control}_{mapping}_per_base_modfrac_10000.pdf",
		modfrac_kmers_table= dirs_dict["PLOTS_DIR"] + "/{genome}/sampleCompare_{mapping}/sampleCompare_{genome}_{sample}_{control}_{mapping}_kmer_modfrac.csv",
		coverage_png= dirs_dict["PLOTS_DIR"] + "/{genome}/sampleCompare_{mapping}/sampleCompare_{genome}_{sample}_{control}_{mapping}_per_base_coverage.pdf",
		dinucleotide= dirs_dict["PLOTS_DIR"] + "/{genome}/sampleCompare_{mapping}/sampleCompare_{genome}_{sample}_{control}_{mapping}_histogram_dinucleotide.pdf",
		trinucleotide= dirs_dict["PLOTS_DIR"] + "/{genome}/sampleCompare_{mapping}/sampleCompare_{genome}_{sample}_{control}_{mapping}_histogram_trinucleotide.pdf",
		tetranucleotide= dirs_dict["PLOTS_DIR"] + "/{genome}/sampleCompare_{mapping}/sampleCompare_{genome}_{sample}_{control}_{mapping}_histogram_tetranucleotide.pdf",
		pentanucleotide= dirs_dict["PLOTS_DIR"] + "/{genome}/sampleCompare_{mapping}/sampleCompare_{genome}_{sample}_{control}_{mapping}_histogram_pentanucleotide.pdf",
		# hexanucleotide= dirs_dict["PLOTS_DIR"] + "/{genome}/sampleCompare_{mapping}/sampleCompare_{genome}_{sample}_{control}_{mapping}_histogram_hexanucleotide.pdf",
	params:
		tombo_dir=dirs_dict["TOMBO"],
		figdir=dirs_dict["PLOTS_DIR"],
		genome="{genome}",
		sample="{sample}",
		control="{control}",
		threshold_modfrac=0.3,
		workdir=OUTPUT_DIR
	log:
		notebook=dirs_dict["NOTEBOOKS_DIR"] + "/04_TOMBO_parsing_sampleCompare_{genome}_{sample}_{control}_{mapping}.ipynb"
	notebook:
		dirs_dict["RAW_NOTEBOOKS"] + "/04_TOMBO_parsing_sampleCompare.py.ipynb"

rule parse_tombo_results_deNovo:
	input:
		stats_deNovo=dirs_dict["TOMBO"] + "/{genome}_{sample}_{mapping}.tombo_denovo/{genome}_{sample}_{mapping}_denovo.tombo.stats",
	output:
		modfrac_png= dirs_dict["PLOTS_DIR"] + "/{genome}/denovo_{mapping}/denovo_{genome}_{sample}_{mapping}_per_base_modfrac_10000.pdf",
		modfrac_kmers_table= dirs_dict["PLOTS_DIR"] + "/{genome}/denovo_{mapping}/denovo_{genome}_{sample}_{mapping}_kmer_modfrac.csv",
		coverage_png= dirs_dict["PLOTS_DIR"] + "/{genome}/denovo_{mapping}/denovo_{genome}_{sample}_{mapping}_per_base_coverage.pdf",
		dinucleotide= dirs_dict["PLOTS_DIR"] + "/{genome}/denovo_{mapping}/denovo_{genome}_{sample}_{mapping}_histogram_dinucleotide.pdf",
		trinucleotide= dirs_dict["PLOTS_DIR"] + "/{genome}/denovo_{mapping}/denovo_{genome}_{sample}_{mapping}_histogram_trinucleotide.pdf",
		tetranucleotide= dirs_dict["PLOTS_DIR"] + "/{genome}/denovo_{mapping}/denovo_{genome}_{sample}_{mapping}_histogram_tetranucleotide.pdf",
		pentanucleotide= dirs_dict["PLOTS_DIR"] + "/{genome}/denovo_{mapping}/denovo_{genome}_{sample}_{mapping}_histogram_pentanucleotide.pdf",
		# hexanucleotide= dirs_dict["PLOTS_DIR"] + "/{genome}/denovo_{mapping}/denovo_{genome}_{sample}_{mapping}_histogram_hexanucleotide.pdf",
	params:
		tombo_dir=dirs_dict["TOMBO"],
		figdir=dirs_dict["PLOTS_DIR"],
		genome="{genome}",
		sample="{sample}",
		threshold_modfrac=0.3,
		workdir=OUTPUT_DIR
	log:
		notebook=dirs_dict["NOTEBOOKS_DIR"] + "/04_TOMBO_parsing_deNovo_{genome}_{sample}_{mapping}.ipynb"
	notebook:
		dirs_dict["RAW_NOTEBOOKS"] + "/04_TOMBO_parsing_deNovo.py.ipynb"

# rule deepsignal:
# 	input:
#  		demultiplexed_dir=dirs_dict["DEMULTIPLEXED"] + "/{barcode}",
# 		genome="{genome}",
# 		resquiggled=(dirs_dict["TOMBO"] + "/resquiggled_checkpoint_{genome}_{barcode}.txt"),
# 	output:
# 		extract=dirs_dict["BASECALLED"] + "/{barcode}_deepsignal-feature.tsv"
# 	conda:
# 		"envs/env2.yaml"
# 	threads: 8
# 	shell:
# 		"""
# 		deepsignal extract --fast5_dir {input.demultiplexed_dir} --reference_path {input.genome} --is_dna true --write_path {output.extract} --nproc {threads}
# 		"""

# rule call_modification_deepsignal:
# 	input:
# 		extract=dirs_dict["BASECALLED"] + "/{barcode}_deepsignal-feature.tsv",
# 		model=(dirs_dict["TOOLS"]+ "/deepsignal/model.CpG.R9.4_1D.human_hx1.bn17.sn360.v0.1.7+"),
# 	output:
# 		dirs_dict["DEEPSIGNAL"] + "{barcode}_deepsignal-prob.tsv"
# 	conda:
# 		"envs/env2.yaml"
# 	shell:
# 		"""
# 		deepsignal call_mods --input_path {input.extract} --is_gpu no --nproc {threads} --model_path {input.model}/bn_17.sn_360.epoch_9.ckpt --result_file {output}
# 		"""

# rule megalodon:
# 	input:
# 		rerio_dir=(dirs_dict["TOOLS"]+ "/rerio"),
# 		demultiplexed_dir=dirs_dict["DEMULTIPLEXED"] + "/{barcode}",
# 		genome="{genome}",
# 	output:
# 		megalodon_dir=directory(dirs_dict["MEGALODON"] + "/{barcode}"),
# 	params:
# 		model_dir=dirs_dict["TOOLS"]+ "/rerio/basecall_models/",
# 		model_name="res_dna_r941_min_modbases_5mC_CpG_v001.cfg",
# 		server=config['guppy_server']
# 	threads: 32
# 	conda:
# 		"envs/env2.yaml"
# 	shell:
# 		"""
# 		megalodon {input.demultiplexed_dir} --guppy-params "-d {params.model_dir}" --guppy-config {params.model_name} \
# 			--outputs  mappings mod_mappings mods per_read_mods  --output-directory {output.megalodon_dir} \
# 			--reference {input.genome} --mod-motif m CG 0 --processes {threads} --guppy-server-path {params.server}
# 		"""
