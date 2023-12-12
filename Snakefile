import os
import re
#import pandas as pd
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

# GENOME=config['genome']
# GENOME_dir=os.path.dirname(os.path.realpath(GENOME))
# GENOME_name=os.path.splitext(os.path.basename(GENOME))[0]

GENOME_dir=config['genome_dir']
GENOME_name,=glob_wildcards(GENOME_dir + "/{GENOME}"+ ".fasta")


SAMPLES=config['samples'].split()
CONTROL=config['control'].split()
BARCODES=SAMPLES+CONTROL
BARCODES = list( dict.fromkeys(BARCODES) )
BARCODES = list(filter(None, BARCODES))

dir_list = ["RULES_DIR","ENVS_DIR","DB", "TOOLS", "SINGLE", "BASECALLED", "DEMULTIPLEXED", "GENOMES", "TOMBO", "MEGALODON", "DEEPSIGNAL", "QC"]
dir_names = ["rules", "../envs", OUTPUT_DIR + "/db", OUTPUT_DIR + "/tools", OUTPUT_DIR + "/03_FAST5_SINGLE", OUTPUT_DIR + "/01_BASECALLED", OUTPUT_DIR + "/02_DEMULTIPLEXED", GENOME_dir , OUTPUT_DIR + "/04_TOMBO", OUTPUT_DIR + "/05_MEGALODON", OUTPUT_DIR + "/06_DEEPSIGNAL", OUTPUT_DIR + "/STATS"]
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

rule all:
	input:
		#plus_corr=expand(dirs_dict["TOMBO"] + "/" + GENOME + "/" + GENOME + "_{sample}_" + CONTROL + "_plusmod_corrected.wig", sample=SAMPLES),
		#genome_oneline=dirs_dict["GENOMES"] + "/" + GENOME + "_one.fasta",
		#stats=expand(dirs_dict["TOMBO"] + "/" + GENOME + "/" + GENOME + "_{sample}_" + CONTROL + ".tombo.stats", sample=SAMPLES) ,
		#directory((dirs_dict["BASECALLED"])),
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
		expand(dirs_dict["QC"] + "/{barcode}_{genome}_nanoQC", barcode=BARCODES, genome=GENOME_name),

rule demultiplex_run:
	input:
		expand(dirs_dict["DEMULTIPLEXED"] + "/{barcode}_checkpoint.txt", barcode=BARCODES),

rule megalodon_run:
	input:
		expand(dirs_dict["MEGALODON"] + "/{barcode}", barcode=BARCODES),

rule deepsignal_run:
	input:
		expand(dirs_dict["DEEPSIGNAL"] + "/{barcode}_deepsignal-prob.tsv", barcode=BARCODES),

rule tombo_run_denovo:
	input:
		#expand(dirs_dict["TOMBO"] + "/{genome}_{barcode}.tombo_denovo", barcode=SAMPLES, genome=GENOME_name),
		#expand(dirs_dict["TOMBO"] + "/{genome}_{barcode}.tombo_denovo/{genome}_{barcode}_tombo_denovo_results.motif_detection.meme.html", barcode=SAMPLES, genome=GENOME_name),
		expand(dirs_dict["TOMBO"] + "/{genome}_{barcode}.tombo_denovo/{genome}_{barcode}_denovo.tombo.stats", barcode=BARCODES, genome=GENOME_name),
		expand(dirs_dict["QC"] + "/{barcode}_{genome}_nanoQC", barcode=BARCODES, genome=GENOME_name),

rule tombo_run_sampleCompare:
	input:
		#expand(dirs_dict["TOMBO"] + "/"+ GENOME_name + "_{sample}.tombo_denovo.stats", sample=SAMPLES),
		#expand(dirs_dict["TOMBO"] + "/{genome}_{sample}_{control}.tombo.stats", sample=SAMPLES, control=CONTROL, genome=GENOME_name),
		#expand(dirs_dict["TOMBO"] + "/{genome}_{sample}_{control}.tombo_sampleCompare/{genome}_{sample}_{control}_tombo_sampleCompare_results.motif_detection.meme.html", sample=SAMPLES, control=CONTROL, genome=GENOME_name),
		expand(dirs_dict["TOMBO"] + "/{genome}_{sample}_{control}.tombo_sampleCompare/{genome}_{sample}_{control}.tombo.stats", sample=SAMPLES, control=CONTROL, genome=GENOME_name),
		expand(dirs_dict["QC"] + "/{barcode}_{genome}_nanoQC", barcode=BARCODES, genome=GENOME_name),

rule tombo_run_alternative:
	input:
		expand(dirs_dict["TOMBO"] + "/{genome}_{sample}.tombo_alternative_{model}/{genome}_{sample}.tombo_alternative_{model}.{model}.tombo.stats", sample=SAMPLES, model=ALTERNATIVE_MODELS, genome=GENOME_name),
		#expand(dirs_dict["TOMBO"] + "/"+ GENOME_name + "_{sample}_{control}.tombo.stats", sample=SAMPLES, control=CONTROL),


# wildcard_constraints:
# 	barcode="barcode..",


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

# else:

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
		guppy_basecaller -i {input.raw_data} -s {params.basecalled_dir} -q 0 -r -x 'cuda:0 cuda:1' -c /opt/ont/guppy/data/{params.model} --barcode_kits {params.kit} --post_out --disable_qscore_filtering --chunks_per_runner 1 --gpu_runners_per_device 1 --num_callers 1
		# --chunks_per_runner 1 --gpu_runners_per_device 1 --num_callers 1
		# --chunks_per_runner 128
		"""


rule map_to_genomes:
	input:
		basecalled_dir=dirs_dict["BASECALLED"] + "/{barcode}",
		genome=GENOME_dir + "/{genome}.fasta",
	output:
		mapped_paf=dirs_dict["BASECALLED"] + "/{barcode}_vs_{genome}.paf",
		mapped_list=dirs_dict["BASECALLED"] + "/{barcode}_vs_{genome}_fast5_list_mapped.txt",
		merged_fastq=temp(dirs_dict["BASECALLED"] + "/{barcode}_vs_{genome}.fastq"),
		mapped_fastq=(dirs_dict["BASECALLED"] + "/{barcode}_vs_{genome}_mapped.fastq"),
	params:
		flowcell=FLOWCELL,
		kit=KIT,
	conda:
		"envs/env1.yaml"
	message:
		"Mapping reads to genomes with minimap2"
	threads: 32
	shell:
		"""
		cat {input.basecalled_dir}/*fastq > {output.merged_fastq}
		# /home/demeter/Storage/lmf/apps/minimap2/minimap2 -ax map-ont {input.genome} {output.merged_fastq} > {output.mapped_paf}
		minimap2 -ax map-ont {input.genome} {output.merged_fastq} > {output.mapped_paf}
		cat {output.mapped_paf} | cut -f1 > {output.mapped_list}
		seqtk subseq {output.merged_fastq} {output.mapped_list} > {output.mapped_fastq}
		"""

rule qualityCheckNanopore:
	input:
		mapped_fastq=(dirs_dict["BASECALLED"] + "/{barcode}_vs_{genome}_mapped.fastq"),
	output:
		nanoqc_dir=directory(dirs_dict["QC"] + "/{barcode}_{genome}_nanoQC"),
	message:
		"Performing nanoQC statistics"
	conda:
		"envs/env3.yaml"
	shell:
		"""
		nanoQC -o {output.nanoqc_dir} {input.mapped_fastq}
		"""

rule demultiplexing:
	input:
		basecalled_summary=dirs_dict["BASECALLED"] + "/sequencing_summary.txt",
		workspace_dir=dirs_dict["BASECALLED"] + "/workspace",
		basecalled_dir=dirs_dict["BASECALLED"] + "/{barcode}",
		mapped_list=dirs_dict["BASECALLED"] + "/{barcode}_vs_{genome}_fast5_list_mapped.txt",
		# annotated=(dirs_dict["BASECALLED"] + "/annotated_checkpoint_{barcode}.txt"),
	output:
		demultiplexed_dir=directory(dirs_dict["DEMULTIPLEXED"] + "/{barcode}_{genome}"),
		length_list=dirs_dict["DEMULTIPLEXED"] + "/{barcode}_{genome}_fast5_list_length.txt",
		demultiplexed_list=dirs_dict["DEMULTIPLEXED"] + "/{barcode}_{genome}_fast5_list_pass.txt",
#		checkpoint=dirs_dict["DEMULTIPLEXED"] + "/{barcode}_checkpoint.txt",
	message:
		"Demultiplexing single fast5 files"
	params:
		min_read_length=config['min_read_length']
	conda:
		"envs/env1.yaml"
	threads: 1
	shell:
		"""
		cat {input.basecalled_dir}/*fastq | sed -n '1~4s/^@/>/p;2~4p' | sed 's/\s.*$//' |
			awk '$0 ~ ">" {{print c; c=0;printf substr($0,2,100) "\t"; }} $0 !~ ">" {{c+=length($0);}} END {{ print c; }}' |
			awk '$2>{params.min_read_length}' | cut -f1 > {output.length_list}
		comm -12  <(sort {input.mapped_list}) <(sort {output.length_list}) > {output.demultiplexed_list}
		fast5_subset -i {input.workspace_dir} -s {output.demultiplexed_dir} -l {output.demultiplexed_list} -n 1000000000
		"""

rule multi_to_single_fast5:
	input:
		demultiplexed_dir=(dirs_dict["DEMULTIPLEXED"] + "/{barcode}_{genome}"),
	output:
		single_data=directory(dirs_dict["SINGLE"]+ "/{barcode}_{genome}"),
	conda:
		"envs/env1.yaml"
	message:
		"Converting multi fast5 to single fast5"
	threads: 16
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
		single_data=(dirs_dict["SINGLE"]+ "/{barcode}_{genome}"),
		genome=GENOME_dir + "/{genome}.fasta",
	output:
		resquiggled=(dirs_dict["TOMBO"] + "/resquiggled_checkpoint_{barcode}_{genome}.txt"),
	threads: 16
	conda:
		"envs/env1.yaml"
	shell:
		"""
		tombo resquiggle --dna {input.single_data} {input.genome} --processes {threads} --overwrite --ignore-read-locks
		touch {output.resquiggled}
		"""


rule tombo_sample_compare:
	input:
		sample=(dirs_dict["SINGLE"] + "/{sample}_{genome}"),
		control=(dirs_dict["SINGLE"] + "/{control}_{genome}"),
		resquiggled=(dirs_dict["TOMBO"] + "/resquiggled_checkpoint_{sample}_{genome}.txt"),
		resquiggled2=(dirs_dict["TOMBO"] + "/resquiggled_checkpoint_{control}_{genome}.txt"),
		genome=GENOME_dir + "/{genome}.fasta",
	output:
		stats=dirs_dict["TOMBO"] + "/{genome}_{sample}_{control}.tombo_sampleCompare/{genome}_{sample}_{control}.tombo.stats" ,
		significant=dirs_dict["TOMBO"] + "/{genome}_{sample}_{control}.tombo_sampleCompare/{genome}_{sample}_{control}_tombo_results.significant_regions.fasta",
	params:
		name="{genome}_{sample}_{control}",
		tombo_results_dir=directory(dirs_dict["TOMBO"] + "/{genome}_{sample}_{control}.tombo_sampleCompare"),
		meme=dirs_dict["TOMBO"] + "/{genome}_{sample}_{control}.tombo_sampleCompare/{genome}_{sample}_{control}_tombo_sampleCompare_results.motif_detection.meme",
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
		tombo filter clear_filters --fast5-basedirs {input.sample}
		tombo filter clear_filters --fast5-basedirs {input.control}
		tombo detect_modifications model_sample_compare --fast5-basedirs {input.sample} --control-fast5-basedirs {input.control} --statistics-file-basename {params.name} --per-read-statistics-basename {params.name} --processes {threads}
		tombo text_output browser_files --fast5-basedirs {input.sample} --control-fast5-basedirs {input.control} --statistics-filename {output.stats} --genome-fasta {input.genome} --browser-file-basename {params.name} --file-types coverage valid_coverage fraction dampened_fraction signal signal_sd
		tombo text_output signif_sequence_context --statistics-filename {output.stats} --genome-fasta {input.genome} --num-regions 100 --num-bases 10 --sequences-filename {output.significant}
		"""

rule tombo_denovo:
	input:
		sample=(dirs_dict["SINGLE"] + "/{barcode}_{genome}"),
		resquiggled=(dirs_dict["TOMBO"] + "/resquiggled_checkpoint_{barcode}_{genome}.txt"),
		genome=GENOME_dir + "/{genome}.fasta",
	output:
		stats=dirs_dict["TOMBO"] + "/{genome}_{barcode}.tombo_denovo/{genome}_{barcode}_denovo.tombo.stats" ,
		significant=dirs_dict["TOMBO"] + "/{genome}_{barcode}.tombo_denovo/{genome}_{barcode}_tombo_denovo_results.significant_regions.fasta",
	params:
		name="{genome}_{barcode}_denovo",
		tombo_results_dir=directory(dirs_dict["TOMBO"] + "/{genome}_{barcode}.tombo_denovo"),
		meme=dirs_dict["TOMBO"] + "/{genome}_{barcode}.tombo_denovo/{genome}_{barcode}_tombo_denovo_results.motif_detection.meme",
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
		tombo filter clear_filters --fast5-basedirs {input.sample}
		tombo detect_modifications de_novo --fast5-basedirs {input.sample} --statistics-file-basename {params.name} --per-read-statistics-basename {params.name} --processes {threads}
		tombo text_output browser_files --fast5-basedirs {input.sample} --statistics-filename {output.stats} --genome-fasta {input.genome} --browser-file-basename {params.name} --file-types coverage valid_coverage fraction dampened_fraction signal signal_sd
		tombo text_output signif_sequence_context --statistics-filename {output.stats} --genome-fasta {input.genome} --num-regions 100 --num-bases 10 --sequences-filename {output.significant}
		"""

rule tombo_alternative:
	input:
		sample=(dirs_dict["SINGLE"] + "/{sample}_{genome}"),
		resquiggled=(dirs_dict["TOMBO"] + "/resquiggled_checkpoint_{sample}_{genome}.txt"),
		genome=GENOME_dir + "/{genome}.fasta",
	output:
		stats=dirs_dict["TOMBO"] + "/{genome}_{sample}.tombo_alternative_{model}/{genome}_{sample}.tombo_alternative_{model}.{model}.tombo.stats" ,
		#minus=dirs_dict["TOMBO"] + "/{genome}_{sample}.tombo_alternative/{genome}_{sample}_alternative_{model}_minusmod.wig",
		#plus=dirs_dict["TOMBO"] + "/{genome}_{sample}.tombo_alternative/{genome}_{sample}_alternative_{model}_plusmod.wig",
		significant=dirs_dict["TOMBO"] + "/{genome}_{sample}.tombo_alternative_{model}/{genome}_{sample}_tombo_alternative_{model}_results.significant_regions.fasta",
	params:
		name="{genome}_{sample}.tombo_alternative_{model}",
		tombo_results_dir=(dirs_dict["TOMBO"] + "/{genome}_{sample}.tombo_alternative_{model}"),
		readstats="{genome}_{sample}.tombo_alternative_{model}_per_read" ,
	# wildcard_constraints:
	# 	control="barcode..",
	# 	sample="barcode..",
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
rule deepsignal:
	input:
 		demultiplexed_dir=dirs_dict["DEMULTIPLEXED"] + "/{barcode}",
		genome="{genome}",
		resquiggled=(dirs_dict["TOMBO"] + "/resquiggled_checkpoint_{barcode}_{genome}.txt"),
	output:
		extract=dirs_dict["BASECALLED"] + "/{barcode}_deepsignal-feature.tsv"
	conda:
		"envs/env2.yaml"
	threads: 8
	shell:
		"""
		deepsignal extract --fast5_dir {input.demultiplexed_dir} --reference_path {input.genome} --is_dna true --write_path {output.extract} --nproc {threads}
		"""

rule call_modification_deepsignal:
	input:
		extract=dirs_dict["BASECALLED"] + "/{barcode}_deepsignal-feature.tsv",
		model=(dirs_dict["TOOLS"]+ "/deepsignal/model.CpG.R9.4_1D.human_hx1.bn17.sn360.v0.1.7+"),
	output:
		dirs_dict["DEEPSIGNAL"] + "{barcode}_deepsignal-prob.tsv"
	conda:
		"envs/env2.yaml"
	shell:
		"""
		deepsignal call_mods --input_path {input.extract} --is_gpu no --nproc {threads} --model_path {input.model}/bn_17.sn_360.epoch_9.ckpt --result_file {output}
		"""

rule megalodon:
	input:
		rerio_dir=(dirs_dict["TOOLS"]+ "/rerio"),
		demultiplexed_dir=dirs_dict["DEMULTIPLEXED"] + "/{barcode}",
		genome="{genome}",
	output:
		megalodon_dir=directory(dirs_dict["MEGALODON"] + "/{barcode}"),
	params:
		model_dir=dirs_dict["TOOLS"]+ "/rerio/basecall_models/",
		model_name="res_dna_r941_min_modbases_5mC_CpG_v001.cfg",
		server=config['guppy_server']
	threads: 32
	conda:
		"envs/env2.yaml"
	shell:
		"""
		megalodon {input.demultiplexed_dir} --guppy-params "-d {params.model_dir}" --guppy-config {params.model_name} \
			--outputs  mappings mod_mappings mods per_read_mods  --output-directory {output.megalodon_dir} \
			--reference {input.genome} --mod-motif m CG 0 --processes {threads} --guppy-server-path {params.server}
		"""
