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
RAW_DATA_DIR =OUTPUT_DIR + "/" + config['input_dir']
BARCODES = config["barcodes"].split()
FLOWCELL=config['flowcell']
KIT=config['kit']

GENOME=config['genome']
SAMPLE=config['sample']
CONTROL=config['control']

dir_list = ["RULES_DIR","ENVS_DIR","DB", "TOOLS", "SINGLE_DATA_DIR", "GUPPY", "DEMULTIPLEXED", "BASECALLED", "GENOMES"]
dir_names = ["rules", "../envs", OUTPUT_DIR + "/db", OUTPUT_DIR + "/tools", OUTPUT_DIR + "/01_GUPPY", OUTPUT_DIR + "/01_GUPPY/02_DEMULTIPLEXED", OUTPUT_DIR + "/01_GUPPY/02_BASECALLED", OUTPUT_DIR + "/01_GUPPY/02_FAST5_SINGLE", OUTPUT_DIR + "/GENOMES"]
dirs_dict = dict(zip(dir_list, dir_names))

#SAMPLES,=glob_wildcards(RAW_DATA_DIR + "/{{input.sample}}_" +".fast5")

#======================================================
# Rules
#======================================================

rule all:
	input:
		plus_corr=dirs_dict["TOMBO"] + "/" + GENOME + "/" + GENOME + "_plusmod_corrected.wig",
		genome_oneline=dirs_dict["GENOMES"] + GENOME + "_one.fasta",
		stats=dirs_dict["TOMBO"] + "/" + GENOME + "/" + GENOME + "_" + SAMPLE + "_" + CONTROL + ".tombo.stats" ,
rule guppy_demultiplexing basecalling:
	input:
		RAW_DATA_DIR,
	output:
		demultiplexed_dir=expand(directory(dirs_dict["DEMULTIPLEXED"] + "/{barcode}"), barcode=BARCODES),
		basecalled_dir=expand(directory(dirs_dict["BASECALLED"] + "/{barcode}"), barcode=BARCODES),
	params:
		guppy_dir=dirs_dict["GUPPY"],
		flowcell=FLOWCELL,
		kit=KIT,
	conda:
		"envs/env1.yaml"
	message:
		"Basecalling single fast5 files with guppy"
	threads: 8
	shell:
		"""
		guppy_basecaller -i {input} -s {output.basecalled_dir} --fast5_out -q 0 -r --trim_barcodes -x 'cuda:0 cuda:1' --flowcell FLO-MIN106 --kit SQK-RBK004
		"""

rule multi_to_single_fast5:
	input:
		basecalled_dir=expand(directory(dirs_dict["BASECALLED"] + "/{barcode}"), barcode=BARCODES),
	output:
		single_data=expand(directory(dirs_dict["BASECALLED"] + "/{barcode}_single"), barcode=BARCODES),
	conda:
		"envs/env1.yaml"
	message:
		"Converting multi fast5 to single fast5"
	threads: 16
	shell:
		"""
		multi_to_single_fast5 --input_path {input} --save_path {output} -t {threads}
		"""
rule tombo:
	input:
		sample=((dirs_dict["BASECALLED"] + "/{sample}_single")),
		control=((dirs_dict["BASECALLED"] + "/{control}_single")),
		genome=dirs_dict["GENOMES"] + "{genome}.fasta",
	output:
		stats=dirs_dict["TOMBO"] + "/{genome}/{genome}_{sample}_{control}.tombo.stats" ,
		significant=dirs_dict["TOMBO"] + "/{genome}/tombo_results.significant_regions.fasta",
		significant_filtered=dirs_dict["TOMBO"] + "/{genome}/{genome}.sig_filtered.fasta",
		minus=dirs_dict["TOMBO"] + "/{genome}/{genome}_minusmod.wig",
		plus=dirs_dict["TOMBO"] + "/{genome}/{genome}_plusmod.wig",
		minus_corr=dirs_dict["TOMBO"] + "/{genome}/{genome}_minusmod_corrected.wig",
		plus_corr=dirs_dict["TOMBO"] + "/{genome}/{genome}_plusmod_corrected.wig",
		genome_oneline=dirs_dict["GENOMES"] + "{genome}_one.fasta",
	params:
		name="{genome}_{sample}_{control}"
	conda:
		"envs/On-rep-seq.yaml"
	message:
		"Demultiplexing step 1"
	threads: 16
	shell:
		"""
		tombo preprocess annotate_raw_with_fastqs --fast5-basedir {input.sample} --fastq-filenames {wildcards.sample}.fastq --overwrite
		tombo preprocess annotate_raw_with_fastqs --fast5-basedir {input.control} --fastq-filenames {wildcards.control}.fastq --overwrite
		tombo resquiggle --processes 40 --overwrite {input.control} {input.genome}
		tombo resquiggle --processes 40 --overwrite {input.sample} {input.genome}
		tombo detect_modifications model_sample_compare --fast5-basedirs {input.sample} --control-fast5-basedirs {input.control} --statistics-file-basename {params.name}
		tombo text_output browser_files --fast5-basedirs {input.sample} --statistics-filename {output.stats} --genome-fasta {input.genome} --browser-file-basename {params.name} --file-types fraction
		tombo text_output signif_sequence_context --statistics-filename {output.stats} --genome-fasta {input.genome} --num-regions 10000 --num-bases 10
		fasta_formatter -i tombo_results.significant_regions.fasta -t | awk '{if ($5 > 0.7) print ">"$1,$2,$3,$4,$5"\n"$6}' > {output.significant_filtered}
		./scripts/format_tombo.py ${name}_plusmod.wig
		./scripts/format_tombo.py ${name}_minusmod.wig
		grep -v ">" {input.genome} | sed 's/./\0\n/g' | sed '/^$/d' > {output.genome_oneline}
		"""

rule cosito:
	input:
		RAW_DATA_DIR
	output:
		temp(expand(OUTPUT_DIR + "/01_porechopped_data/{barcode}.fastq", barcode=BARCODES))
	params:
		output_dir=OUTPUT_DIR + "/01_porechopped_data"
	conda:
		"envs/On-rep-seq.yaml"
	message:
		"Demultiplexing"
	threads: 16
	shell:
		"""
		head -n 25 scripts/logo.txt
		counter=1
		n=$(ls -l {input}/*fastq | wc -l )
		rm -f {params.output_dir}/*fastq
		for filename in {input}/*fastq
		do
			echo "Processing {input.sample} $counter/$n"
			porechop -i $filename -b dir_$filename -t {threads} --discard_unassigned --verbosity 0 > /dev/null 2>&1
			for bar in dir_$filename/*.fastq
			do
				f=$(basename -- $bar)
				cat $bar >> {params.output_dir}/$f
			done
			rm -rf dir_$filename
			counter=$((counter+1))
		done
		line=$(echo {BARCODES})
		for barcode in $line
		do
			touch {params.output_dir}/$barcode.fastq
		done
		"""

rule demultiplexing_Deepbinner:
	input:
		directory(dirs_dict["SINGLE_DATA_DIR"]),
	output:
		demultiplexed_dir=expand(directory(dirs_dict["DEMULTIPLEXED"] + "/{barcode}"), barcode=BARCODES),
		rapid_model=dirs_dict["TOOLS"]+ "/Deepbinner/models/SQK-RBK004_read_starts",
	params:
		tools_dir=dirs_dict["TOOLS"]
	conda:
		"envs/env1.yaml"
	message:
		"Demultiplexing fast5 files with Deepbinner"
	shell:
		"""
		git clone https://github.com/rrwick/Deepbinner.git
		mv Deepbinner {params.tools_dir}
		./tools/Deepbinner/deepbinner-runner.py realtime --in_dir {input} --out_dir {output.demultiplexed_dir} -s {output.rapid_model} --stop
		"""

#include: os.path.join(RULES_DIR, 'resultsParsing.smk')
