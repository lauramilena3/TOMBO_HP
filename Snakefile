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

GENOME=config['genome']
SAMPLES=config['samples'].split()
CONTROL=config['control']
BARCODES=SAMPLES+[CONTROL]
dir_list = ["RULES_DIR","ENVS_DIR","DB", "TOOLS", "GUPPY", "DEMULTIPLEXED", "BASECALLED", "SINGLE", "GENOMES", "TOMBO"]
dir_names = ["rules", "../envs", OUTPUT_DIR + "/db", OUTPUT_DIR + "/tools", OUTPUT_DIR + "/01_GUPPY", OUTPUT_DIR + "/01_GUPPY/02_DEMULTIPLEXED", OUTPUT_DIR + "/01_GUPPY/01_BASECALLED", OUTPUT_DIR + "/01_GUPPY/00_FAST5_SINGLE", OUTPUT_DIR + "/GENOMES", OUTPUT_DIR + "/02_TOMBO"]
dirs_dict = dict(zip(dir_list, dir_names))

#SAMPLES,=glob_wildcards(RAW_DATA_DIR + "/{{input.sample}}_" +".fast5")
print("outdir", OUTPUT_DIR)
print("samples", SAMPLES)
print("control", CONTROL)
print("barcodes", BARCODES)
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
		directory((dirs_dict["DEMULTIPLEXED"])),


rule multi_to_single_fast5:
	input:
		raw_data=RAW_DATA_DIR,
	output:
		single_data=dirs_dict["SINGLE"]
	conda:
		"envs/env1.yaml"
	message:
		"Converting multi fast5 to single fast5"
	threads: 16
	shell:
		"""
		multi_to_single_fast5 --input_path {input} --save_path {output} -t {threads}
		"""

rule guppy_basecalling:
	input:
		single_data=dirs_dict["SINGLE"]
	output:
#		demultiplexed_dir=directory(expand((dirs_dict["DEMULTIPLEXED"] + "/{barcode}"), barcode=BARCODES)),
#		basecalled_dir=directory(expand((dirs_dict["BASECALLED"] + "/{barcode}"), barcode=BARCODES)),
		directory((dirs_dict["BASECALLED"])),
	params:
		guppy_dir=dirs_dict["GUPPY"],
		flowcell=FLOWCELL,
		kit=KIT,
	conda:
		"envs/env1.yaml"
	message:
		"Basecalling single fast5 files with guppy"
	threads: 32
	shell:
		"""
		guppy_basecaller -i {input.single_data} -s {output.basecalled_dir} --fast5_out -q 0 -r --trim_barcodes -x 'cuda:0 cuda:1' --flowcell {params.flowcell} --kit {params.kit} -t {threads}
		"""
rule guppy_demultiplexing:
	input:
		basecalled_dir=directory((dirs_dict["BASECALLED"])),
	output:
#		demultiplexed_dir=directory(expand((dirs_dict["DEMULTIPLEXED"] + "/{barcode}"), barcode=BARCODES)),
#		basecalled_dir=directory(expand((dirs_dict["BASECALLED"] + "/{barcode}"), barcode=BARCODES)),
		demultiplexed_dir=directory((dirs_dict["DEMULTIPLEXED"])),
	params:
		guppy_dir=dirs_dict["GUPPY"],
		flowcell=FLOWCELL,
		kit=KIT,
	conda:
		"envs/env1.yaml"
	message:
		"Demultiplexing single fast5 files with guppy"
	threads: 32
	shell:
		"""
		guppy_barcoder -i {input.basecalled_dir} -s {input.demultiplexed_dir}  -t {threads}
		"""

rule move_fast5_files:
	input:
		raw_data=RAW_DATA_DIR,
	output:
#		demultiplexed_dir=directory(expand((dirs_dict["DEMULTIPLEXED"] + "/{barcode}"), barcode=BARCODES)),
#		basecalled_dir=directory(expand((dirs_dict["BASECALLED"] + "/{barcode}"), barcode=BARCODES)),
		#fastq=dirs_dict["BASECALLED"] + "/{barcode}/fastq/{barcode}.fastq",
		fast5=dirs_dict["GUPPY"] + "/{barcode}/fast5/{barcode}.fast5",
	params:
		dir_fastq=dirs_dict["BASECALLED"]+"/pass/",
		#barcode_number="{barcode}".split("barcode")[0],
	message:
		"Creating folders for fast5 files"
	threads: 1
	shell:
		"""
		nbarcode={wildcards.barcode}
		number=(${{nbarcode//barcode/ }})
		echo "barcode" $number
		cp {input.raw_data}/*_*_${{number}}.fast5 {output.fast5}
		"""

rule guppy_demultiplexing_basecalling:
	input:
		fast5=dirs_dict["GUPPY"] + "/{barcode}/fast5/{barcode}.fast5",
	output:
#		demultiplexed_dir=directory(expand((dirs_dict["DEMULTIPLEXED"] + "/{barcode}"), barcode=BARCODES)),
#		basecalled_dir=directory(expand((dirs_dict["BASECALLED"] + "/{barcode}"), barcode=BARCODES)),
		fastq=dirs_dict["GUPPY"] + "/{barcode}/fastq/{barcode}.fastq",
		basecalled_dir=dirs_dict["GUPPY"] + "/{barcode}/guppy/",
	conda:
		"envs/env1.yaml"
	params:
		fastq_dir=dirs_dict["GUPPY"] + "/{barcode}/guppy/pass",
		fast5_dir=dirs_dict["GUPPY"] + "/{barcode}/fast5",
		flowcell=FLOWCELL,
		kit=KIT,
	message:
		"Basecalling single fast5 files with guppy"
	threads: 8
	shell:
		"""
		nbarcode={wildcards.barcode}
		number=(${{nbarcode//barcode/ }})
		echo "barcode" $number
		guppy_basecaller -i {params.fast5_dir} -s {output.basecalled_dir} --fast5_out -q 0 -r --trim_barcodes -x 'cuda:0 cuda:1' --flowcell {params.flowcell} --kit {params.kit}
		cp {params.fastq_dir}/fastq_runid_*${{number}}_0.fastq {output.fastq}
		"""




rule tombo:
	input:
		sample=((dirs_dict["SINGLE"] + "/{sample}_single/workspace")),
		control=((dirs_dict["SINGLE"] + "/{control}_single/workspace")),
		basecalled_sample=dirs_dict["BASECALLED"] + "/{sample}.fastq",
		basecalled_control=dirs_dict["BASECALLED"] + "/{control}.fastq",
		genome=dirs_dict["GENOMES"] + "/{genome}.fasta",
	output:
		stats=dirs_dict["TOMBO"] + "/{genome}/{genome}_{sample}_{control}.tombo.stats" ,
		significant_filtered=dirs_dict["TOMBO"] + "/{genome}/{genome}_{sample}_{control}.sig_filtered.fasta",
		minus=dirs_dict["TOMBO"] + "/{genome}/{genome}_{sample}_{control}_minusmod.wig",
		plus=dirs_dict["TOMBO"] + "/{genome}/{genome}_{sample}_{control}_plusmod.wig",
		minus_corr=dirs_dict["TOMBO"] + "/{genome}/{genome}_{sample}_{control}_minusmod_corrected.wig",
		plus_corr=dirs_dict["TOMBO"] + "/{genome}/{genome}_{sample}_{control}_plusmod_corrected.wig",
		significant=dirs_dict["TOMBO"] + "/{genome}/{genome}_{sample}_{control}_tombo_results.significant_regions.fasta",
	params:
		name="{genome}_{sample}_{control}",
	conda:
		"envs/env1.yaml"
	message:
		"Detecting modified bases with Tombo"
	threads: 16
	shell:
		"""
		#tombo preprocess annotate_raw_with_fastqs --fast5-basedir {input.sample} --fastq-filenames {input.basecalled_sample} --overwrite
		#tombo preprocess annotate_raw_with_fastqs --fast5-basedir {input.control} --fastq-filenames {input.basecalled_control} --overwrite
		#tombo resquiggle --overwrite {input.control} {input.genome} --processes {threads} --ignore-read-locks
		#tombo resquiggle --overwrite {input.sample} {input.genome} --processes {threads} --ignore-read-locks
		#tombo resquiggle --corrected-group {wildcards.sample} --overwrite {input.sample} {input.genome} --processes {threads} --ignore-read-locks
		#tombo resquiggle --corrected-group {wildcards.control} --overwrite {input.control} {input.genome} --processes {threads} --ignore-read-locks
		tombo detect_modifications model_sample_compare --fast5-basedirs {input.sample} --control-fast5-basedirs {input.control} --statistics-file-basename {params.name}
		mv {params.name}.tombo.stats {output.stats}
		tombo text_output browser_files --fast5-basedirs {input.sample} --statistics-filename {output.stats} --genome-fasta {input.genome} --browser-file-basename {params.name} --file-types fraction
		mv {params.name}.fraction_modified_reads.plus.wig {output.plus}
		mv {params.name}.fraction_modified_reads.minus.wig {output.minus}
		tombo text_output signif_sequence_context --statistics-filename {output.stats} --genome-fasta {input.genome} --num-regions 10000 --num-bases 10
		mv tombo_results.significant_regions.fasta {output.significant}
		fasta_formatter -i {output.significant} -t | awk '{{if ($5 > 0.7) print ">"$1,$2,$3,$4,$5"\\n"$6}}' > {output.significant_filtered}
		python ./scripts/format_tombo.py {output.plus} {output.plus_corr}
		python ./scripts/format_tombo.py {output.minus} {output.minus_corr}
		"""
rule reformat_genome:
	input:
		genome=dirs_dict["GENOMES"] + "/{genome}.fasta",
	output:
		genome_oneline=dirs_dict["GENOMES"] + "/{genome}_one.fasta",
	shell:
		"""
		grep -v ">" {input.genome} | sed 's/./\\0\\n/g' | sed '/^$/d' > {output.genome_oneline}
		"""

rule cosito:
	input:
		RAW_DATA_DIR
	output:
		temp(expand(OUTPUT_DIR + "/01_porechopped_data/{barcode}.fastq", barcode=BARCODES))
	params:
		output_dir=OUTPUT_DIR + "/01_porechopped_data"
	conda:
		"envs/env1.yaml"
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
		(dirs_dict["SINGLE"]),
	output:
		demultiplexed_dir=directory(expand((dirs_dict["DEMULTIPLEXED"] + "/{barcode}"), barcode=BARCODES)),
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
