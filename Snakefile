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
RAW_DATA_DIR =OUTPUT_DIR + config['input_dir']
BARCODES = config["barcodes"].split()
FLOWCELL=config['flowcell']
KIT=config['kit']

dir_list = ["RULES_DIR","ENVS_DIR","DB","SINGLE_DATA_DIR", "DEMULTIPLEXED", "BASECALLED"]
dir_names = ["rules", "../envs", OUTPUT_DIR + "/db", OUTPUT_DIR + "/01_SINGLE_DATA_DIR", OUTPUT_DIR + "/02_DEMULTIPLEXED", OUTPUT_DIR + "/02_BASECALLED"]
dirs_dict = dict(zip(dir_list, dir_names))

SAMPLES,=glob_wildcards(RAW_DATA_DIR + "/{sample}_" +".fast5")

#======================================================
# Rules
#======================================================

rule all:
	input:
		expand(dirs_dict["DEMULTIPLEXED"] + "/{barcode}.fastq", barcode=BARCODES),

rule multi_to_single_fast5:
    input:
        RAW_DATA_DIR,
    output:
        directory(dirs_dict["SINGLE_DATA_DIR"]),
    conda:
        "envs/env1.yaml"
    message:
        "Converting multi fast5 to single fast5"
    threads: 16
    shell:
        """
		multi_to_single_fast5 --input_path {input} --save_path {output} -t {threads}
		"""


rule demultiplexing:
    input:
		directory(dirs_dict["SINGLE_DATA_DIR"]),
	output:
		demultiplexed_dir=directory(dirs_dict["DEMULTIPLEXED"]),
		rapid_model=dirs_dict["DB"]+ "/Deepbinner/RBK004_read_starts",
    params:
        rapid_model= dirs_dict["DB"]+ "/Deepbinner"
    conda:
        "envs/env1.yaml"
    message:
        "Demultiplexing fast5 files with Deepbinner"
    shell:
        """
		git clone https://github.com/rrwick/Deepbinner/
		cp Deepbinner/models/SQK-RBK004_read_starts {params.rapid_model}
		rm -rf Deepbinner
		deepbinner realtime --in_dir {input} --out_dir {output.demultiplexed_dir} -s {output.rapid_model}
		"""

rule basecalling:
    input:
        demultiplexed_dir=directory(dirs_dict["DEMULTIPLEXED"]/{barcode}),
    output:
		basecalled_barcode=directory(dirs_dict["BASECALLED"] + /{barcode})),
    params:
        rapid_model=dirs_dict["DB"]+ "/Deepbinner"
		flowcell=FLOWCELL
		kit=KIT
    conda:
        "envs/env1.yaml"
    message:
        "Basecalling single fast5 files with guppy"
    threads: 16
    shell:
        """
		guppy_basecaller -i {wildcards.barcode} -s {wildcards.barcode} --fast5_out -q 0 -r --trim_barcodes -x 'cuda:0' \
		--flowcell {params.flowcell} --kit {params.kit} –cpu_threads_per_caller {threads}
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
        "Demultiplexing step 1"
    threads: 16
    shell:
        """
        head -n 25 scripts/logo.txt
        counter=1
        n=$(ls -l {input}/*fastq | wc -l )
        rm -f {params.output_dir}/*fastq
        for filename in {input}/*fastq
        do
            echo "Processing sample $counter/$n"
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

#include: os.path.join(RULES_DIR, 'resultsParsing.smk')
