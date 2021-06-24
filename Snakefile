"""
==============================
Snakemake annotation pipeline
==============================
- USAGE
    snakemake --configfile config.yaml --cores 1

dp24
"""

import os
import errno

# Define the common variables
reference = config['reference_data']
reads_1 = config['input_data_1']
reads_2 = config['input_data_2']
mydir = config['working_dir']
analysis = config['analysis']
trim = analysis['s1_trimmomatic']

print(f'Reference: \t{reference}\n'
      f'Reads 1: \t{reads_1}\n'
      f'Reads 2: \t{reads_2}')

# Set the dict to be parsed and stored in a more usable dict format
analysis_dict = {}

for key, value in analysis.items():
    directory = os.path.join(mydir, key)
    for param, vvlaue in value.items():
        analysis_dict[key] = directory
        try:
            print(f'Working dir: {directory}')
            os.makedirs(directory)
        except OSError as e:
            if e.errno != errno.EEXIST:
                raise

rule all:
    input:
        analysis["QC_1_fastqc"]["output_list"],

        analysis["QC_2_fastqc"]["output_list"],

        "s4_sort_done"

rule QC_1_fastqc:
    input:
        config['input_list']
    output:
        analysis["QC_1_fastqc"]["output_list"]

    run:
        for i in input:
            shell("fastqc -t 8 -o QC_1_fastqc/ {i}")

rule s1_trimmomatic:
    input:
        config['input_list']
    output:
        "s1_trimm_done",
        paired = analysis['s1_trimmomatic']['paired_output'],
        unpaired = analysis['s1_trimmomatic']['unpaired_output'],
        list = zip(
            analysis['s1_trimmomatic']['paired_output'],
            analysis['s1_trimmomatic']['unpaired_output']
        ),

    params:
        threads = trim['threads'],
        tlog = trim['tlog'],
        lead = trim['lead'],
        trail = trim['trail'],
        sliding = trim['sliding'],
        mlen = trim['mlen']
    shell:
        "trimmomatic PE {params.threads} -phred33 {input} {output.list} {params};"
        " touch s1_trimmomatic/s1_trimm_done"

rule QC_2_fastqc:
    input:
        rules.s1_trimmomatic.output
    output:
        analysis["QC_2_fastqc"]["output_list"]
    run:
        for i in input[1]:
            shell("fastqc -t 8 -o qc_data/ {i}")

rule s2_ref_index:
    input:
        rules.s1_trimmomatic.output[0],
        reference
    output:
        analysis['s2_ref_index']['output_list'],
        "s2_ref_done"
    shell:
        "bwa index -a {input[1]};"
        " touch s2_ref_index/s2_ref_done"

rule s3_alignment:
    input:
        reference,
        rules.s1_trimmomatic.output.paired,
        rules.s2_ref_index.output[1]
    output:
        analysis['s3_alignment']['output'],
        "s3_alignment_done"
    shell:
        # Possibly switch this out for bowtie
        "bwa mem -R'@RG\\tID:1\\tLB:library\\tPL:Illumina\\tPU:lane1\\tSM:human'"
        " {input[0]} {input[1]} | samtools view -bS - > {output[0]};"
        " touch s3_alignment/s3_alignment_done"

rule s4_sortbam:
    input:
        rules.s3_alignment.output
    output:
        "s4_sort_done",
        analysis['s4_sortbam']['output']
    shell:
        "samtools sort {input} -o {output};"
        " touch s4_sortbam/s4_sort_done"