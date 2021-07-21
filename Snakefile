"""
==============================
Snakemake annotation pipeline
==============================
- USAGE
    snakemake --configfile config.yaml --cores 10

dp24
"""

import os
import errno

singularity: "shub://sschmeier/biotools:latest"

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
        # Branch 1 - Produces alignment statistics
        analysis['s7_align_metrics']['output'],

        # Branch 2 - multiqc brings together FastQC reports
        'QC_3_multiqc/QC_3_done',

        # Main Pipeline
        's13_bgzip_index/s13_done'

rule QC_1_fastqc:
    input:
        reads_1,
        reads_2
    output:
        analysis["QC_1_fastqc"]["output_list"]

    shell:
        "fastqc -t 2 -o QC_1_fastqc/ {input}"

rule s1_trimmomatic:
    input:
        reads_1,
        reads_2
    output:
        "s1_trimmomatic/s1_trimm_done",
        list = zip(
            analysis['s1_trimmomatic']['paired_output'],
            analysis['s1_trimmomatic']['unpaired_output']
        )

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
        rules.s1_trimmomatic.output.list
    output:
        # May need to be changed so that multiqc doesn't take in .zip folder
        analysis["QC_2_fastqc"]["output_list"]
    shell:
        "fastqc -t 8 -o QC_2_fastqc/ {input}"

rule QC_multiqc:
    input:
        rules.QC_1_fastqc.output,
        rules.QC_2_fastqc.output
    output:
        #analysis['QC_3_multiqc']['output'],
        "QC_3_multiqc/QC_3_done"
    shell:
        "multiqc {input};"
        "touch QC_3_multiqc/QC_3_done"

rule s2_ref_index:
    input:
        rules.s1_trimmomatic.output[0],
        reference
    output:
        analysis['s2_ref_index']['output_list'],
        "s2_ref_index/s2_ref_done",
        reference + '.fai'
    shell:
        "bwa index {input[1]};"
        "samtools faidx {input[1]};"
        " touch s2_ref_index/s2_ref_done"

# Alignment stats should to feed into MultiQC, although results don't look right
rule s3_alignment:
    input:
        reference,
        analysis['s1_trimmomatic']['paired_output'],
        "s2_ref_index/s2_ref_done"
    output:
        analysis['s3_alignment']['output'],
        "s3_alignment/s3_alignment_done"
    params:
        threads = 4
    shell:
        # Possibly switch this out for bowtie <--
        "bwa mem -t {params.threads} -R'@RG\\tID:1\\tLB:library\\tPL:Illumina\\tPU:lane1\\tSM:human'"
        " {input[0]} {input[1]} | samtools view -bS - > {output[0]} -@{params.threads};"
        " touch s3_alignment/s3_alignment_done"

rule s4_sortbam:
    input:
        "s3_alignment/s3_alignment_done",
        rules.s3_alignment.output
    output:
        "s4_sortbam/s4_sort_done",
        analysis['s4_sortbam']['output']
    params:
        threads = 4
    shell:
        "samtools sort {input[1]} -o {output[1]} -@{params.threads};"
        " touch s4_sortbam/s4_sort_done"

rule s5_index_bam:
    input:
        rules.s4_sortbam.output[0],
        rules.s4_sortbam.output[1]
    output:
        analysis['s5_index_bam']['output'],
        "s5_index_bam/s5_index_bam_done"
    params:
        threads = 4
    shell:
        "samtools index {input[1]} -@{params.threads}; "
        "touch {output[1]}"

rule s6_mrkdupes:
    input:
        rules.s5_index_bam.output[1],
        rules.s4_sortbam.output[1]
    output:
        analysis['s6_mrkdupes']['output'],
        analysis['s6_mrkdupes']['output'] + '.bai',
        analysis['s6_mrkdupes']['mrkd_output'],
        "s6_mrkdupes/s6_done"
    params:
        threads = 4
    shell:
        "picard MarkDuplicates -I {input[1]} -O {output[0]} -M {output[2]};"
        "samtools index {output[0]} -@{params.threads};"
        "touch {output[3]}"

rule s7_align_metrics:
    input:
        reference,
        rules.s6_mrkdupes.output[0],
        rules.s6_mrkdupes.output[3]
    output:
        analysis['s7_align_metrics']['output']
    shell:
        "picard CollectAlignmentSummaryMetrics -R {reference} -I {input[1]} -O {output}"

rule s8_picard_dictionary:
    input:
        reference,
        "s6_mrkdupes/s6_done"
    output:
        analysis['s8_picard_dict']['output'],
        's8_picard_dict/s8_dict_done'
    shell:
        "picard CreateSequenceDictionary -R {reference} -O {output[0]};"
        "touch {output[1]}"

rule s9_GATK_vcf:
    input:
        "s6_mrkdupes/s6_done",
        's8_picard_dict/s8_dict_done',
        reference,
        rules.s6_mrkdupes.output[0]
    output:
        's9_GATK_vcf/s9_done',
        analysis['s9_GATK_vcf']['output']
    shell:
        "gatk HaplotypeCaller -R {reference} -I {input[3]} -O {output[1]} -ERC GVCF;"
        "touch {output[0]}"

rule s10_excise_chr:
    input:
        "s9_GATK_vcf/s9_done",
        rules.s9_GATK_vcf.output[1],
        reference
    params:
        config['search_chr']
    output:
        analysis['s10_excise_gene']['output'],
        "s10_excise_gene/s10_done"
    shell:
        "bcftools view {input[1]} {params} > {output[0]};"
        "touch {output[1]}"

rule s11_filter_qual:
    input:
        "s10_excise_gene/s10_done",
        rules.s10_excise_chr.output[0]
    params:
        analysis['s11_filter_qual']['qual']
    output:
        analysis['s11_filter_qual']['output'],
        's11_filter_qual/s11_done'
    shell:
        "bcftools view -i '{params}' {input[1]} > {output[0]};"
        "touch {output[1]}"

rule s12_annotate_vcf:
    input:
        's11_filter_qual/s11_done',
        rules.s11_filter_qual.output[0]
    output:
        analysis['s12_annotate_vcf']['output'],
        analysis['s12_annotate_vcf']['output_csv'],
        's12_annotate_vcf/s12_done'
    shell:
        'java -Xmx4g -jar snpEff/snpEff.jar GRCh38.99 {input[1]} > {output[0]} -csvStats {output[1]} -lof;'
        'touch {output[2]}'

rule s13_bgzip_index:
    input:
        rules.s12_annotate_vcf.output[2],
        rules.s12_annotate_vcf.output[0]
    output:
        analysis['s13_bgzip_index']['output_bg'],
        analysis['s13_bgzip_index']['output_tbi'],
        's13_bgzip_index/s13_done'
    shell:
        "bgzip -c {input[1]} > {output[0]};"
        "tabix -p vcf {output[0]};"
        "touch {output[2]}"