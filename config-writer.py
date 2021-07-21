"""
Config-writer
==============
This script will write the config used in executing the
snakemake pipeline

"""

import argparse
import yaml


def parseargs():
    """
    argparse function to collect cli arguments
    :return: args
    """
    parser = argparse.ArgumentParser()

    parser.add_argument("REFERENCE",
                        help="Reference Genome FASTA")
    parser.add_argument("SAMPLE_1",
                        help="Sample FASTQ 1")
    parser.add_argument("SAMPLE_2",
                        help="Sample FASTQ 2")
    parser.add_argument("REGION",
                        help="Chromosome of interest")
    parser.add_argument("-cfile", "--config-name",
                        dest='cfile',
                        help="Configfile naming",
                        default='./config_file.yaml')
    parser.add_argument("-wdir", '--WORKING-DIR',
                        dest='wd',
                        help='Working Directory Loc',
                        default='./')
    parser.add_argument('-q', "--FILTER-QUAL",
                        dest='qual',
                        help="Filter for Quality above this",
                        default='20')
    parser.add_argument("-conda", "--conda",
                        help="Conda Library Location",
                        default='None')
    parser.add_argument("-threads", "--threads",
                        dest='threads',
                        help="Threads used by trimmomatic and samtools",
                        default='4')
    parser.add_argument("-tlog", "--TRIMMER-LOG",
                        dest='tlog',
                        help="The log file for Trimmomatic",
                        default='./trim.log')
    parser.add_argument("-lead", "--LEADING",
                        dest='lead',
                        help="Trim length of leading sequence",
                        default='10')
    parser.add_argument("-trail", "--TRAILING",
                        dest='trail',
                        help="Trim trailing sequence",
                        default='10')
    parser.add_argument("-slide", "--SLIDING-WINDOW",
                        dest='slide',
                        help="Sliding Window",
                        default='4:20')
    parser.add_argument("-minl", "--MINLEN",
                        dest='mlen',
                        help="Minimum Length of sequence",
                        default='36')

    args = parser.parse_args()
    return args


def file_loader():
    with open("sample.config.yaml", 'r') as file:
        base_yaml = yaml.load(file, Loader=yaml.FullLoader)

    return base_yaml


def write_yaml(data, name):
    with open(name, 'w') as stream:
        try:

            yaml.safe_dump(data, stream, default_flow_style=False)
        except yaml.YAMLError as exc:
            print(exc)


def main():
    args = parseargs()
    base_yaml = file_loader()

    if not args.REFERENCE:
        print('Missing REFERENCE positional arg <- you should have this')

    if not args.SAMPLE_1 or not args.SAMPLE_2:
        print('Missing one of the SAMPLE positional arg <- you should have this')

    if not args.REGION:
        print('Missing REGION positional arg <- you should have this')

    base_yaml['input_data_1'] = args.SAMPLE_1
    base_yaml['input_data_2'] = args.SAMPLE_2
    base_yaml['reference_data'] = args.REFERENCE

    base_yaml['samtool_thread'] = args.threads

    base_yaml['working_dir']= args.wd

    base_yaml['search_chr'] = f'-r "chr{args.REGION}"'

    base_yaml['conda_envs'] = args.conda

    # ANALYSIS CHUNK
    anal = base_yaml['analysis']
    sample_a = args.SAMPLE_1.split('/')[-1]
    sample_1 = sample_a.split('.')
    sample_b = args.SAMPLE_2.split('/')[-1]
    sample_2 = sample_b.split('.')

    # ----- QC steps
    anal['QC_1_fastqc']['output_list'] = ['QC_1_fastqc/' + sample_1[0] + '_fastqc.html',
                                          'QC_1_fastqc/' + sample_1[0] + '_fastqc.zip',
                                          'QC_1_fastqc/' + sample_2[0] + '_fastqc.html',
                                          'QC_1_fastqc/' + sample_2[0] + '_fastqc.zip']

    anal['QC_2_fastqc']['output_list'] = ['QC_2_fastqc/' + sample_1[0] + '.paired.trimmed_fastqc.html',
                                          'QC_2_fastqc/' + sample_1[0] + '.paired.trimmed_fastqc.zip',
                                          'QC_2_fastqc/' + sample_2[0] + '.paired.trimmed_fastqc.html',
                                          'QC_2_fastqc/' + sample_2[0] + '.paired.trimmed_fastqc.zip',
                                          'QC_2_fastqc/' + sample_1[0] + '.unpaired.trimmed_fastqc.html',
                                          'QC_2_fastqc/' + sample_1[0] + '.unpaired.trimmed_fastqc.zip',
                                          'QC_2_fastqc/' + sample_2[0] + '.unpaired.trimmed_fastqc.html',
                                          'QC_2_fastqc/' + sample_2[0] + '.unpaired.trimmed_fastqc.zip']

    # Ignore QC_3 as it doesn't need changing

    # ----- S1
    anal['s1_trimmomatic']['paired_output'] = ['s1_trimmomatic/' + sample_1[0] + '.paired.trimmed.fastq',
                                               's1_trimmomatic/' + sample_2[0] + '.paired.trimmed.fastq']

    anal['s1_trimmomatic']['unpaired_output'] = ['s1_trimmomatic/' + sample_1[0] + '.unpaired.trimmed.fastq',
                                                 's1_trimmomatic/' + sample_2[0] + '.unpaired.trimmed.fastq']

    anal['s1_trimmomatic']['threads'] = f"-threads {args.threads}"

    anal['s1_trimmomatic']['tlog'] = f"-trimlog {args.tlog}"

    anal['s1_trimmomatic']['lead'] = f"LEADING:{args.lead}"

    anal['s1_trimmomatic']['trail'] = f"TRAILING:{args.trail}"

    anal['s1_trimmomatic']['sliding'] = f"SLIDINGWINDOW:{args.slide}"

    anal['s1_trimmomatic']['mlen'] = f"MINLEN:{args.mlen}"

    # ----- S2
    anal['s2_ref_index']['output_list'] = [
        f'{args.REFERENCE}.amb',
        f'{args.REFERENCE}.ann',
        f'{args.REFERENCE}.bwt',
        f'{args.REFERENCE}.pac',
        f'{args.REFERENCE}.sa',
    ]

    # ----- S3
    sample = sample_1[0].split('_')[0]
    ref = args.REFERENCE.split('/')[2].split('_')[0]
    anal['s3_alignment']['output'] = f's3_alignment/{sample}_{ref}.aligned.bam'

    # ----- S4
    anal['s4_sortbam']['output'] = f's4_sortbam/{sample}_{ref}.aligned.sort.bam'

    # ----- S5
    # Saved to s4 folder on purpose, indexes should be with the file they index
    anal['s5_index_bam']['output'] = f's4_sortbam/{sample}_{ref}.aligned.sort.bam'

    # ----- S6
    anal['s6_mrkdupes']['output'] = f's6_mrkdupes/{sample}_{ref}.al.s.mkdup.bam'
    anal['s6_mrkdupes']['mrkd_output'] = f's6_mrkdupes/marked_dupes_metrics.txt'

    # ----- S7
    anal['s7_align_metrics']['output'] = f's7_align_metrics/output.txt'

    # ----- S8
    anal['s8_picard_dict']['output'] = f'{args.REFERENCE}.dict'

    # ----- S9
    anal['s9_GATK_vcf']['output'] = f's9_GATK_vcf/{sample}_{ref}.asmrkd.vcf.gz'

    # ----- S10
    anal['s10_excise_gene']['output'] = f's10_excise_gene/chr{args.REGION}.vcf'

    # ----- S11
    anal['s11_filter_qual']['output'] = f's11_filter_qual/chr{args.REGION}.filtered.vcf'

    anal['s11_filter_qual']['qual'] = f'QUAL>{args.qual}'

    # ----- S12
    anal['s12_annotate_vcf']['output'] = f's12_annotate_vcf/chr{args.REGION}.filt.annotated.vcf'
    anal['s12_annotate_vcf']['output_csv'] = f's12_annotate_vcf/snpEff.csv'

    # ----- S13
    anal['s13_bgzip_index']['output_bg'] = f's13_bgzip_index/chr{args.REGION}.filt.annotated.vcf.gz'
    anal['s13_bgzip_index']['output_tbi'] = f's13_bgzip_index/chr{args.REGION}.filt.annotated.vcf.gz.tbi'

    # print(base_yaml)
    write_yaml(base_yaml, args.cfile)


if __name__ == "__main__":
    main()
