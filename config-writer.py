"""
Config-writer
==============
This script will write the config used in executing the snakemake pipeline


"""

import argparse
import yaml


def parseargs():
    parser = argparse.ArgumentParser()
    parser.add_argument("REFERENCE",
                        help="Reference Genome FASTA")
    parser.add_argument("SAMPLE_1",
                        help="Sample FASTQ 1")
    parser.add_argument("SAMPLE_2",
                        help="Sample FASTQ 2")
    parser.add_argument("-C", "--Chromosome",
                        help="Chromosome of interest")
    parser.add_argument("-conda", "--conda",
                        help="Chromosome of interest")
    parser.add_argument("-threads", "--threads",
                        help="Threads used by trimmomatic and samtools")

    args = parser.parse_args()
    return args


def file_loader():
    with open("sample.config.yaml", 'r') as file:
        base_yaml = yaml.load(file, Loader=yaml.FullLoader)

    return base_yaml


def write_yaml(data):
    with open("tested.yaml", 'w') as stream:
        try:
            yaml.dump(data, stream, default_flow_style=False)
        except yaml.YAMLError as exc:
            print(exc)


def main():
    args = parseargs()
    base_yaml = file_loader()

    dir = 'rawdata/sample_data/'
    if not dir:
        dir = None
    base_yaml['input_data_1'] = dir + args.SAMPLE_1
    base_yaml['input_data_2'] = dir + args.SAMPLE_2
    base_yaml['reference_data'] = "rawdata/reference/" + args.REFERENCE

    if not args.Chromosome:
        base_yaml['search_chr'] = None
    else:
        base_yaml['search_chr'] = f'-r "{args.Chromosome}"'

    if not args.conda:
        base_yaml['conda_envs'] = None
    else:
        base_yaml['conda_envs'] = args.conda

    # ANALYSIS CHUNK
    anal = base_yaml['analysis']
    sample_1 = args.SAMPLE_1.split('/')[-1].split('.')
    sample_2 = args.SAMPLE_2.split('/')[-1].split('.')

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

    anal['s1_trimmomatic']['paired_output'] = ['s1_trimmomatic/' + sample_1[0] + '.paired.trimmed.fastq',
                                               's1_trimmomatic/' + sample_2[0] + '.paired.trimmed.fastq']

    anal['s1_trimmomatic']['unpaired_output'] = ['s1_trimmomatic/' + sample_1[0] + '.unpaired.trimmed.fastq',
                                                 's1_trimmomatic/' + sample_2[0] + '.unpaired.trimmed.fastq']

    if not args.threads:
        anal['s1_trimmomatic']['threads'] = f"-threads 4"
    else:
        anal['s1_trimmomatic']['threads'] = f"-threads {args.threads}"

    if not args.tlog:
        anal['s1_trimmomatic']['threads'] = f"-threads 4"
    else:
        anal['s1_trimmomatic']['threads'] = f"-threads {args.threads}"

    if not args.lead:
        anal['s1_trimmomatic']['threads'] = f"-threads 4"
    else:
        anal['s1_trimmomatic']['threads'] = f"-threads {args.threads}"

    if not args.trail:
        anal['s1_trimmomatic']['threads'] = f"-threads 4"
    else:
        anal['s1_trimmomatic']['threads'] = f"-threads {args.threads}"

    if not args.slide:
        anal['s1_trimmomatic']['threads'] = f"-threads 4"
    else:
        anal['s1_trimmomatic']['threads'] = f"-threads {args.threads}"

    if not args.minlen1:
        anal['s1_trimmomatic']['threads'] = f"-threads 4"
    else:
        anal['s1_trimmomatic']['threads'] = f"-threads {args.threads}"

    print(base_yaml)
    write_yaml(base_yaml)


if __name__ == "__main__":
    main()
