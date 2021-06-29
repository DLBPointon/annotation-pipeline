# Annotation Pipeline for BT&S 

This repository will make up the evidence for the Bioinformatics Software and Tools module
of the ARU/Sanger BSc Bioinformatics Year 2 Course.

### Installation of pipeline
This can be achieved through the use of:

`conda create --file environment.yml`

However, this can take quite some time and personally was faster
installing it all separately.

### Usage

Local Machine Usage
1 - Clone repo

2 - `cd annotation-pipeline`

3 - `snakemake --cores 8`

Sanger Farm
1 - Clone repo

2 - `cd annotation-pipeline`

3 - NOT RUN YET - Individual steps are being tested

### Currently in Repo
environment.yml - a list of packages used in this project
multiqc_report.html - a report of all fastqc report, this should soon include stats for samtools flagstats, bwa and bcftools too
testanno.vcf - a vcf before attempted annotation
get_scaff_names.py - used along with "grep -F 'chromosome 10' > components.txt " to return a file of component names used to cut down the original vcf size to focus on region of interest.
annoated.vcf - a vcf from my first attempt using snpeff - chromosomes not found (I'm assuming a naming of scaffolds issue from reference)

### Data
Sample Data comes from the Utah family platinum read set.
```
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR622/SRR622461/SRR622461_1.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR622/SRR622461/SRR622461_2.fastq.gz
```

Mapped against GRCH38.p.15:
```
https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/GCA_000001405.15_GRCh38_genomic.fna.gz
```
Snakefile_old ==> this was the first draft using numerous wildcards to find files, this caused issues.
The new Snakefile although more complex is able to give much more transparency and control over the work involved.

It has been suggested to swap out freebayes for GATK.

SNPeff for variant annotation and filtering.

Move the snakemake pipeline over to a yaml oriented snakemake (yy5 pipes as example) for better control.