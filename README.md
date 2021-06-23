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

3 - NOT RUN YET

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