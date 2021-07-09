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

3 - Download data into the raw_data folders (links below).

4 - `snakemake --configfile config.yaml --cores 10`

This pipeline will require the igv.sh file to have been set to 15gb.
(-Xmx15g compared to -Xmx4g) this is in order to produce the .genome file.

If this can't be used then please use the built-in database for Human (hg38).

snpEff is set to -Xmx4g in order to annotate the vcf.

After completion:

5 - Run IGV

In my case this was `bash {location of installation}/igv.sh`

6 - OPTION A - build .genome file

`Genomes > Create .genome`

The fasta file = the reference genome

The Gene file is the .gff file downloaded earlier

6 - OPTION B - use pre-installed Human (hg38)

7 - Load annotated file

`File > load from file`

Navigate to folder s13 witch should contain something akin to:
`chr10.filt.annotated.vcf.gz`

8 - Navigate to location `chr10:94,760,000-94,860,000` 
which will centre on the gene _CYP2C19_.

### Currently in Repo
environment.yml - a list of packages used in this project

### Data
Files should be downloaded into a `{project dir}/raw_data/{ reference | sample_data }` folder.

Sample Data comes from the Utah family platinum read set.
```
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR622/SRR622461/SRR622461_1.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR622/SRR622461/SRR622461_2.fastq.gz
```

Mapped against GRCH38.p15:
```
https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_full_analysis_set.fna.gz
```

We also used:
```
https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_full_analysis_set.refseq_annotation.gff.gz
```

This was used in conjunction with the reference genome to
build a .genome file to better compare the results of this pipeline when
visualising th end product with IGV.

TO DO:
- Write a python script to write a yaml for the snakemake pipeline.
        
        - will require arg for snpEff jar location