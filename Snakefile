SAMPLES = ["SRR622461_1", "SRR622461_2"]
REF_SAMPLE = ["GRCh38_latest_genomic.fna"]
REF_PREFIX = ["GRCh38"]
TRIM_SAMPLES = ["SRR622461"]
TRIM_OUT  = ["rawdata/sample_data/SRR622461_1.paired.fastq",
             "rawdata/sample_data/SRR622461_2.paired.fastq",
             "rawdata/sample_data/SRR622461_1.unpaired.fastq",
             "rawdata/sample_data/SRR622461_2.unpaired.fastq"]
TRIM_FR = ["1", "2"]
TRIM_PAIR = ["paired","unpaired"]

directories = ["qc_1", "qc_2", "qc_3"]

BWA_OUT = ["amb","ann","bwt","pac","sa"]

# Add config file later
rule all:
    input:
        expand("qc_1/{samples}_fastqc.html", samples = SAMPLES),
        expand("qc_1/{samples}_fastqc.zip", samples = SAMPLES),
        expand("rawdata/sample_data/{read}_1.paired.fastq", read = TRIM_SAMPLES),
        expand("rawdata/sample_data/{read}_2.paired.fastq", read = TRIM_SAMPLES),
        expand("rawdata/sample_data/{read}_1.unpaired.fastq", read = TRIM_SAMPLES),
        expand("rawdata/sample_data/{read}_2.unpaired.fastq", read = TRIM_SAMPLES),
        expand("qc_2/{read}_{no}.{pair}.fastq_fastqc.html", read = TRIM_SAMPLES, no = TRIM_FR, pair = TRIM_PAIR),
        expand("rawdata/reference/{ref_prefix}.amb", ref_prefix = REF_SAMPLE),
        expand("rawdata/reference/{ref_prefix}.ann", ref_prefix = REF_SAMPLE),
        expand("rawdata/reference/{ref_prefix}.bwt", ref_prefix = REF_SAMPLE),
        expand("rawdata/reference/{ref_prefix}.pac", ref_prefix = REF_SAMPLE),
        expand("rawdata/reference/{ref_prefix}.sa", ref_prefix = REF_SAMPLE),
        expand("{read}_{ref_prefix}.vcf", read = TRIM_SAMPLES, ref_prefix = REF_PREFIX)

rule step_1_fastqc:
    input:
        "rawdata/sample_data/{samples}.fastq"
    output:
        "qc_1/{samples}_fastqc.html",
        "qc_1/{samples}_fastqc.zip", # So that it can be looked at by MultiQC
    params:
        q="--quiet"
    threads: 1
    priority: 1
    shell:
        "fastqc -t 8 -o qc_1/ {input}"

rule step_2_trimmomatic: # Testing
    input:
        reads1="rawdata/sample_data/{read}_1.fastq",
        reads2="rawdata/sample_data/{read}_2.fastq"
    output:
        trim="rawdata/sample_data/{read}.log",
        paired1="rawdata/sample_data/{read}_1.paired.fastq",
        paired2="rawdata/sample_data/{read}_2.paired.fastq",
        unpaired1="rawdata/sample_data/{read}_1.unpaired.fastq",
        unpaired2="rawdata/sample_data/{read}_2.unpaired.fastq",
    params:
        threads="-threads 4", # Threads in use for this action
        tlog="-trimlog", # Output a log of trim actions
        lead="LEADING:10",   # Cut leading bases > quality 3
        trail="TRAILING:10", # Cut trailing bases > quality 3
        sliding="SLIDINGWINDOW:4:20", # Cut reads > quality 15 in 4 bp windows
        mlen="MINLEN:36" # Drop Reads > 36bp
    priority: 1
    threads: 4
    shell:
        "trimmomatic PE {params.threads} -phred33 {params.tlog} {output.trim} {input.reads1} {input.reads2}"
        " {output.paired1} {output.unpaired1} {output.paired2} {output.unpaired2}"
        " {params.lead} {params.trail} {params.sliding} {params.mlen}"



rule step_3_fastqc_trimmed:
    input:
        paired1=rules.step_2_trimmomatic.output.paired1,
        paired2=rules.step_2_trimmomatic.output.paired2,
        unpaired1=rules.step_2_trimmomatic.output.unpaired1,
        unpaired2=rules.step_2_trimmomatic.output.unpaired2
    output:
        "qc_2/{read}_{no}.{pair}.fastq_fastqc.html",
        "qc_2/{read}_{no}.{pair}.fastq_fastqc.zip"
    params:
        q="--quiet"
    priority: 2
    threads: 4
    shell:
        "fastqc -t 8 -o qc_2/ {input}"

rule step_4_ref_index:
    input:
        "rawdata/reference/{ref_prefix}"
    output:
        "rawdata/reference/{ref_prefix}.amb",
        "rawdata/reference/{ref_prefix}.ann",
        "rawdata/reference/{ref_prefix}.bwt",
        "rawdata/reference/{ref_prefix}.pac",
        "rawdata/reference/{ref_prefix}.sa",
    priority: 1
    threads: 8
    shell:
        "bwa index -a {input}"

rule step_5_alignment:
    input:
        rules.step_2_trimmomatic.output.paired1,
        rules.step_2_trimmomatic.output.paired2,
        ref_seq="rawdata/reference/{ref_prefix}_latest_genomic.fna"
    output:
        "alignment_1/{read}_{ref_prefix}.aligned.bam",
    priority: 3
    threads: 8
    shell:
        "bwa mem -R '@RG\tID:1\tLB:library\tPL:Illumina\tPU:lane1\tSM:human'"
        " {input.ref_seq} {input[0]} {input[1]} | samtools view -bS - > {output}"

rule step_6_freebayes:
    input:
        rules.step_5_alignment.output,
        "rawdata/sample_data/{read}_1.fastq",
        "rawdata/sample_data/{read}_2.fastq"
    output:
        "{read}_{ref_prefix}.vcf"
    priority: 4
    threads: 8
    shell:
        "touch me"
        #"freebayes -f {input[0]} {input[1]} > something.vcf"
