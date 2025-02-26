# specify config file

configfile: "config.yaml"

# Define a global variable for the sleep duration (in seconds)
sleep_duration = 10  # Adjust as needed

rule all:
    input:
        expand("bam/sort_{sample}.bam", sample=config["samples"]), 
        expand("bam/sort_{sample}.bam.bai", sample=config["samples"]),
        expand("bam/bt2_{sample}.log", sample=config["samples"]),
        expand("bigwig/{sample}.bw", sample=config["samples"])

### Trimm 
rule trim_galore:
    input:
        fq1=lambda wildcards: config["samples"][wildcards.sample]["fq1"],
        fq2=lambda wildcards: config["samples"][wildcards.sample]["fq2"]
    output:
        fq1_trimmed="trimmed/{sample}_R1_trimmed.fastq",
        fq2_trimmed="trimmed/{sample}_R2_trimmed.fastq"
    log:
        "logs/trim_galore_{sample}.log"
    threads: 8
    shell:
        """
        module purge
        module load Bioinformatics
        module load trimgalore/0.6.7-ztb2tpz

        trim_galore --paired --fastqc \
        --output_dir trimmed/ \
        --cores {threads} \
        {input.fq1} {input.fq2} > {log} 2>&1
        """

# Read mapping and converting SAM to BAM
rule bowtie2:
    input:
        fq1=lambda wildcards: config["samples"][wildcards.sample]["fq1"],
        fq2=lambda wildcards: config["samples"][wildcards.sample]["fq2"]
    output:
        bam1 = "bam/ori_{sample}.bam"
    params:
        index = "nfs/turbo/umms-siwase/GENOME_INDICES/HG38/bt2_hg38-index"  # change
    log:
        log = "bam/bt2_{sample}.log"
    resources:
        mem_mb = 40000  # 40 GB
    shell:
        """
        module purge
        module load Bioinformatics
        module load python/3.9.12
        module load bowtie2/2.4.2-3pufpzz
        module load samtools
        bowtie2 -p 16 -q --local \
        -x {params.index} \
        -1 {input.fq1} \
        -2 {input.fq2} | samtools view -h -S -b -q 10 --threads 16 > {output.bam1} 
        samtools stats {output.bam1} > bam/{wildcards.sample}_raw_stats.txt
        sleep {sleep_duration}  # Introduce a delay to avoid overloading the server
        """

# Define the rule for sorting BAM
rule bam_sort:
    input:
        bam = "bam/ori_{sample}.bam"
    output:
        bam_sorted = "bam/sort_{sample}.bam"
    resources:
        mem_mb = 40000  # 40 GB
    shell:
        """
        module load Bioinformatics
        module load python/3.9.12
        module load samtools
        samtools sort -o {output.bam_sorted} -T {wildcards.sample}_temp --threads 16 {input.bam} -m 1G
        sleep {sleep_duration}  # Introduce a delay to avoid overloading the server
        """

# Define the rule for indexing BAM
rule bam_index:
    input:
        bam_sorted = "bam/sort_{sample}.bam"
    output:
        bai1 = "bam/sort_{sample}.bam.bai"
    shell:
        """
        module load Bioinformatics
        module load python/3.9.12
        module load samtools
        samtools index {input.bam_sorted}
        sleep {sleep_duration}  # Introduce a delay to avoid overloading the server
        """

# Define the rule for generating bigwig files
rule bam_to_bigwig:
    input:
        bam_sorted = "bam/sort_{sample}.bam",
        bai1 = "bam/sort_{sample}.bam.bai"
    output:
        bigwig = "bigwig/{sample}.bw"
    shell:
        """
        source ~/miniconda3/etc/profile.d/conda.sh  # change if needed
        conda activate deeptools
        bamCoverage -b {input.bam_sorted} -of bigwig -o {output.bigwig} --normalizeUsing BPM -bs 1 -p 16
        sleep {sleep_duration}  # Introduce a delay to avoid overloading the server
        """

