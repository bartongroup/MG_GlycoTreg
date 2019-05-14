# Snakemake file for processing RNA-seq data with STAR and Salmon
# Using Ensembl genome, annotations and transcriptome
# STAR options are set for GTF annotations

import glob

# Directories
fastqDir = "fastq/"
genomeDir = "genome/"
starIndexDir = "starindex/"
salmonIndexDir = "salmonindex/"

# Genome URLs
genomeURL = "ftp://ftp.ensembl.org/pub/release-96/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna_rm.primary_assembly.fa.gz"
gtfURL = "ftp://ftp.ensembl.org/pub/release-96/gtf/mus_musculus/Mus_musculus.GRCm38.96.gtf.gz"
cdnaURL = "ftp://ftp.ensembl.org/pub/release-96/fasta/mus_musculus/cdna/Mus_musculus.GRCm38.cdna.all.fa.gz"
cdsURL = "ftp://ftp.ensembl.org/pub/release-96/fasta/mus_musculus/cds/Mus_musculus.GRCm38.cds.all.fa.gz"

# Genome files
genomeFile = genomeDir + "Mus_musculus.GRCm38.dna_rm.primary_assembly.fa"
gtfFile = genomeDir + "Mus_musculus.GRCm38.93.gtf"
cdsFile = genomeDir + "Mus_musculus.GRCm38.cds.all.fa"
transcriptomeFile = genomeDir + "Mus_musculus.GRCm38.cdna.all.fa"
genomeIndex = genomeFile + ".fai"
genomeSizeFile = genomeFile + ".txt"

# Files to test if an index was created
starIndexTestFile = starIndexDir + "chrName.txt"
salmonIndexTestFile = salmonIndexDir + "indexing.log"

# for bedgraph
windowSize = "1000000"
windowFile = "bedgraph/window." + windowSize + ".bed"

SAMPLES = ["TregRest_" + str(i) for i in range(1,5)] + ["TregAct_" + str(i) for i in range(1,5)]
PAIRS = ["R1", "R2"]

# Lists for rule all
QCS = expand("qc/{sample}_{pair}_fastqc.html", sample=SAMPLES, pair=PAIRS)
MULTIQC = ["multiqc/report.html"]
NGSREPORT = ["qc/ngsReports_Fastqc.html"]
BAMS = expand("bam/{sample}.bam", sample=SAMPLES)
BAIS = expand("bam/{sample}.bam.bai", sample=SAMPLES)
COUNTS = expand("readcount/{sample}.txt", sample=SAMPLES)
QUANTS = expand("salmon/{sample}/quant.sf", sample=SAMPLES)
BEDGRAPHS = expand("bedgraph/{sample}." + windowSize + ".bedgraph", sample=SAMPLES)

####################################################################

rule all:
    input: [cdsFile] + QCS + MULTIQC + COUNTS + BAIS + BEDGRAPHS


####################################################################
# Merge lanes

rule merge_lanes:
    output:
      R1 = "fastq/{sample}_R1.fastq.gz",
      R2 = "fastq/{sample}_R2.fastq.gz"
    shell:
        """
        cat Samples/{wildcards.sample}/Files/*_R1_001.fastq.gz > {output.R1}
        cat Samples/{wildcards.sample}/Files/*_R2_001.fastq.gz > {output.R2}
        """

####################################################################
# Quality control

rule fastqc:
    input:
        R1 = "fastq/{sample}_R1.fastq.gz",
        R2 = "fastq/{sample}_R2.fastq.gz"
    output:
        "qc/{sample}_R1_fastqc.html",
        "qc/{sample}_R2_fastqc.html"
    threads: 2
    shell:
        "fastqc -o qc --threads {threads} -f fastq {input.R1} {input.R2}"

rule ngsReports:
    input: expand("qc/{sample}_{pair}_fastqc.html", sample=SAMPLES, pair=PAIRS)
    output: "qc/ngsReports_Fastqc.html"
    shell:
        "Rscript Rscript/ngs_reports.R fastq"

rule multiqc:
    input: expand("qc/{sample}_{pair}_fastqc.html", sample=SAMPLES, pair=PAIRS)
    output: "multiqc/report.html"
    shell:
        """
        mkdir -p multiqc
        multiqc -f --filename report --outdir multiqc qc
        """


####################################################################
# Load genome files

rule load_db_files:
    output: genomeFile, gtfFile, transcriptomeFile, cdsFile
    shell:
        """
        wget {genomeURL} -O - | gunzip -c > {genomeFile}
        wget {gtfURL} -O - | gunzip -c > {gtfFile}
        wget {cdnaURL} -O - | gunzip -c > {transcriptomeFile}
        wget {cdsURL} -O - | gunzip -c > {cdsFile}
        """


####################################################################
# Index genome, create chromosome size file

rule index_genome:
    input: genomeFile
    output: genomeIndex
    shell:
        "samtools faidx {input}"

rule size_genome:
    input: genomeIndex
    output: genomeSizeFile
    shell:
        "cut -f 1,2 {input} > {output}"


####################################################################
# Salmon

rule salmon_index:
    input: transcriptomeFile
    output: salmonIndexTestFile
    log: "logs/salmon_index.log"
    shell:
        "salmon index -t {input} -i {salmonIndexDir} &> {log}"


rule salmon_quant:
    input:
        R1 = fastqDir + "{sample}_R1.fastq.gz",
        R2 = fastqDir + "{sample}_R2.fastq.gz",
        testfile = salmonIndexTestFile
    output: "salmon/{sample}/quant.sf"
    params:
        prefix = "salmon/{sample}"
    threads: 12
    log: "logs/salmon_{sample}_quant.log"
    shell:
        """
        salmon quant \
        --index {salmonIndexDir} \
        --libType A \
        --numBootstraps 100 \
        --threads {threads} \
        -1 {input.R1} -2 {input.R2} \
        --output {params.prefix} &> {log}
        """


####################################################################
# STAR

rule star_index:
    input:
      fasta = genomeFile,
      gtf = gtfFile
    output: starIndexTestFile
    threads: 12
    log: "logs/star_index.log"
    shell:
        """
        STAR \
        --runMode genomeGenerate \
        --genomeFastaFiles {input.fasta} \
        --sjdbGTFfile {input.gtf} \
        --genomeDir {starIndexDir} \
        --outFileNamePrefix logs/star_ \
        --runThreadN {threads} &> {log}
        """


rule star_mapping:
    input:
        R1 = fastqDir + "{sample}_R1.fastq.gz",
        R2 = fastqDir + "{sample}_R2.fastq.gz",
        testfile = starIndexTestFile
    output:
        bam = "bam/{sample}.bam",
        readcount = "readcount/{sample}.txt",
        finallog = "starmap/{sample}_Log.final.out"
    threads: 12
    log: "starmap/{sample}_run.log"
    shell:
        """
         STAR \
         --genomeDir {starIndexDir} \
         --sjdbGTFfile {gtfFile} \
         --readFilesIn {input.R1} {input.R2} \
         --outFileNamePrefix starmap/{wildcards.sample}_ \
         --outSAMtype BAM SortedByCoordinate \
         --outFilterMultimapNmax 2 \
         --readFilesCommand zcat \
         --outReadsUnmapped Fastx \
         --quantMode GeneCounts \
         --runThreadN {threads} &> {log}
         mv starmap/{wildcards.sample}_Aligned.sortedByCoord.out.bam {output.bam}
         mv starmap/{wildcards.sample}_ReadsPerGene.out.tab {output.readcount}
        """


####################################################################
# Index BAM files

rule index_bam:
    input: "bam/{sample}.bam"
    output: "bam/{sample}.bam.bai"
    threads: 8
    log: "logs/{sample}.index_bam.log"
    shell:
        "samtools index {input} &> {log}"


####################################################################
# BED and bedGraph files

rule bamtobed:
    input: "bam/{sample}.bam"
    output: "bed/{sample}.bed"
    log: "logs/{sample}_bam2bed.log"
    shell:
        "bedtools bamtobed -i {input} | sort -k1,1 -k2,2n > {output} 2> {log}"

rule bedwindow:
    input: genomeSizeFile
    output: windowFile
    log: "logs/make_window.log"
    shell:
        "bedtools makewindows -w {windowSize} -s {windowSize} -g {input} > {output} 2> {log}"

rule bedgraph:
    input: 
        bed = "bed/{sample}.bed",
        window = windowFile
    output: "bedgraph/{sample}." + windowSize + ".bedgraph"
    log: "logs/{sample}_bedgraph.log"
    shell:
        "bedtools intersect -c -a {input.window} -b {input.bed} 2> {log} > {output}" 


