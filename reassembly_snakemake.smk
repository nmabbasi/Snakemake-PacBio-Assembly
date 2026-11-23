import os

configfile: "config/reassembly_config.yaml"

SAMPLES = config["samples"]
SAMPLE_NAMES = list(SAMPLES.keys())
OUTDIR = config["outdir"]
REF = config["reference"]


rule all:
    input:
        f"{OUTDIR}/spades/contigs.fasta",


rule bwa_index:
    input:
        ref=REF,
    output:
        idx=multiext(REF, ".amb", ".ann", ".bwt", ".pac", ".sa"),
    log:
        f"{OUTDIR}/logs/bwa/index.log",
    conda:
        "envs/reassembly.yaml"
    threads: 1
    resources:
        mem_mb=4096,
        time_min=60,
    shell:
        r"""
        set -euo pipefail
        bwa index {input.ref} > {log} 2>&1
        """


rule bwa_mem:
    input:
        ref=REF,
        idx=multiext(REF, ".amb", ".ann", ".bwt", ".pac", ".sa"),
        r1=lambda wc: SAMPLES[wc.sample]["r1"],
        r2=lambda wc: SAMPLES[wc.sample]["r2"],
    output:
        bam=f"{OUTDIR}/aligned/{{sample}}.sorted.bam",
        bai=f"{OUTDIR}/aligned/{{sample}}.sorted.bam.bai",
    log:
        bwa=f"{OUTDIR}/logs/bwa/{{sample}}_mem.log",
        samtools=f"{OUTDIR}/logs/bwa/{{sample}}_samtools.log",
    conda:
        "envs/reassembly.yaml"
    threads: 8
    resources:
        mem_mb=16384,
        time_min=180,
    shell:
        r"""
        set -euo pipefail
        SORT_THREADS=$(( {threads} > 2 ? 2 : 1 ))
        ALIGN_THREADS=$(( {threads} > 2 ? {threads} - 2 : 1 ))
        bwa mem -t "$ALIGN_THREADS" {input.ref} {input.r1} {input.r2} 2> {log.bwa} \
          | samtools sort -@ "$SORT_THREADS" -o {output.bam} - 2> {log.samtools}
        samtools index -@ "$SORT_THREADS" {output.bam} 2>> {log.samtools}
        """


rule extract_reads:
    input:
        bam=f"{OUTDIR}/aligned/{{sample}}.sorted.bam",
        bai=f"{OUTDIR}/aligned/{{sample}}.sorted.bam.bai",
        read_list=lambda wc: SAMPLES[wc.sample]["read_list"],
    output:
        r1=f"{OUTDIR}/extracted/{{sample}}_R1.fastq.gz",
        r2=f"{OUTDIR}/extracted/{{sample}}_R2.fastq.gz",
    log:
        f"{OUTDIR}/logs/extract/{{sample}}.log",
    conda:
        "envs/reassembly.yaml"
    threads: 4
    resources:
        mem_mb=8192,
        time_min=60,
    shell:
        r"""
        set -euo pipefail
        samtools view -b -N {input.read_list} {input.bam} 2> {log} \
          | samtools sort -n -@ {threads} - 2>> {log} \
          | samtools fastq -1 {output.r1} -2 {output.r2} -0 /dev/null -s /dev/null - 2>> {log}
        """


rule combine_reads:
    input:
        r1=expand(f"{OUTDIR}/extracted/{{sample}}_R1.fastq.gz", sample=SAMPLE_NAMES),
        r2=expand(f"{OUTDIR}/extracted/{{sample}}_R2.fastq.gz", sample=SAMPLE_NAMES),
    output:
        r1=f"{OUTDIR}/combined/all_reads_R1.fastq.gz",
        r2=f"{OUTDIR}/combined/all_reads_R2.fastq.gz",
    log:
        f"{OUTDIR}/logs/combine/combine_reads.log",
    conda:
        "envs/reassembly.yaml"
    threads: 1
    resources:
        mem_mb=2048,
        time_min=30,
    shell:
        r"""
        set -euo pipefail
        cat {input.r1} > {output.r1} 2> {log}
        cat {input.r2} > {output.r2} 2>> {log}
        """


rule spades:
    input:
        r1=f"{OUTDIR}/combined/all_reads_R1.fastq.gz",
        r2=f"{OUTDIR}/combined/all_reads_R2.fastq.gz",
        trusted=config["trusted_contigs"],
    output:
        contigs=f"{OUTDIR}/spades/contigs.fasta",
    params:
        outdir=lambda wc, output: os.path.dirname(output.contigs),
        mem_gb=config["spades"]["memory_gb"],
    log:
        f"{OUTDIR}/logs/spades/spades.log",
    conda:
        "envs/reassembly.yaml"
    threads: 16
    resources:
        mem_mb=config["spades"]["memory_gb"] * 1024,
        time_min=480,
    shell:
        r"""
        set -euo pipefail
        spades.py -o {params.outdir} \
                  -1 {input.r1} -2 {input.r2} \
                  --trusted-contigs {input.trusted} \
                  --only-assembler --careful \
                  -t {threads} -m {params.mem_gb} \
                  > {log} 2>&1
        """
