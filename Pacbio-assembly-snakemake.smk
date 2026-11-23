import os

configfile: "config/assembly_config.yaml"

SAMPLES = config["samples"]
SAMPLE_NAMES = list(SAMPLES.keys())
OUTDIR = config["outdir"]
LINEAGE = config["busco"]["lineage"]


rule all:
    input:
        expand(f"{OUTDIR}/qc/fastqc/{{sample}}_fastqc.html", sample=SAMPLE_NAMES),
        expand(f"{OUTDIR}/contigs/{{sample}}.p_ctg.fasta", sample=SAMPLE_NAMES),
        expand(f"{OUTDIR}/blobtools/{{sample}}/meta.json", sample=SAMPLE_NAMES),


rule fastqc:
    input:
        fastq=lambda wc: SAMPLES[wc.sample],
    output:
        html=f"{OUTDIR}/qc/fastqc/{{sample}}_fastqc.html",
        zip=f"{OUTDIR}/qc/fastqc/{{sample}}_fastqc.zip",
    log:
        f"{OUTDIR}/logs/fastqc/{{sample}}.log",
    conda:
        "envs/assembly.yaml"
    threads: 2
    resources:
        mem_mb=4096,
        time_min=60,
    shell:
        r"""
        set -euo pipefail
        TMPD=$(mktemp -d)
        trap 'rm -rf "$TMPD"' EXIT
        case "{input.fastq}" in
            *.gz) EXT=.fastq.gz ;;
            *)    EXT=.fastq ;;
        esac
        ln -s "$(realpath "{input.fastq}")" "$TMPD/{wildcards.sample}$EXT"
        fastqc -t {threads} --outdir "$TMPD" "$TMPD/{wildcards.sample}$EXT" > {log} 2>&1
        mv "$TMPD/{wildcards.sample}_fastqc.html" {output.html}
        mv "$TMPD/{wildcards.sample}_fastqc.zip" {output.zip}
        """


rule hifiasm:
    input:
        fastq=lambda wc: SAMPLES[wc.sample],
    output:
        gfa=f"{OUTDIR}/hifiasm/{{sample}}.p_ctg.gfa",
    params:
        prefix=lambda wc, output: output.gfa.removesuffix(".p_ctg.gfa"),
        purge_level=config["hifiasm"]["purge_level"],
    log:
        f"{OUTDIR}/logs/hifiasm/{{sample}}.log",
    conda:
        "envs/assembly.yaml"
    threads: 16
    resources:
        mem_mb=65536,
        time_min=480,
    shell:
        r"""
        set -euo pipefail
        hifiasm -l {params.purge_level} --primary -o {params.prefix} -t {threads} {input.fastq} > {log} 2>&1
        if [ ! -f "{output.gfa}" ] && [ -f "{params.prefix}.bp.p_ctg.gfa" ]; then
            mv "{params.prefix}.bp.p_ctg.gfa" "{output.gfa}"
        fi
        """


rule gfa_to_fasta:
    input:
        gfa=f"{OUTDIR}/hifiasm/{{sample}}.p_ctg.gfa",
    output:
        fasta=f"{OUTDIR}/contigs/{{sample}}.p_ctg.fasta",
    log:
        f"{OUTDIR}/logs/gfa_to_fasta/{{sample}}.log",
    conda:
        "envs/assembly.yaml"
    threads: 1
    resources:
        mem_mb=2048,
        time_min=30,
    shell:
        r"""
        set -euo pipefail
        awk '$1 == "S" {{printf ">%s\n%s\n", $2, $3}}' {input.gfa} | fold -w 80 > {output.fasta} 2> {log}
        """


rule map_reads_to_contigs:
    input:
        fasta=f"{OUTDIR}/contigs/{{sample}}.p_ctg.fasta",
        fastq=lambda wc: SAMPLES[wc.sample],
    output:
        bam=f"{OUTDIR}/mapping/{{sample}}.sorted.bam",
        bai=f"{OUTDIR}/mapping/{{sample}}.sorted.bam.bai",
    log:
        minimap2=f"{OUTDIR}/logs/mapping/{{sample}}_minimap2.log",
        samtools=f"{OUTDIR}/logs/mapping/{{sample}}_samtools.log",
    conda:
        "envs/assembly.yaml"
    threads: 16
    resources:
        mem_mb=32768,
        time_min=240,
    shell:
        r"""
        set -euo pipefail
        SORT_THREADS=$(( {threads} > 2 ? 2 : 1 ))
        ALIGN_THREADS=$(( {threads} > 2 ? {threads} - 2 : 1 ))
        minimap2 -ax map-hifi -t "$ALIGN_THREADS" {input.fasta} {input.fastq} 2> {log.minimap2} \
          | samtools sort -@ "$SORT_THREADS" -o {output.bam} - 2> {log.samtools}
        samtools index -@ "$SORT_THREADS" {output.bam} 2>> {log.samtools}
        """


rule busco:
    input:
        fasta=f"{OUTDIR}/contigs/{{sample}}.p_ctg.fasta",
    output:
        full_table=f"{OUTDIR}/busco/{{sample}}/run_{LINEAGE}/full_table.tsv",
    params:
        out_path=lambda wc, output: os.path.dirname(os.path.dirname(os.path.dirname(output.full_table))),
        lineage=LINEAGE,
        downloads_path=config["busco"]["downloads_path"],
    log:
        f"{OUTDIR}/logs/busco/{{sample}}.log",
    conda:
        "envs/assembly.yaml"
    threads: 16
    resources:
        mem_mb=32768,
        time_min=360,
    shell:
        r"""
        set -euo pipefail
        busco -i {input.fasta} \
              -o {wildcards.sample} \
              --out_path {params.out_path} \
              -l {params.lineage} \
              -m genome \
              -c {threads} \
              --download_path {params.downloads_path} \
              --force \
              > {log} 2>&1
        """


rule blastn:
    input:
        fasta=f"{OUTDIR}/contigs/{{sample}}.p_ctg.fasta",
    output:
        hits=f"{OUTDIR}/blast/{{sample}}_blastn.tsv",
    params:
        db=config["blast"]["db_prefix"],
        evalue=config["blast"]["evalue"],
        max_target_seqs=config["blast"]["max_target_seqs"],
    log:
        f"{OUTDIR}/logs/blast/{{sample}}.log",
    conda:
        "envs/assembly.yaml"
    threads: 16
    resources:
        mem_mb=32768,
        time_min=480,
    shell:
        r"""
        set -euo pipefail
        blastn -query {input.fasta} \
               -db {params.db} \
               -outfmt '6 qseqid staxids bitscore std' \
               -max_target_seqs {params.max_target_seqs} \
               -max_hsps 1 \
               -num_threads {threads} \
               -evalue {params.evalue} \
               -out {output.hits} \
               > {log} 2>&1
        """


rule blobtools:
    input:
        fasta=f"{OUTDIR}/contigs/{{sample}}.p_ctg.fasta",
        bam=f"{OUTDIR}/mapping/{{sample}}.sorted.bam",
        bai=f"{OUTDIR}/mapping/{{sample}}.sorted.bam.bai",
        hits=f"{OUTDIR}/blast/{{sample}}_blastn.tsv",
        busco=f"{OUTDIR}/busco/{{sample}}/run_{LINEAGE}/full_table.tsv",
    output:
        meta=f"{OUTDIR}/blobtools/{{sample}}/meta.json",
    params:
        blobdir=lambda wc, output: os.path.dirname(output.meta),
        taxdump=config["blobtools"]["taxdump"],
    log:
        f"{OUTDIR}/logs/blobtools/{{sample}}.log",
    conda:
        "envs/assembly.yaml"
    threads: 8
    resources:
        mem_mb=16384,
        time_min=180,
    shell:
        r"""
        set -euo pipefail
        blobtools create \
                  --fasta {input.fasta} \
                  --cov {input.bam} \
                  --hits {input.hits} \
                  --taxdump {params.taxdump} \
                  --busco {input.busco} \
                  --threads {threads} \
                  {params.blobdir} \
                  > {log} 2>&1
        """
