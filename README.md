# Snakemake PacBio Assembly & Reassembly Pipelines

[![Snakemake](https://img.shields.io/badge/Snakemake-≥7.0-brightgreen.svg?style=flat)](https://snakemake.readthedocs.io)
[![Conda](https://img.shields.io/badge/Conda-Bioconda-blue.svg)](https://bioconda.github.io/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

Modular, configuration-driven Snakemake workflows for *de novo* PacBio HiFi genome assembly, quality/contamination assessment, and targeted paired-end locus reassembly.

---

## Pipelines Overview

### 1. Primary Assembly Pipeline (`Pacbio-assembly-snakemake.smk`)
Configured via `config/assembly_config.yaml` and pinned Conda environment `envs/assembly.yaml`.

**Key Steps:**
1. **Quality Control (`fastqc`):** Evaluates raw HiFi FASTQ reads in isolated temp directories.
2. ***De Novo* Assembly (`hifiasm`):** Assembles PacBio HiFi reads with configurable purge level (`-l`).
3. **Contig Extraction (`gfa_to_fasta`):** Parses `.p_ctg.gfa` primary unitig/contig records cleanly into wrapped FASTA files.
4. **Coverage Mapping (`minimap2` + `samtools`):** Maps HiFi reads back onto assembled contigs (`minimap2 -ax map-hifi`) and produces coordinate-sorted, indexed BAM coverage tracks.
5. **Completeness Assessment (`busco`):** Evaluates gene completeness against the configured OrthoDB lineage (default: `insecta_odb10`).
6. **Taxonomic Profiling (`blastn`):** Queries contigs against a local BLAST database (`nt`).
7. **Contamination Screening (`blobtools`):** Integrates assembly contigs, coverage BAM, BLAST hits, NCBI `taxdump`, and BUSCO `full_table.tsv` into a BlobDir (`meta.json`).

### 2. Targeted Reassembly Pipeline (`reassembly_snakemake.smk`)
Configured via `config/reassembly_config.yaml` and pinned Conda environment `envs/reassembly.yaml`.

**Key Steps:**
1. **Reference Indexing (`bwa_index`):** Indexes the target FASTA reference via `bwa index`.
2. **Paired-End Alignment (`bwa_mem` + `samtools`):** Maps paired-end reads (`R1`, `R2`) using `bwa mem` and streams directly into `samtools sort` and `samtools index`.
3. **Target Read Extraction (`extract_reads`):** Extracts read pairs matching per-sample read ID lists cleanly using native `samtools view -N`, `samtools sort -n`, and `samtools fastq`.
4. **Cohort Read Pooling (`combine_reads`):** Concatenates extracted paired-end `fastq.gz` streams across samples.
5. **Targeted Assembly (`spades`):** Assembles pooled reads with `--trusted-contigs`, `--only-assembler`, and `--careful`.

---

## Repository Structure

```text
Snakemake-PacBio-Assembly/
├── Pacbio-assembly-snakemake.smk   # Primary PacBio HiFi assembly & BlobTools QC workflow
├── reassembly_snakemake.smk        # Targeted locus read extraction & SPAdes reassembly
├── config/
│   ├── assembly_config.yaml        # Sample paths, DBs, and tool parameters for assembly
│   └── reassembly_config.yaml      # Reference, samples, and SPAdes memory config
├── envs/
│   ├── assembly.yaml               # Pinned Bioconda dependencies for assembly workflow
│   └── reassembly.yaml             # Pinned Bioconda dependencies for reassembly workflow
└── .github/workflows/ci.yml        # Automated Snakemake --lint and DAG dry-run CI
```

---

## Usage

### 1. Primary Assembly Pipeline

1. Edit `config/assembly_config.yaml` with your sample FASTQ path(s), `busco.downloads_path`, BLAST `db_prefix`, and `taxdump` directory.
2. Validate the DAG and run with `--use-conda`:
   ```bash
   # Dry-run DAG validation
   snakemake -s Pacbio-assembly-snakemake.smk -n -p --cores 8

   # Execute with pinned Conda environment
   snakemake -s Pacbio-assembly-snakemake.smk --use-conda --cores 16
   ```

### 2. Targeted Reassembly Pipeline

1. Edit `config/reassembly_config.yaml` with your target reference FASTA, trusted contigs FASTA, and per-sample `r1`, `r2`, and `read_list` paths.
2. Validate the DAG and run with `--use-conda`:
   ```bash
   # Dry-run DAG validation
   snakemake -s reassembly_snakemake.smk -n -p --cores 8

   # Execute with pinned Conda environment
   snakemake -s reassembly_snakemake.smk --use-conda --cores 16
   ```

---

## Automated CI Validation

Every push and pull request runs GitHub Actions CI (`.github/workflows/ci.yml`) to verify both workflows against synthetic test fixtures using `snakemake --lint` and `snakemake -n -p --cores 8`.

---

## Author & Contact

**Nasir Mahmood Abbasi, PhD**
Computational Biologist | Open to New Opportunities
- Portfolio: [nmabbasi.github.io](https://nmabbasi.github.io)
- HPC Guide: [nmabbasi.github.io/HPC](https://nmabbasi.github.io/HPC/)
- Email: [nmabbasi@gmail.com](mailto:nmabbasi@gmail.com)
