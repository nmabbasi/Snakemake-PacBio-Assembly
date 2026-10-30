# Snakemake PacBio Assembly & Reassembly Pipelines 🐍

[![Snakemake](https://img.shields.io/badge/Snakemake-≥6.0.0-brightgreen.svg?style=flat)](https://snakemake.github.io)
[![Python](https://img.shields.io/badge/Python-3.x-blue.svg)](https://www.python.org/)
[![License](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

Automated Snakemake workflows for *de novo* genome assembly, quality evaluation, and targeted locus reassembly of PacBio HiFi sequencing reads.

## Pipelines Overview

### 1. Primary Assembly Pipeline (`Pacbio-assembly-snakemake.smk`)
An automated end-to-end workflow for *de novo* genome assembly of PacBio HiFi reads and subsequent taxonomic and completeness evaluation.

**Key Steps:**
1. **Quality Control:** Validates input FastQ reads using [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/).
2. **De Novo Assembly:** Assembles long reads using [Hifiasm](https://github.com/chhylp123/hifiasm), optimized for phased diploid assemblies.
3. **Contig Extraction:** Parses `.gfa` graphs to extract primary contigs (`.fasta`).
4. **Completeness Assessment:** Evaluates assembly completeness using [BUSCO](https://busco.ezlab.org/) against lineage databases (default: `insecta_odb10`).
5. **Taxonomic Profiling:** Queries contigs against NCBI `nt` using `blastn`.
6. **Contamination Screening:** Generates Taxon-Annotated GC-Coverage (TAGC) plots via [BlobTools](https://github.com/DRL/blobtools) by integrating assembly contigs, read depth BAM, BLAST hits, and BUSCO scores.

### 2. Targeted Reassembly Pipeline (`reassembly_snakemake.smk`)
A targeted workflow designed to extract subset reads aligning to a specific genomic region or reference and reassemble them using SPAdes.

**Key Steps:**
1. **Reference Indexing:** Indexes the target reference using [BWA](https://github.com/lh3/bwa).
2. **Alignment:** Maps paired-end reads (`P1`, `P2`) to reference regions with BWA-ALN / BWA-SAMPE.
3. **Read Extraction:** Extracts reads of interest from SAM alignments with a custom Perl script (`extract.fasta.from.sam.using.list.pl`).
4. **Targeted Assembly:** Concatenates extracted reads and reassembles them using [SPAdes](https://github.com/ablab/spades) (`--careful`, `--only-assembler`) incorporating trusted contigs.

---

## Installation & Requirements

Ensure you have Conda or Mamba installed on your HPC cluster or local machine.

### Dependencies
The pipelines require the following bioinformatics tools in your environment or `$PATH`:
- `fastqc`
- `hifiasm`
- `busco`
- `blastn`
- `blobtools`
- `bwa` (v0.6.2+)
- `spades.py` (v3.6.1+)
- `perl`

To install these tools into a dedicated Conda environment:
```bash
conda create -n pacbio_asm -c bioconda -c conda-forge \
  snakemake fastqc hifiasm busco blast blobtools bwa spades perl
conda activate pacbio_asm
```

---

## Usage

### Running the Primary Assembly Pipeline

1. Configure input variables at the top of `Pacbio-assembly-snakemake.smk`:
   ```python
   input_fastq = "your_reads.fastq.gz"
   threads = 92
   taxdump = "/path/to/taxdump"
   busco_lineage = "/path/to/lineage_odb10"
   ```
2. Execute the workflow:
   ```bash
   snakemake -s Pacbio-assembly-snakemake.smk --cores 92
   ```

### Running the Targeted Reassembly Pipeline

1. Define target FASTA references and read subset lists in `reassembly_snakemake.smk`.
2. Execute the workflow:
   ```bash
   snakemake -s reassembly_snakemake.smk --cores 30
   ```

---

## Author & Contact

**Nasir Mahmood Abbasi, PhD**  
Computational Biologist | Open to New Opportunities  
- Portfolio: [nmabbasi.github.io](https://nmabbasi.github.io)  
- HPC Guide: [nmabbasi.github.io/HPC](https://nmabbasi.github.io/HPC/)  
- Email: [nmabbasi@gmail.com](mailto:nmabbasi@gmail.com)
