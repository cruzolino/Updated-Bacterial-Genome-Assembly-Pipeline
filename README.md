# Bacterial Genome Assembly Pipeline

**Version 4.0.0** | Bash · Bioinformatics · Local HPC

A robust, single-script pipeline for bacterial genome assembly from short reads (Illumina PE), long reads (Nanopore / PacBio HiFi), or a combination of both. Covers the full workflow from raw read QC through assembly, polishing, quality assessment, and annotation, producing a structured output directory and a Markdown summary report.

---

## Table of Contents

- [Features](#features)
- [Supported Modes](#supported-modes)
- [Dependencies](#dependencies)
- [Installation](#installation)
- [Usage](#usage)
- [Options Reference](#options-reference)
- [Examples](#examples)
- [Output Structure](#output-structure)
- [Pipeline Stages](#pipeline-stages)
- [Resume Mode](#resume-mode)
- [Changelog (v3 → v4)](#changelog-v3--v4)
- [Suggested Next Steps](#suggested-next-steps)
- [Citation](#citation)

---

## Features

- **Four sequencing modes** — Illumina PE, Nanopore, PacBio HiFi, and Hybrid
- **Six assemblers** — SPAdes, SKESA, Unicycler, Flye, Canu, Raven
- **Resume mode** — checkpoint sentinels allow interrupted runs to restart from the last completed stage
- **Preflight tool check** — all required tools are validated together before the pipeline starts, not mid-run
- **Automatic coverage targeting** — Filtlong is configured dynamically for 20× coverage based on genome size
- **Dual annotation** — Bakta and Prokka both run when available (complementary, not redundant)
- **Structured logging** — every tool writes to a dedicated log file under `logs/`
- **Markdown summary report** — assembly statistics and output file index written to `05_reports/SUMMARY.md`

---

## Supported Modes

| Mode | Input flags | Default assembler |
|:-----|:------------|:------------------|
| Illumina paired-end | `-1 R1.fq.gz -2 R2.fq.gz` | `spades` |
| Nanopore long-read | `-l ont.fq.gz` | `flye` |
| PacBio HiFi | `-l hifi.fq.gz --hifi` | `flye` |
| Hybrid | `-1 R1.fq.gz -2 R2.fq.gz -l ont.fq.gz` | `unicycler` |

---

## Dependencies

### Required (always)

| Tool | Purpose | Install |
|:-----|:--------|:--------|
| `fastqc` | Raw and trimmed read QC | `conda install -c bioconda fastqc` |
| `fastp` | Illumina adapter trimming | `conda install -c bioconda fastp` |
| `quast` | Assembly quality metrics | `conda install -c bioconda quast` |

### Required per mode

| Tool | Mode | Install |
|:-----|:-----|:--------|
| `filtlong` | Long-read / Hybrid | `conda install -c bioconda filtlong` |
| `bowtie2` + `samtools` + `pilon` | Illumina / Hybrid polishing | `conda install -c bioconda bowtie2 samtools pilon` |
| `medaka_consensus` | Nanopore polishing | `pip install medaka` |
| `spades.py` | SPAdes assembly | `conda install -c bioconda spades` |
| `skesa` | SKESA assembly | `conda install -c bioconda skesa` |
| `unicycler` | Unicycler assembly | `conda install -c bioconda unicycler` |
| `flye` | Flye assembly | `conda install -c bioconda flye` |
| `canu` | Canu assembly | `conda install -c bioconda canu` |
| `raven` | Raven assembly | `conda install -c bioconda raven-assembler` |

### Optional (skipped gracefully if absent)

| Tool | Purpose |
|:-----|:--------|
| `NanoPlot` | Long-read QC plots |
| `multiqc` | Aggregated QC report |
| `seqkit` | Fast contig length filtering (recommended; awk fallback used otherwise) |
| `checkm2` / `checkm` | Assembly completeness and contamination |
| `busco` | Single-copy orthologue completeness |
| `plasmidfinder.py` | Plasmid replicon detection |
| `bakta` | Comprehensive bacterial annotation (INSDC-ready) |
| `prokka` | Fast legacy-compatible annotation (.sqn/.gbk) |
| `antismash` | Secondary metabolite biosynthetic gene cluster detection |
| `assembly-stats` | Per-assembly N50 / length statistics in the summary report |

> **Tip:** The pipeline fails fast with a full list of *all* missing required tools before any computation begins. Optional tools are silently skipped with a warning.

---

## Installation

```bash
# Clone or download the script
git clone https://github.com/your-org/bacterial-assembly-pipeline.git
cd bacterial-assembly-pipeline

# Make executable
chmod +x bacterial_assembly_v4.sh

# Create and activate a conda environment with core dependencies
conda create -n bac_assembly -c bioconda -c conda-forge \
    fastqc fastp quast filtlong bowtie2 samtools pilon \
    spades unicycler flye seqkit multiqc nanoplot assembly-stats
conda activate bac_assembly
```

---

## Usage

```bash
./bacterial_assembly_v4.sh [OPTIONS]
```

At least one of `-1` (Illumina R1) or `-l` (long reads) is required.

---

## Options Reference

### Input

| Flag | Description |
|:-----|:------------|
| `-1 FILE` | Illumina forward reads (R1.fastq.gz) |
| `-2 FILE` | Illumina reverse reads (R2.fastq.gz) |
| `-l FILE` | Long reads — Nanopore or PacBio HiFi (.fastq[.gz]) |
| `--hifi` | Treat `-l` reads as PacBio HiFi (activates `flye --pacbio-hifi`) |

### Assembly

| Flag | Description | Default |
|:-----|:------------|:--------|
| `-a STR` | Assembler: `spades` \| `skesa` \| `unicycler` \| `flye` \| `canu` \| `raven` | `spades` |

### General

| Flag | Description | Default |
|:-----|:------------|:--------|
| `-o DIR` | Output directory | `bacterial_assembly` |
| `-t INT` | CPU threads | auto-detected |
| `-m STR` | Memory limit (e.g. `32G`, `64G`) | `32G` |
| `-s STR` | Sample name (used in file prefixes) | `isolate` |
| `-g STR` | Expected genome size (e.g. `5m`, `4.8m`, `0.5g`) | `5m` |
| `-M STR` | Medaka model (Nanopore polishing) | `r941_min_hac_g507` |
| `-d PATH` | Bakta database path | auto-detect |
| `-c INT` | Minimum contig length (bp) | `500` |
| `-P STR` | PlasmidFinder database | `enterobacteriaceae` |
| `-x` | Skip plasmid detection | — |
| `-q` | QC only (stop before assembly) | — |
| `-r` | Resume from last completed checkpoint | — |
| `-h` / `--help` | Show help and exit | — |

---

## Examples

```bash
# Illumina paired-end — basic
./bacterial_assembly_v4.sh \
  -1 R1.fq.gz -2 R2.fq.gz \
  -s Ecoli_K12 -t 16 -m 32G

# Nanopore only — Flye assembler, 4.8 Mb genome
./bacterial_assembly_v4.sh \
  -l ont.fq.gz \
  -s Salmonella -a flye -g 4.8m -t 16

# Hybrid — Unicycler with Illumina + Nanopore
./bacterial_assembly_v4.sh \
  -1 R1.fq.gz -2 R2.fq.gz -l ont.fq.gz \
  -s Klebsiella -a unicycler -g 5.5m -t 24

# PacBio HiFi — Flye, 6.5 Mb genome
./bacterial_assembly_v4.sh \
  -l hifi.fastq.gz --hifi \
  -s Pseudomonas -a flye -g 6.5m -t 16

# QC only (no assembly)
./bacterial_assembly_v4.sh \
  -1 R1.fq.gz -2 R2.fq.gz -s Ecoli -q

# Resume an interrupted run
./bacterial_assembly_v4.sh \
  -1 R1.fq.gz -2 R2.fq.gz -s Ecoli -t 16 -r

# Custom output directory, 750 bp contig filter, skip plasmids
./bacterial_assembly_v4.sh \
  -1 R1.fq.gz -2 R2.fq.gz \
  -s Staph -o /data/results/staph_run1 \
  -c 750 -x -t 32
```

---

## Output Structure

```
bacterial_assembly/
├── 00_raw_data/              Symlinks to input files (no copy)
├── 01_qc/
│   ├── *_fastp.html/json     fastp trimming report
│   ├── *fastqc*.html/zip     FastQC reports (raw and trimmed)
│   ├── post_trim/            Post-trim FastQC output
│   └── nanoplot/             NanoPlot long-read QC (if available)
├── 02_assembly/
│   └── <assembler>/          Raw assembler output
├── 03_polishing/
│   ├── *_pilon.fasta         Pilon-polished assembly (Illumina/Hybrid)
│   └── medaka/consensus.fasta  Medaka-polished assembly (Nanopore)
├── 04_annotation/
│   ├── bakta/                Bakta annotation output
│   ├── prokka/               Prokka annotation output
│   └── antismash/            antiSMASH BGC output
├── 05_reports/
│   ├── quast/                QUAST assembly metrics
│   ├── checkm/               CheckM (legacy) completeness
│   ├── checkm2/              CheckM2 completeness
│   ├── busco/                BUSCO completeness
│   ├── plasmidfinder/        PlasmidFinder results
│   ├── *_multiqc.html        MultiQC aggregated report
│   └── SUMMARY.md            ← Final human-readable report
└── logs/
    ├── fastqc_raw.log
    ├── fastqc_post.log
    ├── fastp.log
    ├── filtlong.log
    ├── <assembler>.log
    ├── bowtie2.log / pilon.log / medaka.log
    ├── quast.log / checkm2.log / busco.log
    ├── bakta.log / prokka.log / antismash.log
    ├── .assembly_path        Internal: persists final assembly path
    └── .done_*               Internal: stage completion sentinels
```

---

## Pipeline Stages

```
[1/6] QC
        Illumina: FastQC (raw) → fastp trimming → FastQC (trimmed) → MultiQC
        Long-read: NanoPlot → Filtlong (target: 20× coverage)

[2/6] Assembly
        Illumina:  SPAdes --isolate | SKESA | Unicycler --conservative
        Long-read: Flye | Canu | Raven
        Hybrid:    Unicycler --normal | SPAdes
        → Contigs filtered to ≥ MIN_CONTIG_LENGTH bp (seqkit or awk fallback)

[3/6] Polishing
        Illumina/Hybrid: bowtie2 + samtools → Pilon --fix all
        Nanopore:        medaka_consensus (threads capped at 4)
        HiFi:            No polishing (HiFi reads are self-correcting)

[4/6] Quality Assessment
        QUAST → CheckM2 (or CheckM legacy) → BUSCO → PlasmidFinder

[5/6] Annotation
        Bakta (INSDC-ready) + Prokka (legacy .sqn/.gbk) + antiSMASH

[6/6] Report
        SUMMARY.md written to 05_reports/
```

---

## Resume Mode

Use `-r` to restart an interrupted pipeline from the last successfully completed stage without re-running earlier steps:

```bash
# Original run fails during annotation
./bacterial_assembly_v4.sh -1 R1.fq.gz -2 R2.fq.gz -s Ecoli -t 16

# Resume — QC, assembly, polishing, and QA are skipped automatically
./bacterial_assembly_v4.sh -1 R1.fq.gz -2 R2.fq.gz -s Ecoli -t 16 -r
```

Checkpoints are stored as hidden sentinel files under `logs/` (`.done_qc`, `.done_assembly`, etc.). To force a stage to re-run, delete its sentinel file or run without `-r`.

> **Note:** Resume mode requires that the output directory (`-o`) and all other flags match the original invocation exactly.

---

## Changelog (v3 → v4)

| # | Category | Change |
|:--|:---------|:-------|
| 1 | **Bug fix** | `awk` contig filter dropped the final FASTA record; fixed. `seqkit` is now the primary filter path; `awk` is an explicit fallback |
| 2 | **Bug fix** | Post-trim FastQC was overwriting `fastqc_raw.log`; separated into `fastqc_post.log` |
| 3 | **Bug fix** | Bakta DB auto-detection parsed formatted table output as a raw path; fixed with `head -1 \| tr -d '[:space:]'` |
| 4 | **Bug fix** | `bc` used for float arithmetic without guaranteed availability; replaced with portable `awk BEGIN` |
| 5 | **Bug fix** | `wtdbg2` accepted by argument parser but had no assembly case block; removed |
| 6 | **Feature** | Resume mode added (`-r` flag + `.done_*` sentinel files) |
| 7 | **Feature** | Preflight `check_tools()` — all required tools reported together at startup |
| 8 | **Feature** | `-P` flag to configure PlasmidFinder database (was hardcoded to `enterobacteriaceae`) |
| 9 | **Feature** | `--help` long-form alias added |
| 10 | **Optimization** | `assembly-stats` called once and cached; was forked three times in report generation |
| 11 | **Optimization** | Medaka threads capped at `min(THREADS, 4)` — tool does not scale beyond 4 CPU threads |
| 12 | **Optimization** | `gzip -1` (fastest) used for Filtlong output; intermediate files don't need archival compression |
| 13 | **Correctness** | Prokka now always runs when available, not only as a Bakta fallback |
| 14 | **Correctness** | Bakta `--db` flag converted from string variable to bash array — safe for paths with spaces |
| 15 | **Correctness** | ERR trap reset before all intentional `exit 0` calls to prevent spurious trap firing |
| 16 | **Correctness** | Pilon RAM guard warns when requested memory exceeds available system RAM |

---

## Suggested Next Steps

After the pipeline completes, consider these downstream analyses:

1. **AMR screening** — `AMRFinderPlus` or `ResFinder` for resistance gene detection
2. **NCBI submission** — use `tbl2asn` with the Prokka `.sqn` file
3. **MLST typing** — `mlst` (Torsten Seemann) for sequence typing
4. **Phylogenetics** — `IQ-TREE2` or `FastTree` for phylogenomic placement
5. **Pan-genome** — `Panaroo` (preferred) or `Roary` for comparative genomics

---

## Citation

If you use this pipeline in a publication, please cite the underlying tools used in your analysis. Key references:

- **SPAdes**: Bankevich et al. *J Comput Biol* 2012
- **Flye**: Kolmogorov et al. *Nature Methods* 2019
- **Unicycler**: Wick et al. *PLOS Computational Biology* 2017
- **Medaka**: Oxford Nanopore Technologies — https://github.com/nanoporetech/medaka
- **Pilon**: Walker et al. *PLOS ONE* 2014
- **Bakta**: Schwengers et al. *Microbial Genomics* 2021
- **Prokka**: Seemann *Bioinformatics* 2014
- **QUAST**: Gurevich et al. *Bioinformatics* 2013
- **CheckM2**: Chklovski et al. *Nature Methods* 2023
- **BUSCO**: Manni et al. *Molecular Biology and Evolution* 2021

---

> Pipeline maintained as a local bash script — no workflow manager required.
> Tested on Ubuntu 22.04 LTS with bash 5.1+.
