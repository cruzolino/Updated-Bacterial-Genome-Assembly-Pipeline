# Bacterial Genome Assembly Pipeline — v5.0.0

A multi-mode bacterial genome assembly pipeline supporting **Illumina paired-end**, **Nanopore**, **Hybrid (Illumina + Nanopore)**, and **PacBio HiFi** sequencing data. The pipeline runs five sequential stages — QC, assembly, polishing, quality assessment, and annotation — with checkpoint-based resume logic, preflight tool validation, and structured output organisation.

---

## Table of Contents

- [Features](#features)
- [Requirements](#requirements)
- [Installation](#installation)
- [Usage](#usage)
- [Pipeline Stages](#pipeline-stages)
- [Output Structure](#output-structure)
- [Examples](#examples)
- [Resume Mode](#resume-mode)
- [Changelog](#changelog)

---

## Features

- Four sequencing modes: Illumina PE, Nanopore, Hybrid, PacBio HiFi
- Six assembler options: SPAdes, SKESA, Unicycler, Flye, Canu, Raven
- Checkpoint/resume logic — interrupted runs can be restarted from the last completed stage
- Preflight tool check — all required tools are validated at startup before any processing begins
- Pilon polishing for Illumina and hybrid runs (long-read assemblers handle internal polishing)
- Bakta annotation for all modes, with optional antiSMASH secondary metabolite detection
- Genome completeness assessment via CheckM2 (or CheckM legacy fallback) and BUSCO
- Plasmid detection via PlasmidFinder (optional, enabled by default)
- Automatic Markdown summary report written to `05_reports/SUMMARY.md`
- Colour-coded, timestamped logging with clear error messages

---

## Requirements

### Mandatory (all modes)

| Tool | Purpose |
|---|---|
| `fastqc` | Raw and post-trim read QC (Illumina) |
| `fastp` | Adapter trimming and read filtering (Illumina) |
| `quast` | Assembly quality statistics |

### Illumina / Hybrid modes

| Tool | Purpose |
|---|---|
| `bowtie2` | Read mapping for Pilon polishing |
| `samtools` | BAM sorting and indexing |
| `pilon` | Illumina-based consensus polishing |

### Long-read / Hybrid modes

| Tool | Purpose |
|---|---|
| `filtlong` | Long-read quality filtering |
| `NanoPlot` *(optional)* | Long-read QC visualisation |

### Assemblers (install the one(s) you need)

| Tool | Modes |
|---|---|
| `spades.py` | Illumina, Hybrid |
| `skesa` | Illumina |
| `unicycler` | Illumina, Hybrid |
| `flye` | Nanopore, PacBio HiFi |
| `canu` | Nanopore, PacBio HiFi |
| `raven` | Nanopore |

### Quality assessment (optional but recommended)

| Tool | Purpose |
|---|---|
| `checkm2` | Genome completeness (preferred) |
| `checkm` | Genome completeness (legacy fallback) |
| `busco` | BUSCO completeness against `bacteria_odb10` |
| `assembly-stats` | Summary N50/length statistics |
| `multiqc` *(optional)* | Aggregate QC report |
| `seqkit` *(optional)* | Contig length filtering (awk fallback if absent) |

### Annotation (optional)

| Tool | Purpose |
|---|---|
| `bakta` | Primary genome annotator — produces INSDC-ready GFF3, GBFF, FAA, FNA, TSV |
| `antismash` *(optional)* | Secondary metabolite biosynthetic gene cluster detection |
| `plasmidfinder.py` *(optional)* | Plasmid replicon typing |

> **Bakta database:** Bakta requires a local sequence database. Download it with `bakta_db download` or specify a custom path with `-d`. The pipeline auto-detects the database if it has been set up via `bakta_db`.

---

## Installation

No installation is required for the pipeline script itself. Clone or download the script and make it executable:

```bash
chmod +x bacterial_assembly_v5.0.0.sh
```

All dependencies must be available on your `PATH`. The recommended approach is to manage them in a conda environment:

```bash
conda create -n bacterial_assembly -c bioconda -c conda-forge \
    fastqc fastp quast bowtie2 samtools pilon filtlong nanoplot \
    spades unicycler flye canu raven checkm2 busco bakta seqkit \
    assembly-stats multiqc plasmidfinder
conda activate bacterial_assembly
```

Install antiSMASH separately if needed, following the [antiSMASH installation guide](https://docs.antismash.secondarymetabolites.org/install/).

---

## Usage

```
Usage:  bacterial_assembly_v5.0.0.sh [OPTIONS]

Input (at least one required):
  -1 FILE   Illumina forward reads (R1.fastq.gz)
  -2 FILE   Illumina reverse reads (R2.fastq.gz)
  -l FILE   Long reads – Nanopore or PacBio HiFi (.fastq[.gz])
  --hifi    Treat -l reads as PacBio HiFi (activates flye --pacbio-hifi)

Assembly:
  -a STR    Assembler: spades|skesa|unicycler (Illumina/hybrid)
                       flye|canu|raven (long-read)        [default: spades]

General:
  -o DIR    Output directory          [default: bacterial_assembly]
  -t INT    Threads                   [default: auto-detected]
  -m STR    Memory limit              [default: 32G]
  -s STR    Sample name               [default: isolate]
  -g STR    Expected genome size      [default: 5m]
  -d PATH   Bakta DB path             [default: auto-detect]
  -c INT    Min contig length (bp)    [default: 500]
  -P STR    PlasmidFinder DB          [default: enterobacteriaceae]
  -x        Skip plasmid detection
  -q        QC only (skip assembly)
  -r        Resume from last checkpoint
  -h|--help Show this help
```

Genome size (`-g`) accepts a number followed by `m` (megabases) or `g` (gigabases), e.g. `5m`, `4.8m`, `0.5g`. This value is used for Filtlong coverage targeting and is passed to Flye/Canu when applicable.

---

## Pipeline Stages

### Stage 1/5 — Quality Control

**Illumina and Hybrid:**
1. FastQC on raw reads → `01_qc/`
2. Adapter trimming and quality filtering with fastp (Phred ≥ 20, length ≥ 50 bp) → `01_qc/`
3. FastQC on trimmed reads → `01_qc/post_trim/`
4. MultiQC aggregate report (if available) → `05_reports/`

**Nanopore and Hybrid:**
1. NanoPlot QC on raw long reads (if available) → `01_qc/nanoplot/`
2. Filtlong filtering: minimum length 1,000 bp, minimum mean Q 7, target bases at 20× coverage → `01_qc/`

---

### Stage 2/5 — Assembly

Assembler is selected with `-a`. Compatibility rules are enforced at startup:

| Assembler | Compatible modes |
|---|---|
| `spades` | Illumina, Hybrid |
| `skesa` | Illumina |
| `unicycler` | Illumina, Hybrid |
| `flye` | Nanopore, PacBio HiFi |
| `canu` | Nanopore, PacBio HiFi |
| `raven` | Nanopore |

After assembly, contigs shorter than the minimum length (`-c`, default 500 bp) are filtered out using `seqkit seq --min-len` (or an `awk` fallback). The filtered assembly FASTA path is persisted to `logs/.assembly_path` for use by downstream stages and resume mode.

---

### Stage 3/5 — Polishing

**Illumina / Hybrid:** Pilon polishing using bowtie2-mapped reads.
- Reads are mapped with `bowtie2 --very-sensitive-local`
- BAM is sorted and indexed with samtools
- Pilon runs with `--fix all` to correct SNPs, small indels, and local misassemblies
- A RAM guard warns if the specified memory limit exceeds available system RAM

**Nanopore / PacBio HiFi:** No external polishing is performed. Flye, Canu, and Raven each include multiple internal polishing rounds during assembly, making a separate polishing step unnecessary.

---

### Stage 4/5 — Quality Assessment

Runs in order; tools that are not installed are skipped with a warning rather than aborting the pipeline.

| Tool | Output |
|---|---|
| QUAST | Assembly contiguity statistics → `05_reports/quast/` |
| CheckM2 *(preferred)* | Completeness and contamination → `05_reports/checkm2/` |
| CheckM *(legacy fallback)* | Completeness and contamination → `05_reports/checkm/` |
| BUSCO | Completeness against `bacteria_odb10` → `05_reports/busco/` |
| PlasmidFinder | Replicon typing → `05_reports/plasmidfinder/` (disable with `-x`) |

---

### Stage 5/5 — Annotation

**Bakta** is the sole annotator. It produces INSDC-ready outputs suitable for direct NCBI/ENA submission.

| Output file | Description |
|---|---|
| `*.gff3` | Genome annotation in GFF3 format |
| `*.gbff` | GenBank flat file |
| `*.faa` | Annotated protein sequences |
| `*.fna` | Annotated nucleotide sequences |
| `*.tsv` | Tabular summary of all features |

If Bakta is not found in `PATH`, the stage issues a warning and continues. If Bakta fails (e.g. missing or misconfigured database), a warning is logged and the pipeline continues to report generation — it does not abort.

**antiSMASH** runs after Bakta if available, using `prodigal` for gene finding (`--minimal` mode). Results are written to `04_annotation/antismash/`.

A final Markdown report is written to `05_reports/SUMMARY.md` with assembly statistics, input file paths, and suggested next steps.

---

## Output Structure

```
<outdir>/
├── 00_raw_data/            # Symlinks to input read files (no copy)
├── 01_qc/
│   ├── post_trim/          # Post-trim FastQC reports
│   ├── nanoplot/           # NanoPlot long-read QC (if applicable)
│   ├── <sample>_R1.trimmed.fastq.gz
│   ├── <sample>_R2.trimmed.fastq.gz
│   ├── <sample>_fastp.html
│   └── <sample>_fastp.json
├── 02_assembly/
│   └── <assembler>/        # Raw assembler output
├── 03_polishing/
│   ├── <sample>_pilon.fasta   # Pilon-polished assembly (Illumina/hybrid)
│   └── mapped.bam             # Bowtie2 alignment used for polishing
├── 04_annotation/
│   ├── bakta/              # Bakta outputs (GFF3, GBFF, FAA, FNA, TSV)
│   └── antismash/          # antiSMASH BGC predictions (if available)
├── 05_reports/
│   ├── quast/              # QUAST assembly statistics
│   ├── checkm2/            # CheckM2 completeness report
│   ├── busco/              # BUSCO completeness report
│   ├── plasmidfinder/      # PlasmidFinder replicon results
│   ├── <sample>_multiqc/   # MultiQC aggregate report (if available)
│   └── SUMMARY.md          # Final pipeline summary report
└── logs/
    ├── fastqc_raw.log
    ├── fastqc_post.log
    ├── fastp.log
    ├── filtlong.log
    ├── <assembler>.log
    ├── bowtie2.log
    ├── pilon.log
    ├── quast.log
    ├── checkm2.log
    ├── busco.log
    ├── bakta.log
    ├── antismash.log
    ├── .assembly_path      # Persisted assembly FASTA path for resume
    └── .done_<stage>       # Checkpoint sentinel files for resume mode
```

---

## Examples

**Illumina paired-end (default SPAdes assembler):**
```bash
./bacterial_assembly_v5.0.0.sh \
    -1 sample_R1.fq.gz -2 sample_R2.fq.gz \
    -s Ecoli -t 16 -m 32G -o results/ecoli
```

**Nanopore with Flye:**
```bash
./bacterial_assembly_v5.0.0.sh \
    -l ont_reads.fq.gz \
    -s Salmonella -a flye -g 4.8m -t 16 -o results/salmonella
```

**Hybrid assembly with Unicycler:**
```bash
./bacterial_assembly_v5.0.0.sh \
    -1 R1.fq.gz -2 R2.fq.gz -l ont.fq.gz \
    -s Klebsiella -a unicycler -t 24 -o results/klebsiella
```

**PacBio HiFi with Flye:**
```bash
./bacterial_assembly_v5.0.0.sh \
    -l hifi_reads.fastq.gz --hifi \
    -s Pseudomonas -a flye -g 6.5m -t 16 -o results/pseudomonas
```

**QC only (skip assembly and downstream steps):**
```bash
./bacterial_assembly_v5.0.0.sh \
    -1 R1.fq.gz -2 R2.fq.gz -s Ecoli -q
```

**Custom Bakta database and SKESA assembler:**
```bash
./bacterial_assembly_v5.0.0.sh \
    -1 R1.fq.gz -2 R2.fq.gz \
    -s Staphylococcus -a skesa -d /db/bakta/db -t 8
```

**Disable plasmid detection:**
```bash
./bacterial_assembly_v5.0.0.sh \
    -1 R1.fq.gz -2 R2.fq.gz -s Ecoli -x
```

---

## Resume Mode

Each stage writes a hidden sentinel file (`logs/.done_<stage>`) upon successful completion. If the pipeline is interrupted (e.g. due to a system timeout or resource limit), it can be restarted from the last completed stage by passing the `-r` flag along with the original arguments:

```bash
# Original run (was interrupted during polishing)
./bacterial_assembly_v5.0.0.sh -1 R1.fq.gz -2 R2.fq.gz -s Ecoli -t 16

# Resume from where it stopped
./bacterial_assembly_v5.0.0.sh -1 R1.fq.gz -2 R2.fq.gz -s Ecoli -t 16 -r
```

Stages already marked as done are skipped entirely. To force a stage to re-run, remove the corresponding `.done_<stage>` sentinel file from `logs/` before resuming, or run without `-r`.

---

## Changelog

### v5.0.0
- **Medaka removed.** External Medaka polishing for Nanopore runs has been removed. Flye, Raven, and Canu perform multiple internal polishing iterations during assembly, making a separate Medaka step redundant. Pilon polishing for Illumina and hybrid runs is unchanged.
- **Prokka removed.** Prokka has been removed as an annotator due to silent failures and environment conflicts. Bakta is now the sole annotator, providing comprehensive INSDC-ready outputs (GFF3, GBFF, FAA, FNA, TSV).
- `-M` flag (Medaka model) removed from CLI.
- `04_annotation/prokka/` removed from directory scaffold.
- `medaka_consensus` removed from preflight tool checks.
- Step counter updated from 6 to 5 throughout all log messages.
- Summary report updated: annotation output path now points specifically to `bakta/`; NCBI submission note updated to reference `table2asn` with Bakta outputs.

### v4.0.0
- Resume logic via step-sentinel files
- Preflight tool check reports all missing tools at once before any processing
- `assembly-stats` output cached to avoid repeated subprocess calls
- `awk` contig-filter bug fixed (final record was previously dropped)
- `seqkit` is now the primary contig-filter path; `awk` is the fallback
- `fastqc_raw.log` renamed to `fastqc_post.log` (was overwriting the raw log)
- Pilon: RAM guard added; warns when heap exceeds available system RAM
- Bakta DB auto-detection fixed (`bakta_db list` output format)
- Strict quoting around all variable expansions
- `GENOME_SIZE` validated as `<number>[m|g]` before numeric arithmetic
- ERR trap reset before clean exit so it does not fire on intentional exits
- `log()` / `ok()` write to stdout; `warn()` / `err()` write to stderr
- PlasmidFinder database flag made configurable (`-P`)
- `--help` added as an alias for `-h`
- `-r` flag added to enable resume mode

