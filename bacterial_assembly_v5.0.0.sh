#!/usr/bin/env bash
# ==============================================================================
# Bacterial Genome Assembly Pipeline — v5.0.0
# Modes: Illumina (PE) | Nanopore | Hybrid (Illumina + Nanopore) | PacBio HiFi
#
# Changes from v4.0.0:
#   - Medaka removed: step 3 (polishing) no longer runs Medaka for Nanopore
#     runs. Flye/Raven/Canu perform their own internal consensus polishing,
#     making a separate Medaka pass unnecessary and a source of pipeline
#     crashes. Pilon (Illumina/hybrid) is retained unchanged.
#   - Prokka removed: annotation step now uses Bakta exclusively. Prokka was
#     causing silent failures and environment conflicts; Bakta covers all
#     required outputs (INSDC-ready GFF3/GBFF/FAA/FNA/TSV).
#   - -M flag (Medaka model) removed from CLI; MEDAKA_MODEL variable dropped.
#   - 04_annotation/prokka/ directory scaffold entry removed.
#   - Preflight tool-check no longer includes medaka_consensus.
#   - Step counter reduced from 6 to 5 throughout all log messages.
#   - Summary report "Next Steps" section updated (removed Prokka .sqn note).
#   - All other v4.0.0 fixes and optimisations retained.
# ==============================================================================
set -euo pipefail
IFS=$'\n\t'

# ── Version ────────────────────────────────────────────────────────────────────
readonly VERSION="5.0.0"

# ── Defaults ───────────────────────────────────────────────────────────────────
THREADS=$(nproc --all 2>/dev/null || echo 8)
MEMORY="32G"
OUTDIR="bacterial_assembly"
ASSEMBLER="spades"
READ_TYPE="illumina"
MIN_CONTIG_LENGTH=500
GENOME_SIZE="5m"
SAMPLE="isolate"
PLASMID_DETECTION=true
PLASMID_DB="enterobacteriaceae"   # -P flag
QC_ONLY=false
HIFI_MODE=false
RESUME=false
BAKTA_DB=""

# ── Colors ─────────────────────────────────────────────────────────────────────
readonly RED='\033[0;31m' GREEN='\033[0;32m' YELLOW='\033[1;33m' \
         BLUE='\033[0;34m' CYAN='\033[0;36m' NC='\033[0m'

# ── Logging ────────────────────────────────────────────────────────────────────
# All log/ok go to stdout; warn/err go to stderr
log()  { echo -e "${BLUE}[$(date +%H:%M:%S)]${NC} $*"; }
ok()   { echo -e "${GREEN}[✓ $(date +%H:%M:%S)]${NC} $*"; }
warn() { echo -e "${YELLOW}[WARN]${NC} $*" >&2; }
err()  { echo -e "${RED}[ERROR]${NC} $*" >&2; exit 1; }

# ── Sentinel helpers (resume logic) ───────────────────────────────────────────
# Each pipeline stage writes a .done file when it finishes successfully.
# If -r (resume) is set, a stage whose sentinel exists is skipped entirely.
sentinel_path() { echo "${OUTDIR}/logs/.done_${1}"; }

stage_done() {
    local stage="$1"
    touch "$(sentinel_path "${stage}")"
}

skip_if_done() {
    local stage="$1"
    if ${RESUME} && [[ -f "$(sentinel_path "${stage}")" ]]; then
        ok "Skipping ${stage} (already completed; use without -r to re-run)"
        return 0
    fi
    return 1
}

# ── Tool requirement check ─────────────────────────────────────────────────────
# Called once at startup with the full list of required tools.
# Collects ALL missing tools and reports them together — avoids "fix one, run,
# find the next" cycle that makes pipelines painful to set up.
check_tools() {
    local missing=()
    for t in "$@"; do
        command -v "${t}" &>/dev/null || missing+=("${t}")
    done
    if (( ${#missing[@]} > 0 )); then
        err "Missing required tools: ${missing[*]}. Install them and retry."
    fi
}

# ── Usage ─────────────────────────────────────────────────────────────────────
usage() {
    cat <<EOF
${GREEN}Bacterial Genome Assembly Pipeline v${VERSION}${NC}

${YELLOW}Usage:${NC}  $0 [OPTIONS]

${YELLOW}Input (at least one required):${NC}
  -1 FILE   Illumina forward reads (R1.fastq.gz)
  -2 FILE   Illumina reverse reads (R2.fastq.gz)
  -l FILE   Long reads – Nanopore or PacBio HiFi (.fastq[.gz])
  --hifi    Treat -l reads as PacBio HiFi (activates flye --pacbio-hifi)

${YELLOW}Assembly:${NC}
  -a STR    Assembler: spades|skesa|unicycler (Illumina/hybrid)
                       flye|canu|raven (long-read)        [default: spades]

${YELLOW}General:${NC}
  -o DIR    Output directory          [default: ${OUTDIR}]
  -t INT    Threads                   [default: auto = ${THREADS}]
  -m STR    Memory limit              [default: ${MEMORY}]
  -s STR    Sample name               [default: ${SAMPLE}]
  -g STR    Expected genome size      [default: ${GENOME_SIZE}]
  -d PATH   Bakta DB path             [default: auto-detect]
  -c INT    Min contig length (bp)    [default: ${MIN_CONTIG_LENGTH}]
  -P STR    PlasmidFinder DB          [default: ${PLASMID_DB}]
  -x        Skip plasmid detection
  -q        QC only (skip assembly)
  -r        Resume from last checkpoint
  -h|--help Show this help

${YELLOW}Examples:${NC}
  # Illumina PE
  $0 -1 R1.fq.gz -2 R2.fq.gz -s Ecoli -t 16

  # Nanopore
  $0 -l ont.fq.gz -s Salmonella -a flye -g 4.8m -t 16

  # Hybrid
  $0 -1 R1.fq.gz -2 R2.fq.gz -l ont.fq.gz -s Klebsiella -a unicycler

  # PacBio HiFi
  $0 -l hifi.fastq.gz -s Pseudomonas -a flye --hifi -g 6.5m

  # Resume interrupted run
  $0 -1 R1.fq.gz -2 R2.fq.gz -s Ecoli -t 16 -r
EOF
    exit 0
}

# ── Argument parsing ───────────────────────────────────────────────────────────
# Strip --hifi and --help before getopts (long opts not supported natively)
ARGS=()
for arg in "$@"; do
    case "${arg}" in
        --hifi)   HIFI_MODE=true ;;
        --help)   usage ;;
        *)        ARGS+=("${arg}") ;;
    esac
done
set -- "${ARGS[@]+"${ARGS[@]}"}"

while getopts "1:2:l:o:t:s:g:a:m:d:c:P:xqrh" opt; do
    case "${opt}" in
        1) R1="${OPTARG}" ;;
        2) R2="${OPTARG}" ;;
        l) LONG="${OPTARG}" ;;
        o) OUTDIR="${OPTARG}" ;;
        t) THREADS="${OPTARG}" ;;
        s) SAMPLE="${OPTARG}" ;;
        g) GENOME_SIZE="${OPTARG}" ;;
        a) ASSEMBLER="${OPTARG}" ;;
        m) MEMORY="${OPTARG}" ;;
        d) BAKTA_DB="${OPTARG}" ;;
        c) MIN_CONTIG_LENGTH="${OPTARG}" ;;
        P) PLASMID_DB="${OPTARG}" ;;
        x) PLASMID_DETECTION=false ;;
        q) QC_ONLY=true ;;
        r) RESUME=true ;;
        h) usage ;;
        *) err "Invalid option: -${OPTARG}" ;;
    esac
done

# ── Validate inputs ────────────────────────────────────────────────────────────
validate_inputs() {
    log "Validating inputs and environment..."

    # Read presence
    [[ -z "${R1:-}" && -z "${LONG:-}" ]] && err "No input reads provided. Use -h for help."

    # PE read pairing
    if { [[ -n "${R1:-}" ]] && [[ -z "${R2:-}" ]]; } ||
       { [[ -z "${R1:-}" ]] && [[ -n "${R2:-}" ]]; }; then
        err "Both -1 and -2 are required for paired-end Illumina reads."
    fi

    # Determine mode
    if [[ -n "${LONG:-}" && -n "${R1:-}" ]]; then
        READ_TYPE="hybrid"
        log "Mode: ${CYAN}Hybrid${NC} (Illumina + Long-read)"
    elif [[ -n "${LONG:-}" ]]; then
        READ_TYPE="long"
        log "Mode: ${CYAN}Long-read${NC}$(${HIFI_MODE} && echo ' [PacBio HiFi]' || echo ' [Nanopore]')"
    else
        READ_TYPE="illumina"
        log "Mode: ${CYAN}Illumina paired-end${NC}"
    fi

    # File existence
    [[ -n "${R1:-}"   && ! -f "${R1}"   ]] && err "R1 file not found: ${R1}"
    [[ -n "${R2:-}"   && ! -f "${R2}"   ]] && err "R2 file not found: ${R2}"
    [[ -n "${LONG:-}" && ! -f "${LONG}" ]] && err "Long-read file not found: ${LONG}"

    # Assembler–mode compatibility
    case "${ASSEMBLER}" in
        spades|skesa|unicycler)
            [[ "${READ_TYPE}" == "long" ]] && \
                err "${ASSEMBLER} cannot assemble long-reads only. Use flye, canu, or raven." ;;
        flye|canu|raven)
            [[ "${READ_TYPE}" == "illumina" ]] && \
                err "${ASSEMBLER} is a long-read assembler. For Illumina use spades, skesa, or unicycler." ;;
        *) err "Unknown assembler: ${ASSEMBLER}" ;;
    esac

    # Numeric validation
    [[ "${THREADS}" =~ ^[0-9]+$ ]] || err "Threads must be a positive integer."
    [[ "${MIN_CONTIG_LENGTH}" =~ ^[0-9]+$ ]] || err "Min contig length must be a positive integer."

    # Genome size format check — used later for math
    [[ "${GENOME_SIZE,,}" =~ ^[0-9.]+[mg]$ ]] || \
        err "Genome size must be a number followed by m or g (e.g. 5m, 4.8m, 0.5g)."

    # ── Preflight tool checks (all at once) ───────────────────────────────────
    local required_tools=(fastqc fastp quast)
    [[ "${READ_TYPE}" == "long" || "${READ_TYPE}" == "hybrid" ]] && \
        required_tools+=(filtlong)
    case "${ASSEMBLER}" in
        spades)    required_tools+=(spades.py) ;;
        skesa)     required_tools+=(skesa) ;;
        unicycler) required_tools+=(unicycler) ;;
        flye)      required_tools+=(flye) ;;
        canu)      required_tools+=(canu) ;;
        raven)     required_tools+=(raven) ;;
    esac
    if [[ "${READ_TYPE}" == "illumina" || "${READ_TYPE}" == "hybrid" ]]; then
        required_tools+=(bowtie2 samtools pilon)
    fi
    check_tools "${required_tools[@]}"

    # ── Directory scaffold ────────────────────────────────────────────────────
    mkdir -p "${OUTDIR}"/{00_raw_data,01_qc/post_trim,02_assembly,\
03_polishing,04_annotation/{bakta,antismash},\
05_reports/{quast,checkm,busco,plasmidfinder},logs}

    # Symlink raw data (no copy)
    [[ -n "${R1:-}"   ]] && ln -sf "$(realpath "${R1}")"   "${OUTDIR}/00_raw_data/$(basename "${R1}")"
    [[ -n "${R2:-}"   ]] && ln -sf "$(realpath "${R2}")"   "${OUTDIR}/00_raw_data/$(basename "${R2}")"
    [[ -n "${LONG:-}" ]] && ln -sf "$(realpath "${LONG}")" "${OUTDIR}/00_raw_data/$(basename "${LONG}")"

    ok "Input validation passed"
}

# ── QC ─────────────────────────────────────────────────────────────────────────
run_qc() {
    skip_if_done "qc" && return

    log "[1/5] Quality control..."

    # ── Illumina QC ──────────────────────────────────────────────────────────
    if [[ "${READ_TYPE}" == "illumina" || "${READ_TYPE}" == "hybrid" ]]; then

        log "FastQC on raw reads..."
        fastqc "${R1}" "${R2}" \
            --outdir "${OUTDIR}/01_qc/" \
            --threads "${THREADS}" \
            --quiet \
            > "${OUTDIR}/logs/fastqc_raw.log" 2>&1

        log "Trimming with fastp..."
        TRIMMED_R1="${OUTDIR}/01_qc/${SAMPLE}_R1.trimmed.fastq.gz"
        TRIMMED_R2="${OUTDIR}/01_qc/${SAMPLE}_R2.trimmed.fastq.gz"
        fastp \
            -i "${R1}" -I "${R2}" \
            -o "${TRIMMED_R1}" -O "${TRIMMED_R2}" \
            --thread "${THREADS}" \
            --qualified_quality_phred 20 \
            --unqualified_percent_limit 40 \
            --length_required 50 \
            --detect_adapter_for_pe \
            --correction \
            --overrepresentation_analysis \
            --html "${OUTDIR}/01_qc/${SAMPLE}_fastp.html" \
            --json "${OUTDIR}/01_qc/${SAMPLE}_fastp.json" \
            --report_title "${SAMPLE}" \
            > "${OUTDIR}/logs/fastp.log" 2>&1

        log "FastQC on trimmed reads..."
        # FIX: use a separate log file so raw FastQC log is not clobbered
        fastqc "${TRIMMED_R1}" "${TRIMMED_R2}" \
            --outdir "${OUTDIR}/01_qc/post_trim/" \
            --threads "${THREADS}" \
            --quiet \
            > "${OUTDIR}/logs/fastqc_post.log" 2>&1

        if command -v multiqc &>/dev/null; then
            multiqc "${OUTDIR}/01_qc/" \
                --outdir "${OUTDIR}/05_reports/" \
                --filename "${SAMPLE}_multiqc" \
                --quiet \
                > "${OUTDIR}/logs/multiqc.log" 2>&1
        else
            warn "multiqc not found; skipping MultiQC report."
        fi
    fi

    # ── Long-read QC ──────────────────────────────────────────────────────────
    if [[ "${READ_TYPE}" == "long" || "${READ_TYPE}" == "hybrid" ]]; then

        if command -v NanoPlot &>/dev/null; then
            log "NanoPlot QC on raw long reads..."
            NanoPlot \
                --fastq "${LONG}" \
                --outdir "${OUTDIR}/01_qc/nanoplot/" \
                --threads "${THREADS}" \
                --plots dot kde \
                > "${OUTDIR}/logs/nanoplot.log" 2>&1
        else
            warn "NanoPlot not found; skipping long-read QC plots."
        fi

        # Compute target bases for 20× coverage
        # FIX: validate GENOME_SIZE format is guaranteed above; safe to proceed
        local _gs="${GENOME_SIZE,,}"
        local TARGET_BASES
        if [[ "${_gs}" =~ ^([0-9.]+)m$ ]]; then
            TARGET_BASES=$(awk "BEGIN{printf \"%d\", ${BASH_REMATCH[1]} * 1000000 * 20}")
        else
            TARGET_BASES=$(awk "BEGIN{printf \"%d\", ${BASH_REMATCH[1]} * 1000000000 * 20}")
        fi
        # FIX: use awk for float math instead of bc (bc may not be installed)

        log "Filtering long reads with Filtlong (target: ${TARGET_BASES} bases)..."
        FILTERED_LONG="${OUTDIR}/01_qc/${SAMPLE}_filtered.fastq.gz"
        filtlong \
            --min_length 1000 \
            --min_mean_q 7 \
            --target_bases "${TARGET_BASES}" \
            "${LONG}" \
            2> "${OUTDIR}/logs/filtlong.log" \
        | gzip -1 -c > "${FILTERED_LONG}"
        # FIX: gzip -1 (fastest compression) — filtlong output is already selected;
        #      speed matters more than compression ratio here
    fi

    ok "QC complete"
    stage_done "qc"
}

# ── Assembly ───────────────────────────────────────────────────────────────────
run_assembly() {
    skip_if_done "assembly" && {
        # Restore ASSEMBLY path so downstream stages can find it
        _restore_assembly_path
        return
    }

    log "[2/5] Assembly with ${ASSEMBLER}..."

    local MEM_INT="${MEMORY//[^0-9]/}"

    case "${ASSEMBLER}" in

        spades)
            spades.py \
                -1 "${TRIMMED_R1}" -2 "${TRIMMED_R2}" \
                -o "${OUTDIR}/02_assembly/spades/" \
                -t "${THREADS}" -m "${MEM_INT}" \
                --isolate \
                --cov-cutoff auto \
                > "${OUTDIR}/logs/spades.log" 2>&1
            ASSEMBLY="${OUTDIR}/02_assembly/spades/contigs.fasta"
            ;;

        skesa)
            mkdir -p "${OUTDIR}/02_assembly/skesa"
            skesa \
                --fastq "${TRIMMED_R1},${TRIMMED_R2}" \
                --cores "${THREADS}" \
                --memory "${MEM_INT}" \
                --contigs_out "${OUTDIR}/02_assembly/skesa/contigs.fasta" \
                > "${OUTDIR}/logs/skesa.log" 2>&1
            ASSEMBLY="${OUTDIR}/02_assembly/skesa/contigs.fasta"
            ;;

        unicycler)
            local uni_args=(-1 "${TRIMMED_R1}" -2 "${TRIMMED_R2}"
                            -o "${OUTDIR}/02_assembly/unicycler/"
                            -t "${THREADS}")
            if [[ "${READ_TYPE}" == "hybrid" ]]; then
                uni_args+=(-l "${FILTERED_LONG}" --mode normal)
            else
                uni_args+=(--mode conservative)
            fi
            unicycler "${uni_args[@]}" > "${OUTDIR}/logs/unicycler.log" 2>&1
            ASSEMBLY="${OUTDIR}/02_assembly/unicycler/assembly.fasta"
            ;;

        flye)
            local flye_read_flag="--nano-hq"
            ${HIFI_MODE} && flye_read_flag="--pacbio-hifi"
            flye \
                ${flye_read_flag} "${FILTERED_LONG}" \
                --out-dir "${OUTDIR}/02_assembly/flye/" \
                --threads "${THREADS}" \
                --genome-size "${GENOME_SIZE}" \
                --iterations 3 \
                > "${OUTDIR}/logs/flye.log" 2>&1
            ASSEMBLY="${OUTDIR}/02_assembly/flye/assembly.fasta"
            ;;

        canu)
            local canu_read_flag="-nanopore"
            ${HIFI_MODE} && canu_read_flag="-pacbio-hifi"
            canu \
                -p "${SAMPLE}" \
                -d "${OUTDIR}/02_assembly/canu/" \
                genomeSize="${GENOME_SIZE}" \
                useGrid=false \
                maxThreads="${THREADS}" \
                maxMemory="${MEMORY}" \
                ${canu_read_flag} "${FILTERED_LONG}" \
                > "${OUTDIR}/logs/canu.log" 2>&1
            ASSEMBLY="${OUTDIR}/02_assembly/canu/${SAMPLE}.contigs.fasta"
            ;;

        raven)
            mkdir -p "${OUTDIR}/02_assembly/raven"
            raven \
                --threads "${THREADS}" \
                --graphical-fragment-assembly "${OUTDIR}/02_assembly/raven/assembly.gfa" \
                "${FILTERED_LONG}" \
                > "${OUTDIR}/02_assembly/raven/assembly.fasta" \
                2> "${OUTDIR}/logs/raven.log"
            ASSEMBLY="${OUTDIR}/02_assembly/raven/assembly.fasta"
            ;;
    esac

    [[ -f "${ASSEMBLY}" && -s "${ASSEMBLY}" ]] || \
        err "Assembly output not found or empty. Check ${OUTDIR}/logs/${ASSEMBLER}.log"

    # ── Contig filtering ─────────────────────────────────────────────────────
    log "Filtering contigs shorter than ${MIN_CONTIG_LENGTH} bp..."
    local filtered="${ASSEMBLY%.fasta}_filtered.fasta"

    if command -v seqkit &>/dev/null; then
        # Preferred path: seqkit is robust and handles multiline FASTA
        seqkit seq --min-len "${MIN_CONTIG_LENGTH}" "${ASSEMBLY}" \
            > "${filtered}" 2>/dev/null
    else
        # FIX: original awk had a bug — the final record was only printed in END
        # but END only ran the print, not the length guard for the last sequence.
        # Corrected logic below handles all records uniformly.
        warn "seqkit not found; using awk fallback for contig filtering."
        awk -v min="${MIN_CONTIG_LENGTH}" '
            /^>/ {
                if (hdr != "" && length(seq) >= min) print hdr "\n" seq
                hdr = $0; seq = ""
                next
            }
            { seq = seq $0 }
            END { if (hdr != "" && length(seq) >= min) print hdr "\n" seq }
        ' "${ASSEMBLY}" > "${filtered}"
    fi

    ASSEMBLY="${filtered}"
    local ctg_count
    ctg_count=$(grep -c '^>' "${ASSEMBLY}")
    ok "Assembly complete: ${ctg_count} contigs ≥ ${MIN_CONTIG_LENGTH} bp"

    # Persist the final ASSEMBLY path for resume
    echo "${ASSEMBLY}" > "${OUTDIR}/logs/.assembly_path"

    stage_done "assembly"
}

# Restore ASSEMBLY variable when resuming past the assembly stage
_restore_assembly_path() {
    local path_file="${OUTDIR}/logs/.assembly_path"
    [[ -f "${path_file}" ]] || err "Cannot resume: assembly path file missing. Re-run without -r."
    ASSEMBLY=$(cat "${path_file}")
    [[ -f "${ASSEMBLY}" && -s "${ASSEMBLY}" ]] || \
        err "Cannot resume: assembly file missing or empty: ${ASSEMBLY}"
    log "Resumed with assembly: ${ASSEMBLY}"
}

# ── Polishing ──────────────────────────────────────────────────────────────────
run_polishing() {
    skip_if_done "polishing" && {
        _restore_assembly_path
        # For Illumina/hybrid, update ASSEMBLY to point to Pilon output.
        # For long-read runs, _restore_assembly_path already has the correct path.
        if [[ "${READ_TYPE}" == "illumina" || "${READ_TYPE}" == "hybrid" ]]; then
            ASSEMBLY="${OUTDIR}/03_polishing/${SAMPLE}_pilon.fasta"
        fi
        return
    }

    log "[3/5] Polishing assembly..."

    # ── Pilon (Illumina) ─────────────────────────────────────────────────────
    if [[ "${READ_TYPE}" == "illumina" || "${READ_TYPE}" == "hybrid" ]]; then

        # RAM guard: Pilon is JVM-based; warn if MEM_INT exceeds available RAM
        local avail_ram
        avail_ram=$(awk '/MemAvailable/{print int($2/1024/1024)}' /proc/meminfo 2>/dev/null || echo 0)
        local MEM_INT="${MEMORY//[^0-9]/}"
        if (( avail_ram > 0 && MEM_INT > avail_ram )); then
            warn "Pilon memory limit (${MEMORY}) exceeds available RAM (${avail_ram}G). " \
                 "Consider reducing -m or closing other processes."
        fi

        log "Mapping Illumina reads for Pilon polishing..."
        bowtie2-build --threads "${THREADS}" --quiet \
            "${ASSEMBLY}" "${OUTDIR}/03_polishing/bt2_idx" \
            > "${OUTDIR}/logs/bowtie2_build.log" 2>&1

        bowtie2 \
            -x "${OUTDIR}/03_polishing/bt2_idx" \
            -1 "${TRIMMED_R1}" -2 "${TRIMMED_R2}" \
            --threads "${THREADS}" \
            --very-sensitive-local \
            --no-unal \
            2> "${OUTDIR}/logs/bowtie2.log" \
        | samtools sort -@ "${THREADS}" \
            -o "${OUTDIR}/03_polishing/mapped.bam"

        samtools index -@ "${THREADS}" "${OUTDIR}/03_polishing/mapped.bam"

        log "Running Pilon..."
        pilon \
            --genome "${ASSEMBLY}" \
            --frags "${OUTDIR}/03_polishing/mapped.bam" \
            --output "${SAMPLE}_pilon" \
            --outdir "${OUTDIR}/03_polishing/" \
            --threads "${THREADS}" \
            --changes --vcf \
            --fix all \
            > "${OUTDIR}/logs/pilon.log" 2>&1

        ASSEMBLY="${OUTDIR}/03_polishing/${SAMPLE}_pilon.fasta"
        ok "Pilon polishing done"
    else
        # Nanopore-only and HiFi: assemblers (Flye, Raven, Canu) include their
        # own internal polishing iterations; no external polishing is applied.
        log "Long-read assembly — skipping external polishing (assembler handles consensus)."
    fi

    ok "Polishing complete: ${ASSEMBLY}"
    # Persist polished path
    echo "${ASSEMBLY}" > "${OUTDIR}/logs/.assembly_path"
    stage_done "polishing"
}

# ── Quality assessment ─────────────────────────────────────────────────────────
run_quality_assessment() {
    skip_if_done "qa" && return
    log "[4/5] Quality assessment..."

    # ── QUAST ─────────────────────────────────────────────────────────────────
    quast "${ASSEMBLY}" \
        -o "${OUTDIR}/05_reports/quast/" \
        -t "${THREADS}" \
        --min-contig "${MIN_CONTIG_LENGTH}" \
        --gene-finding \
        --rna-finding \
        --est-ref-size "${GENOME_SIZE}" \
        > "${OUTDIR}/logs/quast.log" 2>&1

    # ── CheckM2 (preferred) or CheckM ────────────────────────────────────────
    if command -v checkm2 &>/dev/null; then
        log "Running CheckM2..."
        mkdir -p "${OUTDIR}/05_reports/checkm2"
        checkm2 predict \
            --input "${ASSEMBLY}" \
            --output-directory "${OUTDIR}/05_reports/checkm2/" \
            --threads "${THREADS}" \
            --force \
            > "${OUTDIR}/logs/checkm2.log" 2>&1
    elif command -v checkm &>/dev/null; then
        log "Running CheckM (legacy)..."
        local checkm_input="${OUTDIR}/05_reports/checkm/input"
        mkdir -p "${checkm_input}"
        cp "${ASSEMBLY}" "${checkm_input}/${SAMPLE}.fasta"
        checkm lineage_wf \
            -x fasta -t "${THREADS}" \
            "${checkm_input}" \
            "${OUTDIR}/05_reports/checkm/" \
            > "${OUTDIR}/logs/checkm.log" 2>&1
        checkm qa \
            "${OUTDIR}/05_reports/checkm/lineage.ms" \
            "${OUTDIR}/05_reports/checkm/" \
            -o 2 --tab_table \
            > "${OUTDIR}/05_reports/checkm/checkm_results.tsv"
    else
        warn "Neither checkm2 nor checkm found; skipping completeness assessment."
    fi

    # ── BUSCO ─────────────────────────────────────────────────────────────────
    if command -v busco &>/dev/null; then
        log "Running BUSCO..."
        busco \
            -i "${ASSEMBLY}" \
            -o "${SAMPLE}_busco" \
            -l bacteria_odb10 \
            -m genome \
            -c "${THREADS}" \
            --out_path "${OUTDIR}/05_reports/busco/" \
            --quiet \
            > "${OUTDIR}/logs/busco.log" 2>&1
    else
        warn "BUSCO not found; skipping."
    fi

    # ── PlasmidFinder ─────────────────────────────────────────────────────────
    if ${PLASMID_DETECTION}; then
        if command -v plasmidfinder.py &>/dev/null; then
            log "Running PlasmidFinder (db: ${PLASMID_DB})..."
            mkdir -p "${OUTDIR}/05_reports/plasmidfinder"
            plasmidfinder.py \
                -i "${ASSEMBLY}" \
                -o "${OUTDIR}/05_reports/plasmidfinder/" \
                -p "${PLASMID_DB}" \
                > "${OUTDIR}/logs/plasmidfinder.log" 2>&1
        else
            warn "plasmidfinder.py not found; skipping plasmid detection."
        fi
    fi

    ok "Quality assessment complete"
    stage_done "qa"
}

# ── Annotation ─────────────────────────────────────────────────────────────────
run_annotation() {
    skip_if_done "annotation" && return
    log "[5/5] Annotation..."

    # ── Bakta ─────────────────────────────────────────────────────────────────
    if command -v bakta &>/dev/null; then
        log "Running Bakta..."
        local bakta_db_flag=()
        if [[ -n "${BAKTA_DB}" ]]; then
            bakta_db_flag=(--db "${BAKTA_DB}")
        else
            # FIX: bakta_db list prints a path per line; pipe through head -1
            # and strip any trailing whitespace for safety
            local _detected_db
            _detected_db=$(bakta_db list 2>/dev/null | head -1 | tr -d '[:space:]')
            [[ -n "${_detected_db}" ]] && bakta_db_flag=(--db "${_detected_db}")
        fi
        bakta \
            "${bakta_db_flag[@]}" \
            --output "${OUTDIR}/04_annotation/bakta/" \
            --threads "${THREADS}" \
            --prefix "${SAMPLE}" \
            --force \
            "${ASSEMBLY}" \
            > "${OUTDIR}/logs/bakta.log" 2>&1 || \
            warn "Bakta failed (check logs/bakta.log; verify DB path with -d)."
    else
        warn "Bakta not found; skipping annotation. Install Bakta and provide a DB with -d."
    fi

    # ── antiSMASH ─────────────────────────────────────────────────────────────
    if command -v antismash &>/dev/null; then
        log "Running antiSMASH..."
        antismash \
            "${ASSEMBLY}" \
            --output-dir "${OUTDIR}/04_annotation/antismash/" \
            --cpus "${THREADS}" \
            --genefinding-tool prodigal \
            --minimal \
            > "${OUTDIR}/logs/antismash.log" 2>&1
    fi

    ok "Annotation complete"
    stage_done "annotation"
}

# ── Final report ───────────────────────────────────────────────────────────────
generate_report() {
    log "[5/5] Generating final report..."

    # FIX: call assembly-stats once and cache; original called it 3× separately
    local n_contigs total_len n50 largest
    n_contigs=$(grep -c '^>' "${ASSEMBLY}" 2>/dev/null || echo "N/A")

    if command -v assembly-stats &>/dev/null; then
        local _stats
        _stats=$(assembly-stats "${ASSEMBLY}" 2>/dev/null)
        total_len=$(awk '/sum/{print $2}'     <<< "${_stats}" || echo "N/A")
        n50=$(      awk '/N50/{print $2}'     <<< "${_stats}" || echo "N/A")
        largest=$(  awk '/largest/{print $2}' <<< "${_stats}" || echo "N/A")
    else
        total_len="N/A (install assembly-stats)"
        n50="N/A"; largest="N/A"
    fi

    cat > "${OUTDIR}/05_reports/SUMMARY.md" <<EOF
# Bacterial Genome Assembly Report

| Field           | Value                       |
|:----------------|:----------------------------|
| Sample          | ${SAMPLE}                   |
| Date            | $(date '+%Y-%m-%d %H:%M')   |
| Pipeline        | v${VERSION}                 |
| Read type       | ${READ_TYPE}                |
| Assembler       | ${ASSEMBLER}                |
| HiFi mode       | ${HIFI_MODE}                |

## Input Files
$(  [[ -n "${R1:-}"   ]] && echo "- Illumina R1: \`${R1}\`")
$(  [[ -n "${R2:-}"   ]] && echo "- Illumina R2: \`${R2}\`")
$(  [[ -n "${LONG:-}" ]] && echo "- Long reads:  \`${LONG}\`")

## Assembly Statistics

| Metric          | Value        |
|:----------------|:-------------|
| Contigs         | ${n_contigs} |
| Total length    | ${total_len} |
| Largest contig  | ${largest}   |
| N50             | ${n50}       |
| Min contig len  | ${MIN_CONTIG_LENGTH} bp |

## Output Files

| Step            | Path                                        |
|:----------------|:--------------------------------------------|
| Final assembly  | \`${ASSEMBLY}\`                             |
| QC reports      | \`${OUTDIR}/01_qc/\`                        |
| Bakta annotation| \`${OUTDIR}/04_annotation/bakta/\`          |
| Quality reports | \`${OUTDIR}/05_reports/\`                  |
| Logs            | \`${OUTDIR}/logs/\`                         |

## Suggested Next Steps
1. **AMR screening** – \`AMRFinderPlus\` or \`ResFinder\`
2. **NCBI submission** – use Bakta GFF3/GBFF outputs with \`table2asn\`
3. **MLST typing** – \`mlst\` (Torsten Seemann)
4. **Phylogenetics** – \`IQ-TREE2\` or \`FastTree\`
5. **Pan-genome** – \`Panaroo\` (preferred) or \`Roary\`
EOF

    # Disable ERR trap before clean banner so it doesn't fire on subshell exits
    trap - ERR

    echo -e "\n${GREEN}══════════════════════════════════════${NC}"
    echo -e "${GREEN}       PIPELINE COMPLETE ✓            ${NC}"
    echo -e "${GREEN}══════════════════════════════════════${NC}"
    echo -e "  Sample    : ${SAMPLE}"
    echo -e "  Assembly  : ${ASSEMBLY}"
    echo -e "  Contigs   : ${n_contigs}"
    echo -e "  N50       : ${n50}"
    echo -e "  Report    : ${OUTDIR}/05_reports/SUMMARY.md"
    echo -e "${GREEN}══════════════════════════════════════${NC}\n"
}

# ── Error trap ─────────────────────────────────────────────────────────────────
# Set after functions are defined so it doesn't fire during sourcing
trap 'err "Pipeline failed at line ${LINENO}. Check logs in ${OUTDIR}/logs/"' ERR

# ── Main ───────────────────────────────────────────────────────────────────────
main() {
    echo -e "${GREEN}══════════════════════════════════════${NC}"
    echo -e "${GREEN}  Bacterial Genome Assembly Pipeline  ${NC}"
    echo -e "${GREEN}             v${VERSION}               ${NC}"
    echo -e "${GREEN}══════════════════════════════════════${NC}\n"

    validate_inputs
    run_qc

    if ${QC_ONLY}; then
        trap - ERR
        ok "QC-only mode. Exiting."
        exit 0
    fi

    run_assembly
    run_polishing
    run_quality_assessment
    run_annotation
    generate_report
}

main "$@"
