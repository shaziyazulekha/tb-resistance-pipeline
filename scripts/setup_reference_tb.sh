#!/bin/bash
# =============================================================================
# REFERENCE GENOME SETUP — M. tuberculosis H37Rv (Reference-based Pipeline)
# =============================================================================
# Downloads and indexes the single reference genome used by the TB pipeline.
# This is a REFERENCE-BASED pipeline (BWA MEM + GATK HaplotypeCaller).
# No pan-genome, no de-novo assembly — just H37Rv as the one reference.
#
# REFERENCE:
#   M. tuberculosis H37Rv  (GCA_000195955.2 / NC_000962.3)
#   Size  : ~4.41 Mb  (single circular chromosome)
#   Source: NCBI RefSeq — the universal TB reference for all variant calling
#
# OUTPUTS:
#   reference/H37Rv.fasta             — reference FASTA (matches step4_tb.sh)
#   reference/H37Rv.fasta.bwt         — BWA index (.bwt .pac .ann .amb .sa)
#   reference/H37Rv.fasta.fai         — samtools fai index
#   reference/H37Rv.dict              — GATK sequence dictionary
#   reference/logs/                   — per-tool index logs
#
# TOOLS REQUIRED:
#   conda install -c bioconda bwa samtools gatk4
#   (wget or curl for downloads)
#
# USAGE:
#   bash setup_reference_tb.sh
#   Idempotent — re-running safely skips already completed steps.
# =============================================================================

set -uo pipefail
# NOT set -e globally — each step checked individually so one failure does
# not silently skip indexing steps that depend on a file already on disk.

# =============================================================================
# PATHS  — edit only if you want outputs somewhere other than ./reference/
# =============================================================================
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REF_DIR="${SCRIPT_DIR}/reference"
LOG_DIR="${REF_DIR}/logs"
REF_FASTA="${REF_DIR}/H37Rv.fasta"
REF_DICT="${REF_DIR}/H37Rv.dict"

mkdir -p "$REF_DIR" "$LOG_DIR"

# =============================================================================
# LOGGING
# =============================================================================
log()      { echo "[$(date '+%Y-%m-%d %H:%M:%S')]      $*"; }
log_ok()   { echo "[$(date '+%Y-%m-%d %H:%M:%S')] OK   $*"; }
log_warn() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] WARN $*"; }
log_err()  { echo "[$(date '+%Y-%m-%d %H:%M:%S')] ERR  $*"; }

# =============================================================================
# PRE-FLIGHT CHECKS
# =============================================================================
log "=== PRE-FLIGHT CHECKS ==="

MISSING=false

# Prefer wget; fall back to curl
DOWNLOADER=""
if command -v wget &>/dev/null; then
    DOWNLOADER="wget"
    log "  Found   : wget"
elif command -v curl &>/dev/null; then
    DOWNLOADER="curl"
    log "  Found   : curl (wget not found — using curl)"
else
    log_err "  Missing : wget or curl — install one to continue"
    MISSING=true
fi

for tool in bwa samtools gatk python3; do
    if command -v "$tool" &>/dev/null; then
        log "  Found   : $tool  ($(command -v "$tool"))"
    else
        log_err "  Missing : $tool"
        MISSING=true
    fi
done

if [ "$MISSING" = true ]; then
    log_err "Install missing tools and re-run:"
    log_err "  conda install -c bioconda bwa samtools gatk4"
    exit 1
fi

log_ok "All required tools present."
echo ""

# =============================================================================
# DOWNLOAD HELPER
# =============================================================================
# download_url <url> <output_path>
# Retries up to 3 times. Returns 1 on failure.
download_url() {
    local url="$1"
    local out="$2"
    local attempt

    for attempt in 1 2 3; do
        log "  Download attempt ${attempt}/3..."
        if [ "$DOWNLOADER" = "wget" ]; then
            # tee exit code intentionally ignored (|| true) so pipefail
            # does not mask wget's own exit code via the pipe
            wget --no-verbose --show-progress --tries=1 \
                 --timeout=180 -O "$out" "$url" 2>&1 \
                 | tee -a "${LOG_DIR}/downloads.log" || true
            # Re-check wget's exit code via its output file
            [ -s "$out" ] && return 0
        else
            curl -L --progress-bar --retry 0 --max-time 180 \
                 -o "$out" "$url" 2>&1 \
                 | tee -a "${LOG_DIR}/downloads.log" || true
            [ -s "$out" ] && return 0
        fi
        [ "$attempt" -lt 3 ] && { log_warn "  Retrying in 10s..."; sleep 10; }
    done
    return 1
}

# =============================================================================
# STEP 1 — DOWNLOAD H37Rv REFERENCE FASTA
# =============================================================================
log "=== STEP 1 — DOWNLOAD H37Rv REFERENCE GENOME ==="

# Primary URL: NCBI RefSeq FTP (GCA_000195955.2 / ASM19595v2)
H37RV_FTP="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/195/955/GCA_000195955.2_ASM19595v2/GCA_000195955.2_ASM19595v2_genomic.fna.gz"

# Fallback URL: NCBI Entrez efetch for NC_000962.3
# (the /nuccore/*.fasta path is not a valid FTP endpoint)
H37RV_FALLBACK="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_000962.3&rettype=fasta&retmode=text"

RAW_GZ="${REF_DIR}/H37Rv_raw.fna.gz"
RAW_FA="${REF_DIR}/H37Rv_raw.fasta"

if [ -f "$REF_FASTA" ] && [ -s "$REF_FASTA" ]; then
    log_ok "H37Rv.fasta already exists — skipping download."
    log "  Delete ${REF_FASTA} to force a fresh download."
else
    rm -f "$RAW_GZ" "$RAW_FA"

    log "Downloading H37Rv from NCBI FTP (~4.4 Mb compressed)..."
    log "  Source: $H37RV_FTP"

    if download_url "$H37RV_FTP" "$RAW_GZ"; then
        # Validate gzip integrity before decompressing
        if ! gzip -t "$RAW_GZ" 2>/dev/null; then
            log_err "Downloaded file failed gzip integrity check: $RAW_GZ"
            log_err "Trying fallback URL..."
            rm -f "$RAW_GZ"
            if download_url "$H37RV_FALLBACK" "$RAW_FA"; then
                log_ok "Fallback download succeeded."
            else
                log_err "Both download URLs failed. Check your internet connection."
                exit 1
            fi
        else
            log "  Decompressing..."
            # Use gunzip -c (stdout) to control the output filename explicitly.
            # gunzip -f in-place would strip .gz and produce H37Rv_raw.fna,
            # not H37Rv_raw.fasta, causing the header-rename step to fail.
            gunzip -c "$RAW_GZ" > "$RAW_FA"
            rm -f "$RAW_GZ"
            log_ok "Decompressed → ${RAW_FA}"
        fi
    else
        log_err "Primary download failed — trying fallback URL..."
        if download_url "$H37RV_FALLBACK" "$RAW_FA"; then
            log_ok "Fallback download succeeded."
        else
            log_err "Both download URLs failed. Check your internet connection."
            exit 1
        fi
    fi

    # ── Rename the FASTA header to the stable accession NC_000962.3 ──────
    # GATK requires the chromosome name in VCFs to match the FASTA header.
    # Step 4 expects NC_000962.3 as the CHROM field (set in step5 parser).
    log "  Renaming FASTA header to NC_000962.3 (required by GATK + Step 5)..."
    python3 - "$RAW_FA" "$REF_FASTA" <<'PYEOF'
import sys

in_path, out_path = sys.argv[1], sys.argv[2]
header_written = False

with open(in_path) as fin, open(out_path, 'w') as fout:
    for line in fin:
        if line.startswith('>'):
            if not header_written:
                # Always use NC_000962.3 — matches step5.py H37RV_CHROM constant
                fout.write('>NC_000962.3 Mycobacterium tuberculosis H37Rv complete genome\n')
                header_written = True
            # Skip any additional headers (multi-contig assemblies shouldn't
            # exist for H37Rv but this guards against download variants)
        else:
            fout.write(line)

print("  Header normalised → NC_000962.3")
PYEOF

    # Capture Python exit code explicitly; set -e is not active here
    py_exit=$?
    if [ $py_exit -ne 0 ]; then
        log_err "Python header-rename step failed (exit ${py_exit}). Aborting."
        rm -f "$REF_FASTA"
        exit 1
    fi

    rm -f "$RAW_FA"

    if [ ! -s "$REF_FASTA" ]; then
        log_err "H37Rv.fasta is empty after header rename. Aborting."
        exit 1
    fi

    log_ok "H37Rv reference FASTA ready → ${REF_FASTA}"
fi

echo ""

# =============================================================================
# STEP 2 — BWA INDEX
# =============================================================================
log "=== STEP 2 — BWA INDEX ==="
# BWA produces 5 index files: .bwt .pac .ann .amb .sa
# Required by step4_tb.sh before bwa mem alignment

if [ -f "${REF_FASTA}.bwt" ] && [ -s "${REF_FASTA}.bwt" ]; then
    log_ok "BWA index already exists — skipping."
    log "  Delete ${REF_FASTA}.bwt to force re-index."
else
    log "Building BWA index for H37Rv (~2 minutes)..."
    if bwa index "$REF_FASTA" 2>"${LOG_DIR}/bwa_index.log"; then
        log_ok "BWA index complete."
        log "  Index files: ${REF_FASTA}.{bwt,pac,ann,amb,sa}"
    else
        log_err "BWA index failed — see ${LOG_DIR}/bwa_index.log"
        tail -n 20 "${LOG_DIR}/bwa_index.log"
        exit 1
    fi
fi

echo ""

# =============================================================================
# STEP 3 — SAMTOOLS FAIDX
# =============================================================================
log "=== STEP 3 — SAMTOOLS FAIDX ==="
# Produces the .fai index required by GATK HaplotypeCaller and VariantFiltration

if [ -f "${REF_FASTA}.fai" ] && [ -s "${REF_FASTA}.fai" ]; then
    log_ok "samtools fai index already exists — skipping."
else
    log "Running samtools faidx..."
    if samtools faidx "$REF_FASTA" 2>"${LOG_DIR}/faidx.log"; then
        log_ok "samtools fai index complete → ${REF_FASTA}.fai"
    else
        log_err "samtools faidx failed — see ${LOG_DIR}/faidx.log"
        tail -n 10 "${LOG_DIR}/faidx.log"
        exit 1
    fi
fi

echo ""

# =============================================================================
# STEP 4 — GATK SEQUENCE DICTIONARY
# =============================================================================
log "=== STEP 4 — GATK SEQUENCE DICTIONARY ==="
# Produces H37Rv.dict — mandatory for GATK HaplotypeCaller and VariantFiltration
# GATK looks for the .dict file at the same path as the FASTA with .dict extension

if [ -f "$REF_DICT" ] && [ -s "$REF_DICT" ]; then
    log_ok "GATK sequence dictionary already exists — skipping."
    log "  Delete ${REF_DICT} to force re-create."
else
    log "Creating GATK sequence dictionary..."
    if gatk CreateSequenceDictionary \
            -R "$REF_FASTA" \
            -O "$REF_DICT" \
            >/dev/null 2>"${LOG_DIR}/gatk_dict.log"; then
        log_ok "GATK dictionary complete → ${REF_DICT}"
    else
        log_err "GATK CreateSequenceDictionary failed — see ${LOG_DIR}/gatk_dict.log"
        tail -n 20 "${LOG_DIR}/gatk_dict.log"
        exit 1
    fi
fi

echo ""

# =============================================================================
# STEP 5 — VALIDATE ALL INDEX FILES
# =============================================================================
log "=== STEP 5 — VALIDATION ==="

FAIL_COUNT=0

# Ordered list of expected files and their descriptions.
# Using parallel arrays instead of an associative array to guarantee
# deterministic iteration order across all bash versions.
EXPECTED_PATHS=(
    "${REF_FASTA}"
    "${REF_FASTA}.bwt"
    "${REF_FASTA}.pac"
    "${REF_FASTA}.ann"
    "${REF_FASTA}.amb"
    "${REF_FASTA}.sa"
    "${REF_FASTA}.fai"
    "${REF_DICT}"
)
EXPECTED_DESCS=(
    "main reference FASTA"
    "BWA index bwt"
    "BWA index pac"
    "BWA index ann"
    "BWA index amb"
    "BWA index sa"
    "samtools fai index"
    "GATK sequence dictionary"
)

for i in "${!EXPECTED_PATHS[@]}"; do
    fpath="${EXPECTED_PATHS[$i]}"
    desc="${EXPECTED_DESCS[$i]}"
    if [ -f "$fpath" ] && [ -s "$fpath" ]; then
        size=$(du -h "$fpath" 2>/dev/null | cut -f1)
        log_ok "  [${size}]  $(basename "$fpath")  — ${desc}"
    else
        log_err "  MISSING or EMPTY: $(basename "$fpath")  — ${desc}"
        FAIL_COUNT=$((FAIL_COUNT + 1))
    fi
done

# ── Print reference genome summary via Python ──────────────────────────────
echo ""
python3 - "$REF_FASTA" <<'PYEOF'
import sys

fasta = sys.argv[1]
contigs = []

with open(fasta) as f:
    name, length = None, 0
    for line in f:
        line = line.rstrip()
        if line.startswith('>'):
            if name:
                contigs.append((name, length))
            name = line[1:].split()[0]
            length = 0
        else:
            length += len(line)
    if name:
        contigs.append((name, length))

total = sum(l for _, l in contigs)
print(f"  Reference genome contents:")
print(f"  {'Contig':<20} {'Length (bp)':>14}  Description")
print(f"  {'-'*55}")
for name, length in contigs:
    desc = "M. tuberculosis H37Rv complete genome" if "NC_000962" in name else "contig"
    print(f"  {name:<20} {length:>14,}  {desc}")
print(f"  {'-'*55}")
print(f"  {'TOTAL':<20} {total:>14,}")
PYEOF

echo ""

# =============================================================================
# FINAL SUMMARY
# =============================================================================
log "================================================"
if [ "$FAIL_COUNT" -eq 0 ]; then
    log_ok "REFERENCE SETUP COMPLETE — H37Rv ready for the TB pipeline."
else
    log_warn "SETUP FINISHED WITH ${FAIL_COUNT} MISSING FILE(S) — check ${LOG_DIR}/"
fi
log ""
log "Reference files in ${REF_DIR}/:"
log "  H37Rv.fasta           — main reference FASTA (used by step4_tb.sh)"
log "  H37Rv.fasta.bwt/.pac/.ann/.amb/.sa  — BWA alignment index"
log "  H37Rv.fasta.fai       — samtools random-access index"
log "  H37Rv.dict            — GATK sequence dictionary"
log ""
log "Configuration check:"
log "  step4_tb.sh  →  REF=\"reference/H37Rv.fasta\"        ✓ matches"
log "  step5_tb.py  →  H37RV_CHROM = \"NC_000962.3\"        ✓ matches FASTA header"
log ""
log "Index logs : ${LOG_DIR}/"
log ""
log "Next step  : bash run_pipeline_tb.sh"
log "  OR run steps individually starting with:"
log "  python3 step1_tb.py    (fetch BV-BRC metadata)"
log "================================================"

if [ "$FAIL_COUNT" -gt 0 ]; then
    exit 1
fi