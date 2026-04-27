#!/usr/bin/env bash
#
# qc_pca.sh
# ----------
# Quality-control pipeline preparing the merged reference panel for PCA.
#
# Unlike qc_admixture.sh, the SNP-based filters here are computed on the
# FULL merged panel (not a supervised subset), with PCA-specific LD
# pruning parameters and a three-pass kinship strategy that enforces
# the strictest cutoff across the whole cohort.
#
# QC order:
#   1. Genotype missingness filter (--geno)
#   2. Minor allele frequency filter (--maf)
#   3. Exclude long-range LD regions (Price et al. 2008, hg38)
#   4. LD pruning (--indep-pairwise) with PCA-specific window/step/r2
#   5. Kinship — three passes:
#        5a. AMR-only       — KING cutoff KING_CUTOFF_AMR
#                             (AMR = 57 supervised "American" samples)
#        5b. Non-AMR-only   — KING cutoff KING_CUTOFF_NONAMR
#        5c. Cross-group    — merge unrelated sets, re-run KING at
#                             KING_CUTOFF_AMR (catches AMR<->non-AMR
#                             pairs that survive the partition)
#   6. Hardy-Weinberg equilibrium exact test (--hwe)
#      on the final unrelated set, with --make-bed to produce output.
#
# Expected environment (set by pca.sh):
#   MERGE_DIR              — directory containing merged_kg_hgdp_sgdp.{bed,bim,fam}
#   DOWNLOADS_DIR          — directory containing high_ld_regions_hg38.bed
#   PCA_DIR                — output directory
#   PLINK1                 — path to plink1 binary
#   PLINK2                 — path to plink2 binary
#   PYTHON                 — path to venv python
#   PLINK_MEMORY           — memory limit in MB
#   PLINK_THREADS          — number of threads
#   GENO_PCA               — genotype missingness threshold (e.g. 0.01)
#   MAF_PCA                — minor allele frequency threshold (e.g. 0.01)
#   HWE_PVALUE             — Hardy-Weinberg p-value threshold (e.g. 1e-50)
#   LD_WINDOW              — LD pruning window size in SNPs (e.g. 1000)
#   LD_STEP                — LD pruning step size (e.g. 80)
#   LD_R2                  — LD pruning r-squared threshold (e.g. 0.1)
#   KING_CUTOFF_AMR        — kinship cutoff for AMR samples (e.g. 0.088)
#   KING_CUTOFF_NONAMR     — kinship cutoff for non-AMR samples (e.g. 0.05)
#
# Input files (read-only):
#   summary/supervised.csv  — used to identify the AMR sample IDs
#
# Outputs (in PCA_DIR):
#   pca_qc.{bed,bim,fam}    — final QC'd, unrelated, HWE-filtered fileset
#
set -euo pipefail

# ---------------------------------------------------------------------------
# Validate environment
# ---------------------------------------------------------------------------
for var in MERGE_DIR DOWNLOADS_DIR PCA_DIR PLINK1 PLINK2 PYTHON \
           PLINK_MEMORY PLINK_THREADS \
           GENO_PCA MAF_PCA HWE_PVALUE \
           LD_WINDOW LD_STEP LD_R2 \
           KING_CUTOFF_AMR KING_CUTOFF_NONAMR; do
    if [[ -z "${!var:-}" ]]; then
        echo "Error: ${var} is not set." >&2
        exit 1
    fi
done

PROJECT_DIR="$(cd "$(dirname "$0")" && pwd)"
SCRAP="${PCA_DIR}/scrap"
PLINK_FLAGS=(--memory "${PLINK_MEMORY}" --threads "${PLINK_THREADS}" --no-input-missing-phenotype)

mkdir -p "${PCA_DIR}" "${SCRAP}"

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
fmt() { printf "%'d" "$1"; }
count_snps()    { wc -l < "$1.bim"; }
count_samples() { wc -l < "$1.fam"; }

# ---------------------------------------------------------------------------
# Starting point
# ---------------------------------------------------------------------------
INPUT="${MERGE_DIR}/merged_kg_hgdp_sgdp"
SNPS_BEFORE=$(count_snps "${INPUT}")
SAMPLES_BEFORE=$(count_samples "${INPUT}")
echo "  Starting data: $(fmt ${SNPS_BEFORE}) SNPs, $(fmt ${SAMPLES_BEFORE}) samples"
echo ""

# ===================================================================
# STEP 1 — Genotype missingness (--geno)
# ===================================================================
echo "  Step 1/6: Genotype missingness filter (--geno ${GENO_PCA})"

"${PLINK2}" --bfile "${INPUT}" \
    --geno "${GENO_PCA}" \
    --make-bed \
    --out "${SCRAP}/geno_filtered" \
    "${PLINK_FLAGS[@]}"

SNPS_AFTER=$(count_snps "${SCRAP}/geno_filtered")
SNPS_REMOVED=$((SNPS_BEFORE - SNPS_AFTER))
echo "    Removed $(fmt ${SNPS_REMOVED}) SNPs with >${GENO_PCA} missingness"
echo "    Remaining: $(fmt ${SNPS_AFTER}) SNPs"
echo ""

# ===================================================================
# STEP 2 — Minor allele frequency (--maf)
# ===================================================================
echo "  Step 2/6: MAF filter (--maf ${MAF_PCA})"

SNPS_BEFORE_MAF=${SNPS_AFTER}
"${PLINK2}" --bfile "${SCRAP}/geno_filtered" \
    --maf "${MAF_PCA}" \
    --make-bed \
    --out "${SCRAP}/maf_filtered" \
    "${PLINK_FLAGS[@]}"

SNPS_AFTER=$(count_snps "${SCRAP}/maf_filtered")
SNPS_REMOVED=$((SNPS_BEFORE_MAF - SNPS_AFTER))
echo "    Removed $(fmt ${SNPS_REMOVED}) SNPs with MAF < ${MAF_PCA}"
echo "    Remaining: $(fmt ${SNPS_AFTER}) SNPs"
echo ""

# ===================================================================
# STEP 3 — Exclude long-range LD regions
# ===================================================================
echo "  Step 3/6: Exclude long-range LD regions (Price et al. 2008, hg38)"

HIGH_LD_BED="${DOWNLOADS_DIR}/high_ld_regions_hg38.bed"
if [[ ! -f "${HIGH_LD_BED}" ]]; then
    echo "Error: high-LD regions file not found: ${HIGH_LD_BED}" >&2
    echo "Run download_files.sh first." >&2
    exit 1
fi

# Convert BED (0-based) to PLINK range format (1-based)
awk 'BEGIN{OFS="\t"} {
    chr = $1; sub(/^chr/, "", chr)
    print chr, $2+1, $3, $4
}' "${HIGH_LD_BED}" > "${SCRAP}/high_ld_ranges.txt"

SNPS_BEFORE_LD=${SNPS_AFTER}
"${PLINK2}" --bfile "${SCRAP}/maf_filtered" \
    --exclude range "${SCRAP}/high_ld_ranges.txt" \
    --make-bed \
    --out "${SCRAP}/ldregion_excluded" \
    "${PLINK_FLAGS[@]}"

SNPS_AFTER=$(count_snps "${SCRAP}/ldregion_excluded")
SNPS_REMOVED=$((SNPS_BEFORE_LD - SNPS_AFTER))
echo "    Removed $(fmt ${SNPS_REMOVED}) SNPs in long-range LD regions"
echo "    Remaining: $(fmt ${SNPS_AFTER}) SNPs"
echo ""

# ===================================================================
# STEP 4 — LD pruning (--indep-pairwise)
# ===================================================================
echo "  Step 4/6: LD pruning (window=${LD_WINDOW}, step=${LD_STEP}, r2=${LD_R2})"

"${PLINK2}" --bfile "${SCRAP}/ldregion_excluded" \
    --indep-pairwise "${LD_WINDOW}" "${LD_STEP}" "${LD_R2}" \
    --out "${SCRAP}/ld_prune" \
    "${PLINK_FLAGS[@]}"

SNPS_BEFORE_PRUNE=${SNPS_AFTER}

"${PLINK2}" --bfile "${SCRAP}/ldregion_excluded" \
    --extract "${SCRAP}/ld_prune.prune.in" \
    --make-bed \
    --out "${SCRAP}/ld_pruned" \
    "${PLINK_FLAGS[@]}"

SNPS_AFTER=$(count_snps "${SCRAP}/ld_pruned")
SNPS_REMOVED=$((SNPS_BEFORE_PRUNE - SNPS_AFTER))
echo "    Removed $(fmt ${SNPS_REMOVED}) SNPs by LD pruning"
echo "    Remaining: $(fmt ${SNPS_AFTER}) SNPs"
echo ""

# ===================================================================
# STEP 5 — Kinship (three passes)
# ===================================================================
echo "  Step 5/6: Kinship filtering (three-pass: AMR=${KING_CUTOFF_AMR}, non-AMR=${KING_CUTOFF_NONAMR}, cross-group=${KING_CUTOFF_AMR})"

# Build the lenient (0.088) / strict (0.05) keep-lists.
# Lenient group: 57 supervised "American" + the 2 GIAB Ashkenazi parents
# (HG003, HG004). Ashkenazi underwent a severe historical bottleneck, which
# inflates background pairwise kinship — the 0.05 cutoff would otherwise
# drop one of them as falsely related to the other.
"${PYTHON}" -c "
import pandas as pd
sup = pd.read_csv('${PROJECT_DIR}/summary/supervised.csv')
lenient_ids = set(sup.loc[sup['reference_population'] == 'American', 'sample_id'])
lenient_ids |= {'HG003', 'HG004'}
fam = pd.read_csv('${SCRAP}/ld_pruned.fam', sep=r'\s+', header=None,
                  names=['FID','IID','PAT','MAT','SEX','PHENO'])
amr   = fam[ fam['IID'].isin(lenient_ids)][['FID','IID']]
other = fam[~fam['IID'].isin(lenient_ids)][['FID','IID']]
amr.to_csv('${SCRAP}/amr_samples.txt',    sep='\t', header=False, index=False)
other.to_csv('${SCRAP}/nonamr_samples.txt', sep='\t', header=False, index=False)
print(f'    Lenient group (American supervised + GIAB AJ): {len(amr):,}')
print(f'    Strict group  (everyone else):                 {len(other):,}')
"

# --- 5a. Kinship on AMR ---
echo ""
echo "    5a. AMR pass (KING cutoff = ${KING_CUTOFF_AMR})"
AMR_COUNT_BEFORE=$(wc -l < "${SCRAP}/amr_samples.txt")

"${PLINK2}" --bfile "${SCRAP}/ld_pruned" \
    --keep "${SCRAP}/amr_samples.txt" \
    --king-cutoff "${KING_CUTOFF_AMR}" \
    --out "${SCRAP}/amr_king" \
    "${PLINK_FLAGS[@]}"

if [[ -f "${SCRAP}/amr_king.king.cutoff.out.id" ]]; then
    AMR_REMOVED=$(wc -l < "${SCRAP}/amr_king.king.cutoff.out.id")
    if head -1 "${SCRAP}/amr_king.king.cutoff.out.id" | grep -q "^#"; then
        AMR_REMOVED=$((AMR_REMOVED - 1))
    fi
else
    AMR_REMOVED=0
fi
echo "        AMR removed (related): $(fmt ${AMR_REMOVED}) of $(fmt ${AMR_COUNT_BEFORE})"

# --- 5b. Kinship on non-AMR ---
echo ""
echo "    5b. Non-AMR pass (KING cutoff = ${KING_CUTOFF_NONAMR})"
NONAMR_COUNT_BEFORE=$(wc -l < "${SCRAP}/nonamr_samples.txt")

"${PLINK2}" --bfile "${SCRAP}/ld_pruned" \
    --keep "${SCRAP}/nonamr_samples.txt" \
    --king-cutoff "${KING_CUTOFF_NONAMR}" \
    --out "${SCRAP}/nonamr_king" \
    "${PLINK_FLAGS[@]}"

if [[ -f "${SCRAP}/nonamr_king.king.cutoff.out.id" ]]; then
    NONAMR_REMOVED=$(wc -l < "${SCRAP}/nonamr_king.king.cutoff.out.id")
    if head -1 "${SCRAP}/nonamr_king.king.cutoff.out.id" | grep -q "^#"; then
        NONAMR_REMOVED=$((NONAMR_REMOVED - 1))
    fi
else
    NONAMR_REMOVED=0
fi
echo "        Non-AMR removed (related): $(fmt ${NONAMR_REMOVED}) of $(fmt ${NONAMR_COUNT_BEFORE})"

# --- Concatenate the two within-group unrelated sets ---
tail -n +2 "${SCRAP}/amr_king.king.cutoff.in.id"    >  "${SCRAP}/intermediate_unrelated.txt"
tail -n +2 "${SCRAP}/nonamr_king.king.cutoff.in.id" >> "${SCRAP}/intermediate_unrelated.txt"
INTERMEDIATE_COUNT=$(wc -l < "${SCRAP}/intermediate_unrelated.txt")

# --- 5c. Cross-group pass at the AMR cutoff ---
echo ""
echo "    5c. Cross-group pass on merged unrelated set (KING cutoff = ${KING_CUTOFF_AMR})"

"${PLINK2}" --bfile "${SCRAP}/ld_pruned" \
    --keep "${SCRAP}/intermediate_unrelated.txt" \
    --king-cutoff "${KING_CUTOFF_AMR}" \
    --out "${SCRAP}/cross_king" \
    "${PLINK_FLAGS[@]}"

if [[ -f "${SCRAP}/cross_king.king.cutoff.out.id" ]]; then
    CROSS_REMOVED=$(wc -l < "${SCRAP}/cross_king.king.cutoff.out.id")
    if head -1 "${SCRAP}/cross_king.king.cutoff.out.id" | grep -q "^#"; then
        CROSS_REMOVED=$((CROSS_REMOVED - 1))
    fi
else
    CROSS_REMOVED=0
fi
echo "        Cross-group removed (related): $(fmt ${CROSS_REMOVED}) of $(fmt ${INTERMEDIATE_COUNT})"

tail -n +2 "${SCRAP}/cross_king.king.cutoff.in.id" > "${SCRAP}/unrelated_keep.txt"
UNRELATED_FINAL=$(wc -l < "${SCRAP}/unrelated_keep.txt")

TOTAL_REMOVED=$((AMR_REMOVED + NONAMR_REMOVED + CROSS_REMOVED))
echo ""
echo "    Total removed by kinship (all passes): $(fmt ${TOTAL_REMOVED})"
echo "    Unrelated samples after kinship: $(fmt ${UNRELATED_FINAL})"
echo ""

# ===================================================================
# STEP 6 — HWE + final output
# ===================================================================
echo "  Step 6/6: HWE exact test (--hwe ${HWE_PVALUE}) on unrelated samples"

SNPS_BEFORE_HWE=${SNPS_AFTER}

"${PLINK2}" --bfile "${SCRAP}/ld_pruned" \
    --keep "${SCRAP}/unrelated_keep.txt" \
    --hwe "${HWE_PVALUE}" \
    --make-bed \
    --out "${PCA_DIR}/pca_qc" \
    "${PLINK_FLAGS[@]}"

SNPS_FINAL=$(count_snps "${PCA_DIR}/pca_qc")
SAMPLES_FINAL=$(count_samples "${PCA_DIR}/pca_qc")
SNPS_REMOVED_HWE=$((SNPS_BEFORE_HWE - SNPS_FINAL))
echo "    Removed $(fmt ${SNPS_REMOVED_HWE}) SNPs failing HWE (p < ${HWE_PVALUE})"
echo "    Final: $(fmt ${SAMPLES_FINAL}) samples, $(fmt ${SNPS_FINAL}) SNPs"
echo ""

# ===================================================================
# Summary
# ===================================================================
echo "  ================================================"
echo "  PCA QC Summary"
echo "  ================================================"
echo "  Parameters:"
echo "    Geno=${GENO_PCA}, MAF=${MAF_PCA}, HWE=${HWE_PVALUE}"
echo "    LD window=${LD_WINDOW} step=${LD_STEP} r2=${LD_R2}"
echo "    KING AMR=${KING_CUTOFF_AMR}, KING non-AMR=${KING_CUTOFF_NONAMR}, KING cross=${KING_CUTOFF_AMR}"
echo "  SNPs:    $(fmt ${SNPS_BEFORE}) -> $(fmt ${SNPS_FINAL})"
echo "  Samples: $(fmt ${SAMPLES_BEFORE}) -> $(fmt ${SAMPLES_FINAL})"
echo "  Output:  ${PCA_DIR}/pca_qc.{bed,bim,fam}"
echo ""

# Clean up logs and per-step intermediates
rm -f "${SCRAP}"/*.log "${SCRAP}"/*.nosex

echo "  PCA QC complete."
