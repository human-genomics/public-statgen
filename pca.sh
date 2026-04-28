#!/usr/bin/env bash
#
# pca.sh — PCA pipeline orchestrator
# ===================================
# Standalone (not invoked from main.sh). Operates on the merged reference
# panel produced by main.sh (KG + HGDP + SGDP + GIAB).
#
# Usage:
#   bash pca.sh
#
# Steps:
#   1. QC merged panel for PCA (geno/MAF/long-range LD/LD-prune/kinship/HWE)
#   2. Fit PCs and project samples
#   3. Pairwise Hudson FST  (post-MAF, post-long-range-LD, pre-LD-prune SNP set)
#   4. Correlate FST vs squared Euclidean distance between PC centroids
#   5. Calibrate FST against Privé et al. 2022 and fit FST = a·d² + b
#   6. Within-group PC-space variability  (median + RMS distance to centroid)
#
set -euo pipefail

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
PROJECT_DIR="$(cd "$(dirname "$0")" && pwd)"
TOOLS_BIN="${PROJECT_DIR}/tools/bin"

export DOWNLOADS_DIR="${PROJECT_DIR}/downloads"
export MERGE_DIR="${PROJECT_DIR}/merge"
export PCA_DIR="${PROJECT_DIR}/pca"
export PLINK1="${TOOLS_BIN}/plink1"
export PLINK2="${TOOLS_BIN}/plink2"
export PYTHON="${PROJECT_DIR}/tools/venv/bin/python"

# ---------------------------------------------------------------------------
# PLINK runtime
# ---------------------------------------------------------------------------
export PLINK_MEMORY=14000
export PLINK_THREADS=6

# ---------------------------------------------------------------------------
# QC parameters
# ---------------------------------------------------------------------------
# Genotype missingness
export GENO_PCA=0.01

# Minor allele frequency
export MAF_PCA=0.01

# Hardy-Weinberg p-value threshold
export HWE_PVALUE="1e-50"

# LD pruning — wider window than the ADMIXTURE QC (50/10/0.1)
export LD_WINDOW=1000
export LD_STEP=80
export LD_R2=0.1

# Kinship cutoffs (KING coefficient)
export KING_CUTOFF_AMR=0.088          # AMR-only pass + final cross-group pass
export KING_CUTOFF_NONAMR=0.05        # non-AMR-only pass

# PCA fit/projection
export N_PCS=30
export PCA_SEED=0

# Pairwise FST: minimum population sample size
export FST_MIN_POP_SIZE=5

# ---------------------------------------------------------------------------
# Logging — all subsequent output goes to both terminal and log file
# ---------------------------------------------------------------------------
mkdir -p "${PROJECT_DIR}/logs"
LOG_FILE="${PROJECT_DIR}/logs/pca_run_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "${LOG_FILE}") 2>&1

echo "============================================"
echo "PCA Pipeline Run — $(date)"
echo "============================================"
echo "  GENO_PCA=${GENO_PCA}"
echo "  MAF_PCA=${MAF_PCA}"
echo "  HWE_PVALUE=${HWE_PVALUE}"
echo "  LD_WINDOW=${LD_WINDOW}, LD_STEP=${LD_STEP}, LD_R2=${LD_R2}"
echo "  KING_CUTOFF_AMR=${KING_CUTOFF_AMR}"
echo "  KING_CUTOFF_NONAMR=${KING_CUTOFF_NONAMR}"
echo "  N_PCS=${N_PCS}, PCA_SEED=${PCA_SEED}"
echo "  LOG_FILE=${LOG_FILE}"
echo ""

# ---------------------------------------------------------------------------
# Step 1 — QC merged panel for PCA
# ---------------------------------------------------------------------------
echo "============================================"
echo "Step 1: QC merged panel for PCA"
echo "============================================"

if [[ -f "${PCA_DIR}/pca_qc.bed" ]]; then
    echo "  [skip] PCA QC output already exists in ${PCA_DIR}/"
else
    mkdir -p "${PCA_DIR}/scrap"
    bash "${PROJECT_DIR}/qc_pca.sh"
fi
echo ""

# ---------------------------------------------------------------------------
# Step 2 — Fit PCs and project samples
# ---------------------------------------------------------------------------
echo "============================================"
echo "Step 2: Fit PCs (${N_PCS}) and project samples"
echo "============================================"

if [[ -f "${PCA_DIR}/pca_projected.sscore" ]]; then
    echo "  [skip] PCA fit/projection output already exists in ${PCA_DIR}/"
else
    bash "${PROJECT_DIR}/compute_pca.sh"
fi
echo ""

# ---------------------------------------------------------------------------
# Step 3 — Pairwise Hudson FST
# ---------------------------------------------------------------------------
# Same individuals as the PCA panel (post-kinship, post-HWE), but on the
# pre-LD-prune variant set (geno + MAF + long-range-LD-exclusion only).
# LD pruning is appropriate for PCA but shrinks FST estimates by removing
# regions of high between-population differentiation, so we compute FST on
# the larger SNP set to keep values literature-comparable.
echo "============================================"
echo "Step 3: Pairwise Hudson FST"
echo "============================================"

if [[ -f "${PCA_DIR}/fst_pairs/fst_summary.tsv" ]]; then
    echo "  [skip] FST summary already exists in ${PCA_DIR}/fst_pairs/"
else
    "${PLINK2}" --bfile "${PCA_DIR}/scrap/ldregion_excluded" \
        --keep "${PCA_DIR}/pca_qc.fam" \
        --make-bed \
        --out "${PCA_DIR}/fst_input" \
        --memory "${PLINK_MEMORY}" --threads "${PLINK_THREADS}" \
        --no-input-missing-phenotype

    FST_BFILE="${PCA_DIR}/fst_input" \
        "${PYTHON}" "${PROJECT_DIR}/compute_fst.py"
fi
echo ""

# ---------------------------------------------------------------------------
# Step 4 — Correlate FST with squared PC distance
# ---------------------------------------------------------------------------
echo "============================================"
echo "Step 4: Correlate FST vs squared PC distance"
echo "============================================"

if [[ -f "${PCA_DIR}/plots/fst_vs_pcdist_supervised.png" ]]; then
    echo "  [skip] FST-vs-PC-distance plots already exist"
else
    "${PYTHON}" "${PROJECT_DIR}/correlate_fst_pca.py"
fi
echo ""

# ---------------------------------------------------------------------------
# Step 5 — Calibrate FST against Privé et al. 2022 and fit FST = a·d² + b
# ---------------------------------------------------------------------------
echo "============================================"
echo "Step 5: Calibrate FST and fit FST = a·d² + b"
echo "============================================"

if [[ -f "${PCA_DIR}/fst_pcdist_equation.txt" \
   && -f "${PCA_DIR}/plots/fst_vs_pcdist_calibrated.png" \
   && -f "${PCA_DIR}/plots/fst_calibration_vs_reference.png" ]]; then
    echo "  [skip] Calibrated FST equation already exists"
else
    "${PYTHON}" "${PROJECT_DIR}/calibrate_and_fit_fst.py"
fi
echo ""

# ---------------------------------------------------------------------------
# Step 6 — Within-group PC-space variability
# ---------------------------------------------------------------------------
echo "============================================"
echo "Step 6: Within-group PC-space variability"
echo "============================================"

if [[ -f "${PCA_DIR}/within_group_stats_pop.tsv" \
   && -f "${PCA_DIR}/within_group_stats_metadata_superpop.tsv" \
   && -f "${PCA_DIR}/plots/within_group_variability_n5.png" \
   && -f "${PCA_DIR}/plots/within_group_variability_n10.png" \
   && -f "${PCA_DIR}/plots/within_group_variability_n20.png" \
   && -f "${PCA_DIR}/plots/within_group_variability_metadata_superpop.png" ]]; then
    echo "  [skip] Within-group variability outputs already exist"
else
    "${PYTHON}" "${PROJECT_DIR}/within_group_variability.py"
fi
echo ""

echo "============================================"
echo "PCA pipeline steps complete."
echo "Log saved to ${LOG_FILE}"
echo "============================================"
