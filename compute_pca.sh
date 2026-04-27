#!/usr/bin/env bash
#
# compute_pca.sh
# --------------
# Fits PCs on the QC'd panel and projects the same samples back using the
# allele-weights / --score workflow. Projecting the fit set onto its own
# PCs is structurally identical to how external samples would later be
# projected, so the same projection invocation can be reused downstream.
#
# Steps:
#   1. plink2 --pca allele-wts N_PCS approx --seed PCA_SEED
#        Outputs: pca_pcs.{eigenval, eigenvec, eigenvec.allele}
#   2. plink2 --freq counts
#        Outputs: pca_counts.acount
#   3. plink2 --score (project samples onto PCs)
#        Outputs: pca_projected.sscore
#
# Expected environment (set by pca.sh):
#   PCA_DIR, PLINK2, PLINK_MEMORY, PLINK_THREADS, N_PCS, PCA_SEED
#
set -euo pipefail

for var in PCA_DIR PLINK2 PLINK_MEMORY PLINK_THREADS N_PCS PCA_SEED; do
    if [[ -z "${!var:-}" ]]; then
        echo "Error: ${var} is not set." >&2
        exit 1
    fi
done

PLINK_FLAGS=(--memory "${PLINK_MEMORY}" --threads "${PLINK_THREADS}")
INPUT="${PCA_DIR}/pca_qc"
fmt() { printf "%'d" "$1"; }

# ===================================================================
# Step 1 — Fit PCA on the QC'd panel
# ===================================================================
echo "  Step 1/3: Fit PCA (${N_PCS} PCs, exact, seed=${PCA_SEED})"

"${PLINK2}" --bfile "${INPUT}" \
    --pca allele-wts "${N_PCS}" \
    --seed "${PCA_SEED}" \
    --out "${PCA_DIR}/pca_pcs" \
    "${PLINK_FLAGS[@]}"

echo "    Eigenvalues:"
sed 's/^/      /' "${PCA_DIR}/pca_pcs.eigenval"

FIT_SNPS=$(tail -n +2 "${PCA_DIR}/pca_pcs.eigenvec.allele" | wc -l)
FIT_SAMPLES=$(tail -n +2 "${PCA_DIR}/pca_pcs.eigenvec" | wc -l)
echo "    Fit on $(fmt ${FIT_SNPS}) SNP-allele rows, $(fmt ${FIT_SAMPLES}) samples"
echo ""

# ===================================================================
# Step 2 — Allele counts (consumed by --score variance-standardize)
# ===================================================================
echo "  Step 2/3: Compute allele counts"

"${PLINK2}" --bfile "${INPUT}" \
    --freq counts \
    --out "${PCA_DIR}/pca_counts" \
    "${PLINK_FLAGS[@]}"
echo ""

# ===================================================================
# Step 3 — Project samples onto the fitted PCs
# ===================================================================
echo "  Step 3/3: Project samples onto PCs"

# Locate A1 column in eigenvec.allele; PC score columns follow it
A1_COL=$(head -1 "${PCA_DIR}/pca_pcs.eigenvec.allele" | tr '\t' '\n' | grep -n '^A1$' | cut -d: -f1)
FIRST_PC=$((A1_COL + 1))
LAST_PC=$((A1_COL + N_PCS))
echo "    A1 column: ${A1_COL}, score columns: ${FIRST_PC}-${LAST_PC}"

awk 'NR>1 {print $2}' "${PCA_DIR}/pca_pcs.eigenvec.allele" > "${PCA_DIR}/pca_snps.txt"
N_PCA_SNPS=$(wc -l < "${PCA_DIR}/pca_snps.txt")
echo "    Projecting using $(fmt ${N_PCA_SNPS}) PCA SNPs"

"${PLINK2}" --bfile "${INPUT}" \
    --extract "${PCA_DIR}/pca_snps.txt" \
    --read-freq "${PCA_DIR}/pca_counts.acount" \
    --score "${PCA_DIR}/pca_pcs.eigenvec.allele" 2 "${A1_COL}" header-read no-mean-imputation variance-standardize \
    --score-col-nums "${FIRST_PC}-${LAST_PC}" \
    --out "${PCA_DIR}/pca_projected" \
    "${PLINK_FLAGS[@]}"

PROJECTED=$(tail -n +2 "${PCA_DIR}/pca_projected.sscore" | wc -l)
echo "    Projected $(fmt ${PROJECTED}) samples"
echo ""

echo "  ================================================"
echo "  PCA fit + projection complete"
echo "  ================================================"
echo "  Fit outputs:"
echo "    ${PCA_DIR}/pca_pcs.eigenval"
echo "    ${PCA_DIR}/pca_pcs.eigenvec"
echo "    ${PCA_DIR}/pca_pcs.eigenvec.allele"
echo "  Projection outputs:"
echo "    ${PCA_DIR}/pca_counts.acount"
echo "    ${PCA_DIR}/pca_projected.sscore"
