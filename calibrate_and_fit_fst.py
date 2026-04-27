"""
calibrate_and_fit_fst.py — Calibrate Hudson FST against Privé et al. 2022 1KG
reference, then fit a predictive equation FST = a·d² + b where d² is the
squared Euclidean distance between population centroids in PC space.

Pipeline:
  1. Match our unpruned FST against literature_reference/fst_prive_2022.csv
  2. Linear fit: ours = slope * ref + intercept  (we got slope ≈ 1.155)
  3. Calibrate: fst_cal = fst_unpruned / slope
  4. Compute centroids for populations (n ≥ FST_FIT_MIN_POP_SIZE) and
     for the 6 supervised groups, on top N_PCS_FOR_DIST PCs.
  5. Fit a single linear model FST_cal ~ d² across all pop+sup pairs.
  6. Make plot (supervised pairs highlighted in purple, all labeled).
  7. Save the equation as a plain-text file usable by anyone with PC scores.

Outputs:
  pca/fst_pairs/fst_summary_calibrated.tsv
  pca/plots/fst_vs_pcdist_calibrated.png
  pca/centroids_top<N>_pop.tsv
  pca/centroids_top<N>_supervised.tsv
  pca/fst_pcdist_equation.txt
"""

import os
import datetime
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


PROJECT_DIR = Path(os.path.dirname(os.path.abspath(__file__)))
PCA_DIR = Path(os.environ["PCA_DIR"])
N_PCS_FOR_DIST = int(os.environ.get("N_PCS_FOR_DIST", "20"))
MIN_POP_SIZE_FIT = int(os.environ.get("FST_FIT_MIN_POP_SIZE", "10"))

PLOTS_DIR = PCA_DIR / "plots"
PLOTS_DIR.mkdir(exist_ok=True)


# ---------------------------------------------------------------------------
# Load reference + our unpruned FST
# ---------------------------------------------------------------------------
ref = pd.read_csv(PROJECT_DIR / "literature_reference" / "fst_prive_2022.csv")
ref_lookup = {tuple(sorted([r.pop1, r.pop2])): r.fst_hudson for _, r in ref.iterrows()}

ours_all = pd.read_csv(PCA_DIR / "fst_pairs" / "fst_summary.tsv", sep="\t")
ours_pop = ours_all[ours_all.pair_type == "pop_pop"].copy()
ours_sup = ours_all[ours_all.pair_type == "sup_sup"].copy()


# ---------------------------------------------------------------------------
# Calibration: ours = slope * ref + intercept
# ---------------------------------------------------------------------------
matched = [(ref_lookup[k], r.fst_hudson)
           for _, r in ours_pop.iterrows()
           if (k := tuple(sorted([r.entity1, r.entity2]))) in ref_lookup]
m_arr = np.array(matched)
slope, intercept = np.polyfit(m_arr[:, 0], m_arr[:, 1], 1)
calib_r = float(np.corrcoef(m_arr[:, 0], m_arr[:, 1])[0, 1])

print(f"  Calibration vs Privé et al. 2022:")
print(f"    Matched pairs:   {len(matched)}")
print(f"    ours_unpruned = {slope:.4f} * ref + {intercept:.5f}")
print(f"    Pearson r:       {calib_r:.4f}")
print(f"    Calibration:     fst_cal = fst_unpruned / {slope:.4f}")
print()

ours_all = ours_all.copy()
ours_all["fst_calibrated"] = ours_all["fst_hudson"] / slope
ours_pop["fst_calibrated"] = ours_pop["fst_hudson"] / slope
ours_sup["fst_calibrated"] = ours_sup["fst_hudson"] / slope

calib_path = PCA_DIR / "fst_pairs" / "fst_summary_calibrated.tsv"
ours_all.to_csv(calib_path, sep="\t", index=False)
print(f"  Saved {calib_path.relative_to(PROJECT_DIR)}")

# Calibrated symmetric matrix (this is the canonical externally-shareable FST table)
entities = pd.read_csv(PCA_DIR / "fst_pairs" / "entities.tsv", sep="\t")
entity_order = entities["entity"].tolist()
matrix = pd.DataFrame(index=entity_order, columns=entity_order, dtype=float)
for _, r in ours_all.iterrows():
    if r.entity1 in matrix.index and r.entity2 in matrix.columns:
        matrix.loc[r.entity1, r.entity2] = r.fst_calibrated
        matrix.loc[r.entity2, r.entity1] = r.fst_calibrated
matrix_path = PCA_DIR / "fst_pairs" / "fst_matrix.tsv"
matrix.to_csv(matrix_path, sep="\t")
print(f"  Saved {matrix_path.relative_to(PROJECT_DIR)} ({matrix.shape[0]}x{matrix.shape[1]}, calibrated)")


# ---------------------------------------------------------------------------
# Calibration plot: ours vs Privé reference (raw and calibrated)
# ---------------------------------------------------------------------------
ref_x = m_arr[:, 0]
raw_y = m_arr[:, 1]
cal_y = raw_y / slope
ref_pairs_df = pd.DataFrame({
    "ref": ref_x, "raw": raw_y, "calibrated": cal_y,
    "p1": [k[0] for _, r in ours_pop.iterrows()
           if (k := tuple(sorted([r.entity1, r.entity2]))) in ref_lookup],
    "p2": [k[1] for _, r in ours_pop.iterrows()
           if (k := tuple(sorted([r.entity1, r.entity2]))) in ref_lookup],
})

fig, axes = plt.subplots(2, 2, figsize=(13, 11))
lim = max(ref_x.max(), raw_y.max()) * 1.05


def panel_scatter(ax, y, version, title_suffix):
    s, b = np.polyfit(ref_x, y, 1)
    bias = float(np.mean(y - ref_x))
    mae = float(np.mean(np.abs(y - ref_x)))
    rmse = float(np.sqrt(np.mean((y - ref_x) ** 2)))
    pearson_p = float(np.corrcoef(ref_x, y)[0, 1])

    ax.scatter(ref_x, y, s=20, alpha=0.55, c="#0072B2",
               edgecolors="white", linewidths=0.3, label=f"{len(ref_x)} pairs")
    ax.plot([0, lim], [0, lim], "--", color="#888", linewidth=1, alpha=0.6, label="y = x")
    xs = np.linspace(0, lim, 100)
    ax.plot(xs, s * xs + b, "-", color="black", linewidth=1, alpha=0.8,
            label=f"linear fit  slope={s:.3f}")

    diffs = y - ref_x
    top = np.argsort(np.abs(diffs))[-6:]
    ax.scatter(ref_x[top], y[top], s=46, alpha=0.95, c="#D55E00",
               edgecolors="white", linewidths=0.5)
    for i in top:
        ax.annotate(f"{ref_pairs_df.iloc[i].p1} ↔ {ref_pairs_df.iloc[i].p2}",
                    (ref_x[i], y[i]), fontsize=6.5, alpha=0.9,
                    xytext=(4, 4), textcoords="offset points")
    ax.set_xlabel("Reference FST  (Privé et al. 2022)")
    ax.set_ylabel(f"Our FST  ({version})")
    ax.set_title(f"{version} vs reference{title_suffix}\n"
                 f"r={pearson_p:.4f}, slope={s:.3f}, bias={bias:+.4f}, MAE={mae:.4f}, RMSE={rmse:.4f}",
                 fontsize=10)
    ax.legend(fontsize=8, loc="upper left")
    ax.grid(True, alpha=0.2)
    ax.set_xlim(0, lim); ax.set_ylim(0, lim); ax.set_aspect("equal")
    return diffs, top


def panel_bland_altman(ax, diffs, version, top):
    bias = float(np.mean(diffs))
    sd = float(np.std(diffs))
    ax.scatter(ref_x, diffs, s=20, alpha=0.6, c="#0072B2",
               edgecolors="white", linewidths=0.3)
    ax.axhline(0, color="#888", linestyle="--", linewidth=1)
    ax.axhline(bias, color="black", linestyle="-", linewidth=1, alpha=0.8,
               label=f"mean bias = {bias:+.4f}")
    ax.axhline(bias + 1.96 * sd, color="#888", linestyle=":", linewidth=0.8,
               label=f"±1.96σ")
    ax.axhline(bias - 1.96 * sd, color="#888", linestyle=":", linewidth=0.8)
    ax.scatter(ref_x[top], diffs[top], s=46, alpha=0.95, c="#D55E00",
               edgecolors="white", linewidths=0.5)
    for i in top:
        ax.annotate(f"{ref_pairs_df.iloc[i].p1} ↔ {ref_pairs_df.iloc[i].p2}",
                    (ref_x[i], diffs[i]), fontsize=6.5, alpha=0.9,
                    xytext=(4, 4), textcoords="offset points")
    ax.set_xlabel("Reference FST  (Privé et al. 2022)")
    ax.set_ylabel(f"Our FST − Reference   ({version})")
    ax.set_title(f"Bland-Altman: {version} − reference", fontsize=10)
    ax.legend(fontsize=8, loc="upper right")
    ax.grid(True, alpha=0.2)


diffs_raw, top_raw = panel_scatter(axes[0, 0], raw_y, "ours_unpruned (raw)",
                                    "  —  before calibration")
panel_bland_altman(axes[0, 1], diffs_raw, "ours_unpruned (raw)", top_raw)
diffs_cal, top_cal = panel_scatter(axes[1, 0], cal_y, "ours_calibrated",
                                    f"  —  after dividing by {slope:.4f}")
panel_bland_altman(axes[1, 1], diffs_cal, "ours_calibrated", top_cal)

fig.suptitle(f"Calibration of our Hudson FST against Privé et al. 2022  "
             f"({len(ref_x)} matched 1KG pairs)", fontsize=12, y=1.001)
fig.tight_layout()
calib_plot = PLOTS_DIR / "fst_calibration_vs_reference.png"
fig.savefig(calib_plot, dpi=200, bbox_inches="tight")
plt.close(fig)
print(f"  Saved {calib_plot.relative_to(PROJECT_DIR)}")


# ---------------------------------------------------------------------------
# Centroids
# ---------------------------------------------------------------------------
sscore = pd.read_csv(PCA_DIR / "pca_projected.sscore", sep="\t").rename(columns={"#FID": "FID"})
pc_cols = [f"PC{i}_AVG" for i in range(1, N_PCS_FOR_DIST + 1)]

metadata = pd.read_csv(PROJECT_DIR / "summary" / "metadata.csv")
supervised = pd.read_csv(PROJECT_DIR / "summary" / "supervised.csv")

df = sscore.merge(metadata[["sample_id", "population_id", "superpopulation"]],
                  left_on="IID", right_on="sample_id")


def geometric_median(X, eps=1e-8, max_iter=500):
    """Weiszfeld's algorithm for the geometric median.

    Minimizes the sum of Euclidean distances to all points; robust to outliers
    (Privé et al. 2022 use this for population centers in PC space)."""
    X = np.asarray(X, dtype=float)
    if X.ndim == 1:
        return X.copy()
    n = len(X)
    if n == 1:
        return X[0].copy()
    if n == 2:
        return X.mean(axis=0)
    y = np.median(X, axis=0)
    for _ in range(max_iter):
        d = np.linalg.norm(X - y, axis=1)
        nz = d > eps
        if not nz.any():
            return y
        w = 1.0 / d[nz]
        y_new = (w[:, None] * X[nz]).sum(axis=0) / w.sum()
        if np.linalg.norm(y_new - y) < eps:
            return y_new
        y = y_new
    return y


def _gm_centroids(frame, group_col, pcs):
    out = {}
    for g in sorted(frame[group_col].dropna().unique()):
        out[g] = geometric_median(frame.loc[frame[group_col] == g, pcs].values)
    return pd.DataFrame(out, index=pcs).T


pop_centroids = _gm_centroids(df, "population_id", pc_cols)
pop_counts = df["population_id"].value_counts()

sup_groups = sorted(supervised["reference_population"].unique())
sup_data = {sup: geometric_median(
    df.loc[df["IID"].isin(supervised.loc[supervised["reference_population"] == sup, "sample_id"]),
           pc_cols].values)
    for sup in sup_groups}
sup_centroids = pd.DataFrame(sup_data, index=pc_cols).T
sup_counts = supervised["reference_population"].value_counts()

# Save centroids with sample sizes as the leading column
pop_cent_path = PCA_DIR / f"centroids_top{N_PCS_FOR_DIST}_pop.tsv"
sup_cent_path = PCA_DIR / f"centroids_top{N_PCS_FOR_DIST}_supervised.tsv"

pop_out = pop_centroids.copy()
pop_out.insert(0, "n", pop_out.index.map(pop_counts).astype(int))
pop_out.index.name = "population_id"
pop_out.to_csv(pop_cent_path, sep="\t")

sup_out = sup_centroids.copy()
sup_out.insert(0, "n", sup_out.index.map(sup_counts).astype(int))
sup_out.index.name = "supervised_group"
sup_out.to_csv(sup_cent_path, sep="\t")

print(f"  Saved {pop_cent_path.relative_to(PROJECT_DIR)} ({len(pop_out)} populations)")
print(f"  Saved {sup_cent_path.relative_to(PROJECT_DIR)} ({len(sup_out)} supervised groups)")


# ---------------------------------------------------------------------------
# Build pairs (calibrated FST, squared PC distance)
# ---------------------------------------------------------------------------
def build_pairs(centroids, counts, fst_df, label):
    fst_lk = {tuple(sorted([r.entity1, r.entity2])): r.fst_calibrated for _, r in fst_df.iterrows()}
    rows = []
    ents = centroids.index.tolist()
    for i, e1 in enumerate(ents):
        for e2 in ents[i + 1:]:
            k = tuple(sorted([e1, e2]))
            if k not in fst_lk:
                continue
            d2 = float(((centroids.loc[e1] - centroids.loc[e2]) ** 2).sum())
            rows.append({"e1": e1, "e2": e2,
                         "n1": int(counts.get(e1, 0)), "n2": int(counts.get(e2, 0)),
                         "sqdist": d2, "fst": fst_lk[k],
                         "pair_type": label})
    return pd.DataFrame(rows)


qualifying = pop_counts[pop_counts >= MIN_POP_SIZE_FIT].index.tolist()
pop_pairs = build_pairs(pop_centroids.loc[qualifying], pop_counts, ours_pop, "pop_pop")
sup_pairs = build_pairs(sup_centroids, sup_counts, ours_sup, "sup_sup")
combined = pd.concat([pop_pairs, sup_pairs], ignore_index=True)
print(f"\n  Fitting on {len(combined)} pairs "
      f"({len(pop_pairs)} pop with n≥{MIN_POP_SIZE_FIT} + {len(sup_pairs)} supervised)")


# ---------------------------------------------------------------------------
# Linear fit: FST_cal = a * d² + b
# ---------------------------------------------------------------------------
fit_a, fit_b = np.polyfit(combined["sqdist"], combined["fst"], 1)
pearson = float(np.corrcoef(combined["sqdist"], combined["fst"])[0, 1])
pred = fit_a * combined["sqdist"] + fit_b
ss_res = float(np.sum((combined["fst"] - pred) ** 2))
ss_tot = float(np.sum((combined["fst"] - combined["fst"].mean()) ** 2))
r_squared = 1 - ss_res / ss_tot
rmse = float(np.sqrt(np.mean((combined["fst"] - pred) ** 2)))

print(f"\n  predicted FST = {fit_a:.4f} × d² + {fit_b:.5f}")
print(f"  Pearson r: {pearson:.4f}   R²: {r_squared:.4f}   RMSE: {rmse:.5f}")


# ---------------------------------------------------------------------------
# Plot
# ---------------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(11, 8))

# Fit line first (under everything)
xs = np.linspace(0, combined["sqdist"].max() * 1.04, 200)
ax.plot(xs, fit_a * xs + fit_b, "-", color="black", linewidth=1.2, alpha=0.65,
        label=f"fit:  FST = {fit_a:.3f}·d² + {fit_b:.5f}    (R² = {r_squared:.3f})",
        zorder=2)

# Population pairs (blue)
ax.scatter(pop_pairs["sqdist"], pop_pairs["fst"],
           s=20, alpha=0.55, c="#1f77b4",
           edgecolors="white", linewidths=0.3, zorder=3,
           label=f"population pairs (n≥{MIN_POP_SIZE_FIT}):  {len(pop_pairs)}")

# Supervised group pairs (purple, larger, labeled)
PURPLE = "#7B3FB8"
ax.scatter(sup_pairs["sqdist"], sup_pairs["fst"],
           s=110, alpha=0.95, c=PURPLE,
           edgecolors="black", linewidths=0.7, marker="o", zorder=5,
           label=f"supervised super-pop pairs:  {len(sup_pairs)}")

for _, r in sup_pairs.iterrows():
    ax.annotate(f"{r.e1} ↔ {r.e2}",
                (r.sqdist, r.fst),
                fontsize=7.8, color="#3D2570", fontweight="bold",
                xytext=(7, 5), textcoords="offset points", zorder=6)

# Top 3 highest-FST pairs in the combined dataset
TOP_RED = "#C0392B"
top3 = combined.nlargest(3, "fst").reset_index(drop=True)
ax.scatter(top3["sqdist"], top3["fst"],
           s=85, alpha=0.95, c=TOP_RED,
           edgecolors="black", linewidths=0.7, marker="^", zorder=6,
           label=f"top-3 highest FST: {len(top3)}")
for _, r in top3.iterrows():
    ax.annotate(f"{r.e1} ↔ {r.e2}",
                (r.sqdist, r.fst),
                fontsize=7.8, color="#7B1F14", fontweight="bold",
                ha="right", va="bottom",
                xytext=(-8, 6), textcoords="offset points", zorder=7)

ax.set_xlabel(f"squared Euclidean distance between centroids  (top {N_PCS_FOR_DIST} PCs)",
              fontsize=11)
ax.set_ylabel("Hudson FST  (calibrated to Privé et al. 2022)", fontsize=11)
ax.set_title(
    f"FST = a · d²  +  b      —      fit on pop pairs (n ≥ {MIN_POP_SIZE_FIT}) + 6 supervised super-pop pairs\n"
    f"FST calibration: our unpruned FST ÷ {slope:.3f}  (slope vs Privé et al. 2022 1KG reference)",
    fontsize=11.5,
)
ax.legend(loc="upper left", fontsize=10, framealpha=0.94, frameon=True)
ax.grid(True, alpha=0.22, linewidth=0.55)
ax.set_axisbelow(True)
ax.set_xlim(0, combined["sqdist"].max() * 1.04)
ax.set_ylim(min(0, combined["fst"].min() - 0.005), combined["fst"].max() * 1.06)

fig.tight_layout()
fig_path = PLOTS_DIR / "fst_vs_pcdist_calibrated.png"
fig.savefig(fig_path, dpi=200, bbox_inches="tight")
plt.close(fig)
print(f"  Saved {fig_path.relative_to(PROJECT_DIR)}")


# ---------------------------------------------------------------------------
# Equation text file
# ---------------------------------------------------------------------------
eq_path = PCA_DIR / "fst_pcdist_equation.txt"
with open(eq_path, "w") as f:
    f.write("PAIRWISE HUDSON FST  —  predictive equation from this project's PCA\n")
    f.write("=" * 70 + "\n\n")
    f.write("FORMULA\n")
    f.write("-------\n")
    f.write(f"  predicted_FST  =  {fit_a:.6f}  ×  d²  +  {fit_b:.6f}\n\n")
    f.write("where\n")
    f.write(f"  d²  =  squared Euclidean distance between the two populations'\n")
    f.write(f"        geometric-median centroids in the top-{N_PCS_FOR_DIST} PC space\n")
    f.write(f"        (geometric median is the point minimizing the sum of\n")
    f.write(f"        Euclidean distances to all members; robust to outliers,\n")
    f.write(f"        following Privé et al. 2022's approach)\n\n")

    f.write("FIT QUALITY\n")
    f.write("-----------\n")
    f.write(f"  Pearson r            : {pearson:.4f}\n")
    f.write(f"  R²                   : {r_squared:.4f}\n")
    f.write(f"  RMSE                 : {rmse:.5f}\n")
    f.write(f"  N pairs in fit       : {len(combined)}\n")
    f.write(f"    population pairs   : {len(pop_pairs)}  (population n ≥ {MIN_POP_SIZE_FIT})\n")
    f.write(f"    supervised pairs   :  {len(sup_pairs)}  (all 6 super-pop pairs)\n\n")

    f.write("FST CALIBRATION (ours → literature scale)\n")
    f.write("-----------------------------------------\n")
    f.write(f"  Our unpruned Hudson FSTs were calibrated by dividing by the\n")
    f.write(f"  slope of (ours_unpruned vs Privé et al. 2022) = {slope:.4f}.\n")
    f.write(f"  Reference matched on {len(matched)} pairs (Pearson r = {calib_r:.4f}).\n")
    f.write(f"  Reference: literature_reference/fst_prive_2022.csv\n\n")

    f.write("PCA PROVENANCE\n")
    f.write("--------------\n")
    f.write(f"  PC scores from: pca/pca_projected.sscore  (3,640 samples × 30 PCs)\n")
    f.write(f"  PCA fit panel : pca/pca_qc.{{bed,bim,fam}} (~125,457 LD-pruned SNPs)\n")
    f.write(f"  Top {N_PCS_FOR_DIST} PCs used for centroid distance.\n\n")

    f.write("HOW TO USE\n")
    f.write("----------\n")
    f.write("  Given two populations P and Q, each with a set of sample-level\n")
    f.write(f"  PC scores (top {N_PCS_FOR_DIST} PCs in the same coordinate system as this\n")
    f.write("  project's PCA):\n\n")
    f.write("    1. centroid_P = geometric median of P-samples' top-N PC vectors\n")
    f.write("    2. centroid_Q = geometric median of Q-samples' top-N PC vectors\n")
    f.write("    3. d²        = sum_{i=1..N} (centroid_P[i] - centroid_Q[i])²\n")
    f.write("    4. predicted_FST = a × d² + b   (using a, b above)\n\n")
    f.write("  Centroids for known populations and supervised groups are saved as:\n")
    f.write(f"    {pop_cent_path.relative_to(PROJECT_DIR)}\n")
    f.write(f"    {sup_cent_path.relative_to(PROJECT_DIR)}\n\n")

    f.write(f"Generated: {datetime.datetime.now().isoformat(timespec='seconds')}\n")

print(f"  Saved {eq_path.relative_to(PROJECT_DIR)}")
print()
print("Done.")
