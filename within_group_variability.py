"""
within_group_variability.py — Within-group spread analysis in PC space.

For each grouping (population_id, supervised reference_population, and the
metadata.csv superpopulation column) compute the geometric-median centroid
in top-N PC space, then summarize the spread of members around that centroid
via:

  - median Euclidean distance to centroid
  - root mean square (RMS) Euclidean distance to centroid
  - max Euclidean distance to centroid

Outputs:
  pca/within_group_stats_pop.tsv                  (population_id)
  pca/within_group_stats_supervised.tsv           (6 supervised reference pops)
  pca/within_group_stats_metadata_superpop.tsv    (10 metadata superpopulations)
  pca/plots/within_group_variability_n{5,10,20}.png         (per-pop plots)
  pca/plots/within_group_variability_metadata_superpop.png  (per-region plot)
"""

import os
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


PROJECT_DIR = Path(os.path.dirname(os.path.abspath(__file__)))
PCA_DIR = Path(os.environ["PCA_DIR"])
N_PCS_FOR_DIST = int(os.environ.get("N_PCS_FOR_DIST", "20"))

PLOTS_DIR = PCA_DIR / "plots"
PLOTS_DIR.mkdir(parents=True, exist_ok=True)


def geometric_median(X, eps=1e-8, max_iter=500):
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


def spread_stats(samples, centroid):
    if len(samples) == 0:
        return {"n": 0, "median_d": 0.0, "rms_d": 0.0, "max_d": 0.0}
    distances = np.linalg.norm(samples - centroid, axis=1)
    return {
        "n": len(samples),
        "median_d": float(np.median(distances)),
        "rms_d": float(np.sqrt(np.mean(distances ** 2))),
        "max_d": float(np.max(distances)),
    }


# ---------------------------------------------------------------------------
# Load
# ---------------------------------------------------------------------------
sscore = pd.read_csv(PCA_DIR / "pca_projected.sscore", sep="\t").rename(columns={"#FID": "FID"})
pc_cols = [f"PC{i}_AVG" for i in range(1, N_PCS_FOR_DIST + 1)]
metadata = pd.read_csv(PROJECT_DIR / "summary" / "metadata.csv")
supervised = pd.read_csv(PROJECT_DIR / "summary" / "supervised.csv")

df = sscore.merge(metadata[["sample_id", "population_id", "superpopulation"]],
                  left_on="IID", right_on="sample_id")


# ---------------------------------------------------------------------------
# Per-population stats
# ---------------------------------------------------------------------------
pop_rows = []
for pop in sorted(df["population_id"].dropna().unique()):
    sub = df.loc[df["population_id"] == pop, pc_cols].values
    centroid = geometric_median(sub)
    pop_rows.append({"entity": pop, "type": "population", **spread_stats(sub, centroid)})
pop_df = pd.DataFrame(pop_rows)
pop_path = PCA_DIR / "within_group_stats_pop.tsv"
pop_df.to_csv(pop_path, sep="\t", index=False)
print(f"  Wrote {pop_path.relative_to(PROJECT_DIR)} ({len(pop_df)} populations)")


# ---------------------------------------------------------------------------
# Per-supervised stats
# ---------------------------------------------------------------------------
sup_rows = []
for sup in sorted(supervised["reference_population"].unique()):
    iids = supervised.loc[supervised["reference_population"] == sup, "sample_id"]
    sub = df.loc[df["IID"].isin(iids), pc_cols].values
    centroid = geometric_median(sub)
    sup_rows.append({"entity": sup, "type": "supervised", **spread_stats(sub, centroid)})
sup_df = pd.DataFrame(sup_rows)
sup_path = PCA_DIR / "within_group_stats_supervised.tsv"
sup_df.to_csv(sup_path, sep="\t", index=False)
print(f"  Wrote {sup_path.relative_to(PROJECT_DIR)} ({len(sup_df)} supervised super-pops)")


# ---------------------------------------------------------------------------
# Per-superpopulation stats (metadata.csv "superpopulation" column — 10 regions)
# ---------------------------------------------------------------------------
SUPERPOP_DISPLAY_ORDER = [
    "African", "Middle Eastern", "European", "West Eurasian",
    "South Asian", "Central South Asian", "Central Asian Siberian",
    "East Asian", "Oceanian", "American",
]
superpop_rows = []
for sp in sorted(df["superpopulation"].dropna().unique()):
    sub = df.loc[df["superpopulation"] == sp, pc_cols].values
    centroid = geometric_median(sub)
    superpop_rows.append({"entity": sp, "type": "metadata_superpop",
                          **spread_stats(sub, centroid)})
superpop_df = pd.DataFrame(superpop_rows)
superpop_path = PCA_DIR / "within_group_stats_metadata_superpop.tsv"
superpop_df.to_csv(superpop_path, sep="\t", index=False)
print(f"  Wrote {superpop_path.relative_to(PROJECT_DIR)} ({len(superpop_df)} superpopulations)")
print()


# ---------------------------------------------------------------------------
# Plots
# ---------------------------------------------------------------------------
PURPLE = "#7B3FB8"
BLUE = "#1f77b4"
RED = "#C0392B"


def make_plot(min_n, n_label_outliers=10):
    pop = pop_df[pop_df.n >= min_n].copy()

    fig, ax = plt.subplots(figsize=(10, 7.5))

    # y=x reference (RMS == median)
    max_v = max(pop["rms_d"].max(), sup_df["rms_d"].max()) * 1.06
    ax.plot([0, max_v], [0, max_v], "--", color="#888", linewidth=0.8, alpha=0.5,
            label="y = x", zorder=1)

    # Populations (blue)
    ax.scatter(pop["median_d"], pop["rms_d"],
               s=22, alpha=0.6, c=BLUE,
               edgecolors="white", linewidths=0.3, zorder=3,
               label=f"populations (n ≥ {min_n}):  {len(pop)}")

    # Top-N most-variable populations (red triangles, labeled top-left)
    top = pop.nlargest(min(n_label_outliers, len(pop)), "rms_d")
    ax.scatter(top["median_d"], top["rms_d"],
               s=70, alpha=0.95, c=RED,
               edgecolors="black", linewidths=0.6, marker="^", zorder=4,
               label=f"top-{len(top)} most variable")
    for _, r in top.iterrows():
        ax.annotate(r.entity, (r.median_d, r.rms_d),
                    fontsize=7.5, color="#7B1F14", fontweight="bold",
                    ha="right", va="bottom",
                    xytext=(-7, 5), textcoords="offset points", zorder=6)

    # Supervised super-pops (purple, always shown, always labeled)
    ax.scatter(sup_df["median_d"], sup_df["rms_d"],
               s=110, alpha=0.95, c=PURPLE,
               edgecolors="black", linewidths=0.7, marker="o", zorder=5,
               label=f"supervised super-pops:  {len(sup_df)}")
    for _, r in sup_df.iterrows():
        ax.annotate(r.entity, (r.median_d, r.rms_d),
                    fontsize=8, color="#3D2570", fontweight="bold",
                    xytext=(8, 5), textcoords="offset points", zorder=6)

    ax.set_xlabel(f"median Euclidean distance to centroid  (top {N_PCS_FOR_DIST} PCs)", fontsize=11)
    ax.set_ylabel(f"RMS Euclidean distance to centroid  (top {N_PCS_FOR_DIST} PCs)", fontsize=11)
    ax.set_title(
        f"Within-group PC-space variability — populations with n ≥ {min_n}\n"
        f"centroids = geometric median of members; high-variability groups in upper-right",
        fontsize=11,
    )
    ax.legend(loc="upper left", fontsize=9, framealpha=0.93)
    ax.grid(True, alpha=0.22, linewidth=0.55)
    ax.set_axisbelow(True)
    ax.set_xlim(0, max_v); ax.set_ylim(0, max_v)
    ax.set_aspect("equal")

    fig.tight_layout()
    out = PLOTS_DIR / f"within_group_variability_n{min_n}.png"
    fig.savefig(out, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved {out.relative_to(PROJECT_DIR)}  ({len(pop)} pops + {len(sup_df)} super-pops)")


for n in (5, 10, 20):
    make_plot(n)


# ---------------------------------------------------------------------------
# Metadata-superpopulation plot: each of the 10 metadata-defined regions
# ---------------------------------------------------------------------------
def make_superpop_plot():
    fig, ax = plt.subplots(figsize=(9, 7.5))
    max_v = max(superpop_df["rms_d"].max(), superpop_df["median_d"].max()) * 1.10
    ax.plot([0, max_v], [0, max_v], "--", color="#888", linewidth=0.8, alpha=0.5,
            label="y = x", zorder=1)

    # Order points by display order so colors form a region gradient
    order = [s for s in SUPERPOP_DISPLAY_ORDER if s in superpop_df.entity.values]
    cmap = plt.get_cmap("viridis")
    for i, sp in enumerate(order):
        r = superpop_df[superpop_df.entity == sp].iloc[0]
        color = cmap(i / max(1, len(order) - 1))
        ax.scatter(r.median_d, r.rms_d, s=200, alpha=0.93, c=[color],
                   edgecolors="black", linewidths=0.7, zorder=4,
                   label=f"{sp}  (n={int(r.n)})")
        ax.annotate(sp, (r.median_d, r.rms_d),
                    fontsize=9, color="black", fontweight="bold",
                    xytext=(11, 5), textcoords="offset points", zorder=6)

    ax.set_xlabel(f"median Euclidean distance to centroid  (top {N_PCS_FOR_DIST} PCs)", fontsize=11)
    ax.set_ylabel(f"RMS Euclidean distance to centroid  (top {N_PCS_FOR_DIST} PCs)", fontsize=11)
    ax.set_title(
        "Within-region PC-space variability — "
        "10 metadata superpopulations\n"
        "centroids = geometric median of all members in each region",
        fontsize=11,
    )
    ax.legend(loc="upper left", fontsize=8, framealpha=0.93, ncol=1)
    ax.grid(True, alpha=0.22, linewidth=0.55)
    ax.set_axisbelow(True)
    ax.set_xlim(0, max_v); ax.set_ylim(0, max_v)
    ax.set_aspect("equal")

    fig.tight_layout()
    out = PLOTS_DIR / "within_group_variability_metadata_superpop.png"
    fig.savefig(out, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved {out.relative_to(PROJECT_DIR)}")


make_superpop_plot()
print("\nDone.")
