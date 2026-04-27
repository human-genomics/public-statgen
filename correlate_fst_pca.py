"""
correlate_fst_pca.py — Correlate pairwise Hudson FST with squared PC distance.

Theory: pairwise FST and squared Euclidean distance between population
centroids in PC space are proportional under McVean (2009). A linear plot
verifies this; outliers are pairs where one metric disagrees with the other.

Uses unpruned FST (literature-comparable) and the top N PCs from the panel.

Output: pca/plots/fst_vs_pcdist_{pop_n5, pop_n10, pop_n20, supervised}.png
        plus outlier listings to stdout.
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
PLOTS_DIR = PCA_DIR / "plots"
PLOTS_DIR.mkdir(parents=True, exist_ok=True)

N_PCS_FOR_DIST = int(os.environ.get("N_PCS_FOR_DIST", "20"))
OUTLIER_SD = float(os.environ.get("OUTLIER_SD_THRESHOLD", "2.5"))


# ---------------------------------------------------------------------------
# Load
# ---------------------------------------------------------------------------
sscore = pd.read_csv(PCA_DIR / "pca_projected.sscore", sep="\t").rename(columns={"#FID": "FID"})
pc_cols = [f"PC{i}_AVG" for i in range(1, N_PCS_FOR_DIST + 1)]

metadata = pd.read_csv(PROJECT_DIR / "summary" / "metadata.csv")
supervised = pd.read_csv(PROJECT_DIR / "summary" / "supervised.csv")
fst = pd.read_csv(PCA_DIR / "fst_pairs" / "fst_summary.tsv", sep="\t")

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
    y = np.median(X, axis=0)  # robust initialization
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


def group_centroids(frame, group_col, pcs):
    out = {}
    for g in sorted(frame[group_col].dropna().unique()):
        out[g] = geometric_median(frame.loc[frame[group_col] == g, pcs].values)
    return pd.DataFrame(out, index=pcs).T


pop_centroids = group_centroids(df, "population_id", pc_cols)
pop_counts = df["population_id"].value_counts()

sup_groups = sorted(supervised["reference_population"].unique())
sup_data = {sup: geometric_median(
    df.loc[df["IID"].isin(supervised.loc[supervised["reference_population"] == sup, "sample_id"]),
           pc_cols].values)
    for sup in sup_groups}
sup_centroids = pd.DataFrame(sup_data, index=pc_cols).T
sup_centroids.index.name = "group"
sup_counts = supervised["reference_population"].value_counts()


# ---------------------------------------------------------------------------
# Helper: build pairs DataFrame with FST + squared PC distance
# ---------------------------------------------------------------------------
def build_pairs(centroids, counts, fst_subset):
    entities = centroids.index.tolist()
    fst_lookup = {}
    for _, r in fst_subset.iterrows():
        fst_lookup[(r.entity1, r.entity2)] = r.fst_hudson
        fst_lookup[(r.entity2, r.entity1)] = r.fst_hudson
    rows = []
    for i, e1 in enumerate(entities):
        for e2 in entities[i + 1:]:
            f = fst_lookup.get((e1, e2))
            if f is None:
                continue
            d2 = float(((centroids.loc[e1] - centroids.loc[e2]) ** 2).sum())
            rows.append({"e1": e1, "e2": e2,
                         "n1": int(counts.get(e1, 0)),
                         "n2": int(counts.get(e2, 0)),
                         "sqdist": d2, "fst": f})
    return pd.DataFrame(rows)


def make_plot(pairs, title, out_path, label_top_n=12, outlier_sd=OUTLIER_SD):
    """Scatter pairs, fit line, mark outliers, save figure."""
    slope, intercept = np.polyfit(pairs["sqdist"], pairs["fst"], 1)
    pred = slope * pairs["sqdist"] + intercept
    resid = pairs["fst"] - pred
    pairs = pairs.copy()
    pairs["pred"] = pred
    pairs["resid"] = resid
    pairs["resid_sd"] = (resid - resid.mean()) / resid.std()
    pairs["abs_resid_sd"] = pairs["resid_sd"].abs()

    r = float(np.corrcoef(pairs["sqdist"], pairs["fst"])[0, 1])
    is_outlier = pairs["abs_resid_sd"] > outlier_sd

    fig, ax = plt.subplots(figsize=(10, 7))

    # Reference line
    xs = np.linspace(pairs["sqdist"].min(), pairs["sqdist"].max(), 200)
    ax.plot(xs, slope * xs + intercept, "-", color="black", linewidth=1.0, alpha=0.65,
            label=f"linear fit  r={r:.3f}  slope={slope:.3g}", zorder=2)

    # Non-outliers
    nope = pairs[~is_outlier]
    ax.scatter(nope["sqdist"], nope["fst"], s=18, alpha=0.55, c="#0072B2",
               edgecolors="white", linewidths=0.3,
               label=f"{len(pairs)} pairs", zorder=3)

    # Outliers
    yes = pairs[is_outlier]
    if len(yes):
        ax.scatter(yes["sqdist"], yes["fst"], s=32, alpha=0.95,
                   c="#D55E00", edgecolors="white", linewidths=0.4,
                   label=f"outliers (|resid|>{outlier_sd}σ): {len(yes)}", zorder=4)
        # annotate the most extreme N outliers
        labelable = yes.nlargest(min(label_top_n, len(yes)), "abs_resid_sd")
        for _, row in labelable.iterrows():
            ax.annotate(f"{row.e1} ↔ {row.e2}",
                        (row.sqdist, row.fst),
                        fontsize=6.5, alpha=0.9,
                        xytext=(4, 4), textcoords="offset points")

    ax.set_xlabel(f"squared Euclidean distance between centroids  (top {N_PCS_FOR_DIST} PCs)")
    ax.set_ylabel("Hudson FST  (unpruned variants)")
    ax.set_title(title, fontsize=11)
    ax.legend(fontsize=9, loc="upper left", framealpha=0.9)
    ax.grid(True, alpha=0.2, linewidth=0.5)
    ax.set_axisbelow(True)

    fig.tight_layout()
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)

    return pairs, is_outlier, r, slope


def report_outliers(label, pairs, is_outlier, top_n=20):
    outs = pairs[is_outlier].sort_values("abs_resid_sd", ascending=False).head(top_n)
    print(f"\n  ── {label}: {int(is_outlier.sum())} outliers (|resid|>{OUTLIER_SD}σ) ──")
    if len(outs) == 0:
        print("    (none)")
        return
    for _, r in outs.iterrows():
        direction = "above" if r.resid > 0 else "below"
        print(f"    {r.e1:<16} ↔ {r.e2:<16}  Fst={r.fst:.4f}  sqdist={r.sqdist:.4g}  "
              f"|resid|={r.abs_resid_sd:.2f}σ ({direction} fit)  n=({r.n1},{r.n2})")


# ---------------------------------------------------------------------------
# Run
# ---------------------------------------------------------------------------
print(f"  Using top {N_PCS_FOR_DIST} PCs, outlier threshold |resid| > {OUTLIER_SD}σ")
print(f"  FST source: pca/fst_pairs_unpruned/fst_summary.tsv")
print()

fst_pop = fst[fst.pair_type == "pop_pop"]
fst_sup = fst[fst.pair_type == "sup_sup"]

for min_n in [5, 10, 20]:
    pops = pop_counts[pop_counts >= min_n].index.tolist()
    cent = pop_centroids.loc[pops]
    pairs = build_pairs(cent, pop_counts, fst_pop)
    out_path = PLOTS_DIR / f"fst_vs_pcdist_pop_n{min_n}.png"
    p, is_o, r, slope = make_plot(
        pairs,
        title=f"FST vs squared PC distance — populations with n ≥ {min_n}  ({len(pairs)} pairs)",
        out_path=out_path,
    )
    print(f"  Saved {out_path.relative_to(PROJECT_DIR)}  (r={r:.3f}, n_pairs={len(pairs)}, n_outliers={int(is_o.sum())})")
    report_outliers(f"populations with n ≥ {min_n}", p, is_o)
    print()

# Supervised
pairs = build_pairs(sup_centroids, sup_counts, fst_sup)
out_path = PLOTS_DIR / "fst_vs_pcdist_supervised.png"
p, is_o, r, slope = make_plot(
    pairs,
    title=f"FST vs squared PC distance — 6 supervised groups  ({len(pairs)} pairs)",
    out_path=out_path,
    label_top_n=15,
    outlier_sd=1.5,
)
print(f"  Saved {out_path.relative_to(PROJECT_DIR)}  (r={r:.3f}, n_pairs={len(pairs)}, n_outliers={int(is_o.sum())})")
report_outliers("6 supervised groups (using 1.5σ threshold)", p, is_o)
print()
print("Done.")
