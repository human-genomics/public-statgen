"""
compute_fst.py — Pairwise Hudson Fst between populations and supervised groups.

For all pairs of:
  - populations (population_id) with n >= FST_MIN_POP_SIZE in the PCA panel
  - supervised reference populations (from summary/supervised.csv)

Uses plink2's `--fst CATPHENO method=hudson`. Three batches:
  1. Pop × Pop      — single plink2 call (no sample-level overlap)
  2. Sup × Sup      — single plink2 call (no sample-level overlap)
  3. Pop × Sup      — one call per population. Each call's phenotype:
       samples in pop p          → pop_<p>
       samples in sup_X but NOT in p → sup_<X>
     Plink2 returns 21 pairs per call (1 pop + 6 sup → C(7,2)); we keep
     only the 6 pop_<p> × sup_<X> rows (sup×sup duplicates batch 2).

Expected environment (set by pca.sh):
  PCA_DIR, PLINK2, PLINK_MEMORY, PLINK_THREADS, FST_MIN_POP_SIZE

Outputs (in PCA_DIR/fst_pairs/):
  fst_summary.tsv  — long-format master, all pairs, with n1/n2 sample sizes
  fst_matrix.tsv   — symmetric matrix view (NaN diagonal)
  entities.tsv     — list of entities with type and panel sample count
  raw/             — raw plink2 outputs
"""

import os
import subprocess
from pathlib import Path

import pandas as pd


PROJECT_DIR = Path(os.path.dirname(os.path.abspath(__file__)))
PCA_DIR = Path(os.environ["PCA_DIR"])
PLINK2 = os.environ["PLINK2"]
PLINK_MEMORY = os.environ.get("PLINK_MEMORY", "8000")
PLINK_THREADS = os.environ.get("PLINK_THREADS", "4")
MIN_POP_SIZE = int(os.environ.get("FST_MIN_POP_SIZE", "5"))

# Bfile and output subdir are configurable so the same script can compute FST
# on the LD-pruned PCA panel and on a less-pruned variant set.
FST_BFILE = Path(os.environ.get("FST_BFILE", str(PCA_DIR / "pca_qc")))
FST_SUBDIR = os.environ.get("FST_OUTPUT_SUBDIR", "fst_pairs")

FST_DIR = PCA_DIR / FST_SUBDIR
RAW_DIR = FST_DIR / "raw"
FST_DIR.mkdir(parents=True, exist_ok=True)
RAW_DIR.mkdir(parents=True, exist_ok=True)

PLINK_FLAGS = ["--memory", PLINK_MEMORY, "--threads", PLINK_THREADS,
               "--no-input-missing-phenotype"]


def run_plink2(*args, quiet=True):
    cmd = [PLINK2] + list(args) + PLINK_FLAGS
    res = subprocess.run(cmd, capture_output=True, text=True)
    if res.returncode != 0:
        print("plink2 failed:")
        print(res.stdout)
        print(res.stderr)
        raise SystemExit(1)


def safe(s):
    return s.replace(" ", "_") if s else s


def write_pheno(path, mapping):
    """mapping: dict iid -> category. Skips iids with None."""
    rows = [(0, iid, cat) for iid, cat in mapping.items() if cat is not None]
    df = pd.DataFrame(rows, columns=["#FID", "IID", "CAT"])
    df.to_csv(path, sep="\t", index=False)


def parse_fst_summary(path):
    """Parse plink2 .fst.summary into DataFrame[E1, E2, fst_hudson]."""
    df = pd.read_csv(path, sep="\t")
    cols = list(df.columns)
    e1 = cols[0]
    e2 = cols[1]
    fst_col = next((c for c in cols if c.upper().endswith("FST") or "HUDSON" in c.upper()), cols[-1])
    out = df[[e1, e2, fst_col]].rename(columns={e1: "E1", e2: "E2", fst_col: "fst_hudson"})
    return out


# ---------------------------------------------------------------------------
# Load and intersect with PCA panel
# ---------------------------------------------------------------------------
fam = pd.read_csv(f"{FST_BFILE}.fam", sep=r"\s+", header=None,
                  names=["FID", "IID", "PAT", "MAT", "SEX", "PHENO"])
panel_ids = set(fam["IID"])

metadata = pd.read_csv(PROJECT_DIR / "summary" / "metadata.csv")
supervised = pd.read_csv(PROJECT_DIR / "summary" / "supervised.csv")

metadata_p = metadata[metadata["sample_id"].isin(panel_ids)]
supervised_p = supervised[supervised["sample_id"].isin(panel_ids)]

pop_counts = metadata_p["population_id"].value_counts(dropna=False)
qualifying_pops = sorted(pop_counts[pop_counts >= MIN_POP_SIZE].index)
sup_groups = sorted(supervised_p["reference_population"].unique())

n_pp = len(qualifying_pops) * (len(qualifying_pops) - 1) // 2
n_ss = len(sup_groups) * (len(sup_groups) - 1) // 2
n_ps = len(qualifying_pops) * len(sup_groups)

print(f"  Populations with n >= {MIN_POP_SIZE}: {len(qualifying_pops)}")
print(f"  Supervised groups:                {len(sup_groups)}")
print(f"  Total entities:                   {len(qualifying_pops) + len(sup_groups)}")
print(f"  Total pairs:                      {n_pp + n_ss + n_ps:,}")
print(f"    Pop × Pop: {n_pp:,}   Sup × Sup: {n_ss}   Pop × Sup: {n_ps:,}")
print()

sample_to_pop = dict(zip(metadata_p["sample_id"], metadata_p["population_id"]))
sample_to_sup = dict(zip(supervised_p["sample_id"], supervised_p["reference_population"]))
qualifying_set = set(qualifying_pops)

sup_safe = {s: safe(s) for s in sup_groups}
inv_sup_safe = {v: k for k, v in sup_safe.items()}

pop_n = pop_counts.to_dict()
sup_n = supervised_p["reference_population"].value_counts().to_dict()


# ---------------------------------------------------------------------------
# Batch 1 — Pop × Pop
# ---------------------------------------------------------------------------
print(f"Batch 1/3 — Pop × Pop ({n_pp:,} pairs, 1 plink2 call) ...")
mapping = {iid: pop for iid, pop in sample_to_pop.items() if pop in qualifying_set}
pheno1 = RAW_DIR / "pop_pheno.tsv"
write_pheno(pheno1, mapping)
out1 = RAW_DIR / "pop_x_pop"
run_plink2(
    "--bfile", str(FST_BFILE),
    "--pheno", str(pheno1),
    "--pheno-name", "CAT",
    "--fst", "CAT", "method=hudson",
    "--out", str(out1),
)


# ---------------------------------------------------------------------------
# Batch 2 — Sup × Sup
# ---------------------------------------------------------------------------
print(f"Batch 2/3 — Sup × Sup ({n_ss} pairs, 1 plink2 call) ...")
mapping = {iid: sup_safe[sup] for iid, sup in sample_to_sup.items()}
pheno2 = RAW_DIR / "sup_pheno.tsv"
write_pheno(pheno2, mapping)
out2 = RAW_DIR / "sup_x_sup"
run_plink2(
    "--bfile", str(FST_BFILE),
    "--pheno", str(pheno2),
    "--pheno-name", "CAT",
    "--fst", "CAT", "method=hudson",
    "--out", str(out2),
)


# ---------------------------------------------------------------------------
# Batch 3 — Pop × Sup (one call per population)
# ---------------------------------------------------------------------------
print(f"Batch 3/3 — Pop × Sup ({n_ps:,} pairs, {len(qualifying_pops)} plink2 calls) ...")
for i, pop in enumerate(qualifying_pops, 1):
    pop_safe_ = safe(pop)
    mapping = {}
    for iid, p in sample_to_pop.items():
        if p == pop:
            mapping[iid] = f"pop_{pop_safe_}"
    for iid, sup in sample_to_sup.items():
        if iid in mapping:
            continue
        mapping[iid] = f"sup_{sup_safe[sup]}"
    pheno3 = RAW_DIR / f"pop_x_sup_{pop_safe_}.pheno.tsv"
    write_pheno(pheno3, mapping)
    out3 = RAW_DIR / f"pop_x_sup_{pop_safe_}"
    run_plink2(
        "--bfile", str(FST_BFILE),
        "--pheno", str(pheno3),
        "--pheno-name", "CAT",
        "--fst", "CAT", "method=hudson",
        "--out", str(out3),
    )
    if i % 10 == 0 or i == len(qualifying_pops):
        print(f"  [{i}/{len(qualifying_pops)}]")


# ---------------------------------------------------------------------------
# Aggregate
# ---------------------------------------------------------------------------
print("\nAggregating ...")

rows = []

# Batch 1: pop × pop (no overlap)
for _, r in parse_fst_summary(f"{out1}.fst.summary").iterrows():
    rows.append({
        "pair_type":    "pop_pop",
        "entity1":      r["E1"], "entity1_type": "population", "n1": pop_n.get(r["E1"], 0),
        "entity2":      r["E2"], "entity2_type": "population", "n2": pop_n.get(r["E2"], 0),
        "fst_hudson":   r["fst_hudson"],
    })

# Batch 2: sup × sup (full supervised counts)
for _, r in parse_fst_summary(f"{out2}.fst.summary").iterrows():
    e1 = inv_sup_safe.get(r["E1"], r["E1"])
    e2 = inv_sup_safe.get(r["E2"], r["E2"])
    rows.append({
        "pair_type":    "sup_sup",
        "entity1":      e1, "entity1_type": "supervised", "n1": sup_n.get(e1, 0),
        "entity2":      e2, "entity2_type": "supervised", "n2": sup_n.get(e2, 0),
        "fst_hudson":   r["fst_hudson"],
    })

# Batch 3: pop × sup. n_sup_eff = sup_n[X] minus samples in pop_n[P] also in sup_X.
for pop in qualifying_pops:
    pop_safe_ = safe(pop)
    pop_label = f"pop_{pop_safe_}"
    df = parse_fst_summary(f"{RAW_DIR / f'pop_x_sup_{pop_safe_}'}.fst.summary")
    for _, r in df.iterrows():
        e1, e2 = r["E1"], r["E2"]
        if e1 == pop_label and e2.startswith("sup_"):
            sup_name = inv_sup_safe[e2[4:]]
        elif e2 == pop_label and e1.startswith("sup_"):
            sup_name = inv_sup_safe[e1[4:]]
        else:
            continue
        # n_sup excluding overlap with this population
        overlap = sum(1 for iid, ssup in sample_to_sup.items()
                      if ssup == sup_name and sample_to_pop.get(iid) == pop)
        n_sup_eff = sup_n.get(sup_name, 0) - overlap
        rows.append({
            "pair_type":    "pop_sup",
            "entity1":      pop, "entity1_type": "population", "n1": pop_n[pop],
            "entity2":      sup_name, "entity2_type": "supervised", "n2": n_sup_eff,
            "fst_hudson":   r["fst_hudson"],
        })

summary_df = pd.DataFrame(rows)
summary_df = summary_df[["pair_type", "entity1", "entity1_type", "n1",
                          "entity2", "entity2_type", "n2", "fst_hudson"]]
summary_path = FST_DIR / "fst_summary.tsv"
summary_df.to_csv(summary_path, sep="\t", index=False)
print(f"  Wrote {summary_path.relative_to(PROJECT_DIR)} ({len(summary_df):,} rows)")

# Note: we do not write fst_matrix.tsv here. The calibrated matrix is saved
# by calibrate_and_fit_fst.py once the calibration slope (vs Privé 2022) is
# known, so that fst_matrix.tsv always contains literature-scale values.

# Entities table
entities = []
for p in qualifying_pops:
    entities.append({"entity": p, "type": "population", "n": pop_n[p]})
for s in sup_groups:
    entities.append({"entity": s, "type": "supervised", "n": sup_n.get(s, 0)})
entities_df = pd.DataFrame(entities)
ent_path = FST_DIR / "entities.tsv"
entities_df.to_csv(ent_path, sep="\t", index=False)
print(f"  Wrote {ent_path.relative_to(PROJECT_DIR)}")

print("\nDone.")
