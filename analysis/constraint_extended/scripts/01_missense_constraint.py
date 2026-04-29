"""
Extension E: missense and missense-deleterious constraint test for the
CGRP ligand-receptor asymmetry, using gnomAD v4.1.

Probes selection modes beyond LoF tolerance (which Extension A tested
and falsified). Decision rule and seed locked in
analysis/constraint_extended/PRE_REGISTRATION_constraint_extended.md.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT.parents[0] / "data" / "raw" / "gnomad"
PROC = ROOT.parents[0] / "data" / "processed"
RES = ROOT / "results"
RES.mkdir(parents=True, exist_ok=True)

CONSTRAINT_TSV = RAW / "gnomad.v4.1.constraint_metrics.tsv"
BG_TSV = PROC / "background_features.tsv"

LIGAND = ["CALCA", "CALCB"]
RECEPTOR = ["CALCRL", "RAMP1", "CRCP"]
ALL_TEST = LIGAND + RECEPTOR

SEED = 20260503
N_PERM = 10_000

PRIMARY_METRIC = "mis.z_score"
SECONDARY = ["mis.oe", "mis_pphen.oe"]
SANITY = "syn.z_score"


def load_constraint() -> pd.DataFrame:
    df = pd.read_csv(CONSTRAINT_TSV, sep="\t", low_memory=False)
    has_mane = df["mane_select"].fillna(False).astype(bool)
    has_can = df["canonical"].fillna(False).astype(bool)
    df_mane = df[has_mane].drop_duplicates(subset="gene", keep="first")
    df_can_only = df[has_can & ~df["gene"].isin(set(df_mane["gene"]))].drop_duplicates(
        subset="gene", keep="first"
    )
    return pd.concat([df_mane, df_can_only], ignore_index=True)


def percentile_of(value: float, dist: np.ndarray) -> float:
    return float((dist < value).mean() * 100)


def main() -> None:
    print("[1/5] Loading constraint and background...")
    df = load_constraint()
    bg = pd.read_csv(BG_TSV, sep="\t").merge(df, on="gene", how="inner")

    metrics = [PRIMARY_METRIC] + SECONDARY + [SANITY]
    print(f"      {len(bg):,} background genes; metrics: {metrics}")

    # Higher mis.z_score = more constrained; for mis.oe and mis_pphen.oe
    # lower = more constrained.
    direction = {
        "mis.z_score":  "higher_more_constrained",
        "mis.oe":        "lower_more_constrained",
        "mis_pphen.oe":  "lower_more_constrained",
        "syn.z_score":   "higher_more_constrained",
    }

    print("[2/5] Per-gene metrics for the 5 CGRP-pathway genes:")
    test_rows = bg[bg["gene"].isin(ALL_TEST)].set_index("gene")
    show_cols = ["gene"] + metrics
    print(test_rows.reset_index()[show_cols].to_string(index=False))

    results: dict = {
        "spec": "analysis/constraint_extended/PRE_REGISTRATION_constraint_extended.md",
        "seed": SEED,
        "n_permutations": N_PERM,
        "constraint_tsv_sha256":
            "68d8abdb7fc48f570869b02dfaa74b9fecaece7fcc5f301ddca40ec1ce12da00",
        "metrics": {},
        "verdict_primary": None,
        "sanity_flagged": [],
    }

    rng = np.random.default_rng(SEED)
    for metric in metrics:
        col = bg[bg[metric].notna()]
        col_metric = col[metric].to_numpy()
        rec_in = [g for g in RECEPTOR if g in col["gene"].values]
        lig_in = [g for g in LIGAND if g in col["gene"].values]
        per_gene = {g: float(col[col["gene"] == g][metric].iloc[0])
                     for g in ALL_TEST if g in col["gene"].values}

        # Compute "constraint percentile" oriented so higher = more constrained
        if direction[metric] == "higher_more_constrained":
            pct = {g: percentile_of(per_gene[g], col_metric)
                    for g in per_gene}
            bg_pcts = (col[metric].rank(pct=True) * 100).to_numpy()
        else:
            pct = {g: 100 - percentile_of(per_gene[g], col_metric)
                    for g in per_gene}
            bg_pcts = (1 - col[metric].rank(pct=True)) * 100
            bg_pcts = bg_pcts.to_numpy()

        rec_mean = float(np.mean([pct[g] for g in rec_in if g in pct]))
        lig_mean = float(np.mean([pct[g] for g in lig_in if g in pct]))

        # Permutation test: how often does a random size-len(rec_in) draw
        # produce mean constraint-percentile >= the receptor mean?
        idx = np.arange(len(col))
        rec_n = len(rec_in)
        null_means = np.empty(N_PERM)
        for i in range(N_PERM):
            sel = rng.choice(idx, size=rec_n, replace=False)
            null_means[i] = bg_pcts[sel].mean()
        p_perm_receptor = float((null_means >= rec_mean).sum() / N_PERM)

        results["metrics"][metric] = {
            "direction": direction[metric],
            "per_gene_value": per_gene,
            "per_gene_constraint_percentile": pct,
            "receptor_mean_percentile": rec_mean,
            "ligand_mean_percentile": lig_mean,
            "p_perm_receptor": p_perm_receptor,
        }
        print(f"  {metric:14s}  receptor_mean_pct={rec_mean:5.1f}  "
              f"ligand_mean_pct={lig_mean:5.1f}  perm_p={p_perm_receptor:.4f}")

    # Sanity check: synonymous Z should be near 50th percentile for all
    # 5 genes; flag any where |pct - 50| > 25
    syn_pcts = results["metrics"][SANITY]["per_gene_constraint_percentile"]
    for g, p in syn_pcts.items():
        if abs(p - 50) > 25:
            results["sanity_flagged"].append({"gene": g, "syn_percentile": p})
            print(f"  SANITY FLAG: {g} synonymous Z percentile = {p:.1f}")

    # Pre-registered verdict on primary metric (mis.z_score)
    rec_mean_primary = results["metrics"][PRIMARY_METRIC]["receptor_mean_percentile"]
    lig_mean_primary = results["metrics"][PRIMARY_METRIC]["ligand_mean_percentile"]
    p_perm = results["metrics"][PRIMARY_METRIC]["p_perm_receptor"]
    if rec_mean_primary > 75 and lig_mean_primary <= 50 and p_perm < 0.05:
        verdict = "M1-missense supported"
    elif rec_mean_primary > 50:
        verdict = "M1-missense partially supported"
    elif rec_mean_primary <= 50:
        verdict = "M1-missense falsified (Extension A reinforced)"
    else:
        verdict = "ambiguous"
    results["verdict_primary"] = verdict
    print(f"\n[5/5] Verdict on primary metric (pre-registered): {verdict}")

    out_path = RES / "missense_constraint.json"
    out_path.write_text(json.dumps(results, indent=2, default=float))
    print(f"  Wrote {out_path.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
