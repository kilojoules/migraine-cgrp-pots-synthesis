"""
Extension A: gnomAD v4.1 constraint test for the CGRP ligand-receptor asymmetry.

Probes scenario M1 (purifying selection on receptor genes) directly via
LOEUF, the standard gnomAD constraint metric. Decision rule and seed are
locked in analysis/constraint/PRE_REGISTRATION_constraint.md.

Inputs:
  analysis/data/raw/gnomad/gnomad.v4.1.constraint_metrics.tsv
  analysis/data/processed/background_features.tsv

Outputs:
  analysis/constraint/results/gnomad_constraint.json
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

SEED = 20260429
N_PERM = 10_000


def load_constraint() -> pd.DataFrame:
    """Load gnomAD v4.1 constraint, restrict to canonical or MANE-select
    transcript per gene, return per-gene rows."""
    print(f"[1/5] Loading {CONSTRAINT_TSV.name}...")
    df = pd.read_csv(CONSTRAINT_TSV, sep="\t", low_memory=False)
    print(f"      {len(df):,} rows total")
    # Prefer MANE-select; fall back to canonical
    has_mane = df["mane_select"].fillna(False).astype(bool)
    has_can = df["canonical"].fillna(False).astype(bool)
    df_mane = df[has_mane].drop_duplicates(subset="gene", keep="first")
    mane_genes = set(df_mane["gene"])
    df_can_only = df[has_can & ~df["gene"].isin(mane_genes)].drop_duplicates(
        subset="gene", keep="first"
    )
    out = pd.concat([df_mane, df_can_only], ignore_index=True)
    print(f"      {len(out):,} unique genes (MANE-select preferred)")
    return out


def main() -> None:
    df = load_constraint()
    bg = pd.read_csv(BG_TSV, sep="\t")
    print(f"[2/5] Background: {len(bg):,} protein-coding genes")

    # Map background genes to constraint rows
    bg_with_constraint = bg.merge(df, on="gene", how="inner")
    print(f"      {len(bg_with_constraint):,} background genes have gnomAD constraint")

    # Restrict to genes with non-null LOEUF
    has_loeuf = bg_with_constraint["lof.oe_ci.upper"].notna()
    bg_used = bg_with_constraint[has_loeuf].copy()
    print(f"      {len(bg_used):,} background genes have non-null LOEUF")

    # Compute per-gene LOEUF percentile within the background distribution
    loeuf_arr = bg_used["lof.oe_ci.upper"].to_numpy()
    bg_used["loeuf_percentile"] = (
        bg_used["lof.oe_ci.upper"].rank(pct=True) * 100
    )

    print("[3/5] Per-gene metrics for the 5 CGRP-pathway genes:")
    test_rows = bg_used[bg_used["gene"].isin(ALL_TEST)].copy()
    cols = ["gene", "lof.oe_ci.upper", "loeuf_percentile",
             "lof.pLI", "lof.oe", "mis.z_score"]
    cols_present = [c for c in cols if c in test_rows.columns]
    test_rows = test_rows[cols_present].set_index("gene")
    # Reorder to ligand then receptor
    test_rows = test_rows.reindex([g for g in ALL_TEST if g in test_rows.index])
    print(test_rows.to_string())

    # Group means
    rec_in = [g for g in RECEPTOR if g in test_rows.index]
    lig_in = [g for g in LIGAND if g in test_rows.index]
    rec_mean_pct = float(test_rows.loc[rec_in, "loeuf_percentile"].mean()) if rec_in else None
    lig_mean_pct = float(test_rows.loc[lig_in, "loeuf_percentile"].mean()) if lig_in else None
    print(f"\n  Receptor group mean LOEUF percentile (n={len(rec_in)}): {rec_mean_pct:.1f}")
    print(f"  Ligand group mean LOEUF percentile (n={len(lig_in)}):   {lig_mean_pct:.1f}")
    print(f"  Background median LOEUF: {np.median(loeuf_arr):.3f}")

    print(f"\n[4/5] Permutation tests, n={N_PERM:,} draws, seed={SEED}")
    rng = np.random.default_rng(SEED)
    bg_pcts = bg_used["loeuf_percentile"].to_numpy()
    bg_pool_idx = np.arange(len(bg_used))

    # Receptor: random size-3 draws, mean LOEUF percentile
    rec_n = len(rec_in)
    rec_null = np.empty(N_PERM)
    for i in range(N_PERM):
        idx = rng.choice(bg_pool_idx, size=rec_n, replace=False)
        rec_null[i] = bg_pcts[idx].mean()
    p_perm_receptor = float((rec_null <= rec_mean_pct).sum() / N_PERM)
    print(f"  P(random size-{rec_n} mean LOEUF percentile <= {rec_mean_pct:.1f}) = {p_perm_receptor:.4f}")

    # Ligand: random size-2 draws, mean LOEUF percentile (test for >= obs)
    lig_n = len(lig_in)
    lig_null = np.empty(N_PERM)
    for i in range(N_PERM):
        idx = rng.choice(bg_pool_idx, size=lig_n, replace=False)
        lig_null[i] = bg_pcts[idx].mean()
    p_perm_ligand = float((lig_null >= lig_mean_pct).sum() / N_PERM)
    print(f"  P(random size-{lig_n} mean LOEUF percentile >= {lig_mean_pct:.1f}) = {p_perm_ligand:.4f}")

    # Pre-registered decision rule
    if rec_mean_pct < 25.0 and lig_mean_pct > 50.0:
        verdict = "M1 supported"
    elif rec_mean_pct < 50.0:
        verdict = "M1 partially supported"
    elif rec_mean_pct >= 50.0:
        verdict = "M1 falsified (M2 supported)"
    else:
        verdict = "ambiguous"
    print(f"\n  Verdict (pre-registered decision rule): {verdict}")

    print("\n[5/5] Writing JSON output...")
    out: dict = {
        "spec": "analysis/constraint/PRE_REGISTRATION_constraint.md",
        "seed": SEED,
        "n_permutations": N_PERM,
        "constraint_tsv_sha256": "68d8abdb7fc48f570869b02dfaa74b9fecaece7fcc5f301ddca40ec1ce12da00",
        "n_background_genes_with_loeuf": len(bg_used),
        "background_loeuf_quantiles": {
            "p10": float(np.quantile(loeuf_arr, 0.10)),
            "p25": float(np.quantile(loeuf_arr, 0.25)),
            "p50": float(np.median(loeuf_arr)),
            "p75": float(np.quantile(loeuf_arr, 0.75)),
            "p90": float(np.quantile(loeuf_arr, 0.90)),
        },
        "per_gene": test_rows.reset_index().to_dict(orient="records"),
        "receptor_group": {
            "genes": rec_in,
            "mean_loeuf_percentile": rec_mean_pct,
            "p_perm": p_perm_receptor,
        },
        "ligand_group": {
            "genes": lig_in,
            "mean_loeuf_percentile": lig_mean_pct,
            "p_perm": p_perm_ligand,
        },
        "verdict": verdict,
    }
    out_path = RES / "gnomad_constraint.json"
    out_path.write_text(json.dumps(out, indent=2, default=float))
    print(f"  Wrote {out_path.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
