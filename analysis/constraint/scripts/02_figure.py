"""
Figure for the gnomAD constraint test.

Two panels:
  (left)  background histogram of LOEUF with the 5 CGRP-pathway genes
          overlaid as colored vertical lines
  (right) per-gene LOEUF percentile bar chart with thresholds annotated
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT.parents[0] / "data" / "raw" / "gnomad"
PROC = ROOT.parents[0] / "data" / "processed"
RES = ROOT / "results"
FIG = ROOT / "figures"
FIG.mkdir(parents=True, exist_ok=True)


def main() -> None:
    # Re-load constraint to get full background distribution for plotting
    df = pd.read_csv(RAW / "gnomad.v4.1.constraint_metrics.tsv", sep="\t",
                      low_memory=False)
    has_mane = df["mane_select"].fillna(False).astype(bool)
    has_can = df["canonical"].fillna(False).astype(bool)
    df_mane = df[has_mane].drop_duplicates(subset="gene", keep="first")
    mane_genes = set(df_mane["gene"])
    df_can_only = df[has_can & ~df["gene"].isin(mane_genes)].drop_duplicates(
        subset="gene", keep="first"
    )
    constraint = pd.concat([df_mane, df_can_only], ignore_index=True)
    bg = pd.read_csv(PROC / "background_features.tsv", sep="\t")
    bg = bg.merge(constraint, on="gene", how="inner")
    bg = bg[bg["lof.oe_ci.upper"].notna()].copy()
    loeuf = bg["lof.oe_ci.upper"].to_numpy()

    results = json.load((RES / "gnomad_constraint.json").open())
    per_gene = pd.DataFrame(results["per_gene"]).set_index("gene")

    GENE_ORDER = ["CALCA", "CALCB", "CALCRL", "RAMP1", "CRCP"]
    GENE_SIDES = {"CALCA": "ligand", "CALCB": "ligand",
                   "CALCRL": "receptor", "RAMP1": "receptor", "CRCP": "receptor"}
    COLORS = {
        "CALCA":  "#e377c2",
        "CALCB":  "#c0688d",
        "CALCRL": "#d62728",
        "RAMP1":  "#aa3333",
        "CRCP":   "#7f1212",
    }

    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.6),
                              gridspec_kw={"width_ratios": [1.2, 1.0]})

    # ----- LEFT: background histogram with gene overlays -----
    ax = axes[0]
    ax.hist(loeuf, bins=80, color="#cccccc", edgecolor="#555555", lw=0.4,
             label=f"background (n={len(loeuf):,})")
    ax.axvline(np.median(loeuf), color="black", ls=":", lw=1.0,
                label=f"background median = {np.median(loeuf):.2f}")
    ax.axvline(0.6, color="#cc6600", ls="--", lw=1.0,
                label="constrained tier (LOEUF<0.6)")
    for g in GENE_ORDER:
        if g not in per_gene.index:
            continue
        loeuf_g = per_gene.loc[g, "lof.oe_ci.upper"]
        ax.axvline(loeuf_g, color=COLORS[g], lw=2.2,
                    label=f"{g} ({GENE_SIDES[g]}) LOEUF={loeuf_g:.2f}")
    ax.set_xlabel("LOEUF (lower = more constrained)")
    ax.set_ylabel("count of background genes")
    ax.set_title("gnomAD v4.1 LOEUF distribution\nwith CGRP-pathway genes overlaid")
    ax.set_xlim(0, 2.5)
    ax.legend(fontsize=7, loc="upper right", framealpha=0.95)
    ax.grid(True, alpha=0.3)

    # ----- RIGHT: per-gene LOEUF percentile -----
    ax = axes[1]
    pcts = [per_gene.loc[g, "loeuf_percentile"] if g in per_gene.index else np.nan
             for g in GENE_ORDER]
    bars = ax.bar(np.arange(len(GENE_ORDER)), pcts,
                   color=[COLORS[g] for g in GENE_ORDER],
                   edgecolor="black", lw=0.6)
    for b, p in zip(bars, pcts):
        ax.text(b.get_x() + b.get_width() / 2, p + 1.5,
                 f"{p:.1f}", ha="center", va="bottom", fontsize=8,
                 fontweight="bold")
    ax.axhline(50, color="black", ls=":", lw=1.0, label="background median (50th)")
    ax.axhline(25, color="#cc6600", ls="--", lw=1.0, label="constrained tier (<25th)")
    rec_mean = results["receptor_group"]["mean_loeuf_percentile"]
    lig_mean = results["ligand_group"]["mean_loeuf_percentile"]
    ax.axhline(rec_mean, color="#d62728", ls="-.", lw=1.3, alpha=0.6,
                label=f"receptor mean = {rec_mean:.1f}")
    ax.axhline(lig_mean, color="#e377c2", ls="-.", lw=1.3, alpha=0.6,
                label=f"ligand mean = {lig_mean:.1f}")
    ax.set_xticks(np.arange(len(GENE_ORDER)))
    ax.set_xticklabels([
        f"{g}\n({GENE_SIDES[g]})" for g in GENE_ORDER], fontsize=8)
    ax.set_ylim(0, 105)
    ax.set_ylabel("LOEUF percentile (within protein-coding background)")
    ax.set_title(f"Per-gene LOEUF percentile\nverdict (pre-registered): {results['verdict']}")
    ax.legend(fontsize=7, loc="lower left", framealpha=0.95)
    ax.grid(True, axis="y", alpha=0.3)

    fig.tight_layout()
    pdf = FIG / "gnomad_constraint.pdf"
    png = FIG / "gnomad_constraint.png"
    fig.savefig(pdf, bbox_inches="tight")
    fig.savefig(png, bbox_inches="tight", dpi=160)
    print(f"Wrote {pdf.relative_to(ROOT)}")
    print(f"Wrote {png.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
