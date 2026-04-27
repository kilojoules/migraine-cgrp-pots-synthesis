"""
GTEx v8 tissue-expression heatmap for migraine genes + CGRP anchor across
autonomic-relevant tissues, per pre-registration.

Tissues: Brain - Spinal cord (cervical c-1), Brain - Hypothalamus,
Nerve - Tibial, Artery - Tibial, Artery - Aorta, Artery - Coronary,
Whole Blood (control). (Trigeminal/sympathetic ganglia are not in v8.)
"""
from __future__ import annotations

import gzip
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import requests
import seaborn as sns

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw"
PROC = ROOT / "data" / "processed"
RES = ROOT / "results"
FIG = ROOT / "figures"
RES.mkdir(parents=True, exist_ok=True)
FIG.mkdir(parents=True, exist_ok=True)

GTEX_URL = (
    "https://storage.googleapis.com/adult-gtex/bulk-gex/v8/rna-seq/"
    "GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz"
)
GTEX_FILE = RAW / "gtex_v8_median_tpm.gct.gz"

PRIMARY_TISSUES = [
    "Brain - Spinal cord (cervical c-1)",
    "Brain - Hypothalamus",
    "Nerve - Tibial",
    "Artery - Tibial",
    "Artery - Aorta",
    "Artery - Coronary",
    "Whole Blood",
]


def download_gtex() -> None:
    if GTEX_FILE.exists() and GTEX_FILE.stat().st_size > 1_000_000:
        print(f"GTEx already at {GTEX_FILE.relative_to(ROOT)} ({GTEX_FILE.stat().st_size:,} bytes)")
        return
    print(f"Downloading GTEx v8 median TPM matrix...")
    r = requests.get(GTEX_URL, stream=True, timeout=600)
    r.raise_for_status()
    with GTEX_FILE.open("wb") as fh:
        for chunk in r.iter_content(chunk_size=1 << 20):
            fh.write(chunk)
    print(f"  -> {GTEX_FILE.relative_to(ROOT)} ({GTEX_FILE.stat().st_size:,} bytes)")


def load_gtex_median_tpm() -> pd.DataFrame:
    """Load GCT file: skip 2 metadata lines, then header + matrix.

    GCT format: line 1 = '#1.2', line 2 = 'n_rows\\tn_cols', line 3 = header.
    """
    with gzip.open(GTEX_FILE, "rt") as fh:
        next(fh)  # #1.2
        next(fh)  # dims
        df = pd.read_csv(fh, sep="\t")
    # First two columns are 'Name' (Ensembl ID) and 'Description' (gene symbol)
    df = df.rename(columns={"Name": "ensembl_id", "Description": "gene"})
    return df


def main() -> None:
    download_gtex()
    print("Loading GTEx matrix...")
    df = load_gtex_median_tpm()
    print(f"  {len(df):,} genes x {df.shape[1] - 2} tissues")

    missing_tissues = [t for t in PRIMARY_TISSUES if t not in df.columns]
    if missing_tissues:
        print(f"  WARNING: tissues not in matrix: {missing_tissues}")
    avail = [t for t in PRIMARY_TISSUES if t in df.columns]

    # Migraine primary set
    migraine = pd.read_csv(PROC / "migraine_set_primary.tsv", sep="\t")["gene"].tolist()
    cgrp = pd.read_csv(PROC / "cgrp_anchor_strict.tsv", sep="\t")["gene"].tolist()
    cgrp_ext = pd.read_csv(PROC / "cgrp_anchor_extended.tsv", sep="\t")["gene"].tolist()

    target_genes = sorted(set(migraine) | set(cgrp_ext))
    sub = df[df["gene"].isin(target_genes)][["gene"] + avail].copy()
    print(f"  Found expression for {len(sub)}/{len(target_genes)} genes")

    # Average duplicate gene rows (some Ensembl IDs map to same symbol)
    sub = sub.groupby("gene", as_index=False)[avail].mean()

    # Save full matrix
    sub.to_csv(RES / "gtex_median_tpm_subset.tsv", sep="\t", index=False)

    # Log-transform for plotting
    plot = sub.set_index("gene")[avail]
    plot_log = np.log2(plot + 1.0)

    # Mark CGRP genes
    is_cgrp = plot.index.isin(cgrp_ext)
    row_colors = ["#d62728" if c else "#1f77b4" for c in is_cgrp]

    # Heatmap (top 50 by max expression for legibility, plus CGRP genes always shown)
    max_expr = plot_log.max(axis=1)
    top_genes = max_expr.nlargest(50).index.tolist()
    show_genes = sorted(set(top_genes) | set(plot_log.index[is_cgrp]))
    plot_show = plot_log.loc[show_genes]
    is_cgrp_show = plot_show.index.isin(cgrp_ext)
    row_colors_show = ["#d62728" if c else "#1f77b4" for c in is_cgrp_show]

    fig_h = max(8, 0.18 * len(plot_show))
    g = sns.clustermap(
        plot_show,
        cmap="viridis",
        figsize=(8, fig_h),
        row_colors=row_colors_show,
        col_cluster=False,
        linewidths=0.05,
        cbar_kws={"label": "log2(median TPM + 1)"},
    )
    g.ax_heatmap.set_xlabel("")
    g.ax_heatmap.set_ylabel("Gene")
    fig_path = FIG / "gtex_tissue_heatmap.pdf"
    g.savefig(fig_path, bbox_inches="tight")
    plt.close(g.fig)
    fig_png = FIG / "gtex_tissue_heatmap.png"
    g_png = sns.clustermap(
        plot_show,
        cmap="viridis",
        figsize=(8, fig_h),
        row_colors=row_colors_show,
        col_cluster=False,
        linewidths=0.05,
        cbar_kws={"label": "log2(median TPM + 1)"},
    )
    g_png.savefig(fig_png, bbox_inches="tight", dpi=150)
    plt.close(g_png.fig)
    print(f"Heatmap -> {fig_path.relative_to(ROOT)}")

    # Tissue-mean comparison: migraine vs CGRP vs all-genes baseline
    summary = {
        "tissues_used": avail,
        "tissues_missing": missing_tissues,
        "n_migraine_with_expr": int(plot.index.isin(migraine).sum()),
        "n_cgrp_with_expr": int(plot.index.isin(cgrp_ext).sum()),
        "median_log2tpm_per_tissue": {},
    }
    for t in avail:
        mig_med = float(plot_log.loc[plot_log.index.isin(migraine), t].median())
        cgrp_med = float(plot_log.loc[plot_log.index.isin(cgrp_ext), t].median())
        all_med = float(np.log2(df[t].astype(float) + 1.0).median())
        summary["median_log2tpm_per_tissue"][t] = {
            "migraine_set_median": mig_med,
            "cgrp_set_median": cgrp_med,
            "all_genes_median": all_med,
        }

    with (RES / "gtex_tissue_summary.json").open("w") as fh:
        json.dump(summary, fh, indent=2)
    print(f"Summary -> {(RES / 'gtex_tissue_summary.json').relative_to(ROOT)}")


if __name__ == "__main__":
    main()
