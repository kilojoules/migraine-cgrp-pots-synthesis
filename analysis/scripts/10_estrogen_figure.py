"""
Figure for the ESR1 ChIP-seq promoter-occupancy test.

Two-panel figure:
  (left) bar chart of fraction-bound by gene set at the headline threshold
  (right) sensitivity sweep showing fold-enrichment vs stringency
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "results"
FIG = ROOT / "figures"
FIG.mkdir(parents=True, exist_ok=True)


def main() -> None:
    data = json.load((RES / "estrogen_chipseq_results.json").open())
    headline = str(data["data"]["headline_threshold"])
    headline_data = data["by_threshold"][headline]
    thresholds = data["data"]["peak_score_threshold_sweep"]

    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.5))

    set_labels = {
        "cgrp_strict": "CGRP strict\n(CALCRL+RAMP1)",
        "cgrp_extended": "CGRP extended\n(CALCA/B + R/RAMP1 + CRCP)",
        "migraine_primary": "Migraine GWAS\n(123 lead-SNP-nearest)",
        "autonomic": "Autonomic\n(286 GO/KEGG-curated)",
    }
    colors = {
        "cgrp_strict": "#d62728",
        "cgrp_extended": "#e377c2",
        "migraine_primary": "#1f77b4",
        "autonomic": "#2ca02c",
    }
    bg_color = "#999999"

    ax = axes[0]
    keys = ["cgrp_strict", "cgrp_extended", "migraine_primary", "autonomic"]
    bg_frac = headline_data["background_esr1_bound_fraction"]
    fracs = [headline_data[k]["fraction_bound"] for k in keys]
    pvals = [headline_data[k]["length_matched_permutation_p_one_sided_greater"] for k in keys]
    ns = [headline_data[k]["n_in_background"] for k in keys]
    bound = [headline_data[k]["n_with_esr1_peak"] for k in keys]

    xpos = np.arange(len(keys))
    bars = ax.bar(xpos, fracs, color=[colors[k] for k in keys], edgecolor="black", lw=0.6)
    ax.axhline(bg_frac, color=bg_color, ls="--", lw=1.5,
               label=f"background ({100*bg_frac:.0f}%)")
    for i, (b, n_b, n, p) in enumerate(zip(bars, bound, ns, pvals)):
        h = b.get_height()
        ax.text(b.get_x() + b.get_width() / 2, h + 0.015,
                f"{n_b}/{n}\np={p:.2f}",
                ha="center", va="bottom", fontsize=8)
    ax.set_xticks(xpos)
    ax.set_xticklabels([set_labels[k] for k in keys], fontsize=8)
    ax.set_ylim(0, 1.0)
    ax.set_ylabel(f"fraction with ESR1 peak in promoter\n(score $\\geq$ {headline})")
    ax.set_title(f"ESR1 promoter occupancy at headline threshold (score $\\geq$ {headline})\n"
                 f"none of the four gene sets exceeds length-matched expectation")
    ax.legend(loc="lower right", fontsize=8)

    ax = axes[1]
    for k in keys:
        folds = [data["by_threshold"][str(t)][k]["fold_enrichment_vs_background"] for t in thresholds]
        ax.plot(thresholds, folds, "-o", color=colors[k], label=set_labels[k].replace("\n", " "))
    ax.axhline(1.0, color=bg_color, ls="--", lw=1.0, label="no enrichment")
    ax.set_xscale("log")
    ax.set_xticks(thresholds)
    ax.set_xticklabels([str(t) for t in thresholds])
    ax.set_xlabel("ReMap NR peak score threshold (number of supporting datasets)")
    ax.set_ylabel("fold-enrichment of ESR1 promoter occupancy\nvs length-matched background")
    ax.set_title("Sensitivity sweep over peak-score stringency\n"
                 "fold-enrichment stays at or below 1.0 across the range")
    ax.legend(fontsize=7, loc="lower left")
    ax.grid(True, alpha=0.3)

    fig.tight_layout()
    pdf = FIG / "estrogen_chipseq_enrichment.pdf"
    png = FIG / "estrogen_chipseq_enrichment.png"
    fig.savefig(pdf, bbox_inches="tight")
    fig.savefig(png, bbox_inches="tight", dpi=160)
    print(f"Wrote {pdf.relative_to(ROOT)}")
    print(f"Wrote {png.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
