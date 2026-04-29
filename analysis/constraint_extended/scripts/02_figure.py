"""
Figure for Extension E: per-gene constraint percentile across 4 metrics
(missense Z, missense O/E, missense-damaging O/E, synonymous Z control).
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

GENE_ORDER = ["CALCA", "CALCB", "CALCRL", "RAMP1", "CRCP"]
GENE_SIDES = {"CALCA": "ligand", "CALCB": "ligand",
               "CALCRL": "receptor", "RAMP1": "receptor",
               "CRCP": "receptor"}
COLORS = {
    "CALCA":  "#e377c2",
    "CALCB":  "#c0688d",
    "CALCRL": "#d62728",
    "RAMP1":  "#aa3333",
    "CRCP":   "#7f1212",
}

METRIC_LABELS = {
    "mis.z_score":   "missense Z\n(higher = more constrained)",
    "mis.oe":         "missense O/E\n(rescaled to constraint percentile)",
    "mis_pphen.oe":   "PolyPhen-damaging\nmissense O/E (rescaled)",
    "syn.z_score":    "synonymous Z\n(sanity control; expect ~50)",
}


def main() -> None:
    results = json.load((RES / "missense_constraint.json").open())
    metrics = list(METRIC_LABELS.keys())

    fig, axes = plt.subplots(1, len(metrics), figsize=(15, 4.6),
                              sharey=True)

    for idx, metric in enumerate(metrics):
        ax = axes[idx]
        m = results["metrics"][metric]
        pct = m["per_gene_constraint_percentile"]
        rec_mean = m["receptor_mean_percentile"]
        lig_mean = m["ligand_mean_percentile"]
        p_perm = m["p_perm_receptor"]

        ypos = np.arange(len(GENE_ORDER))
        vals = [pct.get(g, np.nan) for g in GENE_ORDER]
        bars = ax.barh(ypos, vals,
                        color=[COLORS[g] for g in GENE_ORDER],
                        edgecolor="black", lw=0.6)
        for b, v in zip(bars, vals):
            if np.isfinite(v):
                ax.text(v + 1, b.get_y() + b.get_height() / 2,
                         f"{v:.0f}", va="center", fontsize=8)
        ax.axvline(50, color="black", ls=":", lw=0.8,
                    label="bg median (50)")
        ax.axvline(75, color="#cc6600", ls="--", lw=0.8,
                    label="constrained (>75)")
        ax.axvline(rec_mean, color="#d62728", ls="-.", lw=1.2, alpha=0.6,
                    label=f"recep mean = {rec_mean:.1f}")
        ax.axvline(lig_mean, color="#e377c2", ls="-.", lw=1.2, alpha=0.6,
                    label=f"lig mean = {lig_mean:.1f}")
        ax.set_xlim(0, 105)
        ax.set_xlabel("constraint percentile")
        ax.set_title(f"{METRIC_LABELS[metric]}\nperm p = {p_perm:.3f}",
                      fontsize=9)
        if idx == 0:
            ax.set_yticks(ypos)
            ax.set_yticklabels([f"{g}\n({GENE_SIDES[g]})" for g in GENE_ORDER],
                                 fontsize=9)
            ax.invert_yaxis()
            ax.legend(fontsize=7, loc="lower right", framealpha=0.95)
        ax.grid(True, axis="x", alpha=0.3)

    fig.suptitle(
        "Constraint metrics beyond LoF tolerance (gnomAD v4.1)\n"
        f"verdict on primary metric (pre-registered): {results['verdict_primary']}",
        fontsize=10, y=1.00)
    fig.tight_layout()
    pdf = FIG / "missense_constraint.pdf"
    png = FIG / "missense_constraint.png"
    fig.savefig(pdf, bbox_inches="tight")
    fig.savefig(png, bbox_inches="tight", dpi=160)
    print(f"Wrote {pdf.relative_to(ROOT)}")
    print(f"Wrote {png.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
