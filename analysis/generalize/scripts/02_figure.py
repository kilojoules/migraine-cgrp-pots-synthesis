"""
Figure for Extension C: cross-pathway ligand vs receptor L2G scatter with
quadrant classifications.
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

CLASS_COLOR = {
    "ligand-asymmetric":   "#d62728",
    "receptor-asymmetric": "#1f77b4",
    "balanced":            "#2ca02c",
    "both-quiet":          "#999999",
}


def main() -> None:
    summary = json.load((RES / "asymmetry_summary.json").open())
    rows = summary["per_pathway_summary"]

    fig, axes = plt.subplots(1, 2, figsize=(13, 5.6),
                              gridspec_kw={"width_ratios": [1.15, 1.0]})

    # ----- LEFT: scatter ligand vs receptor median L2G -----
    ax = axes[0]
    # Quadrant background fills
    ax.fill_between([0.2, 1.05], 0, 0.1, color="#fff0f0", alpha=0.5,
                     zorder=0, label="ligand-asymmetric region")
    ax.fill_betweenx([0.2, 1.05], 0, 0.1, color="#f0f0ff", alpha=0.5,
                      zorder=0, label="receptor-asymmetric region")
    ax.fill_between([0.1, 1.05], 0.1, 1.05, color="#f0fff0", alpha=0.35,
                     zorder=0, label="balanced region (both > 0.1)")
    # Diagonal
    ax.plot([0, 1.05], [0, 1.05], color="#888", ls="--", lw=0.8, zorder=1)

    for pname, row in rows.items():
        x = row["ligand_l2g_median"]
        y = row["receptor_l2g_median"]
        cls = row["class"]
        col = CLASS_COLOR[cls]
        ax.scatter([x], [y], s=140, color=col, edgecolor="black", lw=0.8,
                    zorder=4)
        ax.annotate(f"{pname}\n({row['disease_label'][:25]})",
                     xy=(x, y), xytext=(7, 7), textcoords="offset points",
                     fontsize=8, zorder=5)

    ax.set_xlim(-0.04, 1.05)
    ax.set_ylim(-0.04, 1.05)
    ax.set_xlabel("median L2G genetic-association (ligand genes)")
    ax.set_ylabel("median L2G genetic-association (receptor genes)")
    ax.set_title("Ligand vs receptor genetic signal across 10 peptide-hormone pathways")
    handles_classes = [plt.scatter([], [], s=80, color=c, edgecolor="black",
                                     label=k)
                        for k, c in CLASS_COLOR.items()]
    ax.legend(handles=handles_classes, fontsize=8, loc="lower right",
               framealpha=0.95)
    ax.grid(True, alpha=0.3)

    # ----- RIGHT: per-pathway asymmetry bar chart -----
    ax = axes[1]
    pathways = list(rows.keys())
    asyms = [rows[p]["asymmetry"] for p in pathways]
    cols = [CLASS_COLOR[rows[p]["class"]] for p in pathways]
    # Sort by asymmetry score descending
    order = np.argsort(asyms)[::-1]
    pathways = [pathways[i] for i in order]
    asyms = [asyms[i] for i in order]
    cols = [cols[i] for i in order]

    ypos = np.arange(len(pathways))
    bars = ax.barh(ypos, asyms, color=cols, edgecolor="black", lw=0.6)
    ax.axvline(0, color="black", lw=0.8)
    ax.axvline(0.20, color="#888", ls="--", lw=0.8,
                label="ligand-asym threshold +0.2")
    ax.axvline(-0.20, color="#888", ls="--", lw=0.8,
                label="receptor-asym threshold -0.2")
    for b, a in zip(bars, asyms):
        ax.text(a + (0.02 if a >= 0 else -0.02), b.get_y() + b.get_height() / 2,
                 f"{a:+.2f}", ha="left" if a >= 0 else "right",
                 va="center", fontsize=8)
    ax.set_yticks(ypos)
    ax.set_yticklabels(pathways, fontsize=9)
    ax.invert_yaxis()
    ax.set_xlabel("asymmetry $A$ = median(ligand L2G) - median(receptor L2G)")
    verdict = summary["verdict"]
    n_lig = summary["n_ligand_asymmetric"]
    n_eval = summary["n_evaluable"]
    p_perm = summary["permutation_p"]
    ax.set_title(f"Per-pathway asymmetry score\n"
                  f"verdict (pre-registered): {verdict}  "
                  f"({n_lig} ligand-asym / {n_eval} evaluable, perm $p = {p_perm:.2f}$)")
    ax.legend(fontsize=7, loc="lower right", framealpha=0.95)
    ax.grid(True, axis="x", alpha=0.3)
    ax.set_xlim(-1.0, 1.0)

    fig.tight_layout()
    pdf = FIG / "pathway_asymmetry.pdf"
    png = FIG / "pathway_asymmetry.png"
    fig.savefig(pdf, bbox_inches="tight")
    fig.savefig(png, bbox_inches="tight", dpi=160)
    print(f"Wrote {pdf.relative_to(ROOT)}")
    print(f"Wrote {png.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
