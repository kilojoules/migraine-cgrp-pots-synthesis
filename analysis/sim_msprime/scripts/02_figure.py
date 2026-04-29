"""
Figure for Extension D: side-by-side PRF vs msprime hit-pattern rates.
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
PRF_RES = ROOT.parents[0] / "simulations" / "results" / "asymmetry_simulation.json"
RES = ROOT / "results"
FIG = ROOT / "figures"
FIG.mkdir(parents=True, exist_ok=True)


def main() -> None:
    msp = json.load((RES / "msprime_simulation.json").open())
    prf = json.load(PRF_RES.open())

    models = ["M1_purifying_selection", "M2_allelic_series", "M3_null"]
    labels = ["M1: purifying\nselection", "M2: allelic\nseries", "M3: null"]
    colors = ["#d62728", "#1f77b4", "#7f7f7f"]

    prf_rates = [prf["by_model"][m]["summary"]["hit_pattern_rate"] for m in models]
    msp_rates = [msp["by_model"][m]["hit_pattern_rate"] for m in models]
    prf_lig  = [prf["by_model"][m]["summary"]["ligand_GWS_rate"] for m in models]
    msp_lig  = [msp["by_model"][m]["ligand_GWS_rate"] for m in models]
    prf_rec  = [prf["by_model"][m]["summary"]["receptor_invisible_rate"] for m in models]
    msp_rec  = [msp["by_model"][m]["receptor_invisible_rate"] for m in models]

    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.6))

    # LEFT: hit_pattern_rate, PRF vs msprime, grouped by model
    ax = axes[0]
    xpos = np.arange(len(models))
    w = 0.36
    bars1 = ax.bar(xpos - w/2, prf_rates, w, label="PRF (Wright equilibrium)",
                    color=colors, edgecolor="black", lw=0.6)
    bars2 = ax.bar(xpos + w/2, msp_rates, w, label="msprime + post-hoc selection",
                    color=colors, edgecolor="black", lw=0.6, alpha=0.55,
                    hatch="//")
    for b, r in zip(bars1, prf_rates):
        ax.text(b.get_x() + b.get_width() / 2, r + 0.005,
                 f"{100*r:.1f}%", ha="center", va="bottom", fontsize=8)
    for b, r in zip(bars2, msp_rates):
        ax.text(b.get_x() + b.get_width() / 2, r + 0.005,
                 f"{100*r:.1f}%", ha="center", va="bottom", fontsize=8)
    ax.set_xticks(xpos)
    ax.set_xticklabels(labels)
    ax.set_ylabel("hit-pattern rate")
    ax.set_title("Hit-pattern rate per model: PRF vs msprime\n"
                  "ranking M1 > M2 > M3 preserved across simulation frameworks")
    ax.legend(fontsize=8, loc="upper right", framealpha=0.95)
    ax.grid(True, axis="y", alpha=0.3)
    ax.set_ylim(0, max(prf_rates) * 1.15)

    # RIGHT: scatter of components (ligand_GWS_rate, receptor_invisible_rate)
    ax = axes[1]
    for i, m in enumerate(models):
        ax.plot([prf_lig[i], msp_lig[i]],
                 [prf_rec[i], msp_rec[i]],
                 color=colors[i], lw=1.2, alpha=0.55, zorder=2)
        ax.scatter([prf_lig[i]], [prf_rec[i]], s=130, color=colors[i],
                    edgecolor="black", lw=0.8, marker="o", zorder=3,
                    label=f"{labels[i].splitlines()[0]} PRF" if i == 0 else None)
        ax.scatter([msp_lig[i]], [msp_rec[i]], s=130, color=colors[i],
                    edgecolor="black", lw=0.8, marker="s", zorder=3)
        # Annotate
        ax.annotate(labels[i].split(":")[0], (prf_lig[i], prf_rec[i]),
                     xytext=(8, 6), textcoords="offset points", fontsize=8,
                     color=colors[i], fontweight="bold")
    # Manual legend for marker shape
    ax.scatter([], [], s=80, color="white", edgecolor="black",
                marker="o", label="PRF")
    ax.scatter([], [], s=80, color="white", edgecolor="black",
                marker="s", label="msprime")
    ax.set_xlabel("ligand reaches GWS rate")
    ax.set_ylabel("receptor stays invisible rate")
    ax.set_xlim(0, 1.0)
    ax.set_ylim(0, 0.5)
    ax.set_title("Component decomposition\n"
                  "circles=PRF, squares=msprime; line connects same model")
    ax.legend(fontsize=8, loc="upper right", framealpha=0.95)
    ax.grid(True, alpha=0.3)

    fig.tight_layout()
    pdf = FIG / "msprime_vs_prf.pdf"
    png = FIG / "msprime_vs_prf.png"
    fig.savefig(pdf, bbox_inches="tight")
    fig.savefig(png, bbox_inches="tight", dpi=160)
    print(f"Wrote {pdf.relative_to(ROOT)}")
    print(f"Wrote {png.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
