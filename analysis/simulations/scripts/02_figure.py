"""
Figure for the asymmetry simulation.

Three panels:
  (left)   per-gene max |Z| distribution per model (boxplot)
  (middle) joint distribution of ligand_max_Z vs receptor_max_Z, one
           panel per model, with the observed-pattern region shaded
  (right)  bar chart of hit_pattern_rate per model
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "results"
FIG = ROOT / "figures"
FIG.mkdir(parents=True, exist_ok=True)


MODEL_LABELS = {
    "M1_purifying_selection": "M1: purifying selection\n(strong on receptor)",
    "M2_allelic_series":      "M2: allelic-series\n(small receptor effects)",
    "M3_null":                "M3: null\n(symmetric)",
}
MODEL_COLORS = {
    "M1_purifying_selection": "#d62728",
    "M2_allelic_series":      "#1f77b4",
    "M3_null":                "#7f7f7f",
}
GENE_ORDER = ["CALCA", "CALCB", "CALCRL", "RAMP1", "CRCP"]
GENE_SIDES = {"CALCA": "ligand", "CALCB": "ligand",
              "CALCRL": "receptor", "RAMP1": "receptor", "CRCP": "receptor"}


def main() -> None:
    data = json.load((RES / "asymmetry_simulation.json").open())
    GWS = data["GWS_threshold"]
    INV = data["receptor_invisible_threshold"]

    fig = plt.figure(figsize=(13, 8.5))
    gs = fig.add_gridspec(2, 3, height_ratios=[1.05, 1.0],
                          width_ratios=[1.0, 1.05, 0.85],
                          hspace=0.5, wspace=0.45)

    # Top row, full width: per-gene max |Z| boxplots, one panel per model
    ax_top = fig.add_subplot(gs[0, :])
    n_models = len(MODEL_LABELS)
    width = 0.25
    x_base = np.arange(len(GENE_ORDER))
    for i, mname in enumerate(MODEL_LABELS):
        per_gene = data["by_model"][mname]["per_gene_max_abs_Z"]
        positions = x_base + (i - (n_models - 1) / 2) * width
        bp = ax_top.boxplot(
            [per_gene[g] for g in GENE_ORDER],
            positions=positions, widths=width * 0.85,
            showfliers=False, patch_artist=True,
            medianprops={"color": "black", "lw": 1.2},
            boxprops={"facecolor": MODEL_COLORS[mname], "alpha": 0.7,
                       "edgecolor": "black", "lw": 0.6},
            whiskerprops={"lw": 0.7}, capprops={"lw": 0.7},
        )
    ax_top.axhline(GWS, color="#222", ls="--", lw=1, label=f"GWS |Z|={GWS:.2f}")
    ax_top.axhline(INV, color="#888", ls=":", lw=1, label=f"receptor invisible |Z|={INV:.2f}")
    ax_top.set_xticks(x_base)
    ax_top.set_xticklabels([
        f"{g}\n({GENE_SIDES[g]})" for g in GENE_ORDER])
    ax_top.set_ylabel("per-gene max |Z|  (per replicate)")
    ax_top.set_title("Simulated per-gene max |Z| across $10{,}000$ replicates per model\n"
                      "Box = IQR; whisker = 1.5*IQR; outliers hidden")
    legend_patches = [mpatches.Patch(color=MODEL_COLORS[m], alpha=0.7,
                                       label=MODEL_LABELS[m].replace("\n", " "))
                       for m in MODEL_LABELS]
    legend_patches.append(plt.Line2D([0], [0], color="#222", ls="--",
                                      label=f"GWS |Z|={GWS:.2f}"))
    legend_patches.append(plt.Line2D([0], [0], color="#888", ls=":",
                                      label=f"receptor invisible |Z|={INV:.2f}"))
    ax_top.legend(handles=legend_patches, fontsize=8, loc="upper right",
                   ncol=2, framealpha=0.95)
    ax_top.set_ylim(bottom=0)
    ax_top.grid(True, axis="y", alpha=0.3)

    # Bottom-left and bottom-middle: joint scatter of ligand vs receptor
    # max |Z|, one panel per model
    z_max = max(
        max(max(data["by_model"][m]["ligand_max_Z"]),
             max(data["by_model"][m]["receptor_max_Z"]))
        for m in MODEL_LABELS
    )
    z_max = min(z_max, 30)  # clip for readability
    for i, mname in enumerate(["M1_purifying_selection", "M2_allelic_series", "M3_null"]):
        ax = fig.add_subplot(gs[1, i] if i < 2 else gs[1, 2])
        if i >= 2:
            continue  # last col is bar chart
        lig = np.array(data["by_model"][mname]["ligand_max_Z"])
        rec = np.array(data["by_model"][mname]["receptor_max_Z"])
        # Shade the observed-pattern region: ligand > GWS AND receptor < INV
        ax.axvspan(GWS, z_max + 1, color="#ffe599", alpha=0.45, zorder=0)
        ax.axhspan(0, INV, color="#ffe599", alpha=0.45, zorder=0)
        # The OR-region intersection is the actual hit-pattern square
        ax.add_patch(plt.Rectangle((GWS, 0), z_max + 1 - GWS, INV,
                                     facecolor="#f4cccc", alpha=0.7, zorder=1,
                                     edgecolor="#cc0000", lw=1.0))
        ax.scatter(lig, rec, s=4, alpha=0.18, color=MODEL_COLORS[mname],
                    zorder=2, rasterized=True)
        ax.axvline(GWS, color="#222", ls="--", lw=0.8, zorder=3)
        ax.axhline(INV, color="#888", ls=":", lw=0.8, zorder=3)
        ax.set_xlim(0, z_max + 1)
        ax.set_ylim(0, z_max + 1)
        ax.set_xlabel(r"ligand max $|Z|$")
        if i == 0:
            ax.set_ylabel(r"receptor max $|Z|$")
        rate = data["by_model"][mname]["summary"]["hit_pattern_rate"]
        ax.set_title(f"{MODEL_LABELS[mname]}\nhit pattern: "
                      f"{100*rate:.1f}% of replicates",
                      fontsize=9)
        ax.set_aspect("equal", adjustable="box")
        ax.grid(True, alpha=0.3)

    # Bottom-right: bar chart of hit_pattern_rate per model
    ax_bar = fig.add_subplot(gs[1, 2])
    rates = [data["by_model"][m]["summary"]["hit_pattern_rate"] for m in MODEL_LABELS]
    lig_rates = [data["by_model"][m]["summary"]["ligand_GWS_rate"] for m in MODEL_LABELS]
    rec_rates = [data["by_model"][m]["summary"]["receptor_invisible_rate"] for m in MODEL_LABELS]
    xpos = np.arange(len(MODEL_LABELS))
    barwidth = 0.27
    ax_bar.bar(xpos - barwidth, lig_rates, barwidth, label="ligand reaches GWS",
                color="#666666", edgecolor="black", lw=0.6)
    ax_bar.bar(xpos, rec_rates, barwidth, label="receptor stays invisible",
                color="#bbbbbb", edgecolor="black", lw=0.6)
    bars3 = ax_bar.bar(xpos + barwidth, rates, barwidth, label="both (hit pattern)",
                        color=[MODEL_COLORS[m] for m in MODEL_LABELS],
                        edgecolor="black", lw=0.6)
    for b, r in zip(bars3, rates):
        ax_bar.text(b.get_x() + b.get_width() / 2, r + 0.01,
                     f"{100*r:.1f}%", ha="center", va="bottom", fontsize=8,
                     fontweight="bold")
    ax_bar.set_xticks(xpos)
    ax_bar.set_xticklabels([m.split("_", 1)[0] for m in MODEL_LABELS],
                            fontsize=9)
    ax_bar.set_ylabel("fraction of replicates")
    ax_bar.set_title("Decision rule:\nfraction matching observed pattern",
                      fontsize=9)
    ax_bar.legend(fontsize=8, loc="upper right", framealpha=0.95)
    ax_bar.set_ylim(0, max(max(lig_rates), 1.0) * 1.05)
    ax_bar.grid(True, axis="y", alpha=0.3)

    pdf = FIG / "asymmetry_simulation.pdf"
    png = FIG / "asymmetry_simulation.png"
    fig.savefig(pdf, bbox_inches="tight")
    fig.savefig(png, bbox_inches="tight", dpi=160)
    print(f"Wrote {pdf.relative_to(ROOT)}")
    print(f"Wrote {png.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
