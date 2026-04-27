"""
Headline figure for the network-proximity test.

Uses the saved summary statistics from 04_network_proximity.py to plot a
Gaussian-fit null distribution (the d_c null is approximately Normal for
gene sets of this size) with the observed migraine and control-disease
values overlaid. Avoids recomputing 1000 permutations.
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import norm

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "results"
FIG = ROOT / "figures"
FIG.mkdir(parents=True, exist_ok=True)


def main() -> None:
    results = json.load((RES / "network_proximity_results.json").open())

    primary = results["primary"]
    obs = primary["observed_dc"]
    nm = primary["null_mean"]
    nsd = primary["null_sd"]
    z = primary["z_score"]
    p = primary["empirical_p_two_sided"]
    n_S = primary["n_S_in_graph"]
    n_T = primary["n_T_in_graph"]

    # Plot range covers the full null and the observation
    lo = min(obs, nm - 4 * nsd)
    hi = max(obs, nm + 4 * nsd)
    x = np.linspace(lo, hi, 400)
    pdf = norm.pdf(x, loc=nm, scale=nsd)

    fig, ax = plt.subplots(figsize=(7.2, 4.4))
    ax.fill_between(x, pdf, color="#cccccc", alpha=0.7,
                    label=fr"degree-matched null: $\mathcal{{N}}({nm:.3f}, {nsd:.3f}^2)$")
    ax.plot(x, pdf, color="#888888", lw=1)

    # Overlay null sample dots if available
    sample = primary.get("null_distances_sample") or []
    if sample:
        ys = norm.pdf(sample, loc=nm, scale=nsd) * 0.05  # offset for visibility
        ax.scatter(sample, ys, s=8, alpha=0.4, color="#444444",
                   label=f"null sample (n={len(sample)})")

    ax.axvline(obs, color="#d62728", lw=2.2,
               label=fr"migraine $d_c={obs:.3f}$ ($z={z:+.2f}$, $p={p:.3f}$)")

    # Controls
    palette = {
        "schizophrenia": "#9467bd",
        "psoriasis": "#1f77b4",
        "essential_hypertension": "#2ca02c",
        "type_2_diabetes": "#ff7f0e",
    }
    for name, c in (results.get("specificity_controls") or {}).items():
        if "error" in c:
            continue
        ax.axvline(c["observed_dc"], ls="--", lw=1.4, alpha=0.85,
                   color=palette.get(name, "#777"),
                   label=f"{name.replace('_',' ')}: $d_c={c['observed_dc']:.3f}$, $z={c['z_score']:+.2f}$")

    ax.set_xlabel(r"closest distance $d_c$ (PPI graph; lower = closer)")
    ax.set_ylabel("null density")
    ax.set_title(
        f"Network proximity: migraine GWAS genes (n={n_S}) → "
        f"CGRP-selective receptor (CALCRL+RAMP1, n={n_T})\n"
        "BioGRID ∪ IntAct physical-only; degree-matched random-set null, 1000 permutations"
    )
    ax.legend(fontsize=8, loc="upper left", framealpha=0.95)
    ax.set_yticks([])
    fig.tight_layout()

    pdf_path = FIG / "proximity_null_distribution.pdf"
    png_path = FIG / "proximity_null_distribution.png"
    fig.savefig(pdf_path, bbox_inches="tight")
    fig.savefig(png_path, bbox_inches="tight", dpi=160)
    print(f"Wrote {pdf_path.relative_to(ROOT)}")
    print(f"Wrote {png_path.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
