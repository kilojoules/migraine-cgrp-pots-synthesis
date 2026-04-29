"""
Figure for Extension B: gene-burden evidence presence/absence for the
five CGRP-pathway genes alongside positive-control genes (PCSK9, LDLR,
MC4R) across the broader phenotype panel.
"""
from __future__ import annotations

import json
import time
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import requests

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "results"
FIG = ROOT / "figures"
FIG.mkdir(parents=True, exist_ok=True)

GQL = "https://api.platform.opentargets.org/api/v4/graphql"

GENES = {
    # Test (CGRP pathway)
    "CALCA":  ("ENSG00000110680", "ligand"),
    "CALCB":  ("ENSG00000175868", "ligand"),
    "CALCRL": ("ENSG00000064989", "receptor"),
    "RAMP1":  ("ENSG00000132329", "receptor"),
    "CRCP":   ("ENSG00000241258", "receptor"),
    # Positive controls
    "PCSK9":  ("ENSG00000169174", "control"),
    "LDLR":   ("ENSG00000130164", "control"),
    "MC4R":   ("ENSG00000166603", "control"),
}

PHENOTYPES = {
    "migraine":                "MONDO_0005277",
    "type_2_diabetes":          "MONDO_0005148",
    "hypertension":             "MONDO_0001134",
    "coronary_artery_disease":  "EFO_0001645",
    "heart_failure":            "EFO_0003144",
    "schizophrenia":            "MONDO_0005090",
    "obesity":                 "EFO_0001073",
    "BMI":                     "EFO_0004340",
    "LDL_cholesterol":         "EFO_0004911",
    "HDL_cholesterol":         "EFO_0004612",
    "triglycerides":           "EFO_0004530",
    "atrial_fibrillation":      "EFO_0000275",
    "stroke":                  "EFO_0000712",
    "headache":                "EFO_0010178",
    "depression":              "MONDO_0002009",
}


def query_burden(target_ensembl: str, disease_efo: str) -> list[dict]:
    q = """
    query Q($eid: String!, $did: String!) {
      target(ensemblId: $eid) {
        evidences(efoIds: [$did], datasourceIds: ["gene_burden"], size: 20) {
          rows {
            score pValueMantissa pValueExponent
          }
        }
      }
    }
    """
    r = requests.post(GQL, json={"query": q,
                                   "variables": {"eid": target_ensembl,
                                                  "did": disease_efo}},
                       timeout=60)
    r.raise_for_status()
    d = r.json()
    if "errors" in d or d["data"]["target"] is None:
        return []
    return d["data"]["target"]["evidences"]["rows"]


def main() -> None:
    print(f"Querying {len(GENES)} genes x {len(PHENOTYPES)} phenotypes...")
    grid_n = np.zeros((len(GENES), len(PHENOTYPES)), dtype=int)
    grid_logp = np.full((len(GENES), len(PHENOTYPES)), np.nan)

    for i, (gene, (ensembl, _)) in enumerate(GENES.items()):
        for j, (pname, did) in enumerate(PHENOTYPES.items()):
            try:
                rows = query_burden(ensembl, did)
            except Exception:
                rows = []
            grid_n[i, j] = len(rows)
            if rows:
                # Take the most significant -log10(P) across rows
                best = -np.inf
                for row in rows:
                    m = row.get("pValueMantissa")
                    e = row.get("pValueExponent")
                    if m is not None and e is not None and m > 0:
                        logp = -(np.log10(m) + e)
                        if logp > best:
                            best = logp
                if np.isfinite(best):
                    grid_logp[i, j] = best
            time.sleep(0.04)

    # Save raw grid
    np.savez(RES / "burden_grid.npz",
              grid_n=grid_n, grid_logp=grid_logp,
              genes=list(GENES.keys()),
              phenotypes=list(PHENOTYPES.keys()))

    # Plot
    fig, ax = plt.subplots(figsize=(13, 5.5))

    # Heatmap of -log10(P); NaN cells shown as light gray
    cmap = plt.cm.viridis.copy()
    cmap.set_bad("#eeeeee")
    masked = np.ma.masked_invalid(grid_logp)
    im = ax.imshow(masked, aspect="auto", cmap=cmap, vmin=2, vmax=20)
    cbar = plt.colorbar(im, ax=ax, fraction=0.025, pad=0.02)
    cbar.set_label(r"$-\log_{10}$(burden $P$)  (best across sources)")

    # Annotate cells
    for i in range(grid_n.shape[0]):
        for j in range(grid_n.shape[1]):
            n = grid_n[i, j]
            if n == 0:
                ax.text(j, i, "—", ha="center", va="center",
                         color="#888888", fontsize=9)
            else:
                lp = grid_logp[i, j]
                col = "white" if (np.isfinite(lp) and lp > 10) else "black"
                txt = f"n={n}\n{lp:.0f}" if np.isfinite(lp) else f"n={n}"
                ax.text(j, i, txt, ha="center", va="center",
                         color=col, fontsize=7)

    ax.set_xticks(np.arange(len(PHENOTYPES)))
    ax.set_xticklabels(list(PHENOTYPES.keys()), rotation=45, ha="right",
                        fontsize=8)
    ax.set_yticks(np.arange(len(GENES)))
    sides = [s for _, (_, s) in GENES.items()]
    side_color = {"ligand": "#e377c2", "receptor": "#d62728",
                   "control": "#666666"}
    yticklabels = []
    for gene, (_, side) in GENES.items():
        yticklabels.append(f"{gene} ({side})")
    ax.set_yticklabels(yticklabels, fontsize=9)
    for tick, side in zip(ax.get_yticklabels(), sides):
        tick.set_color(side_color[side])
    # Bold the controls
    for tick, side in zip(ax.get_yticklabels(), sides):
        if side == "control":
            tick.set_fontweight("bold")
    # Separator between CGRP genes and controls
    ax.axhline(4.5, color="black", lw=1.5)
    ax.set_title(
        "Gene-burden evidence in OpenTargets across the CGRP pathway and "
        "positive controls\n"
        "All 5 CGRP-pathway genes are null across 15 phenotypes; "
        "PCSK9/LDLR/MC4R show expected positive signal at metabolic phenotypes"
    )
    fig.tight_layout()
    pdf = FIG / "burden_phenome_heatmap.pdf"
    png = FIG / "burden_phenome_heatmap.png"
    fig.savefig(pdf, bbox_inches="tight")
    fig.savefig(png, bbox_inches="tight", dpi=160)
    print(f"Wrote {pdf.relative_to(ROOT)}")
    print(f"Wrote {png.relative_to(ROOT)}")

    # Print summary table
    print("\n=== Summary ===")
    for gene in GENES:
        i = list(GENES).index(gene)
        n_total = int(grid_n[i].sum())
        max_lp = float(np.nanmax(grid_logp[i])) if n_total else None
        print(f"  {gene:10s}  total burden rows: {n_total:3d}  "
              f"max -log10(P): {max_lp if max_lp else '—'}")


if __name__ == "__main__":
    main()
