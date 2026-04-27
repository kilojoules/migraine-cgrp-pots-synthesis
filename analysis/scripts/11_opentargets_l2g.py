"""
OpenTargets Genetics L2G (locus-to-gene) score check for CGRP-pathway and
selected migraine-relevant genes. Sensitivity layer to complement the
lead-SNP-nearest gene rule used in the primary migraine set.

Reports the highest L2G score across all migraine-credible-set studies
(EFO_0003821 = migraine, MONDO_0005277) for each query gene.
"""
from __future__ import annotations

import json
from pathlib import Path

import requests

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "results"
RES.mkdir(parents=True, exist_ok=True)

GQL = "https://api.platform.opentargets.org/api/v4/graphql"

# OpenTargets disease IDs to query (migraine has multiple ontology mappings)
DISEASES = ["MONDO_0005277", "EFO_0003821"]

QUERY_GENES = [
    # CGRP pathway
    "CALCA", "CALCB", "CALCRL", "RAMP1", "CRCP",
    # Estrogen receptors
    "ESR1", "ESR2", "GPER1",
    # Adrenergic / autonomic anchors
    "ADRA2A", "ADRB1", "ADRB2", "SLC6A2", "DBH", "TH",
    # Cap-style migraine GWAS hits for sanity (known signal)
    "PRDM16", "TRPM8", "FHL5", "MEF2D",
]


def fetch_disease_associations(disease_id: str) -> list[dict]:
    """Fetch top-N gene associations for a disease, including overall and
    geneticAssociation scores (the latter is closest to L2G aggregation)."""
    query = """
    query AssocByDisease($efoId: String!, $size: Int!) {
      disease(efoId: $efoId) {
        id
        name
        associatedTargets(page: {index: 0, size: $size}) {
          count
          rows {
            target {
              id
              approvedSymbol
            }
            score
            datatypeScores {
              id
              score
            }
          }
        }
      }
    }
    """
    r = requests.post(
        GQL,
        json={"query": query, "variables": {"efoId": disease_id, "size": 3000}},
        timeout=120,
    )
    r.raise_for_status()
    payload = r.json()
    if "errors" in payload:
        raise RuntimeError(payload["errors"])
    d = payload["data"]["disease"]
    if d is None:
        return []
    return d["associatedTargets"]["rows"]


def main() -> None:
    out: dict = {"query_genes": QUERY_GENES, "by_disease": {}}
    by_gene: dict[str, dict] = {g: {"per_disease": {}} for g in QUERY_GENES}

    for d in DISEASES:
        try:
            rows = fetch_disease_associations(d)
        except Exception as e:
            print(f"  {d}: ERROR {e}")
            out["by_disease"][d] = {"error": str(e)}
            continue
        out["by_disease"][d] = {"n_rows": len(rows)}
        print(f"  {d}: {len(rows)} associated targets")
        scores_by_symbol = {
            row["target"]["approvedSymbol"]: {
                "overall_score": row["score"],
                "datatype_scores": {ds["id"]: ds["score"] for ds in (row.get("datatypeScores") or [])},
            }
            for row in rows
        }
        for g in QUERY_GENES:
            if g in scores_by_symbol:
                by_gene[g]["per_disease"][d] = scores_by_symbol[g]
            else:
                by_gene[g]["per_disease"][d] = None

    out["by_gene"] = by_gene
    out_path = RES / "opentargets_l2g_results.json"
    out_path.write_text(json.dumps(out, indent=2))
    print(f"\nWrote {out_path.relative_to(ROOT)}")

    print("\nGene | overall (best across diseases) | genetic_association (best)")
    print("-" * 70)
    for g in QUERY_GENES:
        best_overall = None
        best_genetic = None
        for d in DISEASES:
            v = by_gene[g]["per_disease"].get(d)
            if v is None:
                continue
            if best_overall is None or v["overall_score"] > best_overall:
                best_overall = v["overall_score"]
            ga = v["datatype_scores"].get("genetic_association")
            if ga is not None and (best_genetic is None or ga > best_genetic):
                best_genetic = ga
        bo = f"{best_overall:.3f}" if best_overall is not None else "—"
        bg = f"{best_genetic:.3f}" if best_genetic is not None else "—"
        print(f"  {g:8s} | {bo:>8s} | {bg:>8s}")


if __name__ == "__main__":
    main()
