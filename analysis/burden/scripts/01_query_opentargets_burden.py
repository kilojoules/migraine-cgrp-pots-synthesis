"""
Extension B: query OpenTargets GraphQL for gene-burden evidence on the
five CGRP-pathway genes against migraine and four cardiovascular/metabolic
phenotypes plus a schizophrenia negative control.

Decision rule, gene/disease lists, and seed are locked in
analysis/burden/PRE_REGISTRATION_burden.md.
"""
from __future__ import annotations

import json
import time
from pathlib import Path

import requests

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "results"
RES.mkdir(parents=True, exist_ok=True)

GQL = "https://api.platform.opentargets.org/api/v4/graphql"
SEED = 20260430

GENES = {
    "CALCA":  {"ensembl": "ENSG00000110680", "side": "ligand"},
    "CALCB":  {"ensembl": "ENSG00000175868", "side": "ligand"},
    "CALCRL": {"ensembl": "ENSG00000064989", "side": "receptor"},
    "RAMP1":  {"ensembl": "ENSG00000132329", "side": "receptor"},
    "CRCP":   {"ensembl": "ENSG00000241258", "side": "receptor"},
}

# OpenTargets disease IDs. We probe a few alternates per phenotype because
# OpenTargets sometimes routes through MONDO and sometimes EFO.
DISEASES = {
    "migraine":               ["MONDO_0005277", "EFO_0003821"],
    "type_2_diabetes":         ["MONDO_0005148", "EFO_0001360"],
    "hypertension":            ["MONDO_0001134", "EFO_0000537"],
    "coronary_artery_disease": ["EFO_0001645",   "MONDO_0005010"],
    "heart_failure":           ["EFO_0003144",   "MONDO_0005009"],
    "schizophrenia":           ["MONDO_0005090", "EFO_0000692"],
}


def query_evidence_rows(target_ensembl: str, disease_efo: str,
                          datasource: str = "gene_burden",
                          page_size: int = 50) -> list[dict]:
    """Fetch evidence rows for a target-disease pair filtered by datasource."""
    query = """
    query EvidByTargetDisease($eid: String!, $did: String!, $size: Int!,
                                $ds: [String!]) {
      target(ensemblId: $eid) {
        evidences(efoIds: [$did], datasourceIds: $ds,
                    size: $size) {
          count
          rows {
            datasourceId
            datatypeId
            score
            disease { id name }
            target { approvedSymbol }
            studyId
            statisticalMethod
            statisticalMethodOverview
            studySampleSize
            studyCases
            pValueExponent
            pValueMantissa
            beta
            betaConfidenceIntervalLower
            betaConfidenceIntervalUpper
            oddsRatio
            oddsRatioConfidenceIntervalLower
            oddsRatioConfidenceIntervalUpper
            cohortId
            ancestry
            ancestryId
            allelicRequirements
          }
        }
      }
    }
    """
    r = requests.post(
        GQL,
        json={"query": query, "variables": {
            "eid": target_ensembl, "did": disease_efo,
            "size": page_size, "ds": [datasource],
        }},
        timeout=120,
    )
    r.raise_for_status()
    payload = r.json()
    if "errors" in payload:
        raise RuntimeError(payload["errors"])
    target = payload["data"]["target"]
    if target is None:
        return []
    ev = target["evidences"]
    return ev["rows"]


def main() -> None:
    out: dict = {
        "spec": "analysis/burden/PRE_REGISTRATION_burden.md",
        "seed": SEED,
        "datasource_filter": "gene_burden",
        "by_gene": {},
    }
    print(f"Querying OpenTargets gene-burden evidence for "
          f"{len(GENES)} genes x {len(DISEASES)} phenotypes...")
    for gene, ginfo in GENES.items():
        out["by_gene"][gene] = {"side": ginfo["side"], "by_disease": {}}
        for disease_label, ids in DISEASES.items():
            rows: list[dict] = []
            tried = []
            for did in ids:
                try:
                    rows = query_evidence_rows(ginfo["ensembl"], did,
                                                 "gene_burden")
                    tried.append({"id": did, "n_rows": len(rows)})
                    if rows:
                        break  # Use the first hit
                except Exception as e:
                    tried.append({"id": did, "error": str(e)})
                time.sleep(0.05)
            out["by_gene"][gene]["by_disease"][disease_label] = {
                "ids_tried": tried,
                "rows": rows,
            }
            print(f"  {gene:8s}  {disease_label:25s}  rows={len(rows)}")

    out_path = RES / "burden_opentargets.json"
    out_path.write_text(json.dumps(out, indent=2, default=str))
    print(f"\nWrote {out_path.relative_to(ROOT)}  "
          f"({out_path.stat().st_size / 1024:.1f} KB)")

    # Concise summary table
    print("\n=== Summary: gene-burden evidence rows per (gene, disease) ===")
    diseases = list(DISEASES.keys())
    print(f"{'gene':10s}  " + "  ".join(f"{d[:18]:>18s}" for d in diseases))
    for gene in GENES:
        cells = []
        for d in diseases:
            n = len(out["by_gene"][gene]["by_disease"][d]["rows"])
            cells.append(f"{n:>18d}" if n else f"{'-':>18s}")
        print(f"{gene:10s}  " + "  ".join(cells))


if __name__ == "__main__":
    main()
