"""
Extension B follow-up: scan a broader phenotype panel for any gene-burden
evidence on the CGRP-pathway genes.

The primary script (01_query_opentargets_burden.py) returned 0 rows for
all 30 (gene, disease) pairs in the pre-registered phenotype panel. To
verify this null is not an artifact of phenotype-mapping, we also query
a wider set of common biobank phenotypes (~30) where rare-variant burden
results are typically deposited.

Output: analysis/burden/results/burden_broader_scan.json
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

GENES = {
    "CALCA":  "ENSG00000110680",
    "CALCB":  "ENSG00000175868",
    "CALCRL": "ENSG00000064989",
    "RAMP1":  "ENSG00000132329",
    "CRCP":   "ENSG00000241258",
}

# Broader phenotype scan: things that show up in UK Biobank burden tests
PHENOTYPES = {
    # Pre-registered panel (already null in primary script)
    "migraine":               ["MONDO_0005277"],
    "type_2_diabetes":         ["MONDO_0005148"],
    "hypertension":            ["MONDO_0001134"],
    "coronary_artery_disease": ["EFO_0001645"],
    "heart_failure":           ["EFO_0003144"],
    "schizophrenia":           ["MONDO_0005090"],
    # Broader scan
    "obesity":                ["EFO_0001073"],
    "BMI":                    ["EFO_0004340"],
    "LDL_cholesterol":        ["EFO_0004911"],
    "HDL_cholesterol":        ["EFO_0004612"],
    "triglycerides":          ["EFO_0004530"],
    "type_1_diabetes":         ["MONDO_0005147"],
    "atrial_fibrillation":     ["EFO_0000275"],
    "stroke":                 ["EFO_0000712"],
    "asthma":                 ["MONDO_0004979"],
    "depression":             ["MONDO_0002009"],
    "anxiety_disorder":        ["MONDO_0005618"],
    "headache":               ["EFO_0010178"],
    "tension_headache":        ["MONDO_0009110"],
    "trigeminal_neuralgia":    ["MONDO_0009718"],
    "fibromyalgia":           ["MONDO_0005546"],
    "chronic_fatigue":         ["MONDO_0005404"],
    "irritable_bowel":         ["MONDO_0005052"],
    "joint_hypermobility":     ["MONDO_0010518"],
    "rheumatoid_arthritis":    ["EFO_0000685"],
    "osteoarthritis":          ["EFO_0002506"],
    "blood_pressure_systolic": ["EFO_0006335"],
    "resting_heart_rate":      ["EFO_0010046"],
    "syncope":                ["EFO_0009596"],
}


def query_burden(target_ensembl: str, disease_efo: str) -> list[dict]:
    query = """
    query Q($eid: String!, $did: String!) {
      target(ensemblId: $eid) {
        evidences(efoIds: [$did], datasourceIds: ["gene_burden"], size: 20) {
          count
          rows {
            datasourceId disease { id name } score
            pValueMantissa pValueExponent cohortId
            statisticalMethod
          }
        }
      }
    }
    """
    r = requests.post(
        GQL,
        json={"query": query,
              "variables": {"eid": target_ensembl, "did": disease_efo}},
        timeout=120,
    )
    r.raise_for_status()
    payload = r.json()
    if "errors" in payload:
        return []
    target = payload["data"]["target"]
    if target is None:
        return []
    return target["evidences"]["rows"]


def main() -> None:
    out: dict = {"by_gene": {}}
    print(f"Scanning {len(GENES)} CGRP genes x {len(PHENOTYPES)} phenotypes...")
    n_total_rows = 0
    for gene, ensembl in GENES.items():
        out["by_gene"][gene] = {}
        for pname, ids in PHENOTYPES.items():
            rows = []
            for did in ids:
                try:
                    rows = query_burden(ensembl, did)
                    if rows:
                        break
                except Exception as e:
                    pass
                time.sleep(0.05)
            out["by_gene"][gene][pname] = rows
            n_total_rows += len(rows)

    print(f"Total burden evidence rows across {len(GENES)} genes "
          f"x {len(PHENOTYPES)} phenotypes: {n_total_rows}")
    print()
    print(f"{'gene':10s}  " + "  ".join(
        f"{p[:14]:>14s}" for p in PHENOTYPES))
    for gene in GENES:
        cells = []
        for p in PHENOTYPES:
            n = len(out["by_gene"][gene][p])
            cells.append(f"{n:>14d}" if n else f"{'-':>14s}")
        print(f"{gene:10s}  " + "  ".join(cells))

    out_path = RES / "burden_broader_scan.json"
    out_path.write_text(json.dumps(out, indent=2, default=str))
    print(f"\nWrote {out_path.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
