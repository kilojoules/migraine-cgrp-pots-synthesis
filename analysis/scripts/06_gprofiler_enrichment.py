"""
Pathway / ontology enrichment of the migraine gene set using g:Profiler.

Per pre-registration: g:GOSt over GO_BP, GO_MF, KEGG, Reactome with g:SCS correction.
Background: all protein-coding HGNC genes (g:Profiler default 'organism: hsapiens').
"""
import json
from pathlib import Path

import pandas as pd
import requests

ROOT = Path(__file__).resolve().parents[1]
PROC = ROOT / "data" / "processed"
RES = ROOT / "results"
RES.mkdir(parents=True, exist_ok=True)


def gprofiler_enrich(genes: list[str]) -> dict:
    url = "https://biit.cs.ut.ee/gprofiler/api/gost/profile/"
    payload = {
        "organism": "hsapiens",
        "query": genes,
        "sources": ["GO:BP", "GO:MF", "KEGG", "REAC"],
        "user_threshold": 0.05,
        "all_results": False,
        "ordered": False,
        "no_evidences": True,
        "significance_threshold_method": "g_SCS",
        "domain_scope": "annotated",
    }
    r = requests.post(url, json=payload, timeout=120)
    r.raise_for_status()
    return r.json()


def main() -> None:
    primary = pd.read_csv(PROC / "migraine_set_primary.tsv", sep="\t")
    genes = sorted(primary["gene"].dropna().astype(str).str.strip().unique())
    print(f"Submitting {len(genes)} migraine genes to g:Profiler...")
    result = gprofiler_enrich(genes)

    n_terms = len(result.get("result", []))
    print(f"  {n_terms} significant terms returned")

    if n_terms:
        df = pd.DataFrame(result["result"])
        # Sort by adjusted p-value
        df = df.sort_values("p_value")
        keep_cols = [
            c
            for c in [
                "source",
                "native",
                "name",
                "p_value",
                "intersection_size",
                "term_size",
                "query_size",
                "effective_domain_size",
            ]
            if c in df.columns
        ]
        top = df[keep_cols].head(40).copy()
        top_path = RES / "gprofiler_top40.tsv"
        top.to_csv(top_path, sep="\t", index=False)

        full_path = RES / "gprofiler_full.tsv"
        df[keep_cols].to_csv(full_path, sep="\t", index=False)

        # Look for autonomic / adrenergic / CGRP-relevant hits explicitly
        kw_pattern = (
            r"adrener|sympat|parasympat|autonom|catecholamine|nore"
            r"|cgrp|calcit|trigemin|baroreceptor|vascular|pressure"
        )
        flags = df[df["name"].astype(str).str.contains(kw_pattern, case=False, regex=True)]
        flag_path = RES / "gprofiler_autonomic_relevant.tsv"
        flags[keep_cols].to_csv(flag_path, sep="\t", index=False)
        print(f"\n  Top {len(top)} terms -> {top_path.relative_to(ROOT)}")
        print(f"  Full {n_terms} -> {full_path.relative_to(ROOT)}")
        print(f"  Autonomic/CGRP-relevant ({len(flags)} terms) -> {flag_path.relative_to(ROOT)}")
        if len(flags):
            print("  Autonomic/CGRP hits:")
            for _, row in flags.head(20).iterrows():
                print(f"    [{row['source']}] {row['native']} {row['name']}  "
                      f"p={row['p_value']:.2e}  k={row['intersection_size']}/{row['term_size']}")

    out_path = RES / "gprofiler_raw_response.json"
    with out_path.open("w") as fh:
        json.dump(result, fh, indent=2)
    print(f"\nFull response -> {out_path.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
