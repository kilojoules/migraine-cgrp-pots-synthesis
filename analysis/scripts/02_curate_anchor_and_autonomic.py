"""
Build CGRP receptor anchor and autonomic/adrenergic gene sets per pre-registration.

CGRP anchor:
  - strict (primary): CALCRL, RAMP1
  - extended (sensitivity): CALCA, CALCB, CALCRL, RAMP1, CRCP

Autonomic set: union of
  - GO term memberships (queried via mygene.info)
  - KEGG hsa04261 (Adrenergic signaling in cardiomyocytes)
  - KEGG hsa04270 (Vascular smooth muscle contraction)
  - manual seed genes (pre-registered)
"""
import json
import time
from pathlib import Path

import requests

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data" / "processed"
OUT.mkdir(parents=True, exist_ok=True)

# ---------- CGRP anchor ----------
CGRP_STRICT = ["CALCRL", "RAMP1"]
CGRP_EXTENDED = ["CALCA", "CALCB", "CALCRL", "RAMP1", "CRCP"]

# ---------- Autonomic set: pre-registered sources ----------
GO_TERMS = {
    # Pre-registered term name -> GO ID (resolved via QuickGO / OLS lookups)
    "sympathetic nervous system development": "GO:0048485",
    "parasympathetic nervous system development": "GO:0048486",
    "regulation of heart rate by chemical signal": "GO:0003062",
    "adrenergic receptor signaling pathway": "GO:0071875",
    "norepinephrine metabolic process": "GO:0042415",
    "regulation of systemic arterial blood pressure by norepinephrine-epinephrine": "GO:0001994",
    "baroreceptor response to increased systemic arterial blood pressure": "GO:0001978",
}

KEGG_PATHWAYS = {
    "hsa04261": "Adrenergic signaling in cardiomyocytes",
    "hsa04270": "Vascular smooth muscle contraction",
}

MANUAL_SEEDS = [
    "ADRB1", "ADRB2", "ADRB3",
    "ADRA1A", "ADRA1B", "ADRA1D",
    "ADRA2A", "ADRA2B", "ADRA2C",
    "SLC6A2", "COMT", "GNB3",
    "DBH", "TH", "PNMT",
    "CHRM2", "CHRM3", "CHRNB1", "CHRNA3",
    "NPY", "NPY1R", "NPY2R",
]


def query_mygene_for_go(go_id: str) -> set[str]:
    """Use mygene.info to find human gene symbols annotated with a GO term."""
    url = "https://mygene.info/v3/query"
    # mygene `go` field expects the numeric portion with leading zeros; strip 'GO:' prefix only.
    bare = go_id.replace("GO:", "")
    params = {
        "q": f"go:{bare}",
        "species": "human",
        "fields": "symbol,type_of_gene",
        "size": 1000,
    }
    r = requests.get(url, params=params, timeout=30)
    r.raise_for_status()
    data = r.json()
    syms = set()
    for hit in data.get("hits", []):
        s = hit.get("symbol")
        t = hit.get("type_of_gene")
        if s and t == "protein-coding":
            syms.add(s)
    return syms


def query_kegg_pathway_genes(kegg_id: str) -> set[str]:
    """Use KEGG REST to get hsa pathway gene members, then map to HGNC symbols via mygene."""
    # Step 1: KEGG REST - returns 'path:hsaXXXXX\thsa:NNNNN' lines
    r = requests.get(f"https://rest.kegg.jp/link/hsa/pathway:{kegg_id}", timeout=30)
    r.raise_for_status()
    entrez_ids = []
    for line in r.text.strip().splitlines():
        parts = line.split("\t")
        if len(parts) == 2 and parts[1].startswith("hsa:"):
            entrez_ids.append(parts[1].replace("hsa:", ""))
    if not entrez_ids:
        return set()
    # Step 2: mygene batch lookup entrez -> symbol (POST /v3/gene, not /genes)
    url = "https://mygene.info/v3/gene"
    r = requests.post(
        url,
        data={"ids": ",".join(entrez_ids), "fields": "symbol,type_of_gene"},
        timeout=60,
    )
    r.raise_for_status()
    syms = set()
    for hit in r.json():
        s = hit.get("symbol")
        t = hit.get("type_of_gene")
        if s and t == "protein-coding":
            syms.add(s)
    return syms


def main() -> None:
    # CGRP anchor TSVs
    with (OUT / "cgrp_anchor_strict.tsv").open("w") as fh:
        fh.write("gene\tnote\n")
        for g in CGRP_STRICT:
            fh.write(f"{g}\tCGRP-selective receptor heterodimer (CALCRL+RAMP1)\n")
    with (OUT / "cgrp_anchor_extended.tsv").open("w") as fh:
        fh.write("gene\tnote\n")
        for g in CGRP_EXTENDED:
            fh.write(f"{g}\tCALCA/CALCB ligands + CALCRL+RAMP1 receptor + CRCP coreceptor\n")
    print(f"CGRP anchor strict: {len(CGRP_STRICT)} genes ({CGRP_STRICT})")
    print(f"CGRP anchor extended: {len(CGRP_EXTENDED)} genes ({CGRP_EXTENDED})")

    # Autonomic set
    autonomic_provenance: dict[str, list[str]] = {}  # gene -> list of source labels

    print("\nQuerying GO term memberships via mygene.info...")
    for term, go_id in GO_TERMS.items():
        try:
            genes = query_mygene_for_go(go_id)
            print(f"  {go_id} ({term}): {len(genes)} protein-coding genes")
            for g in genes:
                autonomic_provenance.setdefault(g, []).append(f"GO:{go_id}")
            time.sleep(0.5)  # be polite
        except Exception as e:
            print(f"  ! {go_id} failed: {e}")

    print("\nQuerying KEGG pathway memberships...")
    for kegg_id, name in KEGG_PATHWAYS.items():
        try:
            genes = query_kegg_pathway_genes(kegg_id)
            print(f"  {kegg_id} ({name}): {len(genes)} genes")
            for g in genes:
                autonomic_provenance.setdefault(g, []).append(f"KEGG:{kegg_id}")
            time.sleep(0.5)
        except Exception as e:
            print(f"  ! {kegg_id} failed: {e}")

    print("\nAdding manual seed genes...")
    for g in MANUAL_SEEDS:
        autonomic_provenance.setdefault(g, []).append("MANUAL_SEED")

    # Write TSV
    auto_path = OUT / "autonomic_set.tsv"
    with auto_path.open("w") as fh:
        fh.write("gene\tsources\tn_sources\n")
        for g in sorted(autonomic_provenance):
            srcs = autonomic_provenance[g]
            fh.write(f"{g}\t{';'.join(srcs)}\t{len(srcs)}\n")
    print(f"\nAutonomic set: {len(autonomic_provenance)} unique genes -> {auto_path.relative_to(ROOT)}")

    # Provenance JSON
    provenance = {
        "cgrp_anchor_strict": {"genes": CGRP_STRICT, "rationale": "CGRP-selective heterodimer"},
        "cgrp_anchor_extended": {
            "genes": CGRP_EXTENDED,
            "rationale": "Ligands + receptor + coreceptor",
        },
        "autonomic_set": {
            "go_terms": GO_TERMS,
            "kegg_pathways": KEGG_PATHWAYS,
            "manual_seeds": MANUAL_SEEDS,
            "n_unique_genes": len(autonomic_provenance),
            "data_sources": [
                "mygene.info v3 (https://mygene.info/v3/)",
                "KEGG REST API (https://rest.kegg.jp/)",
            ],
        },
    }
    with (OUT / "autonomic_anchor_provenance.json").open("w") as fh:
        json.dump(provenance, fh, indent=2)


if __name__ == "__main__":
    main()
