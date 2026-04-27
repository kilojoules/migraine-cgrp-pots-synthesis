"""
Build the human protein-protein interaction graph for Phase 2.

Per pre-registration: union of BioGRID Homo sapiens **physical-only** interactions
and IntAct human binary interactions. Edges deduplicated, self-loops removed,
nodes mapped to HGNC official symbols.

Output: pickled networkx Graph + edgelist TSV + provenance JSON.
"""
from __future__ import annotations

import csv
import gzip
import hashlib
import json
import pickle
import re
import time
from pathlib import Path

import networkx as nx
import requests

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw"
OUT = ROOT / "data" / "processed"
OUT.mkdir(parents=True, exist_ok=True)

BIOGRID_HUMAN = RAW / "biogrid_homo_sapiens.tab3.txt"
INTACT_ZIP = RAW / "intact.zip"
HUMAN_TAXID = "9606"


def file_sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def parse_biogrid(path: Path) -> set[tuple[str, str]]:
    """Return set of (symA, symB) physical-interaction edges, normalized as sorted tuple, human-only, no self-loops."""
    edges: set[tuple[str, str]] = set()
    n_total = 0
    n_kept = 0
    with path.open(newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            n_total += 1
            if row.get("Experimental System Type") != "physical":
                continue
            if row.get("Organism ID Interactor A") != HUMAN_TAXID:
                continue
            if row.get("Organism ID Interactor B") != HUMAN_TAXID:
                continue
            a = (row.get("Official Symbol Interactor A") or "").strip()
            b = (row.get("Official Symbol Interactor B") or "").strip()
            if not a or not b or a == "-" or b == "-":
                continue
            if a == b:
                continue
            edges.add(tuple(sorted([a, b])))
            n_kept += 1
    print(f"  BioGRID: read {n_total} rows, kept {n_kept} physical human-human, "
          f"{len(edges)} unique undirected edges")
    return edges


# ---------- IntAct ----------
# IntAct PSI-MITAB 2.7 columns (we only need a few). Field indices per
# https://psicquic.github.io/MITAB27Format.html
INTACT_COL_ID_A = 0      # ID(s) interactor A (uniprotkb:Pxxxxx | intact:EBI-...)
INTACT_COL_ID_B = 1      # ID(s) interactor B
INTACT_COL_ALT_A = 2     # alternative ID(s)
INTACT_COL_ALT_B = 3     # alternative ID(s)
INTACT_COL_ALIAS_A = 4   # alias(es) — often contains gene name
INTACT_COL_ALIAS_B = 5   # alias(es)
INTACT_COL_TAXID_A = 9
INTACT_COL_TAXID_B = 10
INTACT_COL_INT_TYPE = 11  # interaction type — we want physical


GENE_NAME_RE = re.compile(r"(?:uniprotkb|psi-mi):([^\s(|]+)\(gene\s*name\)")


def extract_gene_symbol_from_intact_alias(alias_field: str) -> str | None:
    """IntAct alias field example: 'uniprotkb:CALCA(gene name)|psi-mi:"CALC_HUMAN"(display_long)'."""
    if not alias_field or alias_field == "-":
        return None
    m = GENE_NAME_RE.search(alias_field)
    if m:
        return m.group(1).strip().upper()
    return None


def parse_intact_zip(path: Path) -> set[tuple[str, str]]:
    """Stream IntAct PSI-MITAB inside the zip and extract human physical edges keyed on gene names."""
    import zipfile

    edges: set[tuple[str, str]] = set()
    n_total = 0
    n_kept = 0
    n_no_gene = 0
    with zipfile.ZipFile(path) as zf:
        # IntAct main file is usually intact.txt inside
        names = [n for n in zf.namelist() if n.endswith(".txt")]
        if not names:
            raise RuntimeError(f"No .txt found in {path}")
        # Take the largest .txt
        names.sort(key=lambda n: zf.getinfo(n).file_size, reverse=True)
        main_name = names[0]
        print(f"  IntAct main file: {main_name}")
        with zf.open(main_name) as fh:
            header = fh.readline()  # skip header line
            for raw_line in fh:
                n_total += 1
                line = raw_line.decode("utf-8", errors="replace")
                cols = line.rstrip("\n").split("\t")
                if len(cols) < 12:
                    continue
                taxid_a = cols[INTACT_COL_TAXID_A]
                taxid_b = cols[INTACT_COL_TAXID_B]
                # taxid format: 'taxid:9606(...)' — substring check is enough
                if f"taxid:{HUMAN_TAXID}" not in taxid_a:
                    continue
                if f"taxid:{HUMAN_TAXID}" not in taxid_b:
                    continue
                # Filter to physical interactions (interaction type contains 'physical')
                # IntAct uses MI ontology terms; 'physical association' MI:0915 is common.
                # We accept anything containing 'physical' or 'direct interaction' or 'association'.
                int_type = cols[INTACT_COL_INT_TYPE].lower()
                if not (
                    "physical" in int_type
                    or "direct interaction" in int_type
                    or "association" in int_type
                ):
                    continue
                # Extract gene symbols from alias fields
                a = extract_gene_symbol_from_intact_alias(cols[INTACT_COL_ALIAS_A])
                b = extract_gene_symbol_from_intact_alias(cols[INTACT_COL_ALIAS_B])
                if not a or not b:
                    n_no_gene += 1
                    continue
                if a == b:
                    continue
                edges.add(tuple(sorted([a, b])))
                n_kept += 1
                if n_total % 500_000 == 0:
                    print(f"    ...{n_total} rows, {n_kept} kept, {len(edges)} unique edges")

    print(
        f"  IntAct: read {n_total} rows, kept {n_kept} human physical with gene names, "
        f"{len(edges)} unique edges; {n_no_gene} skipped (no gene name)"
    )
    return edges


def main(use_intact: bool = True) -> None:
    print("Parsing BioGRID...")
    biogrid_edges = parse_biogrid(BIOGRID_HUMAN)

    intact_edges: set[tuple[str, str]] = set()
    if use_intact and INTACT_ZIP.exists():
        try:
            print("Parsing IntAct...")
            intact_edges = parse_intact_zip(INTACT_ZIP)
        except Exception as e:
            print(f"!! IntAct parsing failed: {e}; continuing with BioGRID only")
            intact_edges = set()

    union_edges = biogrid_edges | intact_edges
    overlap = biogrid_edges & intact_edges
    print(
        f"\nUnion: {len(union_edges)} edges "
        f"(BioGRID {len(biogrid_edges)}, IntAct {len(intact_edges)}, overlap {len(overlap)})"
    )

    G = nx.Graph()
    G.add_edges_from(union_edges)
    print(
        f"Graph: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges, "
        f"density {nx.density(G):.6f}"
    )

    # Largest connected component
    components = sorted(nx.connected_components(G), key=len, reverse=True)
    lcc = G.subgraph(components[0]).copy()
    print(
        f"Largest connected component: {lcc.number_of_nodes()} nodes "
        f"({100 * lcc.number_of_nodes() / G.number_of_nodes():.1f}%), "
        f"{lcc.number_of_edges()} edges"
    )

    # Save
    edgelist_path = OUT / "ppi_edgelist.tsv.gz"
    with gzip.open(edgelist_path, "wt") as fh:
        fh.write("source\ttarget\n")
        for a, b in sorted(union_edges):
            fh.write(f"{a}\t{b}\n")
    print(f"Wrote edgelist -> {edgelist_path.relative_to(ROOT)}")

    pickle_path = OUT / "ppi_graph.gpickle"
    with pickle_path.open("wb") as fh:
        pickle.dump(G, fh, protocol=pickle.HIGHEST_PROTOCOL)
    print(f"Wrote pickle -> {pickle_path.relative_to(ROOT)}")

    # Provenance
    provenance = {
        "biogrid": {
            "file": str(BIOGRID_HUMAN.relative_to(ROOT)),
            "sha256": file_sha256(BIOGRID_HUMAN) if BIOGRID_HUMAN.exists() else None,
            "filter": "Experimental System Type == physical, both organisms taxid 9606",
            "n_unique_undirected_edges": len(biogrid_edges),
        },
        "intact": {
            "file": str(INTACT_ZIP.relative_to(ROOT)) if INTACT_ZIP.exists() else None,
            "sha256": file_sha256(INTACT_ZIP) if INTACT_ZIP.exists() else None,
            "filter": "human-human, interaction type contains 'physical', 'direct interaction', or 'association'",
            "n_unique_undirected_edges": len(intact_edges),
        },
        "union": {
            "n_edges": len(union_edges),
            "n_nodes": G.number_of_nodes(),
            "n_overlap_with_biogrid_intact": len(overlap),
            "lcc_n_nodes": lcc.number_of_nodes(),
            "lcc_n_edges": lcc.number_of_edges(),
        },
        "build_timestamp_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
    }
    with (OUT / "ppi_provenance.json").open("w") as fh:
        json.dump(provenance, fh, indent=2)


if __name__ == "__main__":
    main()
