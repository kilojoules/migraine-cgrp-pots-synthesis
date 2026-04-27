"""
Network-proximity test: migraine GWAS gene set vs. CGRP receptor anchor.

Statistic: Guney et al. 2016 'closest' distance d_c with degree-binned random
gene-set null. Pre-registered: 1000 permutations, seed 20260426, primary
network = BioGRID U IntAct physical-only, primary anchor = CALCRL+RAMP1
(strict CGRP-selective heterodimer), primary migraine set = lead-SNP nearest
gene from Hautakangas 2022 Supp Table 3a.

Specificity controls (gene sets pulled from DisGeNET / OpenTargets at runtime):
schizophrenia, psoriasis, essential hypertension, type 2 diabetes.
"""
import json
import pickle
import random
import time
from collections import defaultdict
from pathlib import Path

import networkx as nx
import numpy as np
import pandas as pd
import requests

ROOT = Path(__file__).resolve().parents[1]
PROC = ROOT / "data" / "processed"
RES = ROOT / "results"
RES.mkdir(parents=True, exist_ok=True)

SEED = 20260426
N_PERMUTATIONS = 1000

PRIMARY_MIGRAINE = PROC / "migraine_set_primary.tsv"
SENS1_MIGRAINE = PROC / "migraine_set_sensitivity_focus_twas.tsv"
SENS2_MIGRAINE = PROC / "migraine_set_sensitivity_nearest_plus_eqtl.tsv"
PPI_PICKLE = PROC / "ppi_graph.gpickle"


def load_genes(path: Path, col: str = "gene") -> set[str]:
    df = pd.read_csv(path, sep="\t")
    return set(df[col].dropna().astype(str).str.strip().unique())


def degree_bins(G: nx.Graph, n_bins: int = 100) -> tuple[dict[str, int], list[list[str]]]:
    """Bin genes by degree (log scale) so that random draws are degree-matched."""
    degs = dict(G.degree())
    nodes = list(degs.keys())
    log_degs = np.log1p([degs[n] for n in nodes])
    quantiles = np.quantile(log_degs, np.linspace(0, 1, n_bins + 1))
    quantiles[-1] += 1e-9  # ensure highest-degree node falls into last bin
    bin_idx = np.searchsorted(quantiles[1:-1], log_degs, side="right")
    node_to_bin = dict(zip(nodes, bin_idx.tolist()))
    bins: list[list[str]] = [[] for _ in range(n_bins)]
    for n, b in node_to_bin.items():
        bins[b].append(n)
    return node_to_bin, bins


def closest_distance(G: nx.Graph, S: set[str], T: set[str]) -> float:
    """Guney 2016 d_c(S, T) = mean over s in S of min_{t in T} d(s, t)."""
    if not S or not T:
        return float("nan")
    # Single-source shortest paths from each target node restricted to nodes in S
    target_dists: dict[str, dict[str, int]] = {}
    for t in T:
        if t in G:
            target_dists[t] = nx.single_source_shortest_path_length(G, t)
    if not target_dists:
        return float("nan")
    total = 0.0
    n_finite = 0
    for s in S:
        if s not in G:
            continue
        d_min = min((target_dists[t].get(s, np.inf) for t in target_dists), default=np.inf)
        if np.isfinite(d_min):
            total += d_min
            n_finite += 1
    return total / n_finite if n_finite > 0 else float("nan")


def sample_degree_matched(rng: random.Random, gene_set: set[str], node_to_bin: dict[str, int], bins: list[list[str]]) -> set[str]:
    """For each gene in gene_set with bin assignment, draw one replacement from same bin."""
    out: set[str] = set()
    for g in gene_set:
        if g not in node_to_bin:
            continue
        b = node_to_bin[g]
        candidates = bins[b]
        if not candidates:
            continue
        # Sample with replacement; pre-registered protocol allows duplicates within a draw.
        out.add(rng.choice(candidates))
    return out


def proximity_with_null(
    G: nx.Graph,
    S: set[str],
    T: set[str],
    node_to_bin: dict[str, int],
    bins: list[list[str]],
    n_permutations: int,
    seed: int,
) -> dict:
    """Compute observed d_c and degree-matched null distribution."""
    S_in = S & set(G.nodes())
    T_in = T & set(G.nodes())
    obs = closest_distance(G, S_in, T_in)
    rng = random.Random(seed)
    null_distances = np.empty(n_permutations, dtype=float)
    for i in range(n_permutations):
        S_rand = sample_degree_matched(rng, S_in, node_to_bin, bins)
        T_rand = sample_degree_matched(rng, T_in, node_to_bin, bins)
        null_distances[i] = closest_distance(G, S_rand, T_rand)
        if (i + 1) % 100 == 0:
            print(f"    ...{i + 1}/{n_permutations} permutations")
    null_finite = null_distances[np.isfinite(null_distances)]
    null_mean = float(null_finite.mean())
    null_sd = float(null_finite.std(ddof=1))
    z = (obs - null_mean) / null_sd if null_sd > 0 else float("nan")
    p_emp = float(2 * min((null_finite <= obs).mean(), (null_finite >= obs).mean()))
    return {
        "n_S_in_graph": len(S_in),
        "n_T_in_graph": len(T_in),
        "n_S_total": len(S),
        "n_T_total": len(T),
        "observed_dc": obs,
        "null_mean": null_mean,
        "null_sd": null_sd,
        "z_score": z,
        "empirical_p_two_sided": p_emp,
        "n_permutations": int(n_permutations),
        "n_null_finite": int(null_finite.size),
        "null_distances_sample": null_finite[:50].tolist(),
    }


# ---------- Control gene sets ----------
def fetch_disgenet_top(disease_umls: str, top_n: int) -> set[str]:
    """Fetch top-N disease-associated genes from OpenTargets (replacement for DisGeNET).

    OpenTargets GraphQL is open and stable. We query by EFO ID and rank by overall
    association score; return the top-N gene symbols.
    """
    raise NotImplementedError("Replaced by fetch_opentargets_top below")


def fetch_opentargets_top(efo_id: str, top_n: int) -> set[str]:
    query = """
    query AssocByDisease($efoId: String!, $size: Int!) {
      disease(efoId: $efoId) {
        associatedTargets(page: {index: 0, size: $size}) {
          rows {
            target {
              approvedSymbol
            }
            score
          }
        }
      }
    }
    """
    r = requests.post(
        "https://api.platform.opentargets.org/api/v4/graphql",
        json={"query": query, "variables": {"efoId": efo_id, "size": top_n}},
        timeout=60,
    )
    r.raise_for_status()
    data = r.json()
    rows = data["data"]["disease"]["associatedTargets"]["rows"]
    return {row["target"]["approvedSymbol"] for row in rows if row["target"]["approvedSymbol"]}


CONTROL_DISEASES = {
    # OpenTargets has migrated several disease IDs from EFO to MONDO; use current IDs.
    "schizophrenia": "MONDO_0005090",
    "psoriasis": "EFO_0000676",
    "essential_hypertension": "MONDO_0001134",
    "type_2_diabetes": "MONDO_0005148",
}


def main() -> None:
    t0 = time.time()
    print(f"Loading PPI graph from {PPI_PICKLE.relative_to(ROOT)}...")
    with PPI_PICKLE.open("rb") as fh:
        G = pickle.load(fh)
    print(f"  {G.number_of_nodes()} nodes, {G.number_of_edges()} edges")

    # Use largest connected component for proximity (paths must exist)
    components = sorted(nx.connected_components(G), key=len, reverse=True)
    lcc = G.subgraph(components[0]).copy()
    print(f"  LCC: {lcc.number_of_nodes()} nodes, {lcc.number_of_edges()} edges")

    print("Computing degree bins...")
    node_to_bin, bins = degree_bins(lcc, n_bins=100)
    bin_sizes = [len(b) for b in bins]
    print(f"  Bins: median size {np.median(bin_sizes):.0f}, min {min(bin_sizes)}, max {max(bin_sizes)}")

    # Anchors
    cgrp_strict = load_genes(PROC / "cgrp_anchor_strict.tsv")
    cgrp_extended = load_genes(PROC / "cgrp_anchor_extended.tsv")

    # Migraine sets
    migraine_primary = load_genes(PRIMARY_MIGRAINE)
    migraine_sens1 = load_genes(SENS1_MIGRAINE)
    migraine_sens2 = load_genes(SENS2_MIGRAINE)
    print(f"\nMigraine sets: primary={len(migraine_primary)}, "
          f"FOCUS+TWAS={len(migraine_sens1)}, nearest+eQTL={len(migraine_sens2)}")
    print(f"CGRP anchors: strict={len(cgrp_strict)} ({sorted(cgrp_strict)}), "
          f"extended={len(cgrp_extended)} ({sorted(cgrp_extended)})")

    # ---------- Primary analysis ----------
    print("\n=== PRIMARY: migraine_primary -> CGRP_strict ===")
    primary_result = proximity_with_null(
        lcc, migraine_primary, cgrp_strict, node_to_bin, bins, N_PERMUTATIONS, SEED,
    )
    print(json.dumps({k: v for k, v in primary_result.items() if k != "null_distances_sample"}, indent=2))

    # ---------- Sensitivity analyses ----------
    print("\n=== SENSITIVITY: migraine_FOCUS_TWAS -> CGRP_strict ===")
    sens1_result = proximity_with_null(
        lcc, migraine_sens1, cgrp_strict, node_to_bin, bins, N_PERMUTATIONS, SEED + 1,
    )
    print(json.dumps({k: v for k, v in sens1_result.items() if k != "null_distances_sample"}, indent=2))

    print("\n=== SENSITIVITY: migraine_nearest_plus_eQTL -> CGRP_strict ===")
    sens2_result = proximity_with_null(
        lcc, migraine_sens2, cgrp_strict, node_to_bin, bins, N_PERMUTATIONS, SEED + 2,
    )
    print(json.dumps({k: v for k, v in sens2_result.items() if k != "null_distances_sample"}, indent=2))

    print("\n=== SENSITIVITY: migraine_primary -> CGRP_extended ===")
    ext_result = proximity_with_null(
        lcc, migraine_primary, cgrp_extended, node_to_bin, bins, N_PERMUTATIONS, SEED + 3,
    )
    print(json.dumps({k: v for k, v in ext_result.items() if k != "null_distances_sample"}, indent=2))

    # ---------- Specificity controls ----------
    print("\n=== SPECIFICITY CONTROLS: control disease -> CGRP_strict ===")
    n_target = len(migraine_primary)
    control_results = {}
    for name, efo in CONTROL_DISEASES.items():
        try:
            print(f"\n[{name}] fetching top-{n_target} OpenTargets-associated genes for {efo}...")
            ctrl_set = fetch_opentargets_top(efo, n_target)
            print(f"  got {len(ctrl_set)} genes")
            r = proximity_with_null(
                lcc, ctrl_set, cgrp_strict, node_to_bin, bins, N_PERMUTATIONS, SEED + hash(name) % 1000,
            )
            r["disease"] = name
            r["efo_id"] = efo
            r["control_gene_set"] = sorted(ctrl_set)
            control_results[name] = r
            print(f"  {name}: obs={r['observed_dc']:.3f}, "
                  f"null_mean={r['null_mean']:.3f}, z={r['z_score']:.3f}, p={r['empirical_p_two_sided']:.4f}")
            time.sleep(1.0)
        except Exception as e:
            print(f"  ! {name} failed: {e}")
            control_results[name] = {"error": str(e), "disease": name, "efo_id": efo}

    # ---------- Save ----------
    out = {
        "config": {
            "seed": SEED,
            "n_permutations": N_PERMUTATIONS,
            "lcc_nodes": lcc.number_of_nodes(),
            "lcc_edges": lcc.number_of_edges(),
        },
        "primary": primary_result | {
            "migraine_set": "migraine_primary (lead-SNP nearest)",
            "cgrp_anchor": "strict (CALCRL+RAMP1)",
        },
        "sensitivity": {
            "migraine_FOCUS_TWAS_vs_strict": sens1_result,
            "migraine_nearest_plus_eqtl_vs_strict": sens2_result,
            "migraine_primary_vs_extended": ext_result,
        },
        "specificity_controls": control_results,
        "elapsed_seconds": time.time() - t0,
    }
    out_path = RES / "network_proximity_results.json"
    with out_path.open("w") as fh:
        json.dump(out, fh, indent=2)
    print(f"\nWrote results -> {out_path.relative_to(ROOT)}")
    print(f"Elapsed: {out['elapsed_seconds']:.1f}s")


if __name__ == "__main__":
    main()
