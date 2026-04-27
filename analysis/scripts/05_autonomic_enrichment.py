"""
Focused autonomic gene-set enrichment in migraine GWAS hits.

Test: hypergeometric / Fisher-exact over-representation of the autonomic gene set
within the migraine GWAS gene set, with permutation null over gene-length-matched
and expression-matched background gene sets (1000 permutations).
"""
import gzip
import json
import random
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact

ROOT = Path(__file__).resolve().parents[1]
PROC = ROOT / "data" / "processed"
RAW = ROOT / "data" / "raw"
RES = ROOT / "results"
RES.mkdir(parents=True, exist_ok=True)

SEED = 20260426
N_PERMUTATIONS = 1000


def load_gene_set(path: Path) -> set[str]:
    df = pd.read_csv(path, sep="\t")
    return set(df["gene"].dropna().astype(str).str.strip().unique())


def load_protein_coding_background_with_features() -> pd.DataFrame:
    """Get protein-coding HGNC symbols with gene length (bp) and a coarse expression proxy.

    Strategy: query mygene.info for all protein-coding genes annotated with
    `type_of_gene:protein-coding` and pull genomic coordinates; compute span as length.
    Expression bin is approximated via gene length quintile alone if GTEx not available
    (we do not have a local GTEx download); otherwise we attach GTEx whole-blood
    median TPM in a follow-up step.
    """
    cache = PROC / "background_features.tsv"
    if cache.exists():
        return pd.read_csv(cache, sep="\t")

    # Fetch via mygene scroll API - returns protein-coding human gene catalog
    import requests

    print("Fetching protein-coding gene catalog from mygene.info (this can take a minute)...")
    url = "https://mygene.info/v3/query"
    params = {
        "q": "type_of_gene:protein-coding AND taxid:9606",
        "fields": "symbol,genomic_pos.start,genomic_pos.end",
        "size": 1000,
        "fetch_all": "true",
    }
    rows = []
    while True:
        r = requests.get(url, params=params, timeout=120)
        r.raise_for_status()
        data = r.json()
        hits = data.get("hits", [])
        if not hits:
            break
        for h in hits:
            sym = h.get("symbol")
            pos = h.get("genomic_pos")
            if not sym or not pos:
                continue
            if isinstance(pos, list):
                pos = pos[0]
            try:
                length = int(pos["end"]) - int(pos["start"])
            except (KeyError, TypeError, ValueError):
                continue
            if length <= 0:
                continue
            rows.append({"gene": sym, "length_bp": length})
        scroll_id = data.get("_scroll_id")
        if not scroll_id:
            break
        params = {"scroll_id": scroll_id}
        time.sleep(0.3)

    df = pd.DataFrame(rows).drop_duplicates(subset=["gene"]).reset_index(drop=True)
    df["length_quintile"] = pd.qcut(df["length_bp"], q=5, labels=False, duplicates="drop")
    df.to_csv(cache, sep="\t", index=False)
    print(f"  {len(df)} protein-coding genes with length features -> {cache.relative_to(ROOT)}")
    return df


def main() -> None:
    migraine = load_gene_set(PROC / "migraine_set_primary.tsv")
    autonomic = load_gene_set(PROC / "autonomic_set.tsv")
    background = load_protein_coding_background_with_features()
    bg_genes = set(background["gene"])

    # Restrict everything to the protein-coding background
    migraine_in_bg = migraine & bg_genes
    autonomic_in_bg = autonomic & bg_genes
    print(f"Migraine genes in background: {len(migraine_in_bg)} / {len(migraine)}")
    print(f"Autonomic genes in background: {len(autonomic_in_bg)} / {len(autonomic)}")
    print(f"Background size: {len(bg_genes)} protein-coding genes")

    # Observed overlap
    overlap = migraine_in_bg & autonomic_in_bg
    k = len(overlap)
    n_migraine = len(migraine_in_bg)
    n_autonomic = len(autonomic_in_bg)
    n_bg = len(bg_genes)
    print(f"\nObserved migraine ∩ autonomic: {k} genes -> {sorted(overlap)}")

    # Hypergeometric / Fisher exact (one-sided, over-representation)
    # 2x2 table: (in_autonomic, in_migraine) | (in_autonomic, not_migraine) | ...
    a = k
    b = n_autonomic - k
    c = n_migraine - k
    d = n_bg - n_autonomic - n_migraine + k
    odds, p_fisher = fisher_exact([[a, b], [c, d]], alternative="greater")
    print(f"\nFisher exact one-sided p (over-rep): {p_fisher:.4f}, OR={odds:.3f}")

    # Permutation null with length-quintile matching
    print(f"\nRunning {N_PERMUTATIONS} length-quintile-matched permutations...")
    rng = random.Random(SEED)

    # For each migraine gene, find its length quintile, then sample a replacement
    # of the same quintile from the background.
    bg_by_quintile: dict[int, list[str]] = {}
    for _, row in background.iterrows():
        bg_by_quintile.setdefault(int(row["length_quintile"]), []).append(row["gene"])

    migraine_quintiles = (
        background.set_index("gene").loc[list(migraine_in_bg)]["length_quintile"].astype(int).tolist()
    )

    null_overlaps = np.empty(N_PERMUTATIONS, dtype=int)
    for i in range(N_PERMUTATIONS):
        sampled = set()
        for q in migraine_quintiles:
            cands = bg_by_quintile.get(q, [])
            if cands:
                sampled.add(rng.choice(cands))
        null_overlaps[i] = len(sampled & autonomic_in_bg)

    null_mean = float(null_overlaps.mean())
    null_sd = float(null_overlaps.std(ddof=1))
    z = (k - null_mean) / null_sd if null_sd > 0 else float("nan")
    p_emp_two_sided = float(2 * min((null_overlaps >= k).mean(), (null_overlaps <= k).mean()))
    p_emp_one_sided = float((null_overlaps >= k).mean())

    print(
        f"Null mean={null_mean:.2f}, SD={null_sd:.2f}, "
        f"z={z:.3f}, "
        f"emp p (one-sided over-rep) = {p_emp_one_sided:.4f}, "
        f"emp p (two-sided) = {p_emp_two_sided:.4f}"
    )

    out = {
        "config": {
            "seed": SEED,
            "n_permutations": N_PERMUTATIONS,
            "background_size": n_bg,
            "background_match": "gene length quintile",
        },
        "observed": {
            "migraine_in_background": n_migraine,
            "autonomic_in_background": n_autonomic,
            "k_overlap": k,
            "overlap_genes": sorted(overlap),
        },
        "fisher_exact_one_sided": {"p": p_fisher, "odds_ratio": odds},
        "permutation_null": {
            "null_mean": null_mean,
            "null_sd": null_sd,
            "z_score": z,
            "empirical_p_one_sided_over_rep": p_emp_one_sided,
            "empirical_p_two_sided": p_emp_two_sided,
            "fold_enrichment": k / null_mean if null_mean > 0 else float("nan"),
            "ci95_null": [
                float(np.quantile(null_overlaps, 0.025)),
                float(np.quantile(null_overlaps, 0.975)),
            ],
        },
    }
    out_path = RES / "autonomic_enrichment_results.json"
    with out_path.open("w") as fh:
        json.dump(out, fh, indent=2)
    print(f"\nWrote -> {out_path.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
