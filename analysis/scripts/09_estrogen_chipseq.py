"""
Phase 3 post-hoc analysis: ESR1 ChIP-seq enrichment on CGRP-pathway promoters.

Motivated by the Phase 2 synthesis: Qu et al. (2025) recovered
HALLMARK_ESTROGEN_RESPONSE_EARLY in the first systematic POTS genomic study,
and Krause et al. (2021) showed estrogen withdrawal modulates trigeminovascular
CGRP release. If the estrogen->CGRP axis is mechanistically real and not just
parallel literature, the CGRP-pathway gene promoters should be preferentially
bound by the estrogen receptor (ESR1) compared to length-matched background.

Data:
- ReMap 2022 non-redundant ESR1 peaks (866,869 peaks merged across 448 datasets,
  24 biotypes; hg38)
- GENCODE v44 basic annotation (protein-coding gene TSS, hg38)

Definitions (pre-specified before opening the peak file):
- Promoter window: -2000 bp to +500 bp from the most-upstream TSS on the
  gene's strand (the conventional Cistrome/ENCODE promoter window).
- A gene is "ESR1-bound" if any ReMap NR peak overlaps its promoter window.
- Background: 19,567 protein-coding genes used in Phase 2 (length-quintile
  features already curated).
- Statistical tests:
  (1) Two-sided Fisher exact for CGRP-extended (CALCA/CALCB/CALCRL/RAMP1/CRCP)
      vs background.
  (2) Length-quintile-matched permutation null (10,000 permutations,
      seed=20260427) for the same gene set -- consistent with the
      hypergeometric test in Phase 2 (script 05).
  (3) Same two tests on the migraine primary set (123 lead-SNP-nearest
      genes), as a specificity check: estrogen response should be specific
      to the CGRP pathway, not a generic property of migraine GWAS hits.

This is post-hoc relative to the locked Phase 2 pre-registration. We label it
as such; the seed and protocol are fixed before running.
"""
from __future__ import annotations

import bisect
import gzip
import json
import re
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw"
PROC = ROOT / "data" / "processed"
RES = ROOT / "results"
RES.mkdir(parents=True, exist_ok=True)

REMAP_BED = RAW / "remap" / "remap2022_ESR1_nr_macs2_hg38_v1_0.bed.gz"
GTF = RAW / "gencode" / "gencode.v44.basic.annotation.gtf.gz"

PROMOTER_UPSTREAM = 2000
PROMOTER_DOWNSTREAM = 500
N_PERM = 10_000
SEED = 20260427

# ReMap2022 ESR1 NR peak score = number of source datasets supporting the peak.
# Distribution (n=866,869): median=2, p75=6, p90=23, p99=173, max=536.
# 49% of peaks are singletons (score=1), most of which are noise.
# Headline threshold: score>=10 (peaks supported by >=10 of the 448 datasets;
# a regime where binding is reproducible across labs but not extreme).
# Sensitivity sweep run at each level in PEAK_SCORE_THRESHOLDS.
HEADLINE_THRESHOLD = 10
PEAK_SCORE_THRESHOLDS = [1, 5, 10, 30]


def parse_gtf_protein_coding_genes(gtf_path: Path) -> pd.DataFrame:
    """Return DataFrame with gene_name, chrom, strand, tss for each protein-coding gene.

    Uses the `gene` feature line directly; tss is start (1-based) for '+' strand
    and end for '-' strand. When the same gene_name appears multiple times
    (rare, e.g., readthrough loci or PAR genes), keep the longest entry.
    """
    rows = []
    name_re = re.compile(r'gene_name "([^"]+)"')
    type_re = re.compile(r'gene_type "([^"]+)"')
    with gzip.open(gtf_path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[2] != "gene":
                continue
            attrs = parts[8]
            t = type_re.search(attrs)
            if not t or t.group(1) != "protein_coding":
                continue
            n = name_re.search(attrs)
            if not n:
                continue
            chrom = parts[0]
            start = int(parts[3])  # 1-based inclusive
            end = int(parts[4])
            strand = parts[6]
            tss = start if strand == "+" else end
            rows.append((n.group(1), chrom, strand, tss, end - start + 1))
    df = pd.DataFrame(rows, columns=["gene", "chrom", "strand", "tss", "length"])
    df = df.sort_values("length", ascending=False).drop_duplicates("gene", keep="first")
    return df.reset_index(drop=True)


def load_peaks(bed_path: Path, min_score: int) -> dict[str, tuple[np.ndarray, np.ndarray]]:
    """Load ReMap NR peaks at score>=min_score, return {chrom: (starts, ends)} sorted."""
    by_chrom: dict[str, list[tuple[int, int]]] = defaultdict(list)
    with gzip.open(bed_path, "rt") as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 5:
                continue
            score = int(parts[4])
            if score < min_score:
                continue
            chrom = parts[0]
            start = int(parts[1])  # 0-based half-open
            end = int(parts[2])
            by_chrom[chrom].append((start, end))
    out: dict[str, tuple[np.ndarray, np.ndarray]] = {}
    for chrom, items in by_chrom.items():
        items.sort()
        starts = np.fromiter((s for s, _ in items), dtype=np.int64, count=len(items))
        ends = np.fromiter((e for _, e in items), dtype=np.int64, count=len(items))
        out[chrom] = (starts, ends)
    return out


def promoter_has_peak(
    chrom: str,
    strand: str,
    tss: int,
    peaks: dict[str, tuple[np.ndarray, np.ndarray]],
) -> bool:
    if strand == "+":
        win_start = max(0, tss - 1 - PROMOTER_UPSTREAM)  # 0-based
        win_end = tss - 1 + PROMOTER_DOWNSTREAM + 1  # half-open
    else:
        win_start = max(0, tss - 1 - PROMOTER_DOWNSTREAM)
        win_end = tss - 1 + PROMOTER_UPSTREAM + 1
    if chrom not in peaks:
        return False
    starts, ends = peaks[chrom]
    # Find any peak with start < win_end and end > win_start.
    # Binary search for first peak with start >= win_end; consider all earlier peaks.
    idx = bisect.bisect_left(starts, win_end)
    if idx == 0:
        return False
    return bool(np.any(ends[:idx] > win_start))


def annotate_genes(genes_df: pd.DataFrame, peaks: dict) -> pd.DataFrame:
    has_peak = np.zeros(len(genes_df), dtype=bool)
    for i, row in enumerate(genes_df.itertuples(index=False)):
        has_peak[i] = promoter_has_peak(row.chrom, row.strand, row.tss, peaks)
    out = genes_df.copy()
    out["esr1_promoter_peak"] = has_peak
    return out


def length_matched_permutation(
    background_with_peak: pd.DataFrame,
    test_genes: list[str],
    n_perm: int,
    rng: np.random.Generator,
) -> tuple[float, float, int]:
    """Permutation null matched on length quintile (column `length_quintile`).

    Returns (observed_count, empirical_p_one_sided_greater, expected_count_mean).
    """
    bg = background_with_peak
    test_set = set(test_genes)
    test_in_bg = bg[bg["gene"].isin(test_set)]
    obs_count = int(test_in_bg["esr1_promoter_peak"].sum())
    n_test = len(test_in_bg)
    if n_test == 0:
        return 0, 1.0, 0.0

    quintile_to_indices: dict[int, np.ndarray] = {}
    for q, sub in bg.groupby("length_quintile"):
        quintile_to_indices[int(q)] = sub.index.to_numpy()
    test_quintiles = test_in_bg["length_quintile"].astype(int).to_numpy()

    null_counts = np.empty(n_perm, dtype=np.int64)
    has_peak_arr = bg["esr1_promoter_peak"].to_numpy()
    for i in range(n_perm):
        sampled = []
        for q in test_quintiles:
            pool = quintile_to_indices[q]
            sampled.append(pool[rng.integers(0, len(pool))])
        null_counts[i] = int(has_peak_arr[np.array(sampled)].sum())
    p_greater = float(((null_counts >= obs_count).sum() + 1) / (n_perm + 1))
    return obs_count, p_greater, float(null_counts.mean())


def run_one_threshold(
    genes_df: pd.DataFrame,
    bg_features: pd.DataFrame,
    cgrp_ext: list[str],
    cgrp_strict: list[str],
    migraine: list[str],
    autonomic: list[str],
    threshold: int,
    rng: np.random.Generator,
) -> dict:
    print(f"\n  -- threshold score >= {threshold} --")
    peaks = load_peaks(REMAP_BED, min_score=threshold)
    n_peaks = sum(len(s) for s, _ in peaks.values())
    print(f"    peaks retained: {n_peaks:,}")
    annotated = annotate_genes(genes_df, peaks)
    bg = bg_features.merge(
        annotated[["gene", "esr1_promoter_peak", "chrom", "strand", "tss"]],
        on="gene", how="left",
    )
    bg["esr1_promoter_peak"] = bg["esr1_promoter_peak"].fillna(False).astype(bool)
    bg_bound = int(bg["esr1_promoter_peak"].sum())
    bg_total = len(bg)
    print(
        f"    background bound: {bg_bound:,}/{bg_total:,} = "
        f"{100*bg_bound/bg_total:.1f}%"
    )

    out: dict = {
        "threshold_score_min": threshold,
        "n_remap_peaks_retained": n_peaks,
        "n_background_with_esr1_peak": bg_bound,
        "background_esr1_bound_fraction": bg_bound / bg_total,
    }
    for name, gene_list in [
        ("cgrp_extended", cgrp_ext),
        ("cgrp_strict", cgrp_strict),
        ("migraine_primary", migraine),
        ("autonomic", autonomic),
    ]:
        in_bg = bg[bg["gene"].isin(gene_list)]
        n_test = len(in_bg)
        n_test_bound = int(in_bg["esr1_promoter_peak"].sum())
        n_other_bound = bg_bound - n_test_bound
        odds, p_two = fisher_exact(
            [[n_test_bound, n_test - n_test_bound],
             [n_other_bound, (bg_total - n_test) - n_other_bound]],
            alternative="two-sided",
        )
        _, p_greater = fisher_exact(
            [[n_test_bound, n_test - n_test_bound],
             [n_other_bound, (bg_total - n_test) - n_other_bound]],
            alternative="greater",
        )
        obs, perm_p, perm_mean = length_matched_permutation(
            bg, gene_list, N_PERM, rng
        )
        fold = (n_test_bound / n_test) / (bg_bound / bg_total) if n_test and bg_bound else float("nan")
        out[name] = {
            "n_in_set": len(gene_list),
            "n_in_background": n_test,
            "n_with_esr1_peak": n_test_bound,
            "fraction_bound": (n_test_bound / n_test) if n_test else None,
            "fold_enrichment_vs_background": fold,
            "fisher_odds_ratio": float(odds),
            "fisher_p_two_sided": float(p_two),
            "fisher_p_one_sided_greater": float(p_greater),
            "length_matched_permutation_observed": obs,
            "length_matched_permutation_expected_mean": perm_mean,
            "length_matched_permutation_p_one_sided_greater": perm_p,
        }
        print(
            f"    {name}: {n_test_bound}/{n_test} ({100*n_test_bound/max(n_test,1):.0f}%) "
            f"fold={fold:.2f} Fisher_p_two={p_two:.4g} "
            f"perm_p_one>={perm_p:.4g} (E[k]={perm_mean:.2f})"
        )

    if threshold == HEADLINE_THRESHOLD:
        cgrp_table = (
            annotated[annotated["gene"].isin(cgrp_ext)][
                ["gene", "chrom", "strand", "tss", "esr1_promoter_peak"]
            ].reset_index(drop=True)
        )
        out["per_gene_cgrp_extended"] = cgrp_table.to_dict(orient="records")
    return out


def main() -> None:
    print("[1/3] Parsing GENCODE v44 basic protein-coding genes...")
    genes_df = parse_gtf_protein_coding_genes(GTF)
    print(f"      {len(genes_df):,} unique protein-coding gene_names")

    print("[2/3] Loading background features...")
    bg_features = pd.read_csv(PROC / "background_features.tsv", sep="\t")
    print(f"      {len(bg_features):,} background genes")

    cgrp_ext = pd.read_csv(PROC / "cgrp_anchor_extended.tsv", sep="\t")["gene"].tolist()
    cgrp_strict = pd.read_csv(PROC / "cgrp_anchor_strict.tsv", sep="\t")["gene"].tolist()
    migraine = pd.read_csv(PROC / "migraine_set_primary.tsv", sep="\t")["gene"].tolist()
    autonomic = pd.read_csv(PROC / "autonomic_set.tsv", sep="\t")["gene"].tolist()

    rng = np.random.default_rng(SEED)

    print("[3/3] Running enrichment tests at sensitivity-sweep thresholds...")
    by_threshold = {}
    for thr in PEAK_SCORE_THRESHOLDS:
        by_threshold[str(thr)] = run_one_threshold(
            genes_df, bg_features, cgrp_ext, cgrp_strict, migraine, autonomic, thr, rng
        )

    results = {
        "data": {
            "remap_bed": str(REMAP_BED.relative_to(ROOT)),
            "gencode_gtf": str(GTF.relative_to(ROOT)),
            "promoter_window_upstream_bp": PROMOTER_UPSTREAM,
            "promoter_window_downstream_bp": PROMOTER_DOWNSTREAM,
            "n_protein_coding_genes_gencode": len(genes_df),
            "n_background_genes": len(bg_features),
            "headline_threshold": HEADLINE_THRESHOLD,
            "peak_score_threshold_sweep": PEAK_SCORE_THRESHOLDS,
            "seed": SEED,
            "n_permutations": N_PERM,
        },
        "by_threshold": by_threshold,
    }
    out_path = RES / "estrogen_chipseq_results.json"
    out_path.write_text(json.dumps(results, indent=2))
    print(f"\nWrote {out_path.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
