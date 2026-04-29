# Pre-registration: gnomAD constraint test for the CGRP ligand-receptor asymmetry

## Locked before downloading or opening any constraint file. Decision rule and
## seed are fixed before this document is committed; the data fetcher script
## records the SHA-256 of the downloaded TSV against this spec.

## Question
The pre-registered Wright/PRF Monte Carlo (`analysis/simulations/SIM_SPEC.md`)
showed that two biological scenarios produce the observed CGRP ligand-receptor
asymmetry pattern equally well:

  M1: purifying selection on the CGRP receptor side (CALCRL/RAMP1/CRCP)
  M2: intrinsically small per-variant effects on migraine for the receptor side

Both clear the null (M3) at ~8x. The gnomAD population-genetic-constraint
metrics directly probe scenario M1: under purifying selection, receptor-side
genes should be depleted of loss-of-function variants in the population
relative to expectation. M2 makes no such prediction; under M2, receptor genes
should be at baseline tolerance to LoF.

## Data
- Source: gnomAD v4.1 per-gene constraint table.
- URL (canonical): https://gnomad.broadinstitute.org/downloads
  Specific file (per gnomAD docs as of 2026-04): the constraint table accompanying
  release v4.1, distributed at
  https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/constraint/gnomad.v4.1.constraint_metrics.tsv
- The fetcher script will record the SHA-256 of the downloaded TSV in
  `analysis/data/raw/MANIFEST.sha256` so the analysis is replicable from a
  hash-pinned download.

## Metrics
Primary: LOEUF (loss-of-function observed/expected upper-bound fraction).
Lower LOEUF = more constrained. The standard "constrained" tier is LOEUF
< 0.35 (gnomAD docs) or roughly the lowest decile genome-wide.

Secondary (reported but not primary decision metric):
- pLI (probability of LoF intolerance, legacy gnomAD v2 metric retained
  in v4.1 for continuity)
- mis_z (missense Z-score; positive = depleted of missense)
- oe_lof (point estimate of observed/expected LoF; LOEUF is its upper bound)

## Test genes and groups
- Ligand group: CALCA, CALCB
- Receptor group: CALCRL, RAMP1, CRCP

## Background
The 19,567 protein-coding HGNC symbols already curated in
`analysis/data/processed/background_features.tsv`. Genes for which gnomAD v4.1
does not report constraint statistics are excluded (not penalized).

## Decision rule (pre-registered; locked before any data is opened)

**Primary decision metric.** Compare the receptor-group mean LOEUF percentile
against the ligand-group mean LOEUF percentile against the genome-wide
distribution of LOEUF percentiles in the protein-coding background.

  - **M1 supported**: mean(LOEUF percentile of CALCRL, RAMP1, CRCP) <
    25th percentile (i.e., the receptor side is in the most-constrained quartile
    of protein-coding genes) AND mean(LOEUF percentile of CALCA, CALCB) >
    50th percentile (ligand side is at or above background median).
  - **M1 partially supported**: receptor mean percentile < 50th (constrained
    above median) but does not hit the 25th-percentile threshold; ligand at any
    percentile.
  - **M2 supported (M1 falsified)**: mean(LOEUF percentile, receptor) >= 50th
    percentile (no constraint enrichment on receptor side).
  - **Ambiguous**: ligand and receptor both below 25th percentile (both
    constrained), or both above 50th (neither constrained), with no clear
    asymmetry.

**Permutation test for significance.** Under the null that the receptor group
is exchangeable with random size-3 draws of protein-coding genes from
background, what fraction of 10,000 random draws produces a mean LOEUF
percentile <= the observed receptor mean? Reported as `p_perm_receptor`.

A second permutation test for ligand: fraction of 10,000 random size-2 draws
producing mean LOEUF percentile >= the observed ligand mean. Reported as
`p_perm_ligand`.

## Random seed
Master seed: 20260429. All numpy generators derived deterministically.

## Output
- `analysis/constraint/results/gnomad_constraint.json`
  - Per-gene metrics for the 5 test genes
  - Background LOEUF percentile distribution summary
  - Permutation test results
  - Verdict (M1 supported / partial / falsified / ambiguous)
- `analysis/constraint/figures/gnomad_constraint.{pdf,png}`
  - Background histogram of LOEUF with the 5 CGRP genes overlaid

## What this test does NOT do
- Does not test M2 directly. A null receptor-constraint result is consistent
  with M2 but also with several other scenarios (e.g., receptor genes are
  essential at the cellular level but not at the organismal LoF-detection
  level). The test is a one-sided probe of M1; falsification of M1 does not
  prove M2.
- Does not address the lens-(iv) interactome-sparseness explanation directly;
  that is a separate concern about PPI-database completeness.

## Verification before commit
- This file is committed to the repository before
  `analysis/constraint/scripts/01_gnomad_constraint.py` writes any output.
- The data download script records the URL fetched and the SHA-256 hash.
- The analysis script reads the locked seed and decision rule from constants
  matching this spec.
