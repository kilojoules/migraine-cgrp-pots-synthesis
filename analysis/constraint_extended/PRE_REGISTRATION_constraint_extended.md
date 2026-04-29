# Pre-registration: missense and missense-deleterious constraint test for the CGRP asymmetry

## Locked before any output is opened. Decision rule, metric definitions,
## and seed are fixed before this document is committed.

## Question
Extension A tested scenario M1 (purifying selection on receptor genes)
via gnomAD v4.1 LOEUF and falsified it: RAMP1 90.7th LOEUF percentile,
CRCP 80.5th, mean receptor 66.6 (perm p=0.83). The receptor heterodimer
is highly LoF-tolerant.

LoF intolerance is one mode of purifying selection but not the only
one. A gene can be LoF-tolerant yet missense-constrained: small-effect
missense variants disrupting protein function may be selected against
even when complete loss-of-function is tolerated (e.g., where the gene
acts in a redundant pathway and only specific functional domains are
under strong selection). The advisor's recommended next test of M1
specifically called for selection signatures beyond LoF: SDS, iHS, or
the equivalent.

We do a tractable substitute: probe missense constraint metrics from the
same gnomAD v4.1 release we already have, extending beyond LoF to:
  - mis.z_score (missense Z-score; positive = depleted of missense)
  - mis.oe (missense observed/expected; lower = more constrained)
  - mis_pphen.oe (PolyPhen-2-damaging missense observed/expected)
  - syn.z_score (synonymous Z-score; should be near 0 as a sanity check)

If receptor genes show elevated missense Z that LoF tolerance missed,
scenario M1 partially returns at the missense layer. If both LoF and
missense are tolerant, M1 is more strongly falsified.

## Data
Same input as Extension A: gnomAD v4.1 per-gene constraint table at
analysis/data/raw/gnomad/gnomad.v4.1.constraint_metrics.tsv (SHA-256
68d8abdb7fc48f570869b02dfaa74b9fecaece7fcc5f301ddca40ec1ce12da00).

## Test genes
Same as Extension A: ligand={CALCA, CALCB}, receptor={CALCRL, RAMP1,
CRCP}, background=19,567 protein-coding HGNC symbols from
analysis/data/processed/background_features.tsv.

## Metrics
Primary: **mis.z_score** percentile within background. Higher = more
missense-constrained (depleted of missense variants).

Secondary (reported but not primary):
  - mis.oe percentile (lower = more constrained)
  - mis_pphen.oe percentile (lower = more constrained for PolyPhen-2-
    damaging missense; sharper functional signal)

Sanity-check control:
  - syn.z_score percentile (should be near 50 across all genes; if
    receptor genes show extreme synonymous Z, the missense signal could
    be technical and we should reject the analysis)

## Decision rule (pre-registered before any data is opened)

Primary missense test:
  - **M1-missense supported**: mean receptor mis.z_score percentile > 75
    (top quartile of missense constraint) AND ligand mis.z_score
    percentile <= 50, with permutation p < 0.05 over 10,000 random
    size-3 draws.
  - **M1-missense partially supported**: receptor percentile > 50 but
    not above the 75th-percentile threshold.
  - **M1-missense falsified**: receptor mean mis.z_score percentile
    <= 50 (no missense constraint enrichment), reinforcing the
    Extension A conclusion that the receptor heterodimer is not under
    classical purifying selection.
  - **Ambiguous**: ligand and receptor both above 75 (both
    missense-constrained), with no clear asymmetry direction.

Sanity-check rule: if any of the 5 test genes has |syn.z_score
percentile - 50| > 25 (i.e., synonymous Z extreme), flag the result
as potentially technical-artifact-affected.

## Random seed
Master seed: 20260503. Used for the receptor-group permutation test.

## Output
- analysis/constraint_extended/results/missense_constraint.json
- analysis/constraint_extended/figures/missense_constraint.{pdf,png}

## What this test does NOT do
- Does not run a haplotype-based recent-selection scan (iHS, SDS, PBS).
  Those would require phased haplotype VCFs and pre-computed scan
  databases that are not readily accessible at adequate resolution.
  This test is a tractable substitute that probes a different mode of
  purifying selection from the same gnomAD data layer.
- Does not address recent positive selection on common variants
  specifically; gnomAD's metrics summarize cumulative selection over
  longer time scales.
- Does not adjudicate the core M1/M2 question if the result is
  partially supported. A clean falsification or clean support either
  reinforces or rescues M1.

## Verification before commit
- This file is committed before any output JSON is opened.
- The script reads constants matching this spec.
