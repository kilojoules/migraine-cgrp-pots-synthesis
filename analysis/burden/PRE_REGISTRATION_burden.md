# Pre-registration: Rare-variant burden test for the CGRP ligand-receptor asymmetry

## Locked before opening any burden output. Decision rule and seed fixed
## before this document is committed.

## Question
Extension A (gnomAD constraint) falsified scenario M1 for the receptor
heterodimer at the loss-of-function-detection level (RAMP1 LOEUF percentile
90.7, CRCP 80.5; receptor mean 66.6, perm p = 0.83). M1 is therefore not
the explanation for the migraine common-variant ligand-receptor asymmetry.

The remaining live scenario is M2: receptor-side variants reach common
frequency normally but contribute small per-variant effects on migraine
risk specifically. M2 makes a sharper testable prediction at the
rare-variant burden level: receptor-gene rare-variant burden should be
null for migraine across well-powered exome cohorts. A complementary
prediction is that receptor-gene rare-variant burden may be positive for
adjacent cardiovascular and metabolic phenotypes (where receptor
biology is more directly engaged), corroborating the "small migraine
effect, larger CV effect" reading.

## Data
- Source: OpenTargets Genetics Platform (v25.06 or current),
  GraphQL API at https://api.platform.opentargets.org/api/v4/graphql.
  OpenTargets ingests gene-burden test results from three primary sources:
  Genebass (Karczewski et al. 2022; 426,370 UK Biobank exomes), the
  AstraZeneca PheWAS Portal (~450,000 UK Biobank exomes), and the
  Regeneron Genetics Center collection (454,787 individuals). The
  platform standardizes statistical significance thresholds across these
  resources and harmonizes phenotypes to the Experimental Factor Ontology.
  Using OpenTargets here rather than direct Genebass/AZ access is a
  pragmatic choice: the same evidence is reachable through one query
  surface that we already use elsewhere in this project
  (`analysis/scripts/11_opentargets_l2g.py`).

## Test genes
- Ligand group: CALCA (ENSG00000110680), CALCB (ENSG00000175868)
- Receptor group: CALCRL (ENSG00000064989), RAMP1 (ENSG00000132329),
  CRCP (ENSG00000241258)

## Phenotypes
Primary disease panel (in approximate order of relevance):
- Migraine: MONDO_0005277 (or EFO_0003821 if Open Targets re-routes)
- Type 2 diabetes: MONDO_0005148 (CV/metabolic, control)
- Essential hypertension: MONDO_0001134 (CV)
- Coronary artery disease / ischemic heart disease: EFO_0001645 or
  MONDO_0005010 (CV)
- Heart failure: EFO_0003144 or MONDO_0005009 (CV)
- Schizophrenia: MONDO_0005090 (negative control; non-CV, non-metabolic)

Variant masks accepted from OpenTargets: any gene-burden evidence row
the platform ingested. We do not pre-restrict by mask because the three
upstream sources use different mask definitions; we record the mask per
row as reported.

## Decision rule (pre-registered before any query is run)

For each (gene, disease) pair, extract every evidence row in OpenTargets
where the data type is `genetic_association` and the data source is
`gene_burden`. Record P-value, beta/odds-ratio, statistical method, and
the upstream source.

Group-level test:
- For each disease, compute the median -log10(P) across the receptor
  group (CALCRL, RAMP1, CRCP) and the ligand group (CALCA, CALCB),
  using only burden evidence rows (NaN if no row exists for a
  gene-disease pair).

Verdict:
- **M2 confirmed (predominantly null receptor burden)**: median
  -log10(P) for the receptor group is below 2.0 for migraine AND
  below 2.0 for at least 2 of the 4 CV/metabolic phenotypes.
- **M1-like CV signal alongside null migraine**: median -log10(P) for
  the receptor group is above 3.0 (P < 1e-3) for at least one CV
  phenotype AND below 2.0 for migraine. (This would partially rescue
  the M1 reading: selection acting on receptor function for CV
  reasons even though gnomAD LoF-level constraint is null.)
- **Receptor-side migraine signal**: median -log10(P) for receptor
  group above 3.0 for migraine. Would be a striking surprise and
  invalidate the entire ``receptor side is invisible to migraine
  genetics'' framing.
- **Insufficient evidence**: fewer than 50% of (receptor gene, disease)
  pairs have any burden row in OpenTargets. Report as such.

## Random seed
Master seed: 20260430. Used only if any random sampling is required
(e.g., background gene draws); the primary query is deterministic.

## Output
- `analysis/burden/results/burden_opentargets.json`
  - Per (gene, disease, source) burden evidence rows
  - Group-level summary tables
  - Verdict per the decision rule above
- `analysis/burden/figures/burden_phenome_heatmap.{pdf,png}`
  - Heatmap of -log10(P) across genes (rows) x diseases (columns) x
    upstream source

## What this test does NOT do
- Does not run a fresh burden analysis. We reuse OpenTargets-ingested
  results. If a (gene, disease) pair has no burden evidence in OpenTargets,
  it does not necessarily mean no test was run upstream; it may mean the
  upstream test did not pass the source-specific significance threshold.
- Does not test rare-variant gain-of-function or structural variants.
- Does not address the lens-(iv) interactome-completeness explanation.

## Verification before commit
- This file is committed before any output JSON is written.
- The query script reads the locked seed and gene/disease lists from
  constants matching this spec; URLs are hard-coded.
