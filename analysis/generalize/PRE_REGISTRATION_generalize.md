# Pre-registration: ligand-receptor asymmetry across 10 peptide-hormone pathways

## Locked before any output is opened. The pathway list, disease IDs,
## asymmetry metric, decision rule, and seed are fixed before the query
## script runs.

## Question
Extensions A and B (this project) converged on M2 as the explanation for
the CGRP ligand-receptor genetic asymmetry: receptor-pathway genes have
intrinsically small per-variant effects on disease risk at population-
frequency variation levels, even though pharmacological perturbation
works. The natural follow-up question: is this pattern (ligand carries
GWAS signal, receptor does not) specific to CGRP, or a structural feature
of common-variant genetics in peptide-hormone pathways more generally?

We test this by running the same OpenTargets L2G consolidation used in
Section 5.6 of the short note across nine other peptide-hormone pathways
with clean ligand-receptor demarcations and well-powered disease GWAS,
plus the CGRP pathway itself as the positive-control reference. If the
asymmetry pattern recurs in most pathways, the result becomes a methods-
level finding about ligand-receptor asymmetry in GWAS architecture. If
it does not, CGRP is a specific case worth understanding individually.

## Pathways (n=10)
Each entry: pathway name, ligand gene(s), receptor gene(s), Ensembl IDs,
canonical disease GWAS in Open Targets ontology.

  1. CGRP    : ligand={CALCA(ENSG00000110680), CALCB(ENSG00000175868)}
               receptor={CALCRL(ENSG00000064989), RAMP1(ENSG00000132329),
                         CRCP(ENSG00000241258)}
               disease=migraine (MONDO_0005277)

  2. GLP-1   : ligand={GCG(ENSG00000115263)}
               receptor={GLP1R(ENSG00000112164)}
               disease=type 2 diabetes (MONDO_0005148)

  3. GIP     : ligand={GIP(ENSG00000159224)}
               receptor={GIPR(ENSG00000010310)}
               disease=type 2 diabetes (MONDO_0005148)

  4. Kisspeptin: ligand={KISS1(ENSG00000170498)}
                 receptor={KISS1R(ENSG00000116014)}
                 disease=age at menarche (EFO_0004703)

  5. Leptin  : ligand={LEP(ENSG00000174697)}
               receptor={LEPR(ENSG00000116678)}
               disease=body mass index (EFO_0004340) /
                       obesity (EFO_0001073)

  6. Oxytocin: ligand={OXT(ENSG00000101405)}
               receptor={OXTR(ENSG00000180914)}
               disease=major depressive disorder (MONDO_0002009)

  7. Insulin : ligand={INS(ENSG00000254647)}
               receptor={INSR(ENSG00000171105)}
               disease=type 2 diabetes (MONDO_0005148)

  8. Adiponectin: ligand={ADIPOQ(ENSG00000181092)}
                  receptor={ADIPOR1(ENSG00000159346),
                            ADIPOR2(ENSG00000006831)}
                  disease=type 2 diabetes (MONDO_0005148)

  9. Ghrelin : ligand={GHRL(ENSG00000157017)}
               receptor={GHSR(ENSG00000121853)}
               disease=body mass index (EFO_0004340)

 10. Substance P (NK1): ligand={TAC1(ENSG00000006128)}
                        receptor={TACR1(ENSG00000115353)}
                        disease=major depressive disorder (MONDO_0002009)

## Lens
For each pathway, query OpenTargets Genetics for the ligand and receptor
genes against the chosen disease, extract:
- `overall_score`
- `genetic_association` datatype score (this is the L2G aggregate;
  absent value treated as zero)

Per-gene metric: `genetic_association` score, with absent treated as 0.

## Asymmetry score per pathway
$A_{\text{path}} = \mathrm{median}(\text{L2G}_{\text{ligand genes}}) -
                    \mathrm{median}(\text{L2G}_{\text{receptor genes}})$

Positive A means ligand-asymmetric; negative means receptor-asymmetric;
near 0 means balanced.

## Decision rule (pre-registered)
Per-pathway classification:
- ligand-asymmetric: A > +0.2
- receptor-asymmetric: A < -0.2
- balanced: |A| <= 0.2 AND max(ligand or receptor median L2G) >= 0.1
- both-quiet: max(ligand or receptor median L2G) < 0.1 (e.g., disease
  GWAS underpowered, no signal anywhere in the pathway)

Cross-pathway test: under the null that ligand and receptor draws are
exchangeable per pathway, the expected count of ligand-asymmetric
pathways out of 10 is binomial(p=0.5, n=evaluable pathways). Permutation
test: per pathway, randomly relabel ligand and receptor groups; compute
a re-shuffled $A^{*}_{\text{path}}$; report the fraction of 10,000
relabel-permutations producing at least as many ligand-asymmetric
pathways as observed.

Verdict:
- **General asymmetry**: at least 7 of 10 evaluable pathways are
  ligand-asymmetric AND permutation p < 0.05.
- **CGRP-specific**: 0-3 of 10 pathways are ligand-asymmetric (CGRP
  alone or with at most 2 others), permutation p > 0.20.
- **Inconclusive**: anywhere in between, or too many both-quiet
  pathways to evaluate.

## Random seed
20260501 for the relabel-permutation test.

## Output
- `analysis/generalize/results/pathway_lens.json` — per-pathway,
  per-gene OpenTargets evidence
- `analysis/generalize/results/asymmetry_summary.json` — per-pathway
  asymmetry score, classification, permutation p
- `analysis/generalize/figures/pathway_asymmetry_heatmap.{pdf,png}` —
  pathway-by-gene L2G score heatmap with ligand/receptor grouping

## What this test does NOT do
- Does not run a fresh GWAS or fine-mapping per pathway.
- Does not handle the case of multiple disease GWAS per pathway
  (we pick one canonical disease per pathway, locked above).
- Does not address whether the asymmetry direction (ligand-positive)
  has a mechanistic explanation; that is a follow-on question.
- Does not (in this version) re-run the M1/M2/M3 simulation per
  pathway. The simulation per-pathway is queued for a follow-up if
  the consolidation finds general asymmetry.

## Verification before commit
- This file is committed before the query script runs.
- The query script reads gene/disease IDs from constants matching
  this spec.
