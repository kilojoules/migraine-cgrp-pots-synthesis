# CGRP ligand-receptor genetic asymmetry in migraine

**Bottom line.** Anti-CGRP-receptor drugs (erenumab, rimegepant, atogepant) are FDA-approved migraine therapies that reduce monthly migraine days by 3.5–4.7. Yet across **five independent analytical lenses on the largest published migraine GWAS, the genes encoding the CGRP receptor (`CALCRL`, `RAMP1`, `CRCP`) carry zero genetic-association evidence**. The ligand genes (`CALCA`, `CALCB`) carry signal in the canonical GWAS-hit tier. We falsify the obvious explanations — purifying selection, cardiovascular-mediated selection, and "this is what GWAS does to peptide-hormone pathways generally" — and conclude the asymmetry reflects *intrinsic effect-size architecture*: receptor-pathway variants exist at population frequency but contribute too little to migraine risk to surface in current biobanks. Pharmacology at therapeutic doses simply produces far larger effects than common variants ever do.

The paper is in [`paper/main.pdf`](paper/main.pdf), 13 pages. All analyses are pre-registered with locked seeds before data was opened.

---

## The juiciest findings

1. **The CGRP receptor heterodimer is genetically silent for migraine despite being a validated drug target.** OpenTargets Genetics returns `genetic_association = 0.753` for `CALCA` (in the range of canonical migraine GWAS hits *PRDM16* 0.798, *TRPM8* 0.816, *FHL5* 0.855), `0.450` for `CALCB` (graded — intermediate), and **no `genetic_association` datatype recorded at all** for `CALCRL`, `RAMP1`, `CRCP`. The non-zero overall scores those receptor genes do have come entirely from drug-target evidence.

2. **It's not purifying selection.** gnomAD v4.1 LOEUF (730K exomes) places `RAMP1` at the **90.7th percentile** and `CRCP` at the **80.5th percentile** of the protein-coding background — markedly LoF-tolerant (`RAMP1` pLI = 0.0004). Receptor-group mean LOEUF percentile is 66.6, permutation p = 0.83. Only `CALCRL` shows modest constraint (28.7th percentile).

3. **It's not cardiovascular-mediated selection either.** Across `~1.3 million` UK Biobank exome-equivalents (Genebass + Backman 2021 + AstraZeneca PheWAS aggregated via OpenTargets), gene-burden tests return **zero positive evidence rows** for any of the five CGRP genes against any of 15 phenotypes — including all four cardiovascular phenotypes (hypertension, CAD, heart failure, atrial fibrillation). Positive controls (LDLR/LDL at -log10 P=96, MC4R/obesity at 21, PCSK9/LDL at 16) verify the API works.

4. **It's not a generic peptide-hormone-GWAS feature.** Same five-lens consolidation applied to nine other peptide-hormone pathways: only **Kisspeptin/age at menarche reproduces the CGRP pattern**. GLP-1/T2D shows the *opposite* pattern (receptor-asymmetric — *GLP1R* drives the signal). Three pathways are balanced (GIP, Leptin, Insulin); four are both-quiet. Cross-pathway permutation p = 0.51.

5. **There is modest biophysical constraint at the receptor protein level.** `CALCRL` has only 56% of expected PolyPhen-damaging missense (`mis_pphen.oe = 0.556`); `CRCP` 0.615. `RAMP1` is unconstrained at every metric. The receptor proteins are partially constrained at specific functional residues, but the migraine-relevant common-variant signal does not concentrate there.

6. **A pre-registered forward-population-genetic Monte Carlo (Wright/PRF, 10K replicates × 3 models, msprime+demographic confirmation) cannot distinguish "purifying selection" from "intrinsically small effects" on the genetic data alone.** The two scenarios both produce the observed pattern at the same rate. The gnomAD constraint test and rare-variant burden tests are what break the tie — toward "intrinsically small effects."

---

## What this means for migraine genetics

Causal models of CGRP-pathway involvement in migraine genetics must accommodate the asymmetry. The receptor side carries small per-variant effects on disease risk that current biobank-scale GWAS cannot surface, even though pharmacological perturbation at therapeutic doses clearly modulates the disease. **Genetic invisibility ≠ mechanistic uninvolvement**, and pharmaceutical target validation should not over-rely on common-variant GWAS evidence for pathways that look like CGRP.

## What's in the repo

```
paper/                           # Manuscript + bibliography
  main.tex, main.pdf, refs.bib

analysis/
  PRE_REGISTRATION.md            # Original pre-reg (PPI proximity)
  scripts/                       # 11 numbered analysis scripts
  data/                          # raw downloads gitignored; SHA-256 in MANIFEST
  results/                       # JSON results per analysis
  figures/                       # PDFs + PNGs

  simulations/                   # Forward-population-genetic Monte Carlo
    SIM_SPEC.md
    scripts/{01_simulate, 02_figure, 03_sensitivity}_asymmetry.py
    results/asymmetry_simulation.json
    figures/asymmetry_simulation.pdf

  constraint/                    # gnomAD LOEUF test
    PRE_REGISTRATION_constraint.md
    scripts/, results/, figures/

  burden/                        # OpenTargets rare-variant burden
    PRE_REGISTRATION_burden.md
    scripts/, results/, figures/

  constraint_extended/           # gnomAD missense + damaging-missense
    PRE_REGISTRATION_constraint_extended.md
    scripts/, results/, figures/

  generalize/                    # 10 peptide-hormone pathways
    PRE_REGISTRATION_generalize.md
    scripts/, results/, figures/

  sim_msprime/                   # msprime confirmation of the PRF simulation
    PRE_REGISTRATION_msprime.md
    scripts/, results/, figures/
```

## Reproducibility

Every analysis section in the paper has a pre-registration document committed to this repository **before** the corresponding output was opened. Per-analysis random seeds are recorded in those documents (`20260426` PPI proximity, `20260428` PRF simulation, `20260429` constraint, `20260430` burden, `20260501` cross-pathway, `20260502` msprime, `20260503` missense). Raw downloads are not redistributed but are fetched by the scripts from canonical sources, with SHA-256 hashes in `analysis/data/raw/MANIFEST.sha256`.

To rebuild the paper:
```bash
cd paper
make            # builds main.pdf via latexmk -pdf
```

To re-run any analysis (each is self-contained):
```bash
cd analysis/<section>
python3 scripts/01_<analysis_name>.py
```

## History

A long-form 25-page manuscript with the original POTS-comorbidity framing is preserved at git tag `long-form-v1`. The current paper isolates the load-bearing asymmetry observation and its adjudication.
