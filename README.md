# The CGRP migraine paradox: drugs work, the genetics is silent

**Bottom line.** Three FDA-approved migraine drugs (erenumab, rimegepant, atogepant) work by blocking a specific protein on the surface of nerve cells called the CGRP receptor. They cut monthly migraine days by 3.5–4.7. But when you look at the largest human genetic study of migraine — over 100,000 migraine sufferers and 770,000 controls — **the genes that build that receptor (`CALCRL`, `RAMP1`, `CRCP`) show no connection to migraine.** The genes that make the molecule the drugs intercept (CGRP itself, encoded by `CALCA` and `CALCB`) do show up. We checked five different ways of looking at this and the gap was the same every time. Then we ruled out the obvious explanations one by one. The conclusion: the receptor genes carry effects on migraine risk that are too small to detect at population scale, but big enough that drug-strength interference with the same proteins changes the disease.

The paper is in [`paper/main.pdf`](paper/main.pdf), 13 pages. Every analysis was written down in advance, sealed in a git commit, then run.

---

## What we found

1. **The drug target is genetically invisible.** A standard resource (OpenTargets) collects evidence linking each gene to each disease. `CALCA` scores 0.75 on the genetics axis for migraine — in the same range as well-known migraine genes like `TRPM8` (0.82) and `FHL5` (0.86). `CALCB` is moderate at 0.45. The three receptor genes — the actual targets of the marketed drugs — score zero. They have no recorded genetic evidence linking them to migraine at all. Their nonzero overall scores come entirely from the fact that drugs against them work.

2. **It is not because the receptor is too important to mutate.** If breaking the CGRP receptor were lethal or strongly disadvantageous, disabling mutations would be filtered out by natural selection and we would see far fewer than expected in healthy people. The gnomAD project sequenced 730,000 human exomes and counted exactly that. `RAMP1` ranks at the 90th percentile of mutation-tolerance — meaning 90% of all human genes are *less* tolerant of disabling mutations than `RAMP1` is. `CRCP` is at the 80th. The receptor complex is not under strong selection at the loss-of-function level.

3. **It is not because rare mutations are carrying the signal.** Three separate large studies (Genebass, AstraZeneca, Regeneron) have tested whether rare disabling mutations in each gene associate with each of thousands of diseases, across roughly 1.3 million sequenced people. Across migraine, hypertension, heart disease, heart failure, atrial fibrillation, and ten other tested phenotypes, **none of the five CGRP-pathway genes have a single positive hit.** Known positive genes light up as expected in the same lookup (`LDLR` for cholesterol scores higher than 1 in 10⁹⁶, `MC4R` for obesity at 1 in 10²¹), so the lookup works. The CGRP receptor pathway just genuinely has no signal.

4. **It is not how peptide-hormone signaling looks in genetics generally.** We applied the same comparison to nine other signaling systems (insulin, leptin, GLP-1, kisspeptin, oxytocin, ghrelin, GIP, adiponectin, and substance P). Only one — kisspeptin and age at first menstruation — shows the same ligand-strong / receptor-silent pattern. The diabetes drug GLP-1 shows the *opposite*: the receptor (`GLP1R`) carries strong genetic signal and the ligand precursor does not. Three pathways are balanced (insulin, leptin, GIP). Across 10 pathways, no general pattern; the CGRP case is specific.

5. **The receptor protein has functional constraints — but the migraine-relevant variation isn't there.** Two of the three receptor genes (`CALCRL` and `CRCP`) have fewer damaging changes in their protein-coding sequence than expected by chance: about 56% and 62% of expectation respectively. So specific functional residues do matter. But the variants that associate with migraine don't sit on those constrained residues. `RAMP1` is unconstrained at every measure — its protein sequence tolerates damage freely.

6. **A computer simulation cannot tell selection from "small effects" from the genetics alone — but the gnomAD and rare-variant tests above can.** A pre-registered population-genetics simulation (10,000 replicates × 3 competing explanations, plus a more realistic version with population history confirmed in `msprime`) shows that "purifying selection on the receptor" and "intrinsically small effects on the receptor side" produce the observed pattern equally well. The genetics alone is not enough to distinguish them. The constraint and rare-variant lookups break the tie — toward "small effects."

## What this means

A gene that doesn't show up in human genetic studies of a disease isn't necessarily uninvolved in the disease. Pharmaceutical companies sometimes use "genetic evidence supports this target" as a green light for drug development. The CGRP story is a counter-example. Anti-CGRP-receptor drugs were validated by clinical trials and they work — even though five different ways of mining migraine genetics produce no signal at the receptor genes. The reason is mundane: drugs at therapeutic doses produce much larger effects than the small differences carried by common genetic variants. The genetics is silent because population-frequency variation has small effects on disease risk for that side of the pathway, not because the proteins are uninvolved.

## What's in the repo

```
paper/                      Manuscript and bibliography
  main.tex, main.pdf, refs.bib

analysis/
  PRE_REGISTRATION.md       The original sealed prediction (PPI proximity)
  scripts/                  11 numbered analysis scripts
  data/                     Raw downloads ignored by git; SHA-256 hashes recorded
  results/                  JSON output of each analysis
  figures/                  PDF and PNG figures

  simulations/              Population-genetics Monte Carlo simulation
  constraint/               Mutation-tolerance test (gnomAD)
  burden/                   Rare-variant lookup (OpenTargets)
  constraint_extended/      Damaging-missense follow-up (gnomAD)
  generalize/               Same comparison applied to 10 signaling pathways
  sim_msprime/              Simulation re-run under realistic demography
```

Each subdirectory has its own sealed-prediction file (`PRE_REGISTRATION_*.md`) committed before the analysis ran, plus the script that ran it, the JSON of results, and the figure.

## Reproducibility

Every analysis section in the paper has a written-down prediction (decision rule, parameters, random seed) committed to this repository **before** the corresponding output was opened. Per-analysis seeds are recorded in those documents. Raw data files are not redistributed but are fetched by the scripts from their canonical sources, with SHA-256 hashes in `analysis/data/raw/MANIFEST.sha256`.

To rebuild the paper:
```bash
cd paper
make
```

To re-run any analysis (each is self-contained):
```bash
cd analysis/<section>
python3 scripts/01_<analysis_name>.py
```

## History

A longer 25-page manuscript with broader framing about a clinical condition called POTS (Postural Orthostatic Tachycardia Syndrome) is preserved at git tag `long-form-v1`. The current paper isolates the CGRP-genetics observation that survived everything we threw at it.
