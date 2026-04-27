# Phase 2 Pre-Registration

This document locks in the analysis specification **before** any data is opened. It exists so that downstream choices cannot be retroactively tuned to produce a desired result. Once committed to git, every deviation from this spec must be recorded in `DEVIATIONS.md` with a stated reason.

**Date committed**: 2026-04-26
**Repo state at commit**: pre-data-download (Hautakangas xlsx, BioGRID, IntAct downloaded but never opened)
**Headline question**: Do migraine GWAS-prioritized genes sit closer in the protein--protein interaction network to the CGRP receptor complex than degree-matched random gene sets, after controlling for matched-control diseases?

---

## 1. Migraine gene set construction

**Primary rule.** Use the FUMA / positional / eQTL-prioritized gene assignments reported in Hautakangas et al. 2022 (Nat Genet 54:152--160), preferring the FUMA-prioritized gene list when multiple columns disambiguate a locus. Specific source: Supplementary Tables 8--10 (positional, eQTL, chromatin-interaction mapping) of Supplementary Data file MOESM4.

If FUMA-prioritized list is not directly tabulated, fall back order:
1. Lead-SNP nearest gene (Supplementary Table 3).
2. eQTL-mapped gene (Supplementary Table 9).
3. POPS-prioritized gene if available.

Gene symbols mapped to current HGNC symbols via `mygene.info` lookup with timestamp recorded.

**Sensitivity analyses to report** (not primary):
- Lead-SNP nearest-gene only.
- Union of positional + eQTL + chromatin-interaction.

**No subtype-stratified analysis in primary readout.** Subtype splits (with-aura vs without-aura) are reported only as secondary if the all-migraine result is significant.

## 2. CGRP receptor anchor

**Strict definition (primary):** `CALCRL + RAMP1`. This is the CGRP-selective heterodimer; RAMP2 and RAMP3 form heterodimers with CALCRL but bind adrenomedullin (AM) and adrenomedullin 2 (AM2/intermedin) respectively, not CGRP. Including RAMP2/3 would broaden the claim from ``CGRP signaling'' to ``calcitonin-family peptide signaling''---a meaningfully different hypothesis.

**Sensitivity (secondary):** `CALCA + CALCB + CALCRL + RAMP1 + CRCP` (ligands + receptor + RCP coreceptor). Reported as ``calcitonin-family-extended'' set.

**Excluded from anchor:** AMY (amylin) receptors, calcitonin receptor (CALCR), broader GPCR families.

## 3. Autonomic / adrenergic gene set

Set is the **union** of the following sources, deduplicated to HGNC symbols.

**GO terms (queried via `mygene.info` or g:Profiler API at execution time):**
- `sympathetic nervous system development`
- `parasympathetic nervous system development`
- `regulation of heart rate by chemical signal`
- `adrenergic receptor signaling pathway`
- `norepinephrine metabolic process`
- `baroreceptor response to increased systemic arterial blood pressure`
- `regulation of systemic arterial blood pressure by norepinephrine-epinephrine`

(GO IDs resolved at execution and recorded in `analysis/data/processed/autonomic_set_provenance.json`.)

**KEGG pathways:**
- `hsa04261` Adrenergic signaling in cardiomyocytes (precise)
- `hsa04270` Vascular smooth muscle contraction (broader; included)
- KEGG `hsa04080` Neuroactive ligand-receptor interaction is **excluded** as too broad.

**Manual seed genes (fixed list, not tunable):**
`ADRB1, ADRB2, ADRB3, ADRA1A, ADRA1B, ADRA1D, ADRA2A, ADRA2B, ADRA2C, SLC6A2, COMT, GNB3, DBH, TH, PNMT, CHRM2, CHRM3, CHRNB1, CHRNA3, NPY, NPY1R, NPY2R`

Provenance for every gene---which source(s) added it---recorded in TSV.

## 4. Protein--protein interaction network

**Primary network**: union of **BioGRID Homo sapiens physical interactions** and **IntAct human binary interactions** (PSI-MI annotated). Edges deduplicated; nodes mapped to HGNC symbols. Self-loops removed.

**Excluded** from primary: STRING text-mining edges, predicted/inferred edges. STRING high-confidence physical (score $\geq$ 700, text-mining channel disabled) is reported as a sensitivity overlay only.

**Versions and access dates** recorded in `analysis/data/processed/ppi_provenance.json`.

## 5. Network proximity statistic

**Primary statistic**: Guney et al. 2016 ``closest'' distance
$$d_c(S, T) = \frac{1}{|S|}\sum_{s \in S} \min_{t \in T} d(s, t)$$
where $d(s,t)$ is unweighted shortest-path length in the PPI graph, $S$ = migraine gene set restricted to nodes present in the network, and $T$ = CGRP anchor (CALCRL + RAMP1).

**Null model**: degree-binned random gene-set draws (Guney 2016 protocol). For each gene in $S$ and $T$, draw a replacement from the same degree bin (log-spaced bins). $1{,}000$ permutations, fixed RNG seed `20260426`.

**Reported statistics**: observed $d_c$, mean and SD of null distribution, $z$-score, two-sided empirical $p$-value.

**Secondary statistic**: Menche et al. 2015 separation $S_{AB}$.

## 6. Specificity controls

Compute $d_c$ from the migraine gene set to anchor sets for each of the following control diseases. Migraine--CGRP proximity is interpretable as ``CGRP-specific'' only if migraine--CGRP $z < z_{\text{control}}$ for **all four** controls.

**Control disease anchors (from DisGeNET or OpenTargets, accessed at execution date):**
- **Schizophrenia** (psychiatric, vascular-adjacent only via lifestyle).
- **Psoriasis** (immune; no autonomic/CGRP rationale).
- **Essential hypertension** (vascular; could legitimately overlap CGRP biology---a hard control).
- **Type 2 diabetes** (metabolic + vascular; broad; soft control).

Each control gene set is constructed identically to the migraine set rule (top-N GWAS-prioritized genes, $N$ = size of migraine set $\pm 10$).

## 7. Gene-set overlap test (focused autonomic enrichment)

**Test**: hypergeometric over-representation of autonomic gene set in migraine GWAS gene set.

**Null model**: 1{,}000 permutations sampling random gene sets of size $|S_{\text{migraine}}|$ from a background of all protein-coding genes, **matched on**:
- Gene length (bp), 5 bins by quintile.
- Whole-blood mRNA expression decile (GTEx v8 median TPM).

(Per Goeman \& B\"uhlmann 2007: matched-background permutation is the honest null.)

**Reported**: observed overlap $k$, expected from null (mean and 95\% CI), enrichment fold $k / E[k]$, two-sided empirical $p$-value.

## 8. Pathway enrichment

**Primary**: g:Profiler API (`gProfilerAPI` v1.0+ Python client), `g:GOSt` over GO\_BP, GO\_MF, KEGG, Reactome, with `g:SCS` correction. Background: all protein-coding HGNC genes.

**Reported**: top 20 enriched terms with adjusted $p$-values; explicit grep for autonomic/adrenergic/CGRP-relevant terms.

## 9. GTEx tissue expression

**Tissues** (GTEx v8; the trigeminal ganglion and sympathetic ganglion are not in v8, so we use the closest available proxies):
- `Brain - Spinal cord (cervical c-1)`
- `Brain - Hypothalamus`
- `Nerve - Tibial`
- `Artery - Tibial`
- `Artery - Aorta`
- `Artery - Coronary`
- `Whole Blood` (control / baseline)

**Output**: heatmap of median TPM for migraine genes $\cup$ CGRP anchor across these tissues; report whether migraine gene expression is concentrated in nerve/arterial tissues vs. whole-blood baseline.

## 10. Decision rule for the headline

| Outcome | Headline framing |
|---|---|
| migraine--CGRP $z < -2$ AND all four controls fail to separate | ``Migraine genetic signal sits closer to the CGRP receptor complex in the human PPI network than degree-matched chance, with specificity to CGRP signaling over matched control diseases.'' |
| migraine--CGRP $z < -2$ BUT one or more controls also separate | ``Consistent with proximity to CGRP signaling, but the effect is not specific---broader cardiovascular/neurovascular biology cannot be ruled out.'' |
| migraine--CGRP $-2 \leq z \leq 2$ | ``The proposed network-proximity test does not detect significant clustering of migraine GWAS genes near the CGRP receptor complex; the autonomic--CGRP triangulation hypothesis is not supported at the resolution available.'' |
| migraine--CGRP $z > 2$ (separated) | ``Surprising negative finding: migraine GWAS-prioritized genes are further from the CGRP receptor complex than chance. Discuss interactome incompleteness and ligand--receptor asymmetry.'' |

**All four outcomes are publishable.** This is committed in advance.

## 11. Reproducibility

- All scripts in `analysis/scripts/`, runnable as `make all` from `analysis/`.
- All randomness seeded with `20260426`.
- All data downloads dated and SHA-256-hashed in `analysis/data/raw/MANIFEST.sha256`.
- All software pinned in `analysis/requirements.txt` with exact versions.
- Final figures / tables in `analysis/figures/` and `analysis/results/`, regenerable from raw data + code.

## 12. Out of scope

- Fine-mapping or causal-variant analysis at the SNP level.
- Mendelian randomization.
- Trans-ancestry analysis (Hautakangas European-only is the substrate).
- Anything requiring individual-level genotype data.
- Wet-lab follow-up.

---

**Signed off (in commit history) before opening MOESM4 or constructing any analysis.**
