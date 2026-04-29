# Simulation specification: CGRP ligand-receptor asymmetry under three competing models

## Locked before running. Seed and protocol fixed before opening the simulation script's output.

## Question
Three explanations remain compatible with the observed CGRP ligand-receptor
asymmetry in migraine genetics (CALCA L2G genetic_association = 0.753,
CALCB = 0.450, CALCRL/RAMP1/CRCP = 0):

  M1. Purifying selection on receptor genes: variants that affect CGRP
      receptor function (and therefore migraine risk) are also under
      strong purifying selection because the receptor is essential for
      cardiovascular regulation. Large-effect variants are removed before
      reaching common-variant frequency.

  M2. Allelic-series concentration on the ligand: receptor-side variants
      reach common frequency, but their per-variant effect on migraine
      risk is intrinsically smaller than ligand-side variants because the
      migraine-relevant biology is concentrated at the ligand.

  M3. Null: no asymmetry. Both sides equal in selection and effect size.

These cannot be separated by the genetic data alone; we use simulation to
ask which model produces the observed pattern most often.

## Model
Per-gene Poisson Random Field / Wright equilibrium Monte Carlo.

For each of the five genes (CALCA, CALCB, CALCRL, RAMP1, CRCP), in each
replicate:

  1. Draw N_seg ~ Poisson(rate_gene), where rate_gene is the expected
     number of segregating variants from a deleterious-DFE realization
     scaled by gene CDS length L_gene (in bp) and a per-bp mutation
     rate mu = 1.5e-8 per generation, in a population of effective
     size Ne = 10,000. rate_gene = theta_gene * H_(2*N_sample - 1)
     attenuated by the selection efficiency factor described below.

  2. For each segregating variant i:
     - Draw a fitness selection coefficient s_i from gamma DFE with shape
       0.2 and scale specific to the gene's selection regime in each model.
     - Draw a migraine effect size beta_i; coupling to s_i depends on
       the model (see below).
     - Draw a present-day minor allele frequency p_i from the
       selection-drift equilibrium distribution (Wright):
         f(p | s) proportional to (1 - exp(-2*Ne*s*(1-p))) /
                   (p * (1-p) * (1 - exp(-2*Ne*s)))
       sampled by inverse-CDF over a [1/(2Ne), 0.5] grid.

  3. Filter to common variants: MAF >= 0.01.

  4. Simulate a GWAS at sample size N_cases = 102,084, N_controls =
     771,257 (matching Hautakangas et al. 2022). Per-variant test
     statistic under additive model:
       SE_i = sqrt(1 / (2 * N_eff * MAF_i * (1 - MAF_i)))
       Z_i = beta_i / SE_i,  where N_eff = 4 / (1/N_cases + 1/N_controls)

  5. Per-gene summary: max_Z_gene = max_i |Z_i|; n_GWS_gene = count of
     |Z_i| > 5.452 (P < 5e-8 two-sided).

## Model parameterization

Common across all models:
  Ne = 10,000  (effective population size)
  mu = 1.5e-8  (per-bp per-generation mutation rate)
  N_cases = 102,084,  N_controls = 771,257
  Gene CDS lengths (from GENCODE v44 basic, MANE-select transcript or
  longest if MANE-select absent):
    CALCA  ~ 384 bp
    CALCB  ~ 381 bp
    CALCRL ~ 1380 bp
    RAMP1  ~ 444 bp
    CRCP   ~ 444 bp
  (Note: receptor genes are not systematically longer, so any "more
   variants on receptor side" effect from gene-length alone is small.)

M1 (purifying selection on receptor):
  Selection coefficients drawn from gamma DFE.
  - Ligand genes (CALCA, CALCB): mean s = 1e-4 (weakly deleterious;
    variants tolerated at common frequency)
  - Receptor genes (CALCRL, RAMP1, CRCP): mean s = 5e-3 (strongly
    deleterious; variants kept rare)
  Effect sizes: beta = alpha * sign() * sqrt(s), with alpha = 4 (so
  larger-s variants have larger effects, a standard DFE-effect
  coupling). Same alpha for both sides.

M2 (allelic-series concentration):
  Selection coefficients: same on both sides, mean s = 5e-4 (moderate
  selection). All five genes get the same DFE.
  Effect sizes: beta ~ N(0, sigma_gene^2), independent of s.
  - Ligand genes: sigma_lig = 0.20 (large per-variant effects)
  - Receptor genes: sigma_recep = 0.05 (small per-variant effects)
  alpha = 0 (no DFE-effect coupling).

M3 (null):
  Selection coefficients: same on both sides, mean s = 5e-4.
  Effect sizes: beta ~ N(0, sigma^2) with sigma = 0.10 for all five
  genes. alpha = 0.

## Replicates and decision rule
Number of replicates per model: R = 10,000.

For each replicate, compute the asymmetry metric:
  asym = max(max_Z_CALCA, max_Z_CALCB) - max(max_Z_CALCRL, max_Z_RAMP1, max_Z_CRCP)

Observed asymmetry in real data, in Z units:
  ligand_Z = sqrt(-2 * log(P_GWS_ligand))   [P at locus 75 reaches GWS, so |Z| >= 5.45]
  receptor_Z = sqrt(-2 * log(1e-4))         [Hautakangas reports receptor P > 1e-4]
  asym_obs = ligand_Z - receptor_Z ~ 5.45 - 4.26 ~ 1.2

But this lower bound is conservative; the real ligand peak is far above
P=5e-8 (locus 75 in Hautakangas at P ~ 1e-30+). Use as the per-replicate
test:
  hit_pattern = (max_Z_ligand > 5.45) AND (max_Z_receptor < 4.26)

Decision rule: report fraction of replicates per model satisfying
hit_pattern, plus the full distribution of asym across replicates.
Higher fraction means the model is more compatible with the observation.

Pre-registered comparison: rank M1, M2, M3 by hit_pattern fraction.

## Random seed
Master seed: 20260428. All numpy generators derived deterministically.

## Output
- analysis/simulations/results/asymmetry_simulation.json
  per-model summary stats and full per-replicate table
- analysis/simulations/figures/asymmetry_simulation.pdf
  multi-panel figure of per-gene Z distributions and hit-pattern fractions

## Out of scope
- Linkage disequilibrium between sites within a gene (treated as independent)
- Population structure (single panmictic population)
- Variable mutation rate across sites
- Demographic expansion (constant Ne)
These approximations are conservative for the asymmetry question because
they do not differentially affect the five genes; if anything they
homogenize the test across genes.
