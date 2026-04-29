# Pre-registration: msprime confirmation of the asymmetry simulation

## Locked before opening any msprime simulation output. Confirmation
## of the Phase 4 PRF/Wright Monte Carlo under tree-sequence-level
## demographic and recombination realism.

## Question
The Phase 4 simulation (`analysis/simulations/SIM_SPEC.md`) used a
Wright / Poisson Random Field equilibrium Monte Carlo with within-gene
sites treated as independent and a constant effective population size
($N_e = 10{,}000$). Extensions A and B then resolved the M1/M2 confound
empirically (M1 falsified by gnomAD LoF constraint and rare-variant
burden), so the residual question is methodological: would the PRF
ranking M1 > M2 > M3 (or, equivalently, the headline hit-pattern rates)
have come out the same way under a more realistic forward simulation
that incorporates linkage disequilibrium, recombination, and an
out-of-Africa demographic expansion?

We test this by replicating the same three-model comparison using
msprime-generated tree sequences with a realistic European demographic
model and selection applied as a post-hoc allele-frequency filter at
the pre-registered selection coefficients. The msprime side adds:
- LD between sites within a gene (recombination at $1.2 \times 10^{-8}$
  per bp per generation)
- Realistic demographic expansion (Tennessen et al.\ 2012-style three-
  epoch model, parameterized via msprime's StandardCoalescentModel)
- Tree-sequence-level mutation placement

The msprime side does not add within-coalescent selection; selection is
applied post-hoc by filtering allele frequencies through the same
selection-drift equilibrium acceptance probability as the PRF model.
This is a deliberate simplification: the goal is to confirm that the
simpler PRF model is not misleading us, not to introduce a new
mechanism.

## Demographic model
Three-epoch out-of-Africa with European expansion:
- Ancestral diploid Ne = 7,310 (50,000 generations ago to start)
- Expansion 1: at 5,920 generations ago, Ne grows to 14,474 (out of
  Africa, into Eurasia)
- Expansion 2: at 2,040 generations ago, Ne grows to 100,000 (recent
  European expansion)
- Sample 2,000 haplotypes (1,000 diploid individuals) at present

(These are the standard Tennessen 2012 European parameters used as
reference in stdpopsim.)

## Per-gene simulation
For each of CALCA (384 bp), CALCB (381 bp), CALCRL (1380 bp), RAMP1
(444 bp), CRCP (444 bp):
1. Simulate ancestry tree sequence with msprime using the demographic
   model and recombination rate 1.2e-8 per bp per generation.
2. Add neutral mutations at rate 1.5e-8 per bp per generation
   (matching PRF run for direct comparability).
3. Extract present-day variant table.

## Selection filter (post-hoc, applied to msprime output)
For each variant in the simulated tree sequence:
- Draw a per-variant selection coefficient $s \sim \Gamma(\text{shape}=0.2,
  \text{scale}=\bar{s}_{\text{gene-model}} / 0.2)$ with the gene/model-
  specific mean $\bar{s}$ matching SIM_SPEC.md.
- Compute the steady-state acceptance probability of observing the
  variant at its msprime-simulated frequency $p_{\text{obs}}$ under
  Wright equilibrium with selection $s$:
    $\pi_{\text{accept}} = f_{\text{Wright}}(p_{\text{obs}} \mid s) /
                            f_{\text{Wright}}(p_{\text{obs}} \mid s = 0)$
  This reweights the neutral allele-frequency distribution to the
  selected one without resimulating; bounded by 1 for tractability.
- Accept the variant with probability $\pi_{\text{accept}}$.

For accepted variants, assign migraine effect $\beta$ under the same
M1/M2/M3 rules as the PRF run.

Common variants ($\mathrm{MAF} \geq 0.01$) are passed through the same
GWAS test statistic as PRF.

## Models and parameters
Same as `analysis/simulations/SIM_SPEC.md`:
- M1: $\bar{s}_{\text{lig}} = 10^{-4}$, $\bar{s}_{\text{recep}} = 5 \times 10^{-3}$,
       $\beta = 4 \cdot \mathrm{sign}() \cdot \sqrt{|s|}$
- M2: $\bar{s} = 5 \times 10^{-4}$ both sides; $\beta \sim \mathcal{N}(0,
       \sigma^2)$ with $\sigma_{\text{lig}} = 0.20$, $\sigma_{\text{recep}} = 0.05$
- M3: $\bar{s} = 5 \times 10^{-4}$ both sides; $\beta \sim \mathcal{N}(0,
       0.10^2)$ on both

## Replicates
1,000 per model (msprime is slower than PRF; 1k is enough to confirm
the qualitative ranking).

## Decision rule (pre-registered)
Confirmation criterion: the ordering of hit_pattern_rate across the
three models in msprime is the same as in the PRF run. That is:
  msprime: M1 > M2 > M3 (each separated by at least 2x in rate)
We do not require absolute rate equality; the simpler PRF model is
expected to produce slightly different absolute rates due to the
absence of LD attenuation. We do require ranking preservation.

If ranking preserved: PRF result is robust to demographic and LD
realism; the Section~5 simulation conclusion stands.
If ranking flipped: the PRF result is misleading and we have to
revisit; this would be a substantive finding.

## Random seed
Master seed: 20260502.

## Output
- `analysis/sim_msprime/results/msprime_simulation.json`
- `analysis/sim_msprime/figures/msprime_vs_prf.{pdf,png}`

## What this test does NOT do
- Does not add selection within the coalescent (uses post-hoc
  acceptance filtering instead). A future extension via SLiM-coupled
  msprime could test selection-during-coalescent.
- Does not change the biological conclusions (those are settled by
  Extensions A and B).
- Does not test the cross-pathway result (Extension C); those queries
  are deterministic against external resources.

## Verification before commit
- This file is committed before any msprime output is written.
- Selection coefficients and effect-size distributions match SIM_SPEC.md
  exactly.
