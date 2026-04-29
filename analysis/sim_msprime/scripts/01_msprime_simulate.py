"""
Extension D: msprime confirmation of the Phase 4 asymmetry simulation.

Replicates the M1/M2/M3 model comparison using msprime tree sequences
under a three-epoch out-of-Africa European demographic model with
recombination, then applies selection as a post-hoc Wright-equilibrium
acceptance filter.

Locked spec: analysis/sim_msprime/PRE_REGISTRATION_msprime.md.
Output: analysis/sim_msprime/results/msprime_simulation.json
"""
from __future__ import annotations

import json
import time
from pathlib import Path

import msprime
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "results"
RES.mkdir(parents=True, exist_ok=True)

# ---- Parameters matching SIM_SPEC.md ----
NE_FINAL = 100_000  # present-day European Ne (after expansion)
NE_OOA = 14_474     # post-out-of-Africa Ne
NE_ANC = 7_310      # ancestral
T_OOA = 5_920       # generations ago
T_EXP = 2_040       # generations ago
MU = 1.5e-8
RECOMB = 1.2e-8
N_DIPLOID_SAMPLES = 1_000  # 2,000 haplotypes
N_REPLICATES = 1_000
SEED = 20260502

GWS_THRESHOLD = 5.452
RECEPTOR_INVISIBLE_THRESHOLD = 4.26

N_CASES = 102_084
N_CONTROLS = 771_257
N_EFF = 4 / (1 / N_CASES + 1 / N_CONTROLS)

GENES = {
    "CALCA":  {"length_bp":  384, "side": "ligand"},
    "CALCB":  {"length_bp":  381, "side": "ligand"},
    "CALCRL": {"length_bp": 1380, "side": "receptor"},
    "RAMP1":  {"length_bp":  444, "side": "receptor"},
    "CRCP":   {"length_bp":  444, "side": "receptor"},
}

MODELS = {
    "M1_purifying_selection": {
        "ligand":   {"mean_s": 1e-4, "alpha": 4.0, "sigma_beta": None},
        "receptor": {"mean_s": 5e-3, "alpha": 4.0, "sigma_beta": None},
    },
    "M2_allelic_series": {
        "ligand":   {"mean_s": 5e-4, "alpha": 0.0, "sigma_beta": 0.20},
        "receptor": {"mean_s": 5e-4, "alpha": 0.0, "sigma_beta": 0.05},
    },
    "M3_null": {
        "ligand":   {"mean_s": 5e-4, "alpha": 0.0, "sigma_beta": 0.10},
        "receptor": {"mean_s": 5e-4, "alpha": 0.0, "sigma_beta": 0.10},
    },
}

DFE_SHAPE = 0.2

# ---- Demographic model (Tennessen-style three-epoch European) ----
def build_demography() -> msprime.Demography:
    d = msprime.Demography()
    d.add_population(name="EUR", initial_size=NE_FINAL)
    # Going backward in time:
    # At T_EXP, EUR shrinks to NE_OOA
    d.add_population_parameters_change(time=T_EXP, population="EUR",
                                         initial_size=NE_OOA)
    # At T_OOA, shrinks to NE_ANC
    d.add_population_parameters_change(time=T_OOA, population="EUR",
                                         initial_size=NE_ANC)
    return d


def simulate_gene_neutral(rng: np.random.Generator, length_bp: int,
                           seed: int) -> np.ndarray:
    """Run msprime simulation for a single gene; return the array of
    present-day minor allele frequencies for all segregating sites."""
    demog = build_demography()
    ts = msprime.sim_ancestry(
        samples={"EUR": N_DIPLOID_SAMPLES},
        demography=demog,
        sequence_length=length_bp,
        recombination_rate=RECOMB,
        random_seed=seed,
    )
    ts = msprime.sim_mutations(ts, rate=MU, random_seed=seed + 1,
                                model="binary")
    n_haps = ts.num_samples
    if ts.num_sites == 0:
        return np.array([])
    # Per-site allele frequency
    afs = np.zeros(ts.num_sites)
    for j, var in enumerate(ts.variants()):
        # variant.genotypes is binary 0/1 across haplotypes
        derived = var.genotypes.sum()
        af = derived / n_haps
        afs[j] = min(af, 1 - af)  # MAF
    return afs


# ---- Wright equilibrium PDF (same as PRF script) ----
P_GRID = np.linspace(1.0 / (2 * 10_000), 0.5, 1000)


def wright_pdf_at(p: float, s: float) -> float:
    if s <= 0 or s < 1e-8:
        return 1.0 / max(p, 1e-9)
    Ne = 10_000  # we use the standard PRF Ne for the equilibrium PDF
    gamma = 2 * Ne * s
    if gamma > 50:
        return np.exp(-gamma * p) / max(p, 1e-9)
    num = 1.0 - np.exp(-gamma * (1.0 - p))
    den = max(p * (1.0 - p) * (1.0 - np.exp(-gamma)), 1e-12)
    return num / den


def selection_acceptance(p_obs: float, s: float) -> float:
    """Probability that a variant at observed frequency p_obs survives
    selection of coefficient s, relative to the neutral case."""
    if s <= 0 or s < 1e-8:
        return 1.0
    pdf_neutral = wright_pdf_at(p_obs, 0.0)
    pdf_selected = wright_pdf_at(p_obs, s)
    if pdf_neutral <= 0:
        return 0.0
    return min(1.0, pdf_selected / pdf_neutral)


def simulate_one_replicate(rng: np.random.Generator,
                              model_params: dict,
                              base_seed: int) -> dict:
    out = {}
    for gi, (gene, ginfo) in enumerate(GENES.items()):
        side = ginfo["side"]
        length = ginfo["length_bp"]
        params = model_params[side]
        mean_s = params["mean_s"]
        alpha = params["alpha"]
        sigma_beta = params["sigma_beta"]

        # 1. neutral msprime simulation
        gene_seed = base_seed * 7 + gi
        try:
            mafs = simulate_gene_neutral(rng, length, gene_seed)
        except Exception as e:
            out[gene] = {"max_abs_Z": 0.0, "n_GWS": 0, "error": str(e)}
            continue

        if len(mafs) == 0:
            out[gene] = {"max_abs_Z": 0.0, "n_GWS": 0, "n_segregating": 0,
                          "n_common": 0}
            continue

        # 2. assign per-variant selection coefficients from gamma DFE
        s_arr = rng.gamma(shape=DFE_SHAPE,
                            scale=mean_s / DFE_SHAPE,
                            size=len(mafs))
        s_arr = np.clip(s_arr, 1e-6, 0.5)

        # 3. selection acceptance filter
        accept_probs = np.array([selection_acceptance(mafs[i], s_arr[i])
                                   for i in range(len(mafs))])
        accept = rng.random(len(mafs)) < accept_probs
        mafs = mafs[accept]
        s_arr = s_arr[accept]
        n_seg = len(mafs)

        if n_seg == 0:
            out[gene] = {"max_abs_Z": 0.0, "n_GWS": 0, "n_segregating": 0,
                          "n_common": 0}
            continue

        # 4. assign migraine effects
        if alpha > 0:
            sign = rng.choice([-1, 1], size=n_seg)
            beta = alpha * sign * np.sqrt(s_arr)
        else:
            beta = rng.normal(0.0, sigma_beta, size=n_seg)

        # 5. filter to common variants and compute Z
        common = mafs >= 0.01
        if common.sum() == 0:
            out[gene] = {"max_abs_Z": 0.0, "n_GWS": 0,
                          "n_segregating": int(n_seg), "n_common": 0}
            continue
        mafs_c = mafs[common]
        beta_c = beta[common]
        se = np.sqrt(1.0 / (2 * N_EFF * mafs_c * (1.0 - mafs_c)))
        z = beta_c / se
        out[gene] = {
            "max_abs_Z": float(np.max(np.abs(z))),
            "n_GWS": int(np.sum(np.abs(z) > GWS_THRESHOLD)),
            "n_segregating": int(n_seg),
            "n_common": int(common.sum()),
        }
    return out


def run_model(model_name: str, params: dict, n_rep: int, seed: int) -> dict:
    print(f"\n[{model_name}]  running {n_rep:,} replicates...")
    t0 = time.time()
    rng = np.random.default_rng(seed)
    per_gene_maxZ = {g: np.empty(n_rep) for g in GENES}
    for r in range(n_rep):
        rep = simulate_one_replicate(rng, params, base_seed=seed + r * 1000)
        for g in GENES:
            per_gene_maxZ[g][r] = rep[g].get("max_abs_Z", 0.0)
        if (r + 1) % 100 == 0:
            print(f"  replicate {r+1:,}/{n_rep:,}  ({time.time()-t0:.1f}s)")

    lig_maxZ = np.maximum(per_gene_maxZ["CALCA"], per_gene_maxZ["CALCB"])
    rec_maxZ = np.maximum.reduce([per_gene_maxZ["CALCRL"],
                                    per_gene_maxZ["RAMP1"],
                                    per_gene_maxZ["CRCP"]])
    hit_pattern = (lig_maxZ > GWS_THRESHOLD) & (rec_maxZ < RECEPTOR_INVISIBLE_THRESHOLD)
    summary = {
        "model": model_name,
        "n_replicates": n_rep,
        "elapsed_s": time.time() - t0,
        "ligand_GWS_rate": float(np.mean(lig_maxZ > GWS_THRESHOLD)),
        "receptor_invisible_rate": float(np.mean(rec_maxZ < RECEPTOR_INVISIBLE_THRESHOLD)),
        "hit_pattern_rate": float(np.mean(hit_pattern)),
    }
    print(f"  done. hit_pattern_rate = {summary['hit_pattern_rate']:.3f}  "
          f"(elapsed {summary['elapsed_s']:.1f}s)")
    return summary


def main() -> None:
    print(f"Master seed: {SEED}")
    print(f"Replicates per model: {N_REPLICATES:,}")
    out: dict = {
        "spec": "analysis/sim_msprime/PRE_REGISTRATION_msprime.md",
        "seed": SEED,
        "n_replicates": N_REPLICATES,
        "by_model": {},
    }
    for i, (name, params) in enumerate(MODELS.items()):
        out["by_model"][name] = run_model(name, params, N_REPLICATES,
                                            seed=SEED + i * 1_000_000)

    out_path = RES / "msprime_simulation.json"
    out_path.write_text(json.dumps(out, indent=2))
    print(f"\nWrote {out_path.relative_to(ROOT)}")

    print("\n=== Summary ===")
    print(f"{'model':30s}  {'lig_GWS':>10s}  {'rec_inv':>10s}  {'hit_patt':>10s}")
    for name, mout in out["by_model"].items():
        print(f"{name:30s}  "
              f"{mout['ligand_GWS_rate']:10.3f}  "
              f"{mout['receptor_invisible_rate']:10.3f}  "
              f"{mout['hit_pattern_rate']:10.3f}")


if __name__ == "__main__":
    main()
