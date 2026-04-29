"""
Simulate the CGRP ligand-receptor asymmetry under three competing models.

Wright / Poisson Random Field equilibrium Monte Carlo per
analysis/simulations/SIM_SPEC.md (locked before this script's output
was opened). Three models:

  M1: Purifying selection on receptor genes (CALCRL, RAMP1, CRCP);
      ligand variants weakly selected.
  M2: Allelic-series concentration: same selection both sides,
      receptor variants have intrinsically smaller per-variant
      effects on migraine.
  M3: Null: no asymmetry, both sides equal in selection and effect.

Per replicate per gene:
  N_seg     Poisson draw of segregating-variant count
  s_i       gamma DFE; mean s differs by gene/model
  beta_i    effect size; coupling to s differs by model
  MAF_i     selection-drift equilibrium frequency (Wright)
  Z_i      simulated GWAS test statistic at Hautakangas-scale N
"""
from __future__ import annotations

import json
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "results"
RES.mkdir(parents=True, exist_ok=True)

# ----- Locked parameters -----
NE = 10_000                    # effective population size
MU = 1.5e-8                    # per-bp per-generation mutation rate
N_CASES = 102_084              # Hautakangas et al 2022
N_CONTROLS = 771_257
N_EFF = 4 / (1 / N_CASES + 1 / N_CONTROLS)
N_REPLICATES = 10_000
SEED = 20260428

# Gene CDS lengths (bp). MANE-select transcripts (Ensembl).
GENES = {
    "CALCA":  {"length_bp":  384, "side": "ligand"},
    "CALCB":  {"length_bp":  381, "side": "ligand"},
    "CALCRL": {"length_bp": 1380, "side": "receptor"},
    "RAMP1":  {"length_bp":  444, "side": "receptor"},
    "CRCP":   {"length_bp":  444, "side": "receptor"},
}

GWS_THRESHOLD = 5.452           # |Z| corresponding to P = 5e-8 two-sided
RECEPTOR_INVISIBLE_THRESHOLD = 4.26  # |Z| corresponding to P = 1e-4 two-sided

# ----- Model parameterizations -----
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

# DFE shape (gamma); same across models.
DFE_GAMMA_SHAPE = 0.2


# ----- Wright's equilibrium MAF distribution under selection-drift -----
# f(p | s) propto (1 - exp(-2 Ne s (1-p))) / (p (1-p) (1 - exp(-2 Ne s)))
# Sampled via inverse CDF on a fine grid in [1/(2Ne), 0.5].
P_GRID = np.linspace(1.0 / (2 * NE), 0.5, 2000)


def wright_pdf(p_grid: np.ndarray, s: float) -> np.ndarray:
    """Wright equilibrium density at frequency p_grid given selection s.

    s > 0 is deleterious. Limit s -> 0 reduces to f(p) propto 1/p.
    Returns un-normalized density on the grid.
    """
    if s <= 0:
        return 1.0 / p_grid
    gamma = 2 * NE * s
    if gamma > 50:  # numeric overflow; nearly certain to be lost
        # mutation-selection balance: most variants near MAF = mu/s,
        # but at MAF >= 1/(2Ne) only a thin tail; use exponential decay
        return np.exp(-gamma * p_grid) / p_grid
    num = 1.0 - np.exp(-gamma * (1.0 - p_grid))
    den = p_grid * (1.0 - p_grid) * (1.0 - np.exp(-gamma))
    return num / den


def sample_maf(rng: np.random.Generator, s: float, n: int) -> np.ndarray:
    """Sample n MAF values from Wright equilibrium at selection s."""
    pdf = wright_pdf(P_GRID, s)
    cdf = np.cumsum(pdf)
    cdf = cdf / cdf[-1]
    u = rng.random(n)
    idx = np.searchsorted(cdf, u)
    idx = np.clip(idx, 0, len(P_GRID) - 1)
    return P_GRID[idx]


def expected_seg_variants(length_bp: int, mean_s: float) -> float:
    """Expected number of segregating variants at MAF >= 1/(2Ne) for a
    gene of given CDS length under PRF / Wright model with mean s.

    For neutral case: theta * H_(2Ne-1) ~ theta * log(2Ne). With selection,
    integrate Wright density over [1/(2Ne), 0.5].
    """
    theta = 4 * NE * MU * length_bp
    pdf = wright_pdf(P_GRID, mean_s)
    integral = np.trapz(pdf, P_GRID)
    # Normalize relative to neutral (so units of expected variants):
    # under neutral, integral of 1/p from 1/(2Ne) to 0.5 = log(0.5 / (1/(2Ne)))
    # ~ log(Ne)
    return theta * integral


def simulate_one_replicate(
    rng: np.random.Generator,
    model_params: dict,
) -> dict:
    """Run one replicate across all five genes; return per-gene max |Z|
    and number of GWS hits."""
    out = {}
    for gene, ginfo in GENES.items():
        side = ginfo["side"]
        L = ginfo["length_bp"]
        params = model_params[side]
        mean_s = params["mean_s"]
        alpha = params["alpha"]
        sigma_beta = params["sigma_beta"]

        lambda_n = expected_seg_variants(L, mean_s)
        n_seg = rng.poisson(lambda_n)
        if n_seg == 0:
            out[gene] = {"max_abs_Z": 0.0, "n_GWS": 0, "n_segregating": 0,
                          "n_common": 0}
            continue

        # Per-variant selection coefficients from gamma DFE
        s_arr = rng.gamma(shape=DFE_GAMMA_SHAPE,
                          scale=mean_s / DFE_GAMMA_SHAPE,
                          size=n_seg)
        s_arr = np.clip(s_arr, 1e-6, 0.5)

        # Per-variant MAF from Wright equilibrium at the variant's own s
        maf = np.empty(n_seg)
        for i, s in enumerate(s_arr):
            maf[i] = sample_maf(rng, s, 1)[0]

        # Per-variant migraine effect size beta
        if alpha > 0:
            # M1: beta = alpha * sign() * sqrt(s) (DFE-effect coupling)
            sign = rng.choice([-1, 1], size=n_seg)
            beta = alpha * sign * np.sqrt(s_arr)
        else:
            # M2, M3: beta ~ N(0, sigma_beta^2), independent of s
            beta = rng.normal(0.0, sigma_beta, size=n_seg)

        # Filter to common variants (MAF >= 0.01)
        common_mask = maf >= 0.01
        n_common = int(common_mask.sum())
        if n_common == 0:
            out[gene] = {"max_abs_Z": 0.0, "n_GWS": 0,
                          "n_segregating": int(n_seg), "n_common": 0}
            continue

        beta_c = beta[common_mask]
        maf_c = maf[common_mask]
        # GWAS standard error under additive model
        se = np.sqrt(1.0 / (2 * N_EFF * maf_c * (1.0 - maf_c)))
        z = beta_c / se
        max_abs_Z = float(np.max(np.abs(z)))
        n_GWS = int(np.sum(np.abs(z) > GWS_THRESHOLD))
        out[gene] = {
            "max_abs_Z": max_abs_Z,
            "n_GWS": n_GWS,
            "n_segregating": int(n_seg),
            "n_common": n_common,
        }
    return out


def run_model(model_name: str, model_params: dict, n_rep: int,
              master_seed: int) -> dict:
    print(f"\n[{model_name}]  running {n_rep:,} replicates...")
    t0 = time.time()
    rng = np.random.default_rng(master_seed)

    # Storage
    per_gene_maxZ = {g: np.empty(n_rep) for g in GENES}
    per_gene_nGWS = {g: np.empty(n_rep, dtype=int) for g in GENES}

    for r in range(n_rep):
        rep = simulate_one_replicate(rng, model_params)
        for g in GENES:
            per_gene_maxZ[g][r] = rep[g]["max_abs_Z"]
            per_gene_nGWS[g][r] = rep[g]["n_GWS"]
        if (r + 1) % 1000 == 0:
            print(f"  replicate {r+1:,}/{n_rep:,}  ({time.time()-t0:.1f}s)")

    # Aggregate
    lig_maxZ = np.maximum(per_gene_maxZ["CALCA"], per_gene_maxZ["CALCB"])
    rec_maxZ = np.maximum.reduce([per_gene_maxZ["CALCRL"],
                                   per_gene_maxZ["RAMP1"],
                                   per_gene_maxZ["CRCP"]])
    asym = lig_maxZ - rec_maxZ
    hit_pattern = (lig_maxZ > GWS_THRESHOLD) & (rec_maxZ < RECEPTOR_INVISIBLE_THRESHOLD)
    lig_GWS = (lig_maxZ > GWS_THRESHOLD)

    summary = {
        "model": model_name,
        "n_replicates": n_rep,
        "elapsed_s": time.time() - t0,
        "per_gene_max_abs_Z_mean": {g: float(np.mean(per_gene_maxZ[g])) for g in GENES},
        "per_gene_max_abs_Z_median": {g: float(np.median(per_gene_maxZ[g])) for g in GENES},
        "per_gene_max_abs_Z_q95": {g: float(np.quantile(per_gene_maxZ[g], 0.95)) for g in GENES},
        "per_gene_GWS_rate": {g: float(np.mean(per_gene_nGWS[g] > 0)) for g in GENES},
        "ligand_GWS_rate": float(np.mean(lig_GWS)),
        "receptor_invisible_rate": float(np.mean(rec_maxZ < RECEPTOR_INVISIBLE_THRESHOLD)),
        "hit_pattern_rate": float(np.mean(hit_pattern)),
        "asymmetry_mean": float(np.mean(asym)),
        "asymmetry_median": float(np.median(asym)),
        "asymmetry_q05_q95": [float(np.quantile(asym, 0.05)),
                                float(np.quantile(asym, 0.95))],
    }
    print(f"  done. hit_pattern_rate = {summary['hit_pattern_rate']:.3f}")

    return {
        "summary": summary,
        "per_gene_max_abs_Z": {g: per_gene_maxZ[g].tolist() for g in GENES},
        "asymmetry": asym.tolist(),
        "ligand_max_Z": lig_maxZ.tolist(),
        "receptor_max_Z": rec_maxZ.tolist(),
        "hit_pattern": hit_pattern.tolist(),
    }


def main() -> None:
    print(f"Master seed: {SEED}")
    print(f"Replicates per model: {N_REPLICATES:,}")
    print(f"GWAS sample: N_cases={N_CASES:,}  N_controls={N_CONTROLS:,}  N_eff={N_EFF:,.0f}")

    out: dict = {
        "spec": "analysis/simulations/SIM_SPEC.md",
        "seed": SEED,
        "n_replicates": N_REPLICATES,
        "Ne": NE,
        "mu_per_bp_per_gen": MU,
        "N_cases": N_CASES,
        "N_controls": N_CONTROLS,
        "N_effective": N_EFF,
        "GWS_threshold": GWS_THRESHOLD,
        "receptor_invisible_threshold": RECEPTOR_INVISIBLE_THRESHOLD,
        "genes": GENES,
        "models": MODELS,
        "by_model": {},
    }
    for i, (name, params) in enumerate(MODELS.items()):
        # Use a different stream per model to keep results independent
        out["by_model"][name] = run_model(name, params, N_REPLICATES,
                                           master_seed=SEED + i)

    out_path = RES / "asymmetry_simulation.json"
    out_path.write_text(json.dumps(out, indent=2))
    print(f"\nWrote {out_path.relative_to(ROOT)}  ({out_path.stat().st_size / 1e6:.1f} MB)")

    print("\n=== Summary ===")
    print(f"{'model':30s}  {'lig_GWS':>10s}  {'rec_inv':>10s}  {'hit_patt':>10s}  {'asym_med':>10s}")
    for name, mout in out["by_model"].items():
        s = mout["summary"]
        print(f"{name:30s}  "
              f"{s['ligand_GWS_rate']:10.3f}  "
              f"{s['receptor_invisible_rate']:10.3f}  "
              f"{s['hit_pattern_rate']:10.3f}  "
              f"{s['asymmetry_median']:10.2f}")


if __name__ == "__main__":
    main()
