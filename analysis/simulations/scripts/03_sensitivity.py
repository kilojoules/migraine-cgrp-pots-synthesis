"""
Sensitivity analysis on the asymmetry simulation.

Sweep the two parameters most likely to drive the M1 vs M2 ranking:
  M1 receptor mean_s   in {1e-3, 5e-3, 1e-2}
  M2 sigma_beta_recep  in {0.02, 0.05, 0.10}

Use 2,000 replicates per cell (vs 10,000 in the headline run) for speed.
Report the hit_pattern_rate ranking across the grid.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import importlib.util
spec = importlib.util.spec_from_file_location("sim", ROOT / "scripts" / "01_simulate_asymmetry.py")
sim = importlib.util.module_from_spec(spec)
spec.loader.exec_module(sim)

RES = ROOT / "results"
N_REP = 2_000

SWEEP_M1_MEAN_S_RECEP = [1e-3, 5e-3, 1e-2]
SWEEP_M2_SIGMA_RECEP  = [0.02, 0.05, 0.10]

out: dict = {
    "n_replicates_per_cell": N_REP,
    "M1_sweep": [],
    "M2_sweep": [],
    "M3_baseline_2k": None,
}

# M1 sweep (receptor mean_s)
for s_recep in SWEEP_M1_MEAN_S_RECEP:
    params = {
        "ligand":   {"mean_s": 1e-4, "alpha": 4.0, "sigma_beta": None},
        "receptor": {"mean_s": s_recep, "alpha": 4.0, "sigma_beta": None},
    }
    name = f"M1_purifying_recep_s={s_recep:.0e}"
    res = sim.run_model(name, params, N_REP, master_seed=20260428)
    out["M1_sweep"].append({
        "param_mean_s_recep": s_recep,
        "summary": res["summary"],
    })

# M2 sweep (sigma_beta receptor)
for sigma_r in SWEEP_M2_SIGMA_RECEP:
    params = {
        "ligand":   {"mean_s": 5e-4, "alpha": 0.0, "sigma_beta": 0.20},
        "receptor": {"mean_s": 5e-4, "alpha": 0.0, "sigma_beta": sigma_r},
    }
    name = f"M2_allelic_sigma_recep={sigma_r:.2f}"
    res = sim.run_model(name, params, N_REP, master_seed=20260429)
    out["M2_sweep"].append({
        "param_sigma_beta_recep": sigma_r,
        "summary": res["summary"],
    })

# M3 baseline at the same replicate count for comparable noise
m3_params = {
    "ligand":   {"mean_s": 5e-4, "alpha": 0.0, "sigma_beta": 0.10},
    "receptor": {"mean_s": 5e-4, "alpha": 0.0, "sigma_beta": 0.10},
}
res = sim.run_model("M3_null_2k", m3_params, N_REP, master_seed=20260430)
out["M3_baseline_2k"] = {"summary": res["summary"]}

(RES / "asymmetry_sensitivity.json").write_text(json.dumps(out, indent=2))

print("\n=== M1 sweep over receptor mean_s ===")
for cell in out["M1_sweep"]:
    s = cell["summary"]
    print(f"  s_recep={cell['param_mean_s_recep']:.0e}  "
          f"hit_pattern={s['hit_pattern_rate']:.3f}  "
          f"lig_GWS={s['ligand_GWS_rate']:.3f}  "
          f"rec_inv={s['receptor_invisible_rate']:.3f}")

print("\n=== M2 sweep over receptor sigma_beta ===")
for cell in out["M2_sweep"]:
    s = cell["summary"]
    print(f"  sigma_recep={cell['param_sigma_beta_recep']:.2f}  "
          f"hit_pattern={s['hit_pattern_rate']:.3f}  "
          f"lig_GWS={s['ligand_GWS_rate']:.3f}  "
          f"rec_inv={s['receptor_invisible_rate']:.3f}")

print("\n=== M3 baseline (2k) ===")
s = out["M3_baseline_2k"]["summary"]
print(f"  hit_pattern={s['hit_pattern_rate']:.3f}  "
      f"lig_GWS={s['ligand_GWS_rate']:.3f}  "
      f"rec_inv={s['receptor_invisible_rate']:.3f}")
