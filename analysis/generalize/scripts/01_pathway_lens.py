"""
Extension C: query OpenTargets L2G for 10 peptide-hormone pathways
and compute per-pathway ligand-vs-receptor asymmetry score.

Decision rule and seed locked in
analysis/generalize/PRE_REGISTRATION_generalize.md.
"""
from __future__ import annotations

import json
import time
from pathlib import Path

import numpy as np
import requests

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "results"
RES.mkdir(parents=True, exist_ok=True)

GQL = "https://api.platform.opentargets.org/api/v4/graphql"
SEED = 20260501

# Locked pathway definitions (pre-registered)
PATHWAYS = {
    "CGRP": {
        "ligand": {"CALCA": "ENSG00000110680", "CALCB": "ENSG00000175868"},
        "receptor": {"CALCRL": "ENSG00000064989", "RAMP1": "ENSG00000132329",
                      "CRCP": "ENSG00000241258"},
        "disease_id": "MONDO_0005277",
        "disease_label": "migraine",
    },
    "GLP-1": {
        "ligand": {"GCG": "ENSG00000115263"},
        "receptor": {"GLP1R": "ENSG00000112164"},
        "disease_id": "MONDO_0005148",
        "disease_label": "type 2 diabetes",
    },
    "GIP": {
        "ligand": {"GIP": "ENSG00000159224"},
        "receptor": {"GIPR": "ENSG00000010310"},
        "disease_id": "MONDO_0005148",
        "disease_label": "type 2 diabetes",
    },
    "Kisspeptin": {
        "ligand": {"KISS1": "ENSG00000170498"},
        "receptor": {"KISS1R": "ENSG00000116014"},
        "disease_id": "EFO_0004703",
        "disease_label": "age at menarche",
    },
    "Leptin": {
        "ligand": {"LEP": "ENSG00000174697"},
        "receptor": {"LEPR": "ENSG00000116678"},
        "disease_id": "EFO_0004340",
        "disease_label": "body mass index",
    },
    "Oxytocin": {
        "ligand": {"OXT": "ENSG00000101405"},
        "receptor": {"OXTR": "ENSG00000180914"},
        "disease_id": "MONDO_0002009",
        "disease_label": "major depressive disorder",
    },
    "Insulin": {
        "ligand": {"INS": "ENSG00000254647"},
        "receptor": {"INSR": "ENSG00000171105"},
        "disease_id": "MONDO_0005148",
        "disease_label": "type 2 diabetes",
    },
    "Adiponectin": {
        "ligand": {"ADIPOQ": "ENSG00000181092"},
        "receptor": {"ADIPOR1": "ENSG00000159346",
                      "ADIPOR2": "ENSG00000006831"},
        "disease_id": "MONDO_0005148",
        "disease_label": "type 2 diabetes",
    },
    "Ghrelin": {
        "ligand": {"GHRL": "ENSG00000157017"},
        "receptor": {"GHSR": "ENSG00000121853"},
        "disease_id": "EFO_0004340",
        "disease_label": "body mass index",
    },
    "SubstanceP_NK1": {
        "ligand": {"TAC1": "ENSG00000006128"},
        "receptor": {"TACR1": "ENSG00000115353"},
        "disease_id": "MONDO_0002009",
        "disease_label": "major depressive disorder",
    },
}


_DISEASE_CACHE: dict[str, dict[str, dict]] = {}


def fetch_disease_targets(disease_efo: str) -> dict[str, dict]:
    """Fetch all associated targets for a disease, return symbol -> {overall, datatypes}."""
    if disease_efo in _DISEASE_CACHE:
        return _DISEASE_CACHE[disease_efo]
    q = """
    query Q($did: String!, $size: Int!) {
      disease(efoId: $did) {
        id name
        associatedTargets(page: {index: 0, size: $size}) {
          count
          rows {
            target { id approvedSymbol }
            score
            datatypeScores { id score }
          }
        }
      }
    }
    """
    r = requests.post(GQL, json={"query": q,
                                   "variables": {"did": disease_efo,
                                                  "size": 3000}},
                       timeout=180)
    r.raise_for_status()
    payload = r.json()
    if "errors" in payload:
        raise RuntimeError(payload["errors"])
    disease = payload["data"]["disease"]
    if disease is None:
        _DISEASE_CACHE[disease_efo] = {}
        return {}
    out = {}
    for row in disease["associatedTargets"]["rows"]:
        sym = row["target"]["approvedSymbol"]
        out[sym] = {
            "overall": float(row["score"]),
            "datatype_scores": {ds["id"]: float(ds["score"])
                                  for ds in (row.get("datatypeScores") or [])},
        }
    _DISEASE_CACHE[disease_efo] = out
    return out


def query_target_disease(target_ensembl: str, target_symbol: str,
                          disease_efo: str) -> dict:
    """Return overall + datatype scores for a target-disease pair, via disease lookup."""
    disease_targets = fetch_disease_targets(disease_efo)
    if not disease_targets:
        return {"disease_not_found": True}
    if target_symbol not in disease_targets:
        return {"no_evidence": True, "approvedSymbol": target_symbol}
    row = disease_targets[target_symbol]
    return {
        "approvedSymbol": target_symbol,
        "overall": row["overall"],
        "datatype_scores": row["datatype_scores"],
    }


def gene_l2g(score_dict: dict) -> float:
    """Return L2G genetic-association score; 0 if absent or no evidence."""
    if "error" in score_dict or score_dict.get("target_not_found"):
        return 0.0
    if score_dict.get("no_evidence"):
        return 0.0
    return float(score_dict.get("datatype_scores", {})
                  .get("genetic_association", 0.0))


def main() -> None:
    out_per_gene: dict = {}
    out_summary: dict = {}

    print(f"Querying OpenTargets for {len(PATHWAYS)} pathways...")
    for pname, pinfo in PATHWAYS.items():
        out_per_gene[pname] = {
            "disease_id": pinfo["disease_id"],
            "disease_label": pinfo["disease_label"],
            "ligand": {}, "receptor": {},
        }
        for side in ("ligand", "receptor"):
            for sym, eid in pinfo[side].items():
                res = query_target_disease(eid, sym, pinfo["disease_id"])
                out_per_gene[pname][side][sym] = res
                time.sleep(0.05)

        ligand_scores = [gene_l2g(out_per_gene[pname]["ligand"][s])
                          for s in pinfo["ligand"]]
        receptor_scores = [gene_l2g(out_per_gene[pname]["receptor"][s])
                            for s in pinfo["receptor"]]
        med_lig = float(np.median(ligand_scores)) if ligand_scores else 0.0
        med_rec = float(np.median(receptor_scores)) if receptor_scores else 0.0
        A = med_lig - med_rec

        max_side = max(med_lig, med_rec)
        if max_side < 0.10:
            cls = "both-quiet"
        elif A > 0.20:
            cls = "ligand-asymmetric"
        elif A < -0.20:
            cls = "receptor-asymmetric"
        else:
            cls = "balanced"

        out_summary[pname] = {
            "disease_label": pinfo["disease_label"],
            "ligand_genes": list(pinfo["ligand"].keys()),
            "receptor_genes": list(pinfo["receptor"].keys()),
            "ligand_l2g_median": med_lig,
            "receptor_l2g_median": med_rec,
            "asymmetry": A,
            "class": cls,
        }
        print(f"  {pname:18s}  disease={pinfo['disease_label']:25s}  "
              f"lig_med={med_lig:.3f}  rec_med={med_rec:.3f}  "
              f"A={A:+.3f}  class={cls}")

    # Cross-pathway permutation test
    rng = np.random.default_rng(SEED)
    evaluable = [p for p in PATHWAYS
                  if out_summary[p]["class"] != "both-quiet"]
    n_lig_obs = sum(out_summary[p]["class"] == "ligand-asymmetric"
                     for p in evaluable)
    print(f"\nEvaluable pathways: {len(evaluable)}/{len(PATHWAYS)}")
    print(f"Ligand-asymmetric pathways: {n_lig_obs}")

    if evaluable:
        N_PERM = 10_000
        n_relabel_lig = np.zeros(N_PERM, dtype=int)
        for i in range(N_PERM):
            count_lig = 0
            for p in evaluable:
                med_lig = out_summary[p]["ligand_l2g_median"]
                med_rec = out_summary[p]["receptor_l2g_median"]
                if rng.random() < 0.5:
                    A_perm = med_rec - med_lig
                else:
                    A_perm = med_lig - med_rec
                if A_perm > 0.20:
                    count_lig += 1
            n_relabel_lig[i] = count_lig
        p_perm = float((n_relabel_lig >= n_lig_obs).sum() / N_PERM)
        print(f"Permutation p (>= {n_lig_obs} ligand-asym): {p_perm:.4f}")
    else:
        p_perm = float("nan")

    if n_lig_obs >= 7 and p_perm < 0.05:
        verdict = "General asymmetry"
    elif n_lig_obs <= 3 and (p_perm > 0.20 or np.isnan(p_perm)):
        verdict = "CGRP-specific"
    else:
        verdict = "Inconclusive"
    print(f"\nVerdict (pre-registered): {verdict}")

    # Save
    out: dict = {
        "spec": "analysis/generalize/PRE_REGISTRATION_generalize.md",
        "seed": SEED,
        "n_pathways": len(PATHWAYS),
        "n_evaluable": len(evaluable),
        "n_ligand_asymmetric": n_lig_obs,
        "permutation_p": p_perm,
        "verdict": verdict,
        "per_pathway_summary": out_summary,
    }
    (RES / "asymmetry_summary.json").write_text(json.dumps(out, indent=2))
    (RES / "pathway_lens.json").write_text(json.dumps(out_per_gene, indent=2,
                                                          default=str))
    print(f"\nWrote results/asymmetry_summary.json  and  pathway_lens.json")


if __name__ == "__main__":
    main()
