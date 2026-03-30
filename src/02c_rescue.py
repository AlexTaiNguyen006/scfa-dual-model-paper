#!/usr/bin/env python3
"""PPCOACm propionate pathway rescue analysis.

Recon3D lacks functional propionyl-CoA carboxylase (PPCOACm) under the
strict hepatocyte medium.  This script:
  1. Runs unrescued Recon3D (propionate = 0).
  2. Adds PPCOACm *without* the strict medium (unconstrained) — shows
     thermodynamic-loop inflation.
  3. Adds PPCOACm *with* the strict medium (constrained) — converges
     with Human-GEM.
  4. Runs Human-GEM as reference.

Outputs
  results/table_rescue_constrained.csv
  outputs/figs/Fig7.png  /  Fig7.pdf
"""

import os

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import cobra
from cobra import Reaction
from cobra.io import read_sbml_model

from .utils import build_paths, load_config, decompress_gz
from .run_simulation_shared import (
    SCFA_EXCHANGE_IDS, GLUCOSE_IDS, O2_IDS, CO2_IDS, ATPM_IDS,
    _find_rxn, setup_medium,
)

_HGEM_ATPM_METS = {
    "MAM01371c": -1, "MAM02040c": -1, "MAM01285c": 1,
    "MAM02751c": 1, "MAM02039c": 1,
}

DOSES = {
    "Baseline": (0.0, 0.0, 0.0),
    "Low":      (2.0, 0.7, 0.4),
    "Mid":      (4.0, 1.4, 0.8),
    "High":     (8.0, 2.8, 1.6),
}


def _add_ppcoacm(model):
    """Add propionyl-CoA carboxylase if absent."""
    if "PPCOACm" in model.reactions:
        return
    met_stoich = {
        "ppcoa_m": -1, "hco3_m": -1, "atp_m": -1,
        "mmcoa__S_m": +1, "adp_m": +1, "pi_m": +1, "h_m": +1,
    }
    rxn = Reaction("PPCOACm")
    rxn.name = "Propionyl-CoA carboxylase, mitochondrial"
    rxn.bounds = (0.0, 1000.0)
    mets = {}
    for mid, coeff in met_stoich.items():
        if mid in model.metabolites:
            mets[model.metabolites.get_by_id(mid)] = coeff
    rxn.add_metabolites(mets)
    model.add_reactions([rxn])
    print(f"  added PPCOACm: {rxn.reaction}")


def _setup(model, label, cfg):
    """Strict hepatocyte medium + ATPM objective."""
    setup_medium(model, cfg)
    atpm_rxn = _find_rxn(model, ATPM_IDS, silent=True)
    if atpm_rxn is None:
        if "ATPM_added" in model.reactions:
            atpm_rxn = model.reactions.get_by_id("ATPM_added")
        else:
            atpm_rxn = cobra.Reaction("ATPM_added")
            atpm_rxn.name = "ATP maintenance (added)"
            model.add_reactions([atpm_rxn])
            mets = {}
            for mid, coeff in _HGEM_ATPM_METS.items():
                if mid in model.metabolites:
                    mets[model.metabolites.get_by_id(mid)] = coeff
            atpm_rxn.add_metabolites(mets)
    atpm_rxn.bounds = (0, 500)
    model.objective = atpm_rxn
    scfa_rxns = {}
    for name, candidates in SCFA_EXCHANGE_IDS.items():
        scfa_rxns[name] = _find_rxn(model, candidates, silent=True)
    return atpm_rxn, scfa_rxns


def _run_doses(model, label, atpm_rxn, scfa_rxns):
    """Run all dose conditions, return {cond: {atpm, ppa_flux}}."""
    results = {}
    ppa_rxn = scfa_rxns.get("propionate")
    for cond, (ac, ppa, but) in DOSES.items():
        with model:
            for name, dose in [("acetate", ac), ("propionate", ppa),
                               ("butyrate", but)]:
                rxn = scfa_rxns.get(name)
                if rxn and dose > 0:
                    rxn.bounds = (-dose, 0)
            sol = model.optimize()
            atpm_val = float(sol.objective_value or 0)
            ppa_flux = float(sol.fluxes.get(ppa_rxn.id, 0)) if ppa_rxn else 0.0
        results[cond] = {"atpm": atpm_val, "ppa_flux": ppa_flux}
        print(f"    {label:35s} {cond:10s} ATPM={atpm_val:.2f}  PPA={ppa_flux:.4f}")
    return results


def _setup_unconstrained(model, label):
    """ATPM objective only, no strict medium (allows loops)."""
    atpm_rxn = _find_rxn(model, ATPM_IDS, silent=True)
    if atpm_rxn is None:
        if "ATPM_added" in model.reactions:
            atpm_rxn = model.reactions.get_by_id("ATPM_added")
        else:
            atpm_rxn = cobra.Reaction("ATPM_added")
            atpm_rxn.name = "ATP maintenance (added)"
            model.add_reactions([atpm_rxn])
    atpm_rxn.bounds = (0, 500)
    model.objective = atpm_rxn
    scfa_rxns = {}
    for name, candidates in SCFA_EXCHANGE_IDS.items():
        scfa_rxns[name] = _find_rxn(model, candidates, silent=True)
    return atpm_rxn, scfa_rxns


def main():
    paths = build_paths()
    cfg = load_config(paths.config_path)

    print("=" * 50)
    print("PPCOACm Rescue Analysis")
    print("=" * 50)

    # load Recon3D
    sbml_file = decompress_gz(paths.sbml_path)
    print(f"\nloading Recon3D from {sbml_file}")
    recon_base = read_sbml_model(str(sbml_file))
    print(f"  {len(recon_base.reactions)} reactions")

    # 1. original Recon3D (strict medium, no PPCOACm)
    recon_orig = recon_base.copy()
    atpm_o, scfa_o = _setup(recon_orig, "Recon3D-original", cfg)
    results_orig = _run_doses(recon_orig, "Recon3D (original)", atpm_o, scfa_o)

    # 2. unconstrained rescue (PPCOACm added, NO strict medium)
    recon_unc = recon_base.copy()
    _add_ppcoacm(recon_unc)
    atpm_u, scfa_u = _setup_unconstrained(recon_unc, "Recon3D-unconstrained")
    results_unc = _run_doses(recon_unc, "Recon3D+PPCOACm (unconstrained)", atpm_u, scfa_u)

    # 3. constrained rescue (PPCOACm added, strict medium)
    recon_con = recon_base.copy()
    _add_ppcoacm(recon_con)
    atpm_c, scfa_c = _setup(recon_con, "Recon3D-constrained", cfg)
    results_con = _run_doses(recon_con, "Recon3D+PPCOACm (constrained)", atpm_c, scfa_c)

    # 4. Human-GEM reference
    hgem_path = paths.human_gem_path
    print(f"\nloading Human-GEM from {hgem_path}")
    hgem = read_sbml_model(str(hgem_path))
    print(f"  {len(hgem.reactions)} reactions")
    atpm_h, scfa_h = _setup(hgem, "Human-GEM", cfg)
    results_hgem = _run_doses(hgem, "Human-GEM", atpm_h, scfa_h)

    # assemble results table
    conditions = list(DOSES.keys())
    rows = []
    for cond in conditions:
        o, u, c, h = (results_orig[cond], results_unc[cond],
                       results_con[cond], results_hgem[cond])
        rows.append({
            "Condition":                   cond,
            "Recon3D_orig_PPA":            round(o["ppa_flux"], 3),
            "Rescued_unconstrained_PPA":   round(u["ppa_flux"], 3),
            "Rescued_constrained_PPA":     round(c["ppa_flux"], 3),
            "Human_GEM_PPA":              round(h["ppa_flux"], 3),
            "Recon3D_orig_ATPM":           round(o["atpm"], 2),
            "Rescued_unconstrained_ATPM":  round(u["atpm"], 2),
            "Rescued_constrained_ATPM":    round(c["atpm"], 2),
            "Human_GEM_ATPM":             round(h["atpm"], 2),
        })

    df = pd.DataFrame(rows)
    csv_out = paths.results / "table_rescue_constrained.csv"
    df.to_csv(csv_out, index=False)
    print(f"\nwrote {csv_out}")
    print(df.to_string(index=False))

    # Fig 7 (two panels)
    fig, (ax_a, ax_b) = plt.subplots(1, 2, figsize=(11, 5))
    dose_conds = ["Low", "Mid", "High"]
    x = np.arange(len(dose_conds))
    w = 0.25
    COL_ORIG = "#3a86ff"
    COL_RESC = "#06d6a0"
    COL_HGEM = "#ef476f"

    # panel A: propionate flux
    ppa_orig = [-results_orig[c]["ppa_flux"] for c in dose_conds]
    ppa_con  = [-results_con[c]["ppa_flux"] for c in dose_conds]
    ppa_hgem = [-results_hgem[c]["ppa_flux"] for c in dose_conds]

    ax_a.bar(x - w, ppa_orig, w, color=COL_ORIG, label="Recon3D (original)",
             edgecolor="black", linewidth=0.5)
    ax_a.bar(x,     ppa_con, w, color=COL_RESC, label="Recon3D + PPCOACm",
             edgecolor="black", linewidth=0.5)
    ax_a.bar(x + w, ppa_hgem, w, color=COL_HGEM, label="Human-GEM",
             edgecolor="black", linewidth=0.5)
    ax_a.set_xlabel("Dose Condition")
    ax_a.set_ylabel("Propionate uptake (mmol gDW\u207b\u00b9 hr\u207b\u00b9)")
    ax_a.set_title("A   Propionate exchange flux", loc="left", fontweight="bold")
    ax_a.set_xticks(x)
    ax_a.set_xticklabels(dose_conds)
    ax_a.legend(fontsize=8)
    ax_a.spines["top"].set_visible(False)
    ax_a.spines["right"].set_visible(False)
    ax_a.set_ylim(0, max(ppa_hgem) * 1.3)

    # panel B: ATPM (constrained rescue vs originals)
    atpm_orig_all = [results_orig[c]["atpm"] for c in conditions]
    atpm_con_all  = [results_con[c]["atpm"] for c in conditions]
    atpm_hgem_all = [results_hgem[c]["atpm"] for c in conditions]
    x_b = np.arange(len(conditions))

    ax_b.bar(x_b - w, atpm_orig_all, w, color=COL_ORIG,
             label="Recon3D (original)", edgecolor="black", linewidth=0.5)
    ax_b.bar(x_b,     atpm_con_all, w, color=COL_RESC,
             label="Recon3D + PPCOACm (constrained)", edgecolor="black", linewidth=0.5)
    ax_b.bar(x_b + w, atpm_hgem_all, w, color=COL_HGEM,
             label="Human-GEM", edgecolor="black", linewidth=0.5)
    ax_b.set_xlabel("Dose Condition")
    ax_b.set_ylabel("Predicted ATPM (mmol gDW\u207b\u00b9 hr\u207b\u00b9)")
    ax_b.set_title("B   Predicted ATPM flux", loc="left", fontweight="bold")
    ax_b.set_xticks(x_b)
    ax_b.set_xticklabels(conditions)
    ax_b.legend(fontsize=8)
    ax_b.spines["top"].set_visible(False)
    ax_b.spines["right"].set_visible(False)

    plt.tight_layout()
    fig.savefig(paths.figs_dir / "Fig7.png", dpi=300, bbox_inches="tight")
    fig.savefig(paths.figs_dir / "Fig7.pdf", dpi=300, bbox_inches="tight")
    plt.close()
    print(f"saved Fig7.png / Fig7.pdf")


if __name__ == "__main__":
    main()
