#!/usr/bin/env python3


import os
import sys
import copy

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import cobra
from cobra import Reaction
from cobra.io import read_sbml_model

#paths 
REPO   os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  #repo root
FIGS   os.path.join(REPO, "outputs", "figs")
TABS   os.path.join(REPO, "outputs", "tables")
os.makedirs(FIGS, exist_okTrue)
os.makedirs(TABS, exist_okTrue)

sys.path.insert(0, REPO)
from src.run_simulation_shared import (
    SCFA_EXCHANGE_IDS, GLUCOSE_IDS, O2_IDS, CO2_IDS, ATPM_IDS,
    _find_rxn, setup_medium,
)

#human-gem atpm metabolites (for adding demand reaction) 
_HGEM_ATPM_METS  {
    "MAM01371c": -1, "MAM02040c": -1, "MAM01285c": 1,
    "MAM02751c": 1, "MAM02039c": 1,
}

#scfa dose conditions (mirrors 02_run_simulation.py) 
DOSES  {
    "Baseline": (0.0, 0.0, 0.0),
    "Low":      (2.0, 0.7, 0.4),
    "Mid":      (4.0, 1.4, 0.8),
    "High":     (8.0, 2.8, 1.6),
}


def _add_ppcoacm(model):
    
    if "PPCOACm" in model.reactions:
        print("  PPCOACm already present, skipping.")
        return

    met_stoich  {
        "ppcoa_m":    -1,
        "hco3_m":     -1,
        "atp_m":      -1,
        "mmcoa__S_m": +1,
        "adp_m":      +1,
        "pi_m":       +1,
        "h_m":        +1,
    }

    rxn  Reaction("PPCOACm")
    rxn.name  "Propionyl-CoA carboxylase (mitochondrial) [added for rescue]"
    rxn.bounds  (0.0, 1000.0)

    mets  {}
    missing  []
    for mid, coeff in met_stoich.items():
        if mid in model.metabolites:
            mets[model.metabolites.get_by_id(mid)]  coeff
        else:
            missing.append(mid)

    if missing:
        print(f"  WARNING: metabolites not found in model: {missing}")
        print("  Attempting partial rescue with available metabolites.")

    rxn.add_metabolites(mets)
    model.add_reactions([rxn])
    print(f"  Added PPCOACm: {rxn.reaction}  bounds{rxn.bounds}")
    return rxn


def _load_cfg():
    import yaml
    cfg_path  os.path.join(REPO, "data", "inputs", "project_config.yml")
    with open(cfg_path) as fh:
        return yaml.safe_load(fh)


def _setup(model, model_label, cfgNone):
    
    if cfg is None:
        cfg  _load_cfg()
    setup_medium(model, cfg)

    atpm_rxn  _find_rxn(model, ATPM_IDS, silentTrue)
    if atpm_rxn is None:
        if "ATPM_added" in model.reactions:
            atpm_rxn  model.reactions.get_by_id("ATPM_added")
        else:
            atpm_rxn  cobra.Reaction("ATPM_added")
            atpm_rxn.name  "ATP maintenance (added)"
            model.add_reactions([atpm_rxn])
            mets  {}
            for mid, coeff in _HGEM_ATPM_METS.items():
                if mid in model.metabolites:
                    mets[model.metabolites.get_by_id(mid)]  coeff
            atpm_rxn.add_metabolites(mets)

    atpm_rxn.bounds  (0, 500)
    model.objective  atpm_rxn

    scfa_rxns  {}
    for name, candidates in SCFA_EXCHANGE_IDS.items():
        scfa_rxns[name]  _find_rxn(model, candidates, silentTrue)
    return atpm_rxn, scfa_rxns


def run_doses(model, model_label, atpm_rxn, scfa_rxns):
    
    results  {}
    ppa_rxn  scfa_rxns.get("propionate")
    for cond, (ac, ppa, but) in DOSES.items():
        with model:
            for name, dose in [("acetate", ac), ("propionate", ppa), ("butyrate", but)]:
                rxn  scfa_rxns.get(name)
                if rxn and dose > 0:
                    rxn.bounds  (-dose, 0)
            sol  model.optimize()
            atpm_val   float(sol.objective_value or 0)
            ppa_flux   float(sol.fluxes.get(ppa_rxn.id, 0)) if ppa_rxn else 0.0
        results[cond]  {"atpm": atpm_val, "ppa_flux": ppa_flux}
        print(f"    {model_label:30s}  {cond:10s}  "
              f"ATPM  {atpm_val:.4f}   PPA_flux  {ppa_flux:.4f}")
    return results


def load_hgem():
    
    hgem_path  os.path.join(REPO, "data", "models", "Human-GEM.xml")
    print(f"\nLoading Human-GEM from {hgem_path} ...")
    model  read_sbml_model(hgem_path)
    print(f"  Loaded: {len(model.reactions)} reactions")
    return model


#main 
print("" * 60)
print("PPCOACm Rescue Analysis")
print("" * 60)

#1. load recon3d
recon_path  os.path.join(REPO, "data", "models", "cache", "Recon3D.xml")
print(f"\nLoading Recon3D from {recon_path} …")
recon_orig  read_sbml_model(recon_path)
recon_rescued  recon_orig.copy()
print(f"  Loaded: {len(recon_orig.reactions)} reactions")

#2. setup unrescued recon3d
print("\n Setting up unrescued Recon3D ")
cfg  _load_cfg()
atpm_orig, scfa_orig  _setup(recon_orig, "Recon3D (original)", cfg)
results_orig  run_doses(recon_orig, "Recon3D (original)", atpm_orig, scfa_orig)

#3. add ppcoacm and rerun
print("\n Adding PPCOACm to Recon3D ")
_add_ppcoacm(recon_rescued)
print("\n Setting up rescued Recon3D ")
atpm_resc, scfa_resc  _setup(recon_rescued, "Recon3D + PPCOACm", cfg)
results_resc  run_doses(recon_rescued, "Recon3D + PPCOACm", atpm_resc, scfa_resc)

#4. load human-gem
hgem  load_hgem()
print("\n Setting up Human-GEM ")
atpm_hgem, scfa_hgem  _setup(hgem, "Human-GEM", cfg)
results_hgem  run_doses(hgem, "Human-GEM", atpm_hgem, scfa_hgem)

#5. assemble results table
conditions  list(DOSES.keys())
rows  []
for cond in conditions:
    o    results_orig[cond]
    r    results_resc[cond]
    h    results_hgem[cond]
    rows.append({
        "Condition":               cond,
        "Recon3D_orig_PPA":        round(o["ppa_flux"], 3),
        "Recon3D_rescued_PPA":     round(r["ppa_flux"], 3),
        "Human_GEM_PPA":           round(h["ppa_flux"], 3),
        "Recon3D_orig_ATPM":       round(o["atpm"], 2),
        "Recon3D_rescued_ATPM":    round(r["atpm"], 2),
        "Human_GEM_ATPM":          round(h["atpm"], 2),
    })

df  pd.DataFrame(rows)
csv_out  os.path.join(TABS, "table_rescue.csv")
df.to_csv(csv_out, indexFalse)
print(f"\nSaved rescue table to {csv_out}")
print(df.to_string(indexFalse))

#figure 7  (two panels) 
#panel a:  propionate exchange flux magnitude  (primary rescue metric)
#panel b:  atpm comparison  (shows overshoot / secondary finding)
fig, (ax_a, ax_b)  plt.subplots(1, 2, figsize(11, 5))

#dose conditions excluding baseline (no scfa in baseline)
dose_conds  ["Low", "Mid", "High"]
x  np.arange(len(dose_conds))
w  0.25

COL_ORIG   "#3a86ff"   #blue
COL_RESC   "#06d6a0"   #teal
COL_HGEM   "#ef476f"   #red/pink

#panel a: propionate flux magnitude 
ppa_orig  [-results_orig[c]["ppa_flux"] for c in dose_conds]   #magnitude (≥0)
ppa_resc  [-results_resc[c]["ppa_flux"] for c in dose_conds]
ppa_hgem  [-results_hgem[c]["ppa_flux"] for c in dose_conds]

ax_a.bar(x - w, ppa_orig, w, colorCOL_ORIG, label"Recon3D (original)",
         edgecolor"black", linewidth0.5)
ax_a.bar(x,     ppa_resc, w, colorCOL_RESC, label"Recon3D + PPCOACm",
         edgecolor"black", linewidth0.5)
ax_a.bar(x + w, ppa_hgem, w, colorCOL_HGEM, label"Human-GEM",
         edgecolor"black", linewidth0.5)

ax_a.set_xlabel("Dose Condition", fontsize11)
ax_a.set_ylabel("Propionate uptake (mmol gDW⁻¹ hr⁻¹)", fontsize11)
ax_a.set_title("A   Propionate exchange flux", fontsize11, loc"left", fontweight"bold")
ax_a.set_xticks(x)
ax_a.set_xticklabels(dose_conds, fontsize10)
ax_a.legend(framealpha0.9, fontsize8)
ax_a.spines["top"].set_visible(False)
ax_a.spines["right"].set_visible(False)
ax_a.set_ylim(0, max(ppa_hgem) * 1.3)

#panel b: atpm 
atpm_orig_all  [results_orig[c]["atpm"] for c in conditions]
atpm_resc_all  [results_resc[c]["atpm"] for c in conditions]
atpm_hgem_all  [results_hgem[c]["atpm"] for c in conditions]
x_b  np.arange(len(conditions))

ax_b.bar(x_b - w, atpm_orig_all, w, colorCOL_ORIG, label"Recon3D (original)",
         edgecolor"black", linewidth0.5)
ax_b.bar(x_b,     atpm_resc_all, w, colorCOL_RESC, label"Recon3D + PPCOACm",
         edgecolor"black", linewidth0.5)
ax_b.bar(x_b + w, atpm_hgem_all, w, colorCOL_HGEM, label"Human-GEM",
         edgecolor"black", linewidth0.5)

#annotate overshoot arrow at high dose
hi_idx  conditions.index("High")
ax_b.annotate(
    "Overshoot:\ninternal propionyl-CoA\nhub unlocked",
    xy     (x_b[hi_idx], atpm_resc_all[hi_idx]),
    xytext (x_b[hi_idx] + 0.55, atpm_resc_all[hi_idx] * 0.85),
    fontsize7.5,
    arrowpropsdict(arrowstyle"->", color"gray", lw1.0),
    ha"left", va"center", color"gray",
)

ax_b.set_xlabel("Dose Condition", fontsize11)
ax_b.set_ylabel("Predicted ATPM (mmol gDW⁻¹ hr⁻¹)", fontsize11)
ax_b.set_title("B   Predicted ATPM flux", fontsize11, loc"left", fontweight"bold")
ax_b.set_xticks(x_b)
ax_b.set_xticklabels(conditions, fontsize10)
ax_b.legend(framealpha0.9, fontsize8)
ax_b.spines["top"].set_visible(False)
ax_b.spines["right"].set_visible(False)

plt.tight_layout()
fig_out  os.path.join(FIGS, "Fig7.png")
plt.savefig(fig_out, dpi150, bbox_inches"tight")
plt.close()
print(f"Saved rescue figure to {fig_out}")
print("\nDone.")
