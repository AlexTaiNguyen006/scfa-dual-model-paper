#!/usr/bin/env python3


import pandas as pd
import numpy as np
import cobra
from cobra.io import read_sbml_model

from .utils import build_paths, load_config, decompress_gz
from .run_simulation_shared import (
    SCFA_EXCHANGE_IDS, GLUCOSE_IDS, O2_IDS, ATPM_IDS,
    _find_rxn, setup_medium,
)

#human-gem atpm metabolite ids (same as 02_run_simulation.py)
_HGEM_ATPM_METS = {
    "MAM01371c": -1,   #atp
    "MAM02040c": -1,   #h2o
    "MAM01285c": 1,    #adp
    "MAM02751c": 1,    #pi
    "MAM02039c": 1,    #h+
}


def _run_single(model, atpm_rxn, scfa_rxns, glc_rxn, o2_rxn,
                ac, ppa, but, glc_ub, o2_ub):
    
    with model:
        if glc_rxn:
            glc_rxn.bounds = (-glc_ub, 0)
        if o2_rxn:
            o2_rxn.bounds = (-o2_ub, 0)

        for name, dose in [("acetate", ac), ("propionate", ppa), ("butyrate", but)]:
            rxn = scfa_rxns.get(name)
            if rxn and dose > 0:
                rxn.bounds = (-dose, 0)

        sol = model.optimize()
        return float(sol.objective_value or 0)


def _setup_model(model, model_label, cfg):
    
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
            for met_id, coeff in _HGEM_ATPM_METS.items():
                if met_id in model.metabolites:
                    mets[model.metabolites.get_by_id(met_id)] = coeff
            atpm_rxn.add_metabolites(mets)

    atpm_rxn.bounds = (0, 500)
    model.objective = atpm_rxn

    scfa_rxns = {}
    for name, candidates in SCFA_EXCHANGE_IDS.items():
        scfa_rxns[name] = _find_rxn(model, candidates)

    glc_rxn = _find_rxn(model, GLUCOSE_IDS)
    o2_rxn = _find_rxn(model, O2_IDS)
    return atpm_rxn, scfa_rxns, glc_rxn, o2_rxn


def _run_sweeps(model, model_label, cfg, sens_cfg):
    
    sim_cfg = cfg["host_simulation"]
    default_glc = float(sim_cfg["glucose_uptake"])
    default_o2 = float(sim_cfg["oxygen_uptake"])
    default_ac = 4.0
    default_ppa = 1.4
    default_but = 0.8

    ac_range = sens_cfg.get("acetate_range", [0.5, 1.0, 2.0, 4.0, 8.0, 12.0])
    ppa_range = sens_cfg.get("propionate_range", [0.0, 0.7, 1.4, 2.8, 5.0])
    but_range = sens_cfg.get("butyrate_range", [0.0, 0.4, 0.8, 1.6, 3.0])
    glc_range = sens_cfg.get("glucose_range", [0.1, 0.5, 1.0, 2.0, 5.0])
    o2_range = sens_cfg.get("oxygen_range", [10.0, 50.0, 100.0, 200.0])

    atpm_rxn, scfa_rxns, glc_rxn, o2_rxn = _setup_model(model, model_label, cfg)

    rows = []

    def _sweep(sweep_name, values, ac_f, ppa_f, but_f, glc_f, o2_f):
        print(f"  [{model_label}] sweeping {sweep_name}...")
        for val in values:
            ac, ppa, but, glc, o2 = ac_f(val), ppa_f(val), but_f(val), glc_f(val), o2_f(val)
            v = _run_single(model, atpm_rxn, scfa_rxns, glc_rxn, o2_rxn,
                             ac, ppa, but, glc, o2)
            rows.append({"model": model_label, "sweep": sweep_name,
                          "sweep_value": val,
                          "acetate": ac, "propionate": ppa, "butyrate": but,
                          "glucose": glc, "oxygen": o2, "ATPM": v})

    _sweep("acetate", ac_range,
           lambda v: v, lambda v: default_ppa, lambda v: default_but,
           lambda v: default_glc, lambda v: default_o2)
    _sweep("propionate", ppa_range,
           lambda v: default_ac, lambda v: v, lambda v: default_but,
           lambda v: default_glc, lambda v: default_o2)
    _sweep("butyrate", but_range,
           lambda v: default_ac, lambda v: default_ppa, lambda v: v,
           lambda v: default_glc, lambda v: default_o2)
    _sweep("glucose", glc_range,
           lambda v: default_ac, lambda v: default_ppa, lambda v: default_but,
           lambda v: v, lambda v: default_o2)
    _sweep("oxygen", o2_range,
           lambda v: default_ac, lambda v: default_ppa, lambda v: default_but,
           lambda v: default_glc, lambda v: v)

    #ratio sweep
    total_scfa = default_ac + default_but
    ratio_vals = list(np.linspace(0.1, 0.9, 9))
    print(f"  [{model_label}] sweeping acetate:butyrate ratio...")
    for frac_ac in ratio_vals:
        ac = total_scfa * frac_ac
        but = total_scfa * (1 - frac_ac)
        v = _run_single(model, atpm_rxn, scfa_rxns, glc_rxn, o2_rxn,
                         ac, default_ppa, but, default_glc, default_o2)
        rows.append({"model": model_label, "sweep": "ac_but_ratio",
                      "sweep_value": round(frac_ac, 2),
                      "acetate": round(ac, 3), "propionate": default_ppa,
                      "butyrate": round(but, 3),
                      "glucose": default_glc, "oxygen": default_o2, "ATPM": v})
    return rows


def main():
    paths = build_paths()
    cfg = load_config(paths.config_path)
    sens_cfg = cfg.get("sensitivity", {})

    all_rows = []

    #recon3d
    sbml_file = decompress_gz(paths.sbml_path)
    print(f"loading Recon3D from {sbml_file}...")
    recon = read_sbml_model(str(sbml_file))
    print(f"  {len(recon.reactions)} rxns")
    all_rows.extend(_run_sweeps(recon, "Recon3D", cfg, sens_cfg))

    #human-gem
    hgem_path = paths.human_gem_path
    if hgem_path and hgem_path.exists():
        print(f"loading Human-GEM from {hgem_path}...")
        hgem = read_sbml_model(str(hgem_path))
        print(f"  {len(hgem.reactions)} rxns")
        all_rows.extend(_run_sweeps(hgem, "Human-GEM", cfg, sens_cfg))

    df = pd.DataFrame(all_rows)
    outfile = paths.results / "sensitivity_analysis.csv"
    df.to_csv(outfile, index=False)
    print(f"\nwrote {len(df)} conditions to {outfile}")
    for model_name in df["model"].unique():
        sub = df[df["model"] == model_name]
        print(f"\n  {model_name}: {len(sub)} conditions")
        print(sub.groupby("sweep")[["sweep_value", "ATPM"]].apply(
            lambda g: g.to_string(index=False)).to_string())


if __name__ == "__main__":
    main()
