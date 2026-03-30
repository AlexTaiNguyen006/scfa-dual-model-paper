#!/usr/bin/env python3


import pandas as pd
import numpy as np
import cobra
from cobra.io import read_sbml_model
from cobra.flux_analysis import flux_variability_analysis, pfba

from .utils import build_paths, load_config, decompress_gz
from .run_simulation_shared import (
    SCFA_EXCHANGE_IDS, GLUCOSE_IDS, O2_IDS, CO2_IDS, ATPM_IDS,
    PATHWAY_RXNS, _find_rxn, setup_medium,
)

_HGEM_ATPM_METS = {
    "MAM01371c": -1, "MAM02040c": -1,
    "MAM01285c": 1, "MAM02751c": 1, "MAM02039c": 1,
}


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


def _apply_condition(scfa_rxns, ac, ppa, but):
    for name, dose in [("acetate", ac), ("propionate", ppa), ("butyrate", but)]:
        rxn = scfa_rxns.get(name)
        if rxn and dose > 0:
            rxn.bounds = (-dose, 0)


#1. scfa ratio sensitivity
RATIO_PROFILES = {
    "Stachyose-like (65:23:13)": (0.645, 0.226, 0.129),
    "Inulin-like (65:20:15)":    (0.650, 0.200, 0.150),
    "High-butyrate (55:15:30)":  (0.550, 0.150, 0.300),
    "High-propionate (50:35:15)":(0.500, 0.350, 0.150),
    "Acetate-dominated (80:10:10)":(0.800, 0.100, 0.100),
    "Equimolar (33:33:33)":      (0.333, 0.333, 0.333),
}


def run_ratio_sensitivity(model, model_label, cfg):
    
    atpm_rxn, scfa_rxns, glc_rxn, o2_rxn = _setup_model(model, model_label, cfg)
    total_scfa = 6.2  #mid dose total
    rows = []
    for profile_name, (fac, fppa, fbut) in RATIO_PROFILES.items():
        ac = total_scfa * fac
        ppa = total_scfa * fppa
        but = total_scfa * fbut
        with model:
            _apply_condition(scfa_rxns, ac, ppa, but)
            sol = model.optimize()
            atpm_val = float(sol.objective_value or 0)
        rows.append({
            "model": model_label,
            "ratio_profile": profile_name,
            "acetate_frac": fac,
            "propionate_frac": fppa,
            "butyrate_frac": fbut,
            "acetate": round(ac, 3),
            "propionate": round(ppa, 3),
            "butyrate": round(but, 3),
            "ATPM": round(atpm_val, 4),
        })
        print(f"    {profile_name}: ATPM={atpm_val:.2f}")
    return rows


#2. alternative objective (pfba)
COND_MAP = {
    "StachysDose_Low":  (2.0, 0.7, 0.4),
    "StachysDose_Mid":  (4.0, 1.4, 0.8),
    "StachysDose_High": (8.0, 2.8, 1.6),
}


def run_pfba_comparison(model, model_label, cfg):
    
    atpm_rxn, scfa_rxns, glc_rxn, o2_rxn = _setup_model(model, model_label, cfg)
    rows = []
    for cond, (ac, ppa, but) in COND_MAP.items():
        with model:
            _apply_condition(scfa_rxns, ac, ppa, but)
            #standard fba
            sol_fba = model.optimize()
            atpm_fba = float(sol_fba.objective_value or 0)
            total_flux_fba = sol_fba.fluxes.abs().sum()

            #pfba (minimize total flux while maintaining optimal atpm)
            sol_pfba = pfba(model)
            atpm_pfba = float(sol_pfba.fluxes[atpm_rxn.id])
            total_flux_pfba = sol_pfba.fluxes.abs().sum()

            #pathway fluxes under pfba
            pathway_fluxes = {}
            for canonical_id, label, group, candidates in PATHWAY_RXNS:
                rxn = _find_rxn(model, candidates, silent=True)
                if rxn:
                    pathway_fluxes[f"pfba_{canonical_id}"] = float(
                        sol_pfba.fluxes.get(rxn.id, 0))
                    pathway_fluxes[f"fba_{canonical_id}"] = float(
                        sol_fba.fluxes.get(rxn.id, 0))

        row = {
            "model": model_label,
            "condition": cond,
            "ATPM_FBA": round(atpm_fba, 4),
            "ATPM_pFBA": round(atpm_pfba, 4),
            "total_flux_FBA": round(total_flux_fba, 2),
            "total_flux_pFBA": round(total_flux_pfba, 2),
            "flux_reduction_pct": round(
                100 * (total_flux_fba - total_flux_pfba) / total_flux_fba, 1)
                if total_flux_fba > 0 else 0,
        }
        row.update(pathway_fluxes)
        rows.append(row)
        print(f"    {cond}: FBA={atpm_fba:.2f} pFBA={atpm_pfba:.2f}  "
              f"flux reduction={row['flux_reduction_pct']:.1f}%")
    return rows


#3. multi-threshold fva
def run_multi_fva(model, model_label, cfg):
    
    atpm_rxn, scfa_rxns, glc_rxn, o2_rxn = _setup_model(model, model_label, cfg)

    co2_rxn = _find_rxn(model, CO2_IDS)
    fva_targets = [atpm_rxn]
    for name in ["acetate", "propionate", "butyrate"]:
        if scfa_rxns.get(name):
            fva_targets.append(scfa_rxns[name])
    if co2_rxn:
        fva_targets.append(co2_rxn)
    if _find_rxn(model, O2_IDS):
        fva_targets.append(_find_rxn(model, O2_IDS))

    rows = []
    for cond, (ac, ppa, but) in COND_MAP.items():
        for frac in [0.90, 0.95, 1.00]:
            with model:
                _apply_condition(scfa_rxns, ac, ppa, but)
                fva_result = flux_variability_analysis(
                    model, reaction_list=fva_targets,
                    fraction_of_optimum=frac)
                for rxn_id in fva_result.index:
                    fmin = float(fva_result.loc[rxn_id, "minimum"])
                    fmax = float(fva_result.loc[rxn_id, "maximum"])
                    rows.append({
                        "model": model_label,
                        "condition": cond,
                        "fraction_optimum": frac,
                        "reaction": rxn_id,
                        "fva_minimum": round(fmin, 4),
                        "fva_maximum": round(fmax, 4),
                        "range": round(fmax - fmin, 4),
                    })
            print(f"    {cond} @ {int(frac*100)}%: done")
    return rows


#4. flux cap sensitivity
def run_flux_cap_sensitivity(model, model_label, cfg):
    
    rows = []
    for cap in [100, 200, 500]:
        with model:  #preserve original bounds across caps
            setup_medium(model, cfg)
            for rxn in model.reactions:
                if not rxn.id.startswith("EX_") and not rxn.id.startswith("MAR09"):
                    lb = max(rxn.lower_bound, -cap)
                    ub = min(rxn.upper_bound, cap)
                    rxn.bounds = (lb, ub)

            atpm_rxn = _find_rxn(model, ATPM_IDS, silent=True)
            if atpm_rxn is None:
                if "ATPM_added" in model.reactions:
                    atpm_rxn = model.reactions.get_by_id("ATPM_added")
            if atpm_rxn:
                atpm_rxn.bounds = (0, cap)
                model.objective = atpm_rxn

            scfa_rxns = {}
            for name, candidates in SCFA_EXCHANGE_IDS.items():
                scfa_rxns[name] = _find_rxn(model, candidates)

            for cond, (ac, ppa, but) in COND_MAP.items():
                with model:
                    _apply_condition(scfa_rxns, ac, ppa, but)
                    sol = model.optimize()
                    atpm_val = float(sol.objective_value or 0)
                rows.append({
                    "model": model_label,
                    "condition": cond,
                    "flux_cap": cap,
                    "ATPM": round(atpm_val, 4),
                })
                print(f"    cap={cap} {cond}: ATPM={atpm_val:.2f}")
    return rows


#main
def main():
    paths = build_paths()
    cfg = load_config(paths.config_path)

    ratio_rows, pfba_rows, fva_rows, cap_rows = [], [], [], []

    for model_path, label in [
        (paths.sbml_path, "Recon3D"),
        (paths.human_gem_path, "Human-GEM"),
    ]:
        if model_path is None or not model_path.exists():
            print(f"  {label} not found, skipping")
            continue

        if label == "Recon3D":
            model_path = decompress_gz(model_path)
        print(f"\nLoading {label} from {model_path}...")
        model = read_sbml_model(str(model_path))
        print(f"  {len(model.reactions)} reactions")

        #with model: so each analysis starts clean
        print(f"\n  [{label}] SCFA ratio sensitivity...")
        with model:
            ratio_rows.extend(run_ratio_sensitivity(model, label, cfg))

        print(f"\n  [{label}] pFBA comparison...")
        with model:
            pfba_rows.extend(run_pfba_comparison(model, label, cfg))

        print(f"\n  [{label}] Multi-threshold FVA...")
        with model:
            fva_rows.extend(run_multi_fva(model, label, cfg))

        print(f"\n  [{label}] Flux cap sensitivity...")
        with model:
            cap_rows.extend(run_flux_cap_sensitivity(model, label, cfg))

    #save results
    pd.DataFrame(ratio_rows).to_csv(
        paths.results / "ratio_sensitivity.csv", index=False)
    pd.DataFrame(pfba_rows).to_csv(
        paths.results / "pfba_comparison.csv", index=False)
    pd.DataFrame(fva_rows).to_csv(
        paths.results / "fva_multi_threshold.csv", index=False)
    pd.DataFrame(cap_rows).to_csv(
        paths.results / "flux_cap_sensitivity.csv", index=False)

    print("\n Robustness analyses complete ")
    print(f"  ratio_sensitivity.csv:     {len(ratio_rows)} rows")
    print(f"  pfba_comparison.csv:        {len(pfba_rows)} rows")
    print(f"  fva_multi_threshold.csv:    {len(fva_rows)} rows")
    print(f"  flux_cap_sensitivity.csv:   {len(cap_rows)} rows")


if __name__ == "__main__":
    main()

