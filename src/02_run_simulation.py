#!/usr/bin/env python3
"""Run FBA on Recon3D and Human-GEM with SCFA dose conditions.

Outputs
  results/host_fluxes_by_condition.csv  (both models)
  results/merged_dose_scfa_host.csv     (Recon3D only, for backward compat)
  results/model_comparison.csv          (side-by-side ATPM + propionate)
  results/fva_ranges.csv                (99% optimality FVA)
"""

import pandas as pd
import numpy as np
import cobra
from cobra.io import read_sbml_model
from cobra.flux_analysis import flux_variability_analysis

from .utils import build_paths, load_config, decompress_gz, read_scfa_inputs
from .run_simulation_shared import (
    SCFA_EXCHANGE_IDS, GLUCOSE_IDS, O2_IDS, CO2_IDS, ATPM_IDS,
    PATHWAY_RXNS, _find_rxn, setup_medium,
)

# Human-GEM ATPM stoichiometry (for adding demand reaction)
_HGEM_ATPM_METS = {
    "MAM01371c": -1, "MAM02040c": -1,
    "MAM01285c": 1,  "MAM02751c": 1,  "MAM02039c": 1,
}


def collect_pathway_fluxes(model, solution):
    result = {}
    for canonical_id, label, group, candidates in PATHWAY_RXNS:
        rxn = _find_rxn(model, candidates, silent=True)
        if rxn:
            val = solution.fluxes.get(rxn.id, float("nan"))
            result[f"pathway_{canonical_id}"] = float(val)
        else:
            result[f"pathway_{canonical_id}"] = float("nan")
    return result


def _setup_model(model, model_label, cfg):
    """Apply hepatocyte medium, set ATPM objective, locate SCFA exchanges."""
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
        found_id = scfa_rxns[name].id if scfa_rxns[name] else "NOT FOUND"
        print(f"  {name}: {found_id}")

    glc_rxn = _find_rxn(model, GLUCOSE_IDS)
    o2_rxn = _find_rxn(model, O2_IDS)
    co2_rxn = _find_rxn(model, CO2_IDS, silent=True)
    return atpm_rxn, scfa_rxns, glc_rxn, o2_rxn, co2_rxn


def _run_model(model, model_label, scfa_df, cfg):
    """Run baseline + dose conditions for one model; return list of dicts."""
    atpm_rxn, scfa_rxns, glc_rxn, o2_rxn, co2_rxn = _setup_model(
        model, model_label, cfg)

    # baseline (no SCFAs)
    print(f"\n  [{model_label}] baseline...")
    baseline = model.optimize()
    baseline_atpm = float(baseline.objective_value or 0)
    print(f"    baseline ATPM = {baseline_atpm:.4f}")

    rows = []
    for _, row in scfa_df.iterrows():
        cond = row["condition"]
        ac  = float(row["acetate_mmol_gDW_hr"])
        ppa = float(row["propionate_mmol_gDW_hr"])
        but = float(row["butyrate_mmol_gDW_hr"])

        with model:
            for name, dose in [("acetate", ac), ("propionate", ppa),
                               ("butyrate", but)]:
                rxn = scfa_rxns.get(name)
                if rxn and dose > 0:
                    rxn.bounds = (-dose, 0)

            sol = model.optimize()
            atpm_val = float(sol.objective_value or 0)
            delta = atpm_val - baseline_atpm
            pct = 100.0 * delta / baseline_atpm if baseline_atpm else float("nan")

            def flux(r):
                return float(sol.fluxes.get(r.id, 0)) if r else None

            result = {
                "model": model_label,
                "condition": cond,
                "objective_id": atpm_rxn.id,
                "objective_value": atpm_val,
                "baseline_objective": baseline_atpm,
                "objective_delta": delta,
                "objective_pct_change": pct,
                "glucose_flux": flux(glc_rxn),
                "oxygen_flux": flux(o2_rxn),
                "co2_flux": flux(co2_rxn),
                "acetate_flux": flux(scfa_rxns.get("acetate")),
                "propionate_flux": flux(scfa_rxns.get("propionate")),
                "butyrate_flux": flux(scfa_rxns.get("butyrate")),
                "solver": str(sol.solver),
            }
            result.update(collect_pathway_fluxes(model, sol))
            rows.append(result)
            print(f"    {cond}: ATPM={atpm_val:.2f} ({pct:+.1f}%)")

    return rows, baseline_atpm, atpm_rxn, scfa_rxns, co2_rxn, o2_rxn


def _run_fva(model, model_label, scfa_df, atpm_rxn, scfa_rxns, co2_rxn, o2_rxn):
    """Run FVA at 99% optimality for ATPM and key exchanges."""
    targets = [atpm_rxn]
    for name in ["acetate", "propionate", "butyrate"]:
        if scfa_rxns.get(name):
            targets.append(scfa_rxns[name])
    if co2_rxn:
        targets.append(co2_rxn)
    if o2_rxn:
        targets.append(o2_rxn)

    rows = []
    for _, row in scfa_df.iterrows():
        cond = row["condition"]
        ac  = float(row["acetate_mmol_gDW_hr"])
        ppa = float(row["propionate_mmol_gDW_hr"])
        but = float(row["butyrate_mmol_gDW_hr"])
        with model:
            for name, dose in [("acetate", ac), ("propionate", ppa),
                               ("butyrate", but)]:
                rxn = scfa_rxns.get(name)
                if rxn and dose > 0:
                    rxn.bounds = (-dose, 0)
            fva = flux_variability_analysis(
                model, reaction_list=targets, fraction_of_optimum=0.99)
            for rxn_id in fva.index:
                fmin = float(fva.loc[rxn_id, "minimum"])
                fmax = float(fva.loc[rxn_id, "maximum"])
                rows.append({
                    "model": model_label,
                    "condition": cond,
                    "reaction": rxn_id,
                    "fva_minimum": round(fmin, 4),
                    "fva_maximum": round(fmax, 4),
                })
        print(f"    [{model_label}] FVA {cond}: done")
    return rows


def main():
    paths = build_paths()
    cfg = load_config(paths.config_path)
    conditions = cfg["project"]["conditions"]
    scfa_df = read_scfa_inputs(
        paths.results / "scfa_inputs_canonical.csv", conditions)

    all_rows = []
    fva_rows = []

    # Recon3D
    sbml_file = decompress_gz(paths.sbml_path)
    print(f"loading Recon3D from {sbml_file}")
    recon = read_sbml_model(str(sbml_file))
    print(f"  {len(recon.reactions)} rxns, {len(recon.metabolites)} mets")
    rows_r, bl_r, atpm_r, scfa_r, co2_r, o2_r = _run_model(
        recon, "Recon3D", scfa_df, cfg)
    all_rows.extend(rows_r)
    fva_rows.extend(_run_fva(recon, "Recon3D", scfa_df,
                             atpm_r, scfa_r, co2_r, o2_r))

    # Human-GEM
    hgem_path = paths.human_gem_path
    if hgem_path and hgem_path.exists():
        print(f"\nloading Human-GEM from {hgem_path}")
        hgem = read_sbml_model(str(hgem_path))
        print(f"  {len(hgem.reactions)} rxns, {len(hgem.metabolites)} mets")
        rows_h, bl_h, atpm_h, scfa_h, co2_h, o2_h = _run_model(
            hgem, "Human-GEM", scfa_df, cfg)
        all_rows.extend(rows_h)
        fva_rows.extend(_run_fva(hgem, "Human-GEM", scfa_df,
                                 atpm_h, scfa_h, co2_h, o2_h))

    # save host_fluxes_by_condition.csv (both models)
    host_df = pd.DataFrame(all_rows)
    host_df.to_csv(paths.results / "host_fluxes_by_condition.csv", index=False)

    # merged: Recon3D only (backward compat for 03_figures.py input data)
    recon_df = host_df[host_df["model"] == "Recon3D"].drop(columns="model")
    merged = scfa_df.merge(recon_df, on="condition", how="left")
    merged.to_csv(paths.results / "merged_dose_scfa_host.csv", index=False)

    # model_comparison.csv
    mc = host_df[["model", "condition", "objective_value",
                   "objective_pct_change", "propionate_flux"]].copy()
    mc.to_csv(paths.results / "model_comparison.csv", index=False)

    # fva_ranges.csv
    pd.DataFrame(fva_rows).to_csv(
        paths.results / "fva_ranges.csv", index=False)

    print("\ndone")
    for _, r in host_df.iterrows():
        print(f"  {r['model']:12s} {r['condition']:25s} "
              f"ATPM={r['objective_value']:.2f}")


if __name__ == "__main__":
    main()
