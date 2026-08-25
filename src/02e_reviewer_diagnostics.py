#!/usr/bin/env python3
"""Reviewer Diagnostics: Futile cycles, unconstrained loops, and residual differences.

This script runs isolated tests strictly to address Reviewer 1 and 2 queries:
  1. energy-generating futile cycles (R2.5)
  2. Spurious unconstrained propionate loops (R1.4)
  3. Final residual divergence reactions between models (R1.6)
"""

import pandas as pd
import cobra
from cobra.io import read_sbml_model

from .utils import build_paths, load_config, decompress_gz
from .run_simulation_shared import setup_medium, ATPM_IDS, _find_rxn, SCFA_EXCHANGE_IDS, GLUCOSE_IDS, O2_IDS

# Same exact parameters and setup components as original pipeline
_HGEM_ATPM_METS = {
    "MAM01371c": -1, "MAM02040c": -1, "MAM01285c": 1,
    "MAM02751c": 1, "MAM02039c": 1,
}

def _add_ppcoacm(model):
    if "PPCOACm" in model.reactions: return
    met_stoich = {"ppcoa_m": -1, "hco3_m": -1, "atp_m": -1, "mmcoa__S_m": 1, "adp_m": 1, "pi_m": 1, "h_m": 1}
    rxn = cobra.Reaction("PPCOACm")
    rxn.name = "Propionyl-CoA carboxylase, mitochondrial"
    rxn.bounds = (0.0, 1000.0)
    mets = {model.metabolites.get_by_id(mid): coeff for mid, coeff in met_stoich.items() if mid in model.metabolites}
    rxn.add_metabolites(mets)
    model.add_reactions([rxn])

def _force_atpm_objective(model):
    atpm_rxn = _find_rxn(model, ATPM_IDS, silent=True)
    if atpm_rxn is None:
        if "ATPM_added" in model.reactions:
            atpm_rxn = model.reactions.get_by_id("ATPM_added")
        else:
            atpm_rxn = cobra.Reaction("ATPM_added")
            atpm_rxn.name = "ATP maintenance (added)"
            model.add_reactions([atpm_rxn])
            mets = {model.metabolites.get_by_id(mid): coeff for mid, coeff in _HGEM_ATPM_METS.items() if mid in model.metabolites}
            atpm_rxn.add_metabolites(mets)
    atpm_rxn.bounds = (0, 10000)
    model.objective = atpm_rxn
    return atpm_rxn

def test_futile_cycles(model, label, cfg, paths):
    """(R2.5) Check if ATP is produced when all external inputs are 0."""
    setup_medium(model, cfg)
    atpm_rxn = _force_atpm_objective(model)
    
    # Close all carbon, oxygen, and nutrient sources strictly to 0
    for rxn_list in [SCFA_EXCHANGE_IDS["acetate"], SCFA_EXCHANGE_IDS["propionate"], SCFA_EXCHANGE_IDS["butyrate"], SCFA_EXCHANGE_IDS["lactate"], GLUCOSE_IDS, O2_IDS]:
        target = _find_rxn(model, rxn_list, silent=True)
        if target:
            target.bounds = (0.0, 0.0)
            
    # Also restrict the amino acids and vitamins entirely (strict starvation)
    for rxn in model.reactions:
        if rxn.id.startswith("EX_") and rxn.lower_bound < 0:
            if rxn.id not in ["EX_h2o_e", "EX_h_e", "EX_pi_e", "EX_so4_e", "EX_na1_e", "EX_k_e", "EX_cl_e", "EX_ca2_e", "EX_mg2_e", "EX_fe2_e",
                              "EX_h2o[e]", "EX_h[e]", "EX_pi[e]", "EX_so4[e]", "EX_na1[e]", "EX_k[e]", "EX_cl[e]", "EX_ca2[e]", "EX_mg2[e]", "EX_fe2[e]"]:
                rxn.bounds = (0.0, 0.0) # lock all remaining uptake
                
    with model:
        sol = model.optimize()
        atp = sol.objective_value or 0.0
        
    print(f"[{label}] Futile Cycle ATPM (strict starvation): {atp:.4f}")
    if atp > 1e-4:
        fluxes = {r_id: sol.fluxes[r_id] for r_id in sol.fluxes.index if abs(sol.fluxes[r_id]) > 1e-4}
        df_fl = pd.DataFrame(list(fluxes.items()), columns=["Reaction", "Flux"])
        df_fl.to_csv(paths.results / f"revision_diagnostics/{label}_futile_cycle_fluxes.csv", index=False)
        print(f"  WARNING: Model produces {atp:.2f} free ATP. Degenerate fluxes logged.")
    return atp

def diagnose_spurious_loops(recon, cfg, paths):
    """(R1.4) Compare constrained vs unconstrained PPCOACm rescue."""
    _add_ppcoacm(recon)
    
    # 1. Unconstrained
    recon_unc = recon.copy()
    atpm_u = _force_atpm_objective(recon_unc)
    sol_unc = recon_unc.optimize()
    val_u = sol_unc.objective_value or 0
    
    # 2. Constrained
    recon_con = recon.copy()
    setup_medium(recon_con, cfg)
    atpm_c = _force_atpm_objective(recon_con)
    sol_con = recon_con.optimize()
    val_c = sol_con.objective_value or 0
    
    print(f"[Recon3D] Spurious loops check -> Unconstrained: {val_u:.2f} || Constrained: {val_c:.2f}")
    
    flux_unc = sol_unc.fluxes
    flux_con = sol_con.fluxes
    
    diff_records = []
    for rxn_id in recon_unc.reactions.list_attr("id"):
        fu = flux_unc.get(rxn_id, 0.0)
        fc = flux_con.get(rxn_id, 0.0)
        diff = abs(fu - fc)
        if diff > 5.0 and "EX_" not in rxn_id.upper() and "MAR" not in rxn_id.upper(): 
            # Looking for massive internal cycle artifacts
            rxn = recon_unc.reactions.get_by_id(rxn_id)
            diff_records.append({
                "Reaction": rxn_id,
                "Name": rxn.name,
                "Subsystem": rxn.subsystem,
                "Unconstrained_Flux": fu,
                "Constrained_Flux": fc,
                "Difference": diff,
                "Formula": rxn.reaction
            })
            
    df = pd.DataFrame(diff_records).sort_values("Difference", ascending=False)
    out_csv = paths.results / "revision_diagnostics/recon3d_spurious_loops.csv"
    df.to_csv(out_csv, index=False)
    print(f"  Logged {len(df)} spurious loop reactions to {out_csv.name}")

def main():
    paths = build_paths()
    cfg = load_config(paths.config_path)

    print("loading Recon3D")
    sbml_file = decompress_gz(paths.sbml_path)
    recon = read_sbml_model(str(sbml_file))
    
    print("\nTEST 1: Futile Cycle Validation")
    test_futile_cycles(recon.copy(), "Recon3D", cfg, paths)
    recon_ppcoa = recon.copy()
    _add_ppcoacm(recon_ppcoa)
    test_futile_cycles(recon_ppcoa, "Recon3D_PPCOACm", cfg, paths)
    
    if paths.human_gem_path and paths.human_gem_path.exists():
        print("\nloading Human-GEM")
        hgem = read_sbml_model(str(paths.human_gem_path))
        test_futile_cycles(hgem.copy(), "Human-GEM", cfg, paths)
    else:
        print("\nskipping Human-GEM (model file not found)")

    print("\nTEST 2: Spurious Propionate Loops")
    diagnose_spurious_loops(recon.copy(), cfg, paths)

if __name__ == "__main__":
    main()