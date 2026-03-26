#!/usr/bin/env python3


import pandas as pd
import numpy as np
import cobra
from cobra.io import read_sbml_model


#reaction id candidates (recon3d + human-gem ids)
SCFA_EXCHANGE_IDS  {
    "acetate":    ["EX_ac_e", "EX_ac[e]", "EX_ac(e)", "MAR09086"],
    "propionate": ["EX_ppa_e", "EX_propn_e", "EX_ppn_e", "EX_ppa[e]", "MAR09808"],
    "butyrate":   ["EX_but_e", "EX_btn_e", "EX_but[e]", "MAR09809"],
}
GLUCOSE_IDS  ["EX_glc__D_e", "EX_glc_D_e", "EX_glc[e]", "MAR09034"]
O2_IDS       ["EX_o2_e", "EX_o2[e]", "MAR09048"]
CO2_IDS      ["EX_co2_e", "EX_co2[e]", "MAR09058"]
ATPM_IDS     ["ATPM", "DM_atp_c_"]

#medium components (recon3d ex_ + human-gem mar ids)
FREE_EXCHANGE  [
    "EX_h2o_e", "EX_h_e", "EX_h2o[e]", "EX_h[e]",
    "MAR09047", "MAR09079",  #human-gem: h2o, h+
]

IONS  {
    "EX_pi_e": 10, "EX_so4_e": 10, "EX_k_e": 10, "EX_na1_e": 10,
    "EX_ca2_e": 10, "EX_cl_e": 10, "EX_mg2_e": 10, "EX_fe2_e": 10,
    "EX_pi[e]": 10, "EX_so4[e]": 10, "EX_k[e]": 10, "EX_na1[e]": 10,
    "EX_ca2[e]": 10, "EX_cl[e]": 10, "EX_mg2[e]": 10, "EX_fe2[e]": 10,
    #human-gem ions
    "MAR09072": 10, "MAR09074": 10, "MAR09076": 10, "MAR09077": 10,
    "MAR09081": 10, "MAR09082": 10, "MAR09096": 10, "MAR09150": 10,
    "MAR13072": 10,
}

SECRETION_RXNS  {
    "EX_co2_e": 1000, "EX_nh4_e": 100,
    "EX_co2[e]": 1000, "EX_nh4[e]": 100,
    "MAR09058": 1000, "MAR11420": 100,  #human-gem: co2, nh4+
}

ESSENTIAL_AA  [
    "EX_his__L_e", "EX_ile__L_e", "EX_leu__L_e", "EX_lys__L_e",
    "EX_met__L_e", "EX_phe__L_e", "EX_thr__L_e", "EX_trp__L_e",
    "EX_val__L_e",
    "EX_his__L[e]", "EX_ile__L[e]", "EX_leu__L[e]", "EX_lys__L[e]",
    "EX_met__L[e]", "EX_phe__L[e]", "EX_thr__L[e]", "EX_trp__L[e]",
    "EX_val__L[e]",
    #human-gem
    "MAR09038", "MAR09039", "MAR09040", "MAR09041",  #his, ile, leu, lys
    "MAR09042", "MAR09043", "MAR09044", "MAR09045",  #met, phe, thr, trp
    "MAR09046",  #val
]

VITAMINS  [
    "EX_thm_e", "EX_ribflv_e", "EX_ncam_e", "EX_pnto__R_e",
    "EX_pydxn_e", "EX_fol_e", "EX_cbl1_e", "EX_chol_e", "EX_inost_e",
    "EX_thm[e]", "EX_ribflv[e]", "EX_ncam[e]", "EX_pnto__R[e]",
    "EX_pydxn[e]", "EX_fol[e]", "EX_cbl1[e]", "EX_chol[e]", "EX_inost[e]",
    #human-gem vitamins
    "MAR09159",  #thiamin
    "MAR09143",  #riboflavin
    "MAR09378",  #nicotinamide
    "MAR09145",  #pantothenate
    "MAR09144",  #pyridoxine
    "MAR09146",  #folate
    "MAR09083",  #choline
    "MAR09361",  #inositol
    "MAR09109",  #cobalamin (b12) — human-gem
]

#(canonical_id, label, group, [candidate_reaction_ids])
#candidate list: first match found in the model is used
PATHWAY_RXNS  [
    ("PYK",         "Pyruvate kinase",        "Glycolysis",
     ["PYK", "MAR04358"]),
    ("PDHm",        "Pyruvate dehydrogenase", "Glycolysis",
     ["PDHm", "MAR06412"]),
    ("CSm",         "Citrate synthase",       "TCA",
     ["CSm", "MAR04145"]),
    ("ACONTm",      "Aconitase",              "TCA",
     ["ACONTm", "MAR04454"]),
    ("ICDHxm",      "Isocitrate DH",          "TCA",
     ["ICDHxm", "MAR03957"]),
    ("AKGDm",       "α-Ketoglutarate DH",     "TCA",
     ["AKGDm", "MAR06414"]),
    ("SUCDi",       "Succinate DH",           "TCA",
     ["SUCDi", "MAR04652"]),
    ("FUMm",        "Fumarase",               "TCA",
     ["FUMm", "MAR04408"]),
    ("MDHm",        "Malate DH",              "TCA",
     ["MDHm", "MAR04139"]),
    ("PCm",         "Pyruvate carboxylase",   "Gluconeogenesis",
     ["PCm", "MAR04143"]),
    ("PEPCK",       "PEPCK",                  "Gluconeogenesis",
     ["PEPCK", "MAR04101"]),
    ("G6PDH2r",     "G6PD",                   "Pentose phosphate",
     ["G6PDH2r", "MAR04304"]),
    ("FBA",         "Aldolase",               "Glycolysis",
     ["FBA", "MAR04375"]),
    ("PFK",         "PFK",                    "Glycolysis",
     ["PFK", "MAR04379"]),
    ("ATPS4mi",     "ATP synthase",           "ETC",
     ["ATPS4mi", "MAR06918"]),
    ("NADH2_u10mi", "Complex I",              "ETC",
     ["NADH2_u10mi", "MAR06921"]),
    ("CYOOm3i",     "Complex IV",             "ETC",
     ["CYOOm3i", "MAR06914"]),
]


def _find_rxn(model, id_list, silentFalse):
    
    for rid in id_list:
        if rid in model.reactions:
            return model.reactions.get_by_id(rid)
    if not silent:
        import sys
        print(f"  WARNING: none of {id_list} found in model", filesys.stderr)
    return None


def setup_medium(model, cfg):
    
    sim_cfg  cfg["host_simulation"]

    #identify boundary (exchange/demand/sink) reactions model-agnostically
    boundary_ids  set()
    for rxn in model.boundary:
        boundary_ids.add(rxn.id)
    #also include anything with ex_/dm_/sink_/sk_ prefix
    for rxn in model.reactions:
        if rxn.id.startswith(("EX_", "DM_", "sink_", "SK_")):
            boundary_ids.add(rxn.id)

    #cap internal reaction bounds
    for rxn in model.reactions:
        if rxn.id in boundary_ids:
            continue
        if rxn.lower_bound < -500:
            rxn.lower_bound  -500
        if rxn.upper_bound > 500:
            rxn.upper_bound  500

    #close all boundary reactions
    for rid in boundary_ids:
        model.reactions.get_by_id(rid).bounds  (0, 0)

    for rid in FREE_EXCHANGE:
        if rid in model.reactions:
            model.reactions.get_by_id(rid).bounds  (-1000, 1000)

    for rid, ub in IONS.items():
        if rid in model.reactions:
            model.reactions.get_by_id(rid).bounds  (-ub, 100)

    for rid, ub in SECRETION_RXNS.items():
        if rid in model.reactions:
            model.reactions.get_by_id(rid).bounds  (0, ub)

    for rid in ["EX_nh4_e", "EX_nh4[e]", "MAR11420"]:
        if rid in model.reactions:
            model.reactions.get_by_id(rid).lower_bound  -0.5

    o2_uptake  float(sim_cfg["oxygen_uptake"])
    o2_rxn  _find_rxn(model, O2_IDS)
    if o2_rxn:
        o2_rxn.bounds  (-o2_uptake, 0)

    glc_uptake  float(sim_cfg["glucose_uptake"])
    glc_rxn  _find_rxn(model, GLUCOSE_IDS)
    if glc_rxn:
        glc_rxn.bounds  (-glc_uptake, 0)

    aa_rate  float(sim_cfg.get("amino_acid_uptake", 0.01))
    for rid in ESSENTIAL_AA:
        if rid in model.reactions:
            model.reactions.get_by_id(rid).bounds  (-aa_rate, 0)

    vit_rate  float(sim_cfg.get("vitamin_uptake", 0.01))
    for rid in VITAMINS:
        if rid in model.reactions:
            model.reactions.get_by_id(rid).bounds  (-vit_rate, 0)

