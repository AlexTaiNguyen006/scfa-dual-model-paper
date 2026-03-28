#!/usr/bin/env python3
"""Generate all manuscript tables (Tables 1–6) as formatted CSVs."""

import pandas as pd
from .utils import build_paths

DOSE_ORDER = {
    "StachysDose_High": 0,
    "StachysDose_Mid": 1,
    "StachysDose_Low": 2,
}


def main():
    paths = build_paths()
    df = pd.read_csv(paths.results / "merged_dose_scfa_host.csv")

    # sort
    df["_s"] = df["condition"].map(DOSE_ORDER)
    df = df.sort_values("_s").drop(columns="_s").reset_index(drop=True)

    # ── Table 1: SCFA input vectors ──────────────────────────────
    t1 = df[["condition", "acetate_mmol_gDW_hr", "propionate_mmol_gDW_hr",
             "butyrate_mmol_gDW_hr"]].copy()
    t1.columns = ["Condition", "Acetate (mmol/gDW/hr)",
                   "Propionate (mmol/gDW/hr)", "Butyrate (mmol/gDW/hr)"]
    t1.to_csv(paths.tables_dir / "table1_scfa_vectors.csv", index=False)
    print("  table1_scfa_vectors.csv")

    # ── Table 2a: ATPM values ────────────────────────────────────
    atpm_cols = ["condition", "objective_value", "baseline_objective",
                 "objective_delta", "objective_pct_change"]
    atpm_cols = [c for c in atpm_cols if c in df.columns]
    t2a = df[atpm_cols].copy()
    t2a.to_csv(paths.tables_dir / "table2a_atpm_values.csv", index=False)
    print("  table2a_atpm_values.csv")

    # ── Table 2b: Exchange fluxes ────────────────────────────────
    flux_cols = ["condition", "glucose_flux", "oxygen_flux", "co2_flux",
                 "acetate_flux", "propionate_flux", "butyrate_flux"]
    flux_cols = [c for c in flux_cols if c in df.columns]
    t2b = df[flux_cols].copy()
    t2b.to_csv(paths.tables_dir / "table2b_exchange_fluxes.csv", index=False)
    print("  table2b_exchange_fluxes.csv")

    # ── Table 3: Summary statistics ──────────────────────────────
    summary_cols = ["condition", "objective_value"]
    for scfa in ["acetate", "propionate", "butyrate"]:
        for suf in ["_mmol_gDW_hr", "_flux"]:
            c = scfa + suf
            if c in df.columns:
                summary_cols.append(c)
    t3 = df[[c for c in summary_cols if c in df.columns]]
    t3.to_csv(paths.tables_dir / "table3_summary.csv", index=False)
    print("  table3_summary.csv")

    # ── Table 4: FVA ranges at 99% optimality ────────────────────
    fva_path = paths.results / "fva_ranges.csv"
    if fva_path.exists():
        fva = pd.read_csv(fva_path)
        fva["range"] = fva["fva_maximum"] - fva["fva_minimum"]
        fva = fva.round(4)
        fva.to_csv(paths.tables_dir / "table4_fva_ranges.csv", index=False)
        print("  table4_fva_ranges.csv")

    # ── Table 5: Propionate rescue results ───────────────────────
    rescue_path = paths.results / "table_rescue_constrained.csv"
    if rescue_path.exists():
        rescue = pd.read_csv(rescue_path)
        rescue.to_csv(paths.tables_dir / "table5_propionate_rescue.csv",
                       index=False)
        print("  table5_propionate_rescue.csv")

    # ── Table 6: Predicted vs published ATP yields ───────────────
    # Assemble from model_comparison.csv + literature reference values
    mc_path = paths.results / "model_comparison.csv"
    if mc_path.exists():
        mc = pd.read_csv(mc_path)
        # Literature ATP yields (Rolfe & Brown 1997; Clausen & Mortensen 1994)
        lit = pd.DataFrame([
            {"Source": "Rolfe & Brown (1997)", "Metric": "Basal hepatocyte ATPM",
             "Value_mmol_gDW_hr": 34.5, "Notes": "Isolated rat hepatocytes"},
            {"Source": "Clausen & Mortensen (1994)", "Metric": "Butyrate ATP yield",
             "Value_mmol_gDW_hr": None,
             "Notes": "Rank: butyrate > acetate per mole (qualitative)"},
        ])
        # Model predictions summary
        pred = mc.groupby("model").agg(
            baseline_ATPM=("objective_value", "min"),
            max_ATPM=("objective_value", "max"),
        ).reset_index()
        pred.columns = ["Model", "Min ATPM (Low dose)", "Max ATPM (High dose)"]

        pred.to_csv(paths.tables_dir / "table6_atp_yield_comparison.csv",
                     index=False)
        lit.to_csv(paths.tables_dir / "table6_literature_reference.csv",
                    index=False)
        print("  table6_atp_yield_comparison.csv")
        print("  table6_literature_reference.csv")

    print("All manuscript tables generated.")


if __name__ == "__main__":
    main()

