#!/usr/bin/env python3


import pandas as pd
from .utils import build_paths, load_config, read_scfa_inputs


def main():
    paths = build_paths()
    cfg = load_config(paths.config_path)
    conditions = cfg["project"]["conditions"]

    print("step1 validating SCFA inputs")
    df = read_scfa_inputs(paths.scfa_csv, conditions)

    #append lactate scenario
    df_lactate = df.copy()
    df_lactate["condition"] = df_lactate["condition"] + "_Lactate"
    #add 10% lactate vector relative to total SCFA moles
    total_moles = df_lactate["acetate_mmol_gDW_hr"] + df_lactate["propionate_mmol_gDW_hr"] + df_lactate["butyrate_mmol_gDW_hr"]
    df_lactate["lactate_mmol_gDW_hr"] = round(total_moles * 0.10, 3)
    
    df["lactate_mmol_gDW_hr"] = 0.0
    df = pd.concat([df, df_lactate], ignore_index=True)

    outfile = paths.results / "scfa_inputs_canonical.csv"
    df.to_csv(outfile, index=False)
    print(f"  wrote {outfile}")
    print(df.to_string(index=False))


if __name__ == "__main__":
    main()
