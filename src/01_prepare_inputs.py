#!/usr/bin/env python3


import pandas as pd
from .utils import build_paths, load_config, read_scfa_inputs


def main():




    paths = build_paths()
    cfg = load_config(paths.config_path)
    conditions = cfg["project"]["conditions"]

    print("step1 validating SCFA inputs")
    df = read_scfa_inputs(paths.scfa_csv, conditions)

    outfile = paths.results / "scfa_inputs_canonical.csv"
    df.to_csv(outfile, index=False)
    print(f"  wrote {outfile}")
    print(df.to_string(index=False))


if __name__ == "__main__":
    main()
