# Stachys affinis SCFA modeling -- reproducibility code

Code and data for the constraint-based modeling part of the manuscript. We use
COBRApy + Recon3D to simulate how SCFAs from stachyose fermentation affect
hepatic ATP maintenance.

## Setup

You'll need conda (or mamba) and Python 3.11.

```bash
conda env create -f environment.yml
conda activate stachys-scfa
```

You also need Recon3D. Download the SBML from BiGG:
https://bigg.ucsd.edu/models/Recon3D

Put `Recon3D.xml.gz` in `data/models/`. The script will decompress it
automatically the first time you run it (and cache the result so it
doesn't have to do it again).

## Running

```bash
make all
```

This runs 4 steps in order:
1. `01_prepare_inputs.py` -- validates the SCFA dose csv
2. `02_run_simulation.py` -- loads Recon3D, sets up medium, runs FBA
   (takes a few min because loading the model is slow)
3. `03_figures.py` -- generates the 6 figures
4. `04_tables.py` -- generates the 3 csv tables

You can also run them individually:
```bash
python -m src.01_prepare_inputs
python -m src.02_run_simulation
# etc.
```

## Outputs

- `results/` -- intermediate csvs (fluxes, merged data)
- `outputs/figs/` -- the figures as PNGs
- `outputs/tables/` -- tables as csv

## SCFA dose scenarios

Three levels based on estimated S. affinis tuber intake:

| Condition | Acetate | Propionate | Butyrate |
|-----------|---------|------------|----------|
| Low (~25g)  | 2.0 | 0.7 | 0.4 |
| Mid (~50g)  | 4.0 | 1.4 | 0.8 |
| High (~100g) | 8.0 | 2.8 | 1.6 |

All in mmol/gDW/hr.

## Known issues

- Propionate shows 0 flux in every condition. I think this is a Recon3D
  thing -- the methylmalonyl-CoA pathway doesn't seem to carry flux with
  the strict boundary constraints we use. Might work better with Human-GEM
  but haven't tried yet.
- The medium setup is pretty aggressive (close all 1800+ boundary rxns,
  only reopen a curated set). This is necessary to avoid the thermodynamic
  loops that otherwise give absurd ATP yields.
- We use ATPM as the objective, not biomass, since hepatocytes don't divide.
