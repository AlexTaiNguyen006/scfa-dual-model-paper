# Reproducibility Submission Checklist

This checklist is intended for peer reviewers and editors.

Archived submission snapshot DOI: https://doi.org/10.5281/zenodo.19293353

## 1. Environment

- Python: 3.11
- Solver: GLPK (via COBRApy/optlang)
- Environment lock: `environment.yml`

Create environment:

```bash
conda env create -f environment.yml
conda activate stachys-scfa
```

## 2. Required model files

Place the following files before execution:

- `data/models/Recon3D.xml.gz`
- `data/models/Human-GEM.xml`

## 3. Pipeline execution

Core manuscript pipeline:

```bash
make all
```

Optional microbiome cross-validation (Step 2a):

```bash
make microbiome
```

Single-command full run including Step 2a:

```bash
make all-with-microbiome
```

## 4. Expected core outputs

Tables:

- `outputs/tables/table1_scfa_vectors.csv`
- `outputs/tables/table2_dual_model_atpm.csv`
- `outputs/tables/table3_exchange_fluxes.csv`
- `outputs/tables/table4_fva_ranges.csv`
- `outputs/tables/table5_propionate_rescue.csv`
- `outputs/tables/table6_atp_yield_comparison.csv`

Figures:

- `outputs/figs/Fig1.png` ... `outputs/figs/Fig10.png`
- `outputs/figs/FigS1.png`

## 5. Quick numerical sanity checks

`table2_dual_model_atpm.csv` should include:

- Recon3D High ATPM: `133.91`
- Human-GEM High ATPM: `178.77`
- Recon3D Propionate Flux (all dose conditions): `0.0`
- Human-GEM Propionate Flux High/Mid/Low: `-2.8`, `-1.4`, `-0.7`

`table5_propionate_rescue.csv` should include convergence:

- Low: `+0.3%`
- Mid: `+0.6%`
- High: `+0.7%`

## 6. Runtime expectations

On Apple Silicon (16 GB RAM), the full core pipeline typically completes in a few minutes.

## 7. Notes for manuscript alignment

- Manuscript Tables 1-6 map directly to `outputs/tables/table1` to `table6` files listed above.
- Step 2a (microbiome integration) is optional by design and does not affect core host-model results.
