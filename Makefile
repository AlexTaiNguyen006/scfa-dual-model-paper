SHELL := /bin/bash

.PHONY: all clean

all:
	python -m src.01_prepare_inputs
	python -m src.02_run_simulation
	python -m src.02a_microbiome_integration
	python -m src.02b_sensitivity
	python src/02c_rescue.py
	python -m src.02d_robustness
	python src/03_figures.py
	python -m src.04_tables

clean:
	rm -rf results/*.csv outputs/figs/*.png outputs/tables/*.csv
