.PHONY: all all-with-microbiome clean inputs simulate sensitivity rescue robustness diagnostics figures tables microbiome
PYTHON = python -m
all: inputs simulate sensitivity rescue robustness diagnostics figures tables
	@echo "Pipeline complete."
all-with-microbiome: all microbiome
	@echo "Pipeline complete (including microbiome cross-validation)."
inputs:
	$(PYTHON) src.01_prepare_inputs
simulate:
	$(PYTHON) src.02_run_simulation
microbiome:
	$(PYTHON) src.02a_microbiome_integration
sensitivity:
	$(PYTHON) src.02b_sensitivity
rescue:
	$(PYTHON) src.02c_rescue
robustness:
	$(PYTHON) src.02d_robustness
diagnostics:
	$(PYTHON) src.02e_reviewer_diagnostics
figures:
	$(PYTHON) src.03_figures
tables:
	$(PYTHON) src.04_tables
pipeline-fig:
	$(PYTHON) src.make_pipeline_fig
clean:
	rm -f results/*.csv
	rm -rf results/revision_diagnostics/
	rm -f outputs/figs/*.png outputs/figs/*.pdf
	rm -f outputs/tables/*.csv
