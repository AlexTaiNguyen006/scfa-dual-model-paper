.PHONY: all clean inputs simulate sensitivity rescue robustness figures tables

PYTHON = python -m

all: inputs simulate sensitivity rescue robustness figures tables
	@echo "Pipeline complete."

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

figures:
	$(PYTHON) src.03_figures

tables:
	$(PYTHON) src.04_tables

pipeline-fig:
	$(PYTHON) src.make_pipeline_fig

clean:
	rm -f results/*.csv
	rm -f outputs/figs/*.png outputs/figs/*.pdf
	rm -f outputs/tables/*.csv
