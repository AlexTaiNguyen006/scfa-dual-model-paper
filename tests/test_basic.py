"""Minimal unit tests for the stachys-affinis SCFA pipeline."""

import os
import sys
from pathlib import Path

import pytest
import pandas as pd

# Ensure src is importable
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))


# ---------- Config & paths ---------------------------------------------------

def test_config_loads():
    from src.utils import load_config
    cfg = load_config(ROOT / "data" / "inputs" / "project_config.yml")
    assert "human_model" in cfg
    assert "host_simulation" in cfg
    assert cfg["host_simulation"]["objective"]["prefer_atp_maintenance"] is True


def test_build_paths():
    from src.utils import build_paths
    paths = build_paths()
    assert paths.root == ROOT
    assert paths.scfa_csv.exists()
    assert paths.sbml_path.suffix in (".xml", ".gz")


# ---------- SCFA inputs ------------------------------------------------------

def test_scfa_inputs_valid():
    df = pd.read_csv(ROOT / "data" / "inputs" / "scfa_inputs.csv")
    assert "condition" in df.columns
    for col in ["acetate_mmol_gDW_hr", "propionate_mmol_gDW_hr", "butyrate_mmol_gDW_hr"]:
        assert col in df.columns
        assert (df[col] >= 0).all(), f"Negative values in {col}"


# ---------- Shared constants --------------------------------------------------

def test_pathway_rxns_format():
    from src.run_simulation_shared import PATHWAY_RXNS
    for entry in PATHWAY_RXNS:
        assert len(entry) == 4, f"Expected 4-tuple, got {len(entry)}: {entry}"
        canonical, label, group, candidates = entry
        assert isinstance(candidates, list) and len(candidates) >= 1


def test_find_rxn_returns_none_for_missing():
    from src.run_simulation_shared import _find_rxn
    # Use a minimal mock-like model with no reactions
    class FakeModel:
        class _rxns:
            @staticmethod
            def __contains__(x):
                return False
        reactions = _rxns()
    assert _find_rxn(FakeModel(), ["NONEXISTENT"], silent=True) is None


# ---------- AGORA2 data -------------------------------------------------------

def test_agora2_csv_loads():
    path = ROOT / "data" / "inputs" / "agora2_community_scfa.csv"
    if not path.exists():
        pytest.skip("agora2_community_scfa.csv not present")
    df = pd.read_csv(path, comment="#")
    assert "total_scfa" in df.columns
    assert len(df) >= 5


# ---------- Results (only when pipeline has run) ------------------------------

@pytest.fixture
def results_dir():
    d = ROOT / "results"
    if not d.exists() or not any(d.iterdir()):
        pytest.skip("No results yet — run the pipeline first")
    return d


def test_host_fluxes_csv(results_dir):
    csv = results_dir / "host_fluxes_by_condition.csv"
    if not csv.exists():
        pytest.skip("host_fluxes_by_condition.csv not generated yet")
    df = pd.read_csv(csv)
    assert "model" in df.columns
    assert "condition" in df.columns
    assert "objective_value" in df.columns
    assert (df["objective_value"] > 0).all()


def test_sensitivity_csv(results_dir):
    csv = results_dir / "sensitivity_analysis.csv"
    if not csv.exists():
        pytest.skip("sensitivity_analysis.csv not generated yet")
    df = pd.read_csv(csv)
    assert "sweep" in df.columns
    assert "ATPM" in df.columns
