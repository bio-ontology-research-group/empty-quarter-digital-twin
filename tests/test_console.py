"""API-level tests for the console endpoints."""
import json
from pathlib import Path

import pytest
from fastapi.testclient import TestClient

from eq_biodt.v03.api import app

client = TestClient(app)
ROOT = Path(__file__).resolve().parents[1]
CANDIDATE = "11Dr3__single_SemiBin_120"


def test_console_served_and_current():
    html = client.get("/").text
    assert "EQ-BioDT v0.6" in html and "eq-biodt-v06-scenario.json" in html
    assert "v04" not in html and "/v3/mags?limit=5438" not in html


def test_mags_search_is_limited():
    d = client.get("/v3/mags", params={"candidates_only": True, "limit": 5, "q": "Actinomycetota"}).json()
    assert len(d["mags"]) <= 5 and d["matching_mags"] >= len(d["mags"])
    assert all(m["h2_cbb_all_steps"] and m["phylum"] == "Actinomycetota" for m in d["mags"])
    assert client.get("/v3/mags", params={"q": "no-such-taxon-xyz"}).json()["matching_mags"] == 0


def test_mag_detail_uses_reasoner_entailments():
    d = client.get(f"/v3/mags/{CANDIDATE}").json()
    classes = {c["class"] for c in d["entailed_classes"]}
    assert {"rak:H2AutotrophicCandidateMAG", "rak:CBBCompleteMAG"} <= classes
    assert len(d["cbb_steps"]) == 11 and all(s["present"] for s in d["cbb_steps"])
    just = d["justifications"]["H2AutotrophicCandidateMAG"]
    assert just["mag_axioms"] and just["tbox_axioms"]
    assert client.get("/v3/mags/not-a-mag").status_code == 404


def test_transfer_numbers_match_manuscript():
    t = client.get("/v3/transfer").json()
    fold = {c["site"]: c for c in t["namib"]["contrasts_quantified"]}
    assert [round(fold["C14-2"][k], 1) for k in ("fold", "fold_lcl", "fold_ucl")] == [41.3, 31.1, 56.0]
    assert [round(fold["C14-4"][k], 1) for k in ("fold", "fold_lcl", "fold_ucl")] == [19.4, 15.6, 24.0]
    rmse = {m["model"]: round(m["rmse"], 2) for m in t["sonoran"]["field_metrics"]}
    assert rmse == {"frozen_reference": 1.65, "calibrated_kinetics": 1.72, "calibrated_diffusion_only": 1.72,
                    "training_mean": 0.79, "linear_moisture_temperature": 1.35, "calibrated_retention": 0.69}
    tex = (ROOT/"paper"/"main.tex").read_text()
    for token in ("41.3", "19.4", "1.65", "0.69", "0.79", "1.35"):
        assert token in tex
    assert {r["soil"] for r in t["temperate"]["rmse"]} == {"Forest", "Meadow", "PB_sand"}


@pytest.mark.parametrize("node,value,version,admissible,stage", [
    ("K", None, None, False, "zone"), ("q_h2", None, None, False, "zone"), ("J_h2", None, None, False, "zone"),
    ("W", 999, None, False, "feasibility"), ("W", 5, None, True, "feasible"),
    ("psi", -1, None, False, "zone"), ("psi", -0.5, "matric", True, "feasible")])
def test_admissibility(node, value, version, admissible, stage):
    d = client.post("/v3/admissibility", json={"node": node, "value": value, "version": version}).json()
    assert d["admissible"] is admissible and d["stage"] == stage and d["message"]


def test_daily_12mm_and_gate_zero_co():
    ex = json.loads((ROOT/"data"/"v03"/"example_intervention.json").read_text())
    days = {str(int(i*14/14)): {"water_mm": 12/14} for i in range(14)}
    body = {**ex, "intervention": {"forcing": days}, "community": {"sample": "51PRr3"}, "coupling": "gate"}
    d = client.post("/v3/intervene", json=body).json()
    assert round(d["gas_uptake"]["h2_mol_m2"]["intervened"]*1e6, 2) == 2.87
    assert d["gas_uptake"]["co_mol_m2"]["intervened"] == 0
    assert d["deduced"][0]["effect"] == "zero_by_deduction"
    assert all(d["checks"].values())
