"""Numerical sanity checks of the twelve v0.7 composition experiments (run into a temporary directory)."""
import importlib.util
import json
from pathlib import Path

import pytest

from eq_biodt.v03 import external_models as xm
from eq_biodt.v03.composition import load_constants

ROOT = Path(__file__).resolve().parents[1]
_spec = importlib.util.spec_from_file_location("v07_composition", ROOT / "reproducibility/v07/composition.py")
EXP = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(EXP)


@pytest.fixture(scope="module")
def docs(tmp_path_factory):
    EXP.OUT = tmp_path_factory.mktemp("composition")
    return {fn.__name__: fn() for fn in EXP.EXPERIMENTS}


def test_all_experiments_complete_and_follow_the_spec(docs):
    assert len(docs) == 12
    for d in docs.values():
        for k in ("question", "chain", "levels", "refusals", "deduced_nulls", "plain_language", "caveats", "headline", "baseline"):
            assert k in d
        assert d["question"].endswith("?") and len(d["caveats"]) >= 2 and len(d["plain_language"].split(". ")) <= 8
        for h in d["headline"]:
            assert h["provenance"].split(":")[0] in ("twin", "external", "deduced")


def test_committed_results_are_reproduced_exactly(docs):
    for d in docs.values():
        txt = json.dumps(EXP.rnd(d), indent=2, sort_keys=True, allow_nan=False) + "\n"
        assert txt == (ROOT / "reviews/v07/composition" / f"{d['id']}.json").read_text()


def test_registry_mirror_of_the_digitised_ratio():
    reg = load_constants()
    for t in ("fresh", "dry"):
        o = xm.observed_mock_ratio(t)
        for k in ("value", "low", "high"):
            assert reg[f"sa187_mock_ratio_{t}"][k] == pytest.approx(o[k], abs=6e-5)
    assert not any(k.startswith("sa187_mock_") and k.endswith("t_ha") for k in reg)


def test_exp1_salinity_lowers_yield_and_raises_soil_salt(docs):
    lv = docs["exp01_brackish_irrigation"]["levels"]
    salt = [l["twin"]["mean_nacl_mol_l"]["value"] for l in lv]
    assert salt == sorted(salt) and salt[0] > 0
    ry = [l["downstream"]["maas_hoffman_alfalfa"]["relative_yield_alfalfa"]["value"] for l in lv]
    assert ry == sorted(ry, reverse=True) and ry[-1] < 1 and ry[0] == 1


def test_exp2_daily_pulses_have_lowest_uptake(docs):
    d = docs["exp02_drip_vs_flood"]
    u = [l["twin"]["h2_uptake_mol_m2"]["value"] for l in d["levels"]]
    assert u[0] < u[2] and d["contrast"]["v06_finding_reproduced_daily_lowest"]


def test_exp3_mulch_helps_shade_cools(docs):
    lv = {l["label"]: l["difference_vs_baseline"]["h2_uptake_mol_m2"]["value"] for l in docs["exp03_mulch_shade"]["levels"]}
    assert lv["mulch 1, shade 0"] > 0 > lv["mulch 0, shade 1"] and lv["mulch 0, shade 0"] == 0


def test_exp4_biochar_zero_dose_is_baseline_and_effect_small(docs):
    lv = docs["exp04_biochar"]["levels"]
    assert lv[0]["difference_vs_baseline"]["h2_uptake_mol_m2"]["value"] == 0
    assert abs(lv[-1]["difference_vs_baseline"]["h2_uptake_mol_m2"]["percent"]["value"]) < 10


def test_exp5_gate_opening_negative_control_and_modes(docs):
    d = docs["exp05_inoculation"]
    g = d["gate"]
    assert g["co_gate_opened_by_carrier"] and g["co_uptake_before_umol_m2"] == 0 and g["co_uptake_after_umol_m2"] > 0
    assert g["gate_already_open"] and g["open_sample_co_change_gate_umol_m2"] == 0
    assert g["negative_control_co_change_gate_umol_m2"] == 0 and g["negative_control_h2_change_gate_umol_m2"] == 0
    prop = [l for l in d["levels"] if l["case"].startswith("51PRr3 (CO") and l["coupling"] == "proportional"]
    co = [l["twin"]["co_uptake_mol_m2"]["value"] for l in prop]
    assert co == sorted(co) and co[0] > 0                # graded dose-response in proportional mode
    gate = [l["twin"]["co_uptake_mol_m2"]["value"] for l in d["levels"] if l["coupling"] == "gate" and l["case"].startswith("51PRr3 (CO")]
    assert len(set(gate)) == 1                           # gate: no dose dependence


def test_exp6_bound_equals_h2_times_yield_and_zero_is_deduced(docs):
    d = docs["exp06_carbon_fixation_bound"]
    reg = load_constants()
    for l in d["levels"]:
        k = l["downstream"]["knallgas_carbon_fixation"]["co2_fixed_upper_bound_mol_m2"]
        if l["sample"] == d["gate"]["zero_share_sample"]:
            assert k["value"] == 0 and k["provenance"] == "deduced"
        else:
            assert k["value"] == pytest.approx(l["twin"]["h2_uptake_mol_m2"]["value"] * reg["knallgas_co2_per_h2"]["value"])
            c = l["c_fixed_upper_bound"]
            assert c["per_year_ug_c_m2"]["value"] == pytest.approx(c["per_period_ug_c_m2"]["value"] * 365.25 / 14)
        assert l["co_oxidation"]["net_atmospheric_c_change_mol_m2"] == 0 and l["context"]["low"] < l["context"]["high"]
    assert d["gate"]["positive_sample_share"] > 0


def test_exp7_linear_and_not_oxygen_limited(docs):
    d = docs["exp07_h2_leakage"]
    lv = d["levels"]
    u = [l["twin"]["h2_uptake_mol_m2"]["value"] for l in lv]
    assert u[2] / u[0] == pytest.approx(2.0) and u[3] / u[0] == pytest.approx(10.0)
    assert d["linearity"]["linear"] and not d["linearity"]["oxygen_limited"]


def test_exp8_to_11_signs_and_zero_levels(docs):
    sp = docs["exp08_solar_park"]["levels"]
    assert sp[0]["difference_vs_baseline"]["h2_uptake_mol_m2"]["value"] == 0 and sp[-1]["difference_vs_baseline"]["h2_uptake_mol_m2"]["value"] > 0
    cs = [l["twin"]["h2_uptake_mol_m2"]["value"] for l in docs["exp09_cloud_seeding"]["levels"]]
    assert cs == sorted(cs) and cs[0] == docs["exp09_cloud_seeding"]["levels"][0]["twin_baseline"]["h2_uptake_mol_m2"]
    co = docs["exp10_compaction"]["levels"]
    por = [l["porosity_after"]["value"] for l in co]
    assert por == sorted(por, reverse=True) and por[0] == pytest.approx(0.4)
    w = {(l["warming_c"], l["rain_change_frac"]): l["twin"]["h2_uptake_mol_m2"]["value"] for l in docs["exp11_warming"]["levels"]}
    assert w[(4, 0)] > w[(0, 0)] and w[(0, 0.2)] > w[(0, 0)] > w[(0, -0.2)]


def test_exp12_fit_reproduces_observed_ratio_and_validation(docs):
    d = docs["exp12_inoculated_alfalfa_salinity"]
    f = d["fitted_leaching"]
    for t in ("fresh", "dry"):
        assert f[t]["observed_ratio"] == pytest.approx({"fresh": 0.7958, "dry": 0.8282}[t], abs=1e-4)
        assert f[t]["f_range_from_digitisation_error"][0] < f[t]["f_central"] < f[t]["f_range_from_digitisation_error"][1]
    ry = {(l["model"], l["ec_w_dS_m"]): l for l in d["levels"]}
    for m, t in (("C_fitted_fresh", "fresh"), ("C_fitted_dry", "dry")):
        hi = ry[(m, 7.81)]["relative_yield_mock"]["value"]; lo = ry[(m, 3.12)]["relative_yield_mock"]["value"]
        assert hi / lo == pytest.approx(f[t]["observed_ratio"], abs=1e-4)
        assert ry[(m, 3.12)]["relative_yield_sa187"]["value"] == pytest.approx(lo * 1.15)
    assert abs(d["validation_mock_ratio"]["A_fao_fixed_leaching"]["gap_to_observed_fresh"]) > abs(d["validation_mock_ratio"]["B_twin_physics"]["gap_to_observed_fresh"])
    assert d["mock_yields_digitised"]["fresh"]["low_3.12_t_ha"] == pytest.approx(5.411)
