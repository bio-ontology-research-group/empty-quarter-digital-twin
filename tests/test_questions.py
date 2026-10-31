from __future__ import annotations

import json
from pathlib import Path

import pytest

from eq_biodt.v03.composition import run_chain, twin_request
from eq_biodt.v03.external_models import MODELS
from eq_biodt.v03.questions import QUESTIONS, evaluate

ROOT = Path(__file__).resolve().parents[1]
FORCING = [{"water_mm": 12.0 if d == 0 else 0.0, "nacl_mol_l": 0.0, "evaporation_mm": 1.0, "temp_c": 25.0} for d in range(14)]


def _rnd6(x: float) -> float:
    return float(f"{x:.6g}")


def stored_central_ry_alfalfa(stored: dict, ec: float) -> float:
    lvl = next(l for l in stored["levels"] if l["action"]["ec_w_dS_m"] == ec)
    stored_6g = lvl["downstream"]["maas_hoffman_alfalfa"]["relative_yield_alfalfa"]["value"]
    ups = [MODELS[m] for m in stored["chain"] if m in MODELS and MODELS[m].kind == "upstream"]
    downs = [MODELS[m] for m in stored["chain"] if m in MODELS and MODELS[m].kind == "downstream"]
    r = run_chain(lvl["action"], ups, twin_request(FORCING), downs, envelope=False)
    exact = r["outputs"]["maas_hoffman_alfalfa:relative_yield_alfalfa"]["value"]
    assert _rnd6(exact) == pytest.approx(stored_6g, rel=1e-12)
    return exact


def test_fourteen_questions_with_themes():
    assert set(QUESTIONS) == {"rain-pattern","drip-or-flood","cloud-seeding","salty-water","alfalfa-bacterium",
                              "hotter","mulch-shade","solar-park","compaction","biochar","missing-co-genes",
                              "adding-microbes","air-hydrogen","carbon-fixation"}
    assert {q.theme for q in QUESTIONS.values()} == {"water","salt","heat","soil","microbes"}
    assert all(1 <= len(s.levels) <= 7 for q in QUESTIONS.values() for s in q.series)


def test_v07_levels_reproduce_stored_central_values():
    # each series sourced from a v0.7 experiment reproduces reviews/v07 central values exactly
    stored = json.loads((ROOT/"reviews/v07/composition/brackish_irrigation.json").read_text())
    spec = next(l for l in QUESTIONS["salty-water"].series[1].levels if l.x == 8.0)
    out = evaluate(spec)
    assert out["relative_yield_alfalfa"] == pytest.approx(stored_central_ry_alfalfa(stored, 8.0), rel=1e-12)

    checks = [
        ("drip-or-flood", 0, 10.0, "drip_vs_flood.json",
         lambda s: next(l for l in s["levels"] if l["action"]["irrigation_interval_days"] == 10)["twin"]["h2_uptake_mol_m2"]["value"] * 1e6,
         "h2_uptake_umol_m2"),
        ("cloud-seeding", 0, 1.0, "cloud_seeding.json",
         lambda s: next(l for l in s["levels"] if l["action"]["seeding_intensity"] == 1)["twin"]["h2_uptake_mol_m2"]["value"] * 1e6,
         "h2_uptake_umol_m2"),
        ("alfalfa-bacterium", 2, 7.81, "inoculated_alfalfa_salinity.json",
         lambda s: next(l for l in s["levels"] if l["model"] == "C_fitted_fresh" and l["ec_w_dS_m"] == 7.81)["relative_yield_sa187"]["value"],
         "relative_yield_alfalfa_sa187"),
        ("hotter", 1, 4.0, "warming.json",
         lambda s: next(l for l in s["levels"] if l["warming_c"] == 4 and l["rain_change_frac"] == 0)["twin"]["h2_uptake_mol_m2"]["value"] * 1e6,
         "h2_uptake_umol_m2"),
        ("mulch-shade", 0, 1.0, "mulch_shade.json",
         lambda s: next(l for l in s["levels"] if l["action"] == {"mulch_cover_fraction": 1, "shade_net_cover_fraction": 0})["twin"]["h2_uptake_mol_m2"]["value"] * 1e6,
         "h2_uptake_umol_m2"),
        ("solar-park", 0, 1.0, "solar_park.json",
         lambda s: next(l for l in s["levels"] if l["action"]["panel_cover_fraction"] == 1)["twin"]["h2_uptake_mol_m2"]["value"] * 1e6,
         "h2_uptake_umol_m2"),
        ("compaction", 0, 50.0, "compaction.json",
         lambda s: next(l for l in s["levels"] if l["action"]["offroad_passes"] == 50)["twin"]["h2_uptake_mol_m2"]["value"] * 1e6,
         "h2_uptake_umol_m2"),
        ("biochar", 0, 5.0, "biochar.json",
         lambda s: next(l for l in s["levels"] if l["action"]["biochar_dose_pct_w_w"] == 5)["twin"]["h2_uptake_mol_m2"]["value"] * 1e6,
         "h2_uptake_umol_m2"),
        ("adding-microbes", 0, 1.0, "inoculation.json",
         lambda s: next(l for l in s["levels"] if l["sample"] == "51PRr3" and l["coupling"] == "gate" and l["share"] == 0.01)["twin"]["co_uptake_mol_m2"]["value"] * 1e6,
         "co_uptake_umol_m2"),
        ("air-hydrogen", 0, 2.0, "h2_leakage.json",
         lambda s: next(l for l in s["levels"] if l["action"]["h2_multiplier"] == 2)["twin"]["h2_uptake_mol_m2"]["value"] * 1e6,
         "h2_uptake_umol_m2"),
        ("carbon-fixation", 0, "F14Dr1", "carbon_fixation_bound.json",
         lambda s: next(l for l in s["levels"] if l["sample"] == "F14Dr1" and l["scenario"].startswith("baseline"))["c_fixed_upper_bound"]["per_year_ug_c_m2"]["value"] / 1e6,
         "c_fixed_upper_bound_g_m2_yr"),
    ]
    for slug, s_idx, x_val, fname, extract_fn, out_key in checks:
        doc = json.loads((ROOT / "reviews/v07/composition" / fname).read_text())
        lvl_spec = next(l for l in QUESTIONS[slug].series[s_idx].levels if l.x == x_val)
        res = evaluate(lvl_spec)
        assert _rnd6(res[out_key]) == pytest.approx(_rnd6(extract_fn(doc)), rel=1e-12), f"mismatch on {slug}"


def test_v06_points_placed_on_axis():
    hot = QUESTIONS["hotter"]; assert any(l.x == 10.0 for l in hot.series[0].levels)          # 35 C fortnight = +10 C
    comp = QUESTIONS["compaction"]; assert any(l.twin.get("params",{}).get("porosity") == 0.3 for s in comp.series for l in s.levels)
