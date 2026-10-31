"""v0.7 composition of the twin with external models: contract, refusals, gates, DAG, determinism, twin extensions,
and one unit test per external model."""
from dataclasses import replace
import json
from math import exp

import pytest

from eq_biodt.v03 import community as comm
from eq_biodt.v03 import external_models as xm
from eq_biodt.v03 import identification, theory
from eq_biodt.v03.composition import (CompositionError, ExternalModel, Port, check_dag, constant_values, constants_sha256,
                                      load_constants, run_chain, twin_request)
from eq_biodt.v03.interventions import FORCING_TARGETS, intervene
from eq_biodt.v03.process import Parameters, State, simulate

M = xm.MODELS
REG = load_constants()
P = Parameters(k_h2_s=.001, k_co_s=.001)
DATA = twin_request([{"temp_c": 25., "evaporation_mm": 1.}])["data"]
DAYS = 14
FLAT = [{"water_mm": 0., "nacl_mol_l": 0., "evaporation_mm": 1., "temp_c": 25.} for _ in range(DAYS)]
WET = [{**f, "water_mm": 12. if d == 0 else 0.} for d, f in enumerate(FLAT)]
CO_MAG = "11Dr3__single_SemiBin_57"               # entailed COCarrierMAG
C = constant_values(sorted(REG))                  # central values of every registry entry


# ---------------------------------------------------------------- registry

def test_registry_complete_and_unverified():
    assert len(REG) == 44 and constants_sha256()
    for e in REG.values():
        assert e["verified"] is False and e["low"] <= e["value"] <= e["high"]
        assert set(e["source"]) == {"authors", "year", "title", "doi_or_url", "page_or_table"}
    used = {c for m in M.values() for c in m.constants}
    assert used <= set(REG)
    for k in ("sa187_effect_fresh_low_s1", "sa187_benefit_b", "mh_alfalfa_threshold_dS_m", "fao29_ecsw_per_ece"):
        assert k in REG


def test_no_literature_constant_is_hidden_in_external_models():
    # every constant a model function reads is declared in model.constants (fn gets exactly those)
    for m in M.values():
        used = {c: C[c] for c in m.constants}
        inputs = {p.name: 0.5 for p in m.inputs}
        # a model that reads an undeclared id raises KeyError in fn
        sample = {"baseline_water_mm": [0., 3.], "baseline_nacl_mol_l": [0., 0.], "baseline_temp_c": [25., 26.],
                  "baseline_evaporation_mm": [1., 1.], "baseline_porosity": .4, "baseline_h2_mol_m3": 2.13e-5,
                  "irrigation_interval_days": 1, "irrigation_mm_per_day": 4., "inoculum_mag_id": "x", "inoculum_share": .01,
                  "h2_multiplier": 2., "offroad_passes": 3, "warming_c": 2., "rain_change_frac": .1, "ece_dS_m": 4.,
                  "ec_w_dS_m": 3.12, "mean_nacl_mol_l": .02, "relative_yield_alfalfa": .8, "h2_uptake_mol_m2": .001,
                  "co_uptake_mol_m2": .0001, "irrigation_nacl_mol_l": .03}
        inputs = {p.name: sample.get(p.name, .5) for p in m.inputs}
        out = m.fn(inputs, used)
        assert set(out) == {o.name for o in m.outputs}, m.id


# ---------------------------------------------------------------- contract

def test_port_mismatch_raises():
    bad = ExternalModel("bad_yield", "wrong unit", "downstream", (Port("ece_dS_m", "rak:SaturationExtractElectricalConductivity", "mS/cm"),),
                        (Port("bad_out", "rak:X", "-"),), {}, (), lambda i, c: {"bad_out": 1.})
    with pytest.raises(CompositionError, match="expects"):
        run_chain({"ec_w_dS_m": 3.}, [], twin_request(FLAT), [M["fao29_leaching"], bad])


def test_upstream_to_twin_port_mismatch_raises():
    bad = ExternalModel("bad_temp", "wrong unit", "upstream", (), (Port("temp_c", "rak:SoilTemperature", "K"),), {}, (),
                        lambda i, c: {"temp_c": 300.})
    with pytest.raises(CompositionError, match="twin input"):
        run_chain({}, [bad], twin_request(FLAT), [])


def test_missing_input_and_duplicate_producer_raise():
    with pytest.raises(CompositionError, match="not produced"):
        run_chain({}, [M["ec_to_nacl"]], twin_request(FLAT), [])
    with pytest.raises(CompositionError, match="more than one producer"):
        run_chain({"ec_w_dS_m": 1.}, [M["ec_to_nacl"]], twin_request(FLAT), [M["fao29_leaching"], replace(M["fao29_leaching"], id="fao29_again")])


def test_wrong_kind_and_gate_on_upstream_rejected():
    with pytest.raises(CompositionError):
        run_chain({"ec_w_dS_m": 1.}, [M["fao29_leaching"]], twin_request(FLAT), [])
    with pytest.raises(CompositionError, match="gate"):
        ExternalModel("g", "g", "upstream", (), (), {}, (), lambda i, c: {}, gate="H2CarrierMAG")


# ---------------------------------------------------------------- validity as admissibility

def test_validity_refusal_is_a_result_with_interval():
    r = run_chain({"ec_w_dS_m": 20., "irrigation_interval_days": 3},
                  [M["ec_to_nacl"], M["crop_water_demand_generic"], M["irrigation_schedule"]], twin_request(FLAT), [M["fao29_leaching"]])
    assert r["refused"]["refused"] is True and r["refused"]["model"] == "ec_to_nacl"
    assert r["refused"]["port"] == "ec_w_dS_m" and r["refused"]["value"] == 20. and r["refused"]["interval"] == [0.0, 16.0]
    assert r["twin"] is None and r["outputs"] == {}


def test_downstream_validity_refusal_and_twin_refusal():
    wet = twin_request(FLAT)
    r = run_chain({"warming_c": 6., "rain_change_frac": 0.}, [M["climate_warming"]], twin_request(
        [{**f, "temp_c": 48.} for f in FLAT]), [])
    assert r["refused"]["model"] == "twin" and "Invalid forcing" in r["refused"]["reason"]    # 54 degC leaves the twin's 0-50 degC domain
    r = run_chain({"warming_c": 9., "rain_change_frac": 0.}, [M["climate_warming"]], wet, [])
    assert r["refused"]["model"] == "climate_warming" and r["refused"]["interval"] == [0.0, 6.0]


# ---------------------------------------------------------------- deductive gate

def _spy(calls):
    inner = M["knallgas_carbon_fixation"]

    def fn(i, c):
        calls.append(1)
        return inner.fn(i, c)
    return replace(inner, fn=fn)


def test_closed_gate_gives_deduced_zero_and_does_not_call_fn():
    calls = []
    req = twin_request(WET, community={"sample": "51PRr3"}, coupling="gate")      # H2 autotroph share is 0 in this sample
    assert comm.carrier_shares("51PRr3")["shares"]["H2AutotrophicCandidateMAG"] == 0
    r = run_chain({}, [], req, [_spy(calls)])
    o = r["outputs"]["knallgas_carbon_fixation:co2_fixed_upper_bound_mol_m2"]
    assert o["value"] == 0.0 and o["provenance"] == "deduced" and not calls
    assert "absence of entailment is not read as biological absence" in o["derivation"]
    assert r["deduced_nulls"][0]["model"] == "knallgas_carbon_fixation"


def test_open_gate_calls_fn_and_is_bounded_by_h2_uptake():
    calls = []
    req = twin_request(WET, community={"sample": "10Dr2"}, coupling="gate")
    assert comm.carrier_shares("10Dr2")["shares"]["H2AutotrophicCandidateMAG"] > 0
    r = run_chain({}, [], req, [_spy(calls)])
    o = r["outputs"]["knallgas_carbon_fixation:co2_fixed_upper_bound_mol_m2"]
    h2 = r["outputs"]["twin:h2_uptake_mol_m2"]["value"]
    assert calls and o["provenance"].startswith("external:") and o["value"] == pytest.approx(h2 * C["knallgas_co2_per_h2"])
    assert o["low"] <= o["value"] <= o["high"] and o["low"] < o["high"]


def test_gate_without_community_is_a_structural_error():
    with pytest.raises(CompositionError, match="no community"):
        run_chain({}, [], twin_request(WET), [M["knallgas_carbon_fixation"]])


# ---------------------------------------------------------------- DAG

def test_feedback_is_refused():
    fb = ExternalModel("feedback", "downstream output feeds an upstream input", "downstream", (xm.H2_UPTAKE,), (xm.EC_W,), {}, (),
                       lambda i, c: {"ec_w_dS_m": 1.})
    with pytest.raises(CompositionError, match="cycle"):
        run_chain({}, [M["ec_to_nacl"], M["crop_water_demand_generic"], M["irrigation_schedule"]], twin_request(FLAT), [fb])
    with pytest.raises(CompositionError, match="cycle"):
        check_dag([M["ec_to_nacl"], M["crop_water_demand_generic"], M["irrigation_schedule"]], [fb])


def test_acyclic_chain_passes():
    check_dag([M["ec_to_nacl"], M["crop_water_demand_generic"], M["irrigation_schedule"]], [M["fao29_leaching"], M["maas_hoffman_alfalfa"]])


# ---------------------------------------------------------------- envelope, provenance, determinism

def _brackish(ec=4.0):
    return run_chain({"ec_w_dS_m": ec, "irrigation_interval_days": 3},
                     [M["ec_to_nacl"], M["crop_water_demand_generic"], M["irrigation_schedule"]], twin_request(FLAT),
                     [M["fao29_leaching"], M["maas_hoffman_alfalfa"]])


def test_envelope_is_one_at_a_time_and_wraps_every_value():
    r = _brackish()
    assert r["refused"] is None and r["twin"]["checks_passed"]
    ry = r["outputs"]["maas_hoffman_alfalfa:relative_yield_alfalfa"]
    assert ry["value"] == pytest.approx(1 - 0.073 * (1.5 * 4.0 - 2.0))
    # one constant at a time: the widest single effect is the leaching factor (ECe = 2.0 x 4 or 1.0 x 4); slope and
    # threshold (+-20%) move the value less, and the envelope is not the joint worst case
    assert ry["low"] == pytest.approx(1 - 0.073 * (2.0 * 4.0 - 2.0))
    assert ry["high"] == pytest.approx(1 - 0.073 * (1.0 * 4.0 - 2.0))
    for w in r["outputs"].values():
        assert set(w) >= {"value", "low", "high", "unit", "provenance", "constants"}
    assert r["outputs"]["twin:h2_uptake_mol_m2"]["provenance"] == "twin"
    assert r["twin"]["trace"]["graph_sha256"] and r["registry_sha256"] == constants_sha256()
    assert r["constants_used"]["fao29_ece_per_ecw"]["verified"] is False


def test_scenario_low_high_and_no_envelope():
    lo = run_chain({"ec_w_dS_m": 4.}, [], twin_request(FLAT), [M["fao29_leaching"]], scenario="low")
    hi = run_chain({"ec_w_dS_m": 4.}, [], twin_request(FLAT), [M["fao29_leaching"]], scenario="high")
    assert lo["outputs"]["fao29_leaching:ece_dS_m"]["value"] == 4.0 and hi["outputs"]["fao29_leaching:ece_dS_m"]["value"] == 8.0
    assert lo["outputs"]["fao29_leaching:ece_dS_m"]["low"] == lo["outputs"]["fao29_leaching:ece_dS_m"]["high"] == 4.0


def test_determinism():
    a, b = json.dumps(_brackish(), sort_keys=True), json.dumps(_brackish(), sort_keys=True)
    assert a == b


def test_overrides_are_reported():
    r = run_chain({"ec_w_dS_m": 4.}, [], twin_request(FLAT), [M["fao29_leaching"]], constants={"fao29_ece_per_ecw": 1.2}, envelope=False)
    assert r["outputs"]["fao29_leaching:ece_dS_m"]["value"] == pytest.approx(4.8) and r["constants_used"]["fao29_ece_per_ecw"]["value"] == 1.2


# ---------------------------------------------------------------- twin extension 1: atmospheric H2 / CO forcing

def test_atmospheric_forcing_is_admissible_and_in_the_graph():
    assert theory.check_admissible(["h2_mol_m3", "co_mol_m3"]) == ["C_h2", "C_co"]
    assert {"h2_mol_m3", "co_mol_m3"} <= FORCING_TARGETS
    g = theory.graph()
    edges = {(e["source"], e["target"]) for e in g["edges"]}
    assert {("C_h2", "J_h2"), ("C_co", "J_co"), ("C_a", "C_h2"), ("C_a", "C_co")} <= edges and ("C_a", "J_h2") not in edges
    assert theory.slice_is_acyclic()
    G = identification.unrolled_graph(2)
    assert ("C_h2@0", "J_h2@0") in G.edges and G.is_acyclic()


def test_forcing_h2_scales_uptake_and_defaults_to_parameter():
    f2 = [{**f, "h2_mol_m3": 2 * P.h2_mol_m3} for f in WET]
    base, dbl = simulate(State(1.8), WET, P), simulate(State(1.8), f2, P)
    assert dbl["summary"]["h2_uptake_mol_m2"] / base["summary"]["h2_uptake_mol_m2"] == pytest.approx(2.0, rel=1e-6)    # linear regime
    same = simulate(State(1.8), [{**f, "h2_mol_m3": P.h2_mol_m3, "co_mol_m3": P.co_mol_m3} for f in WET], P)
    assert same["trajectory"] == base["trajectory"] and same["summary"] == base["summary"]
    with pytest.raises(ValueError):
        simulate(State(1.8), [{**WET[0], "h2_mol_m3": -1.}], P)


def test_v06_numbers_unchanged_without_the_new_keys():
    ex = json.loads((DATA / "example_scenario.json").read_text())
    s = simulate(State(**ex["initial"]), ex["forcing"], Parameters(**ex["parameters"]))
    assert s["summary"]["h2_uptake_mol_m2"] == 0.00028263236248340614
    ref = json.loads((DATA.parents[1] / "reviews" / "v06" / "paired_regimes_summary.json").read_text())["reference_uptake_mol_m2"]
    n = 4
    forcing = [{"water_mm": 12 / n if day in {int(i * 14 / n) for i in range(n)} else 0, "nacl_mol_l": 0., "evaporation_mm": 1., "temp_c": 25.}
               for day in range(14)]
    rows = simulate(State(1.8), forcing, P)["trajectory"]
    assert sum(r["gas"]["h2_mol_m2_s"] * 86400 for r in rows) == ref["four_pulses"]["h2"]
    assert sum(r["gas"]["co_mol_m2_s"] * 86400 for r in rows) == ref["four_pulses"]["co"]
    r = intervene(State(1.8), WET, P, {"forcing": {"0": {"h2_mol_m3": 2 * P.h2_mol_m3}}}, DATA)
    assert all(r["checks"].values()) and r["gas_uptake"]["h2_mol_m2"]["difference"] > 0


def test_do_c_h2_through_the_chain_reports_sink_per_m2():
    r = run_chain({"h2_multiplier": 2.}, [M["h2_atmosphere"]], twin_request(WET), [], envelope=False)
    assert r["upstream"]["h2_atmosphere"]["h2_mixing_ratio_ppm"]["value"] == pytest.approx(2 * REG["h2_background_ppm"]["value"])
    base = simulate(State(1.8), WET, P)["summary"]["h2_uptake_mol_m2"]
    assert r["outputs"]["twin:h2_uptake_mol_m2"]["value"] == pytest.approx(2 * base, rel=1e-6)


# ---------------------------------------------------------------- twin extension 2: inoculation

def test_inoculate_opens_the_co_gate_for_51PRr3():
    s = comm.carrier_shares("51PRr3")["shares"]
    assert s["COCarrierMAG"] == 0 and comm.gamma(s, "gate")["co"] == 0           # closed in v0.6
    assert CO_MAG in comm._load()["classes"] and "COCarrierMAG" in comm._load()["classes"][CO_MAG]
    r = comm.inoculate("51PRr3", CO_MAG, 0.01)
    assert r["gamma"]["co"] > 0 and r["baseline_gamma"]["co"] == 0 and r["shares"]["COCarrierMAG"] == pytest.approx(0.01)
    assert r["shares"]["H2CarrierMAG"] == pytest.approx(0.99 * s["H2CarrierMAG"])      # others rescaled
    assert all(r["consistency"].values())
    assert r["evidence_sha256"] == comm.carrier_shares("51PRr3")["evidence_sha256"]    # E and K untouched


def test_inoculate_validation_and_zero_share():
    with pytest.raises(ValueError):
        comm.inoculate("51PRr3", "no_such_mag", 0.01)
    with pytest.raises(ValueError):
        comm.inoculate("51PRr3", CO_MAG, 1.5)
    assert comm.inoculate("51PRr3", CO_MAG, 0.0)["shares"] == comm.carrier_shares("51PRr3")["shares"]
    prop = comm.inoculate("51PRr3", CO_MAG, 0.05, mode="proportional")
    assert prop["gamma"]["co"] > 0


def test_intervene_with_inoculation_shows_gate_opening():
    r = intervene(State(1.8), WET, P, {"community": {"sample": "51PRr3"}, "coupling": "gate",
                                       "inoculation": {"mag_id": CO_MAG, "share": 0.01}}, DATA)
    assert r["inoculation"]["gate_opened"] == {"h2": False, "co": True}
    assert r["gas_uptake"]["co_mol_m2"]["baseline"] == 0.0 and r["gas_uptake"]["co_mol_m2"]["intervened"] > 0
    assert r["gas_uptake"]["h2_mol_m2"]["difference"] == 0.0
    assert all(r["checks"].values()) and not any(d["gas"] == "co" for d in r["deduced"])
    with pytest.raises(ValueError, match="needs a community"):
        intervene(State(1.8), WET, P, {"inoculation": {"mag_id": CO_MAG, "share": 0.01}}, DATA)


def test_inoculation_chain_and_already_open_sample():
    r = run_chain({"inoculum_mag_id": CO_MAG, "inoculum_share": 0.01}, [M["inoculant"]],
                  twin_request(WET, community={"sample": "51PRr3"}, coupling="gate"), [M["co_oxidation_carbon_balance"]])
    assert r["twin"]["inoculation"]["gate_opened"]["co"] is True
    assert r["outputs"]["twin:co_uptake_mol_m2"]["value"] > 0
    assert r["outputs"]["co_oxidation_carbon_balance:net_atmospheric_c_change_mol_m2"]["value"] == 0.0
    open_sample = max(comm._load()["table"].index, key=lambda s: comm.carrier_shares(s)["shares"]["COCarrierMAG"])
    r2 = run_chain({"inoculum_mag_id": CO_MAG, "inoculum_share": 0.01}, [M["inoculant"]],
                   twin_request(WET, community={"sample": open_sample}, coupling="gate"), [], envelope=False)
    assert r2["twin"]["inoculation"]["gate_opened"]["co"] is False       # gate was already open
    assert r2["twin"]["checks_passed"]


# ---------------------------------------------------------------- unit tests of external model functions

def test_ec_to_nacl():
    assert M["ec_to_nacl"].fn({"ec_w_dS_m": 4.}, C) == {"irrigation_nacl_mol_l": 4. * C["nacl_mol_l_per_dS_m"]}


def test_irrigation_schedule_equal_totals_and_salt_mass():
    inp = {"baseline_water_mm": [1.] + [0.] * 9, "baseline_nacl_mol_l": [0.] * 10, "irrigation_mm_per_day": 2., "irrigation_interval_days": 5,
           "irrigation_nacl_mol_l": 0.05}
    out = M["irrigation_schedule"].fn(inp, {})
    assert sum(out["water_mm"]) == pytest.approx(1. + 20.) and out["water_mm"][0] == pytest.approx(11.) and out["water_mm"][5] == pytest.approx(10.)
    assert sum(w * c for w, c in zip(out["water_mm"], out["nacl_mol_l"])) == pytest.approx(20. * 0.05)
    daily = M["irrigation_schedule"].fn({**inp, "irrigation_interval_days": 1}, {})
    assert sum(daily["water_mm"]) == pytest.approx(sum(out["water_mm"]))


def test_water_demand_models():
    assert M["crop_water_demand_generic"].fn({}, C) == {"irrigation_mm_per_day": C["crop_etc_mm_per_day"]}
    assert M["crop_water_demand_alfalfa"].fn({}, C) == {"irrigation_mm_per_day": C["alfalfa_etc_mm_per_day"]}


def test_fao29_and_maas_hoffman():
    assert M["fao29_leaching"].fn({"ec_w_dS_m": 4.}, C)["ece_dS_m"] == pytest.approx(6.0)
    for crop, thr, slope in (("date_palm", 4.0, 3.6), ("barley", 8.0, 5.0), ("alfalfa", 2.0, 7.3), ("tomato", 2.5, 9.9)):
        m = M[f"maas_hoffman_{crop}"]
        assert m.fn({"ece_dS_m": thr}, C)[f"relative_yield_{crop}"] == 1.0
        assert m.fn({"ece_dS_m": thr + 1}, C)[f"relative_yield_{crop}"] == pytest.approx(1 - slope / 100)
        assert m.fn({"ece_dS_m": 40.}, C)[f"relative_yield_{crop}"] == 0.0


def test_twin_solution_to_ece():
    out = M["twin_solution_to_ece"].fn({"mean_nacl_mol_l": 0.04}, C)["ece_dS_m"]
    assert out == pytest.approx(0.04 / C["nacl_mol_l_per_dS_m"] / 2.0)


@pytest.mark.parametrize("trait", ["fresh", "dry"])
def test_fitted_leaching_reproduces_observed_ratio(trait):
    f = xm.fit_leaching_factor(C, trait)
    obs = C[f"sa187_mock_{trait}_high_s1_t_ha"] / C[f"sa187_mock_{trait}_low_s1_t_ha"]
    ry = lambda ec: xm.maas_hoffman(f * ec, C["mh_alfalfa_threshold_dS_m"], C["mh_alfalfa_slope_pct_per_dS_m"])
    assert ry(7.81) / ry(3.12) == pytest.approx(obs, rel=1e-9)
    out = M[f"fitted_leaching_{trait}"].fn({"ec_w_dS_m": 3.12}, C)
    assert out["ece_dS_m"] == pytest.approx(f * 3.12) and out["fitted_leaching_factor"] == f
    r = xm.observed_mock_ratio(trait)
    assert r["low"] < r["value"] < r["high"] and r["value"] == pytest.approx(obs)


def test_sa187_benefit():
    assert M["sa187_benefit"].fn({"relative_yield_alfalfa": 0.8}, C)["relative_yield_alfalfa_sa187"] == pytest.approx(0.8 * 1.15)


def test_mulch_shade_biochar():
    assert M["mulch_evaporation"].fn({"mulch_cover_fraction": 1., "baseline_evaporation_mm": [2., 4.]}, C)["evaporation_mm"] == \
        pytest.approx([2 * (1 - C["mulch_evap_reduction_frac"]), 4 * (1 - C["mulch_evap_reduction_frac"])])
    assert M["shade_net_temperature"].fn({"shade_net_cover_fraction": 1., "baseline_temp_c": [30., 2.]}, C)["temp_c"] == \
        pytest.approx([30 - C["shade_net_soil_temp_reduction_c"], 0.0])
    assert M["biochar_porosity"].fn({"biochar_dose_pct_w_w": 5., "baseline_porosity": .4}, C)["porosity"] == \
        pytest.approx(.4 + 5 * C["biochar_porosity_increase_per_pct"])
    assert M["biochar_porosity"].fn({"biochar_dose_pct_w_w": 0., "baseline_porosity": .4}, C)["porosity"] == .4


def test_inoculant():
    assert M["inoculant"].fn({"inoculum_mag_id": "m", "inoculum_share": .05}, C) == {"inoculation": {"mag_id": "m", "share": .05}}


def test_knallgas_and_co_balance():
    out = M["knallgas_carbon_fixation"].fn({"h2_uptake_mol_m2": 6.0}, C)
    assert out["co2_fixed_upper_bound_mol_m2"] == pytest.approx(6.0 * C["knallgas_co2_per_h2"])
    assert out["c_fixed_upper_bound_ug_m2"] == pytest.approx(out["co2_fixed_upper_bound_mol_m2"] * 12.011e6)
    assert M["co_oxidation_carbon_balance"].fn({"co_uptake_mol_m2": .3}, C) == {"co2_from_co_mol_m2": .3, "net_atmospheric_c_change_mol_m2": 0.0}
    assert xm.annualise(1.0, 14.0) == pytest.approx(C["days_per_year"] / 14)


def test_h2_atmosphere():
    out = M["h2_atmosphere"].fn({"h2_multiplier": 10., "baseline_h2_mol_m3": 2e-5}, C)
    assert out["h2_mol_m3"] == pytest.approx(2e-4) and out["h2_mixing_ratio_ppm"] == pytest.approx(5.3)


def test_solar_park_models():
    out = M["solar_park_microclimate"].fn({"panel_cover_fraction": .5, "baseline_temp_c": [30.], "baseline_evaporation_mm": [2.]}, C)
    assert out["temp_c"] == [pytest.approx(30 - .5 * C["solar_panel_soil_temp_reduction_c"])]
    assert out["evaporation_mm"] == [pytest.approx(2 * (1 - .5 * C["solar_panel_evap_reduction_frac"]))]
    w = M["panel_washing"].fn({"panel_cover_fraction": 1., "baseline_water_mm": [0.] * 30}, C)["water_mm"]
    k = int(C["panel_wash_interval_days"])
    assert [d for d, x in enumerate(w) if x > 0] == [d for d in range(30) if (d + 1) % k == 0]
    assert max(w) == pytest.approx(C["panel_wash_l_m2_per_wash"])


def test_cloud_seeding_only_changes_rain_days():
    w = M["cloud_seeding"].fn({"seeding_intensity": 2., "baseline_water_mm": [0., 5., 0.]}, C)["water_mm"]
    assert w == [0., pytest.approx(5 * (1 + 2 * C["cloud_seeding_rain_increase_frac"])), 0.]
    assert M["cloud_seeding"].fn({"seeding_intensity": 0., "baseline_water_mm": [3.]}, C)["water_mm"] == [3.]


def test_offroad_compaction():
    f = M["offroad_compaction"].fn
    assert f({"offroad_passes": 0., "baseline_porosity": .4}, C)["porosity"] == pytest.approx(.4)
    rs = C["particle_density_g_cm3"]
    expect = 1 - ((1 - .4) * rs + C["offroad_bd_max_increase_g_cm3"] * (1 - exp(-10 / C["offroad_bd_passes_efold"]))) / rs
    assert f({"offroad_passes": 10., "baseline_porosity": .4}, C)["porosity"] == pytest.approx(expect)
    assert f({"offroad_passes": 50., "baseline_porosity": .4}, C)["porosity"] < f({"offroad_passes": 1., "baseline_porosity": .4}, C)["porosity"] < .4


def test_climate_warming():
    out = M["climate_warming"].fn({"warming_c": 2., "rain_change_frac": -.2, "baseline_temp_c": [25.], "baseline_evaporation_mm": [1.],
                                   "baseline_water_mm": [10.]}, C)
    assert out["temp_c"] == [pytest.approx(25 + 2 * C["soil_temp_per_air_temp_warming"])]
    assert out["evaporation_mm"] == [pytest.approx(1 + 2 * C["evap_scaling_per_c"])] and out["water_mm"] == [pytest.approx(8.)]


def test_every_model_is_covered_by_a_unit_test():
    # guard: adding a model requires adding a test; ids tested above
    tested = {"ec_to_nacl", "irrigation_schedule", "crop_water_demand_generic", "crop_water_demand_alfalfa", "fao29_leaching",
              "twin_solution_to_ece", "fitted_leaching_fresh", "fitted_leaching_dry", "sa187_benefit", "mulch_evaporation",
              "shade_net_temperature", "biochar_porosity", "inoculant", "knallgas_carbon_fixation", "co_oxidation_carbon_balance",
              "h2_atmosphere", "solar_park_microclimate", "panel_washing", "cloud_seeding", "offroad_compaction", "climate_warming",
              "maas_hoffman_date_palm", "maas_hoffman_barley", "maas_hoffman_alfalfa", "maas_hoffman_tomato"}
    assert tested == set(M)


def test_alfalfa_chain_b_and_twin_abstention():
    ups = [M["ec_to_nacl"], M["crop_water_demand_alfalfa"], M["irrigation_schedule"]]
    downs = [M["twin_solution_to_ece"], M["maas_hoffman_alfalfa"], M["sa187_benefit"]]
    r = run_chain({"ec_w_dS_m": 3.12, "irrigation_interval_days": 1}, ups, twin_request(FLAT), downs)
    assert r["refused"] is None
    ry, ryi = (r["outputs"][k] for k in ("maas_hoffman_alfalfa:relative_yield_alfalfa", "sa187_benefit:relative_yield_alfalfa_sa187"))
    assert ryi["value"] == pytest.approx(ry["value"] * 1.15) and 0 <= ry["value"] <= 1
    # very saline water drives the soil solution beyond the twin's ideal-dilute domain: the twin abstains and the chain refuses downstream
    r = run_chain({"ec_w_dS_m": 16., "irrigation_interval_days": 1}, ups, twin_request(FLAT), downs, envelope=False)
    assert r["outputs"]["twin:h2_uptake_mol_m2"]["value"] is None and r["outputs"]["twin:h2_uptake_mol_m2"]["low"] is None
    h2 = run_chain({"ec_w_dS_m": 16., "irrigation_interval_days": 1}, ups, twin_request(FLAT, community={"sample": "10Dr2"}, coupling="gate"),
                   [M["knallgas_carbon_fixation"]], envelope=False)
    assert h2["refused"]["model"] == "knallgas_carbon_fixation" and "unavailable" in h2["refused"]["reason"]
