"""Scientific contracts for the current public model, including falsification cases."""
from dataclasses import replace
from pathlib import Path
import json

import numpy as np
import pandas as pd
import pytest
from fastapi.testclient import TestClient

from eq_biodt.v03.api import app, DATA
from eq_biodt.v03.contract import contract
from eq_biodt.v03.evidence import CBB_STEPS, sha256
from eq_biodt.v03.forecast import backtest, freeze, score
from eq_biodt.v03.ontology import evaluate_ontology_reasoning
from eq_biodt.v03.planning import compare, experiment, parameter_ensemble
from eq_biodt.v03.process import Parameters, State, deposition_flux, gas_exchange, simulate
from eq_biodt.v03.resources import resource_ceiling, carbon_account
from eq_biodt.v03.semantics import Observation, additive_total, carrier_partition, validate_composition
from eq_biodt.models.causal_ontological_twin import OntologicalPolytope


P = Parameters(k_h2_s=.001, k_co_s=.001)


def test_dry_activation_and_wet_diffusion_limit():
    theta = np.linspace(0, P.porosity, 401)
    uptake = np.array([gas_exchange(x, 25, P)["h2_mol_m2_s"] for x in theta])
    assert uptake[0] == uptake[-1] == 0
    assert 0 < uptake.argmax() < len(uptake)-1
    assert uptake.max() > uptake[10] and uptake.max() > uptake[-2]
    assert np.any(np.diff(uptake) > 0) and np.any(np.diff(uptake) < 0)


def test_gas_concentration_oxygen_and_maintenance_are_connected():
    base = gas_exchange(.08, 25, P)
    absent = gas_exchange(.08, 25, replace(P, h2_mol_m3=0, co_mol_m3=0))
    doubled = gas_exchange(.08, 25, replace(P, h2_mol_m3=2*P.h2_mol_m3))
    assert absent["chemical_power_w_m2"] == 0
    assert doubled["h2_mol_m2_s"] == pytest.approx(2*base["h2_mol_m2_s"])
    assert gas_exchange(.08, 25, replace(P, oxygen_supply_mol_m2_s=0))["chemical_power_w_m2"] == 0
    demand = replace(P, maintenance_w_m2=1e-5)
    assert gas_exchange(.08, 35, demand)["maintenance_coverage"] != gas_exchange(.08, 25, demand)["maintenance_coverage"]
    assert base["maintenance_coverage"] is None
    assert gas_exchange(.08, 25, replace(P, d_h2_m2_s=P.d_h2_m2_s/100))["h2_mol_m2_s"] < base["h2_mol_m2_s"]


def test_diffusion_reaction_has_correct_limiting_solution():
    # For a shallow layer the analytical flux tends to concentration*k*depth.
    assert deposition_flux(2e-5, 1e-4, .001, 1e-6) == pytest.approx(2e-14, rel=1e-7)
    assert deposition_flux(0, 1e-4, .001, .1) == 0
    with pytest.raises(ValueError):
        deposition_flux(1, -1, 1, 1)


@pytest.mark.parametrize("kwargs", [{"porosity": 1}, {"k_h2_s": -1}, {"k_h2_s": float("nan")},
                                  {"energy_capture_fraction": 2}, {"maintenance_w_m2": 0}])
def test_invalid_parameters_rejected(kwargs):
    with pytest.raises(ValueError):
        replace(P, **kwargs)


def test_water_and_salt_conservation_with_shared_memory():
    rng = np.random.default_rng(7)
    forcing = [{"water_mm": float(rng.uniform(0, 40)), "nacl_mol_l": .03,
                "evaporation_mm": float(rng.uniform(0, 6)), "temp_c": 30} for _ in range(50)]
    initial = State(1.8, .1)
    out = simulate(initial, forcing, P)
    rows = out["trajectory"]
    assert rows[-1]["water_mm"] + sum(r["outflow_mm"]+r["evaporation_mm"] for r in rows) == pytest.approx(initial.water_mm+sum(f["water_mm"] for f in forcing))
    assert rows[-1]["salt_mol_m2"]+sum(r["salt_outflow_mol_m2"] for r in rows) == pytest.approx(initial.salt_mol_m2+sum(f["water_mm"]*f["nacl_mol_l"] for f in forcing))
    assert all(0 <= r["water_mm"] <= 40 for r in rows)


def test_salt_abstention_and_unsupported_species():
    out = simulate(State(2, 2), [{"evaporation_mm": 1, "temp_c": 25}], P)
    assert out["summary"]["h2_uptake_mol_m2"] is None
    assert out["trajectory"][0]["salt_mol_m2"] == 2
    assert out["summary"]["net_co2_removal_g_m2"] is None
    with pytest.raises(ValueError, match="Unsupported causal"):
        simulate(State(2), [{"temp_c": 25, "evaporation_mm": 1, "shrub_species": "Haloxylon"}], P)


def test_ontological_maps_preserve_endpoints_missingness_and_overlaps():
    a = Observation(1.2, "GE", "SCG", "reads", "gene_copy", "1l")
    b = Observation(1.2, "GE", "SCG", "reads", "gene_copy", "1h")
    missing = Observation(None, "GE", "SCG", "reads", "gene_copy", "2a")
    assert additive_total([a, b], disjoint=True).value == 2.4
    with pytest.raises(ValueError):
        additive_total([a, a], disjoint=True)
    assert additive_total([a, missing], disjoint=True).value is None
    with pytest.raises(ValueError):
        additive_total([a, Observation(.2, "fraction", "mapped_reads", "coverm", "composition")], disjoint=True)
    with pytest.raises(ValueError):
        additive_total([a, a], disjoint=False)
    w = {"a": .5, "b": .2, "c": .3}
    assert carrier_partition(w, {"a", "b"})["share"] == .7
    assert carrier_partition(w, {"a", "c"})["share"] == .8
    assert validate_composition({"one_phylum": 1.})["unassigned"] == 0
    partial = validate_composition({"known": .8, "unknown": None})
    assert partial["parts"]["unknown"] is None and partial["status"] == "partial"


def test_legacy_projection_missing_and_invalid_metric_fail():
    p = OntologicalPolytope(["x", "y"])
    with pytest.raises(ValueError):
        p.project({"x": 1})
    with pytest.raises(ValueError):
        p.project({"x": 1, "y": 2}, np.array([[1, 2], [2, 1]]))


def test_resource_network_balances_atoms_and_requires_maintenance():
    args = dict(h2=1e-8, co=1e-9, oxygen=1e-8, external_co2=1e-8,
                capture_fraction=.25, maintenance_w_m2=1e-5, synthesis_j_mol_c=500000,
                pathway_supported=True)
    r = resource_ceiling(**args)
    x = r["fluxes"]
    assert x["h2_respiration"]+2*x["carbohydrate_c"] <= args["h2"]+1e-16
    assert .5*(x["h2_respiration"]+x["co_respiration"]) <= args["oxygen"]+1e-16
    e = r["net_species_exchange_mol_m2_s"]
    assert e["CO"]+e["CO2"]+e["CH2O"] == pytest.approx(0, abs=1e-20)
    assert 2*e["H2"]+2*e["H2O"]+2*e["CH2O"] == pytest.approx(0, abs=1e-20)
    assert r["energy_slack_w_m2"] >= -1e-12
    assert resource_ceiling(**(args|{"maintenance_w_m2": 1.}))["status"] == "infeasible_resource_or_maintenance_demand"
    assert resource_ceiling(**(args|{"pathway_supported": False}))["status"] == "insufficient_pathway_evidence"


def test_net_carbon_needs_sources_and_counts_imports_once():
    args = dict(stock_change_g_c_m2=20, imported_carbon_g_c_m2=10, geological_carbon_g_c_m2=2,
                lifecycle_g_co2e_m2=5, counterfactual_stock_change_g_c_m2=3,
                non_co2_g_co2e_m2=1, boundary="plot A, all ecosystem pools", interval="2026-2027")
    with pytest.raises(ValueError):
        carbon_account(**args, evidence={})
    r = carbon_account(**args, evidence={k:"test measurement" for k in ["stock_change","imports","geological_sources","lifecycle","counterfactual","non_co2"]})
    assert r["net_balance_g_co2e_m2"] == pytest.approx(5*44/12-6)


def test_complete_metagenomic_pathway_screen_and_artifact_identity():
    manifest = json.loads((DATA/"manifest.json").read_text())
    for filename, h in manifest["artifacts"].items():
        assert sha256(DATA/filename) == h
    mag = pd.read_csv(DATA/"mag_evidence.tsv", sep="\t")
    assert len(mag) == 5438
    assert (mag.cbb_all_steps_detected == mag[list(CBB_STEPS)].all(axis=1)).all()
    assert ((mag.cbb_steps_detected == 11) == mag.cbb_all_steps_detected).all()
    assert (mag.rbcl_prk_detected & ~mag.cbb_all_steps_detected).any()
    assert not (mag.h2_cbb_all_steps & ~mag.h2_carrier_union).any()
    panel = pd.read_csv(DATA/"sample_evidence.tsv", sep="\t")
    assert panel.loc[panel["sample"].eq("F5Dr1"), "mapped_mag_shannon"].isna().all()


def test_forecast_no_future_target_leakage(tmp_path):
    panel = pd.read_csv(DATA/"sample_evidence.tsv", sep="\t")
    result = backtest(panel)
    changed = panel.copy()
    changed.loc[changed.trip_num.eq(5), "total_h2_marker_ge"] *= 100
    after = backtest(changed)
    before_values = [r["predicted_log1p_ge"] for r in result["predictions"]]
    assert before_values == [r["predicted_log1p_ge"] for r in after["predictions"]]
    path = tmp_path/"fixture-only.json"
    frozen = freeze(DATA/"sample_evidence.tsv", "2027-01-01", path, issued_on="2026-10-27")
    assert frozen["status"] == "test_fixture"
    assert score(frozen, panel)["status"] == "awaiting_matching_future_observations"
    with pytest.raises(FileExistsError):
        freeze(DATA/"sample_evidence.tsv", "2027-01-01", path, issued_on="2026-10-27")
    future = pd.DataFrame([{ "site": r["site"], "compartment": r["compartment"],
                            "collection_date": "2027-01-01", "total_h2_marker_ge": 1.}
                           for r in frozen["predictions"]])
    assert score(frozen, future)["status"] == "scored"
    frozen["target_date"] = "2028-01-01"
    with pytest.raises(ValueError, match="checksum"):
        score(frozen, future)


def test_design_holdout_and_robust_conservation_baseline():
    d = experiment(P)
    assert len(d["schedule"]) == 54
    assert all(r["role"] == "held_out_condition" for r in d["schedule"] if r["target_water_potential_mpa"] == -3)
    force = [{"temp_c":25, "evaporation_mm":1, "water_mm":0} for _ in range(14)]
    scenarios = [{"name":"intact", "forcing":force, "disturbance_m2":0},
                 {"name":"water", "forcing":[f|{"water_mm":1} for f in force], "disturbance_m2":1}]
    r = compare(State(1.8), scenarios, parameter_ensemble(P))
    assert len(r["parameters"]) == 24
    assert "intact" in r["pareto_in_every_member"]


def test_public_api_and_explicit_retirement():
    client = TestClient(app)
    assert client.get("/health").json()["version"] == "0.5.0"
    assert client.get("/v3/evidence?limit=1").json()["manifest"]["reference_MAGs"] == 5438
    assert client.post("/v1/intervene", json={}).status_code == 410
    assert "What happens after the rain?" in client.get("/").text
    request = {"parameters":{"k_h2_s":.001,"k_co_s":.001}, "initial":{"water_mm":1.8},
               "forcing":[{"water_mm":12,"evaporation_mm":1,"temp_c":25}]}
    assert client.post("/v3/simulate", json=request).status_code == 200
    request["parameters"]["shrub_species"] = "x"
    assert client.post("/v3/simulate", json=request).status_code == 422
    ids = [c["id"] for c in contract()["clauses"]]
    assert len(ids) == len(set(ids))


def test_resource_api_uses_actual_mag_evidence_and_comparison_validates_inputs():
    client = TestClient(app)
    mag = pd.read_csv(DATA/"mag_evidence.tsv", sep="\t")
    args = {"h2":1e-8, "co":1e-9, "oxygen":1e-8, "external_co2":1e-8,
            "capture_fraction":.25, "maintenance_w_m2":1e-5, "synthesis_j_mol_c":500000}
    for support in (True, False):
        ref = mag.loc[mag.h2_cbb_all_steps.eq(support)].iloc[0].rep
        response = client.post("/v3/resource-ceiling", json=args|{"mag_id":ref})
        assert response.status_code == 200
        assert response.json()["genomic_evidence"]["h2_and_all_cbb_steps"] == support
        assert (response.json()["fluxes"] is not None) == support
    assert client.post("/v3/resource-ceiling", json=args|{"pathway_supported":True}).status_code == 422
    assert client.post("/v3/compare", json={"initial":{}, "parameters":{}, "scenarios":[{}]}).status_code == 422


def test_owl2_reasoner_classification_justifications_and_frechet_bounds(tmp_path):
    ofn = tmp_path / "eq_biodt_v03.ofn"
    mag = pd.read_csv(DATA / "mag_evidence.tsv", sep="\t")
    sample = pd.read_csv(DATA / "sample_evidence.tsv", sep="\t")
    res = evaluate_ontology_reasoning(mag, sample, ofn)
    assert res["consistent"] is True
    assert res["unsatisfiable_classes"] == 0
    assert res["dropped_clauses"] == 0
    assert res["exact_agreement_all_5438_mags"] is True
    assert res["entailed_class_counts"]["CBBCompleteMAG"] == 232
    assert res["entailed_class_counts"]["ArchaealRuBisCOMAG"] == 241
    assert res["entailed_class_counts"]["H2AutotrophicCandidateMAG"] == 136
    assert res["entailed_class_counts"]["HighQualityH2AutotrophicCandidateMAG"] == 100
    assert res["entailed_class_counts"]["DualH2COAutotrophicCandidateMAG"] == int(mag.dual_h2_co_cbb_all_steps.sum())
    assert res["compartment_niche_partitioning"]["PlantRoot"]["samples_exceeding_1pct_h2_cbb"] == 151
    assert res["compartment_niche_partitioning"]["Deep"]["samples_exceeding_1pct_h2_cbb"] == 139
    assert res["compartment_niche_partitioning"]["Surface"]["samples_exceeding_1pct_h2_cbb"] == 76
    just = res["example_candidate_justification"]
    assert just["verified"] is True
    assert just["subset_minimal"] is True
    assert just["axiom_count"] > 0
    disj = res["endpoint_disjointness_justification"]
    assert disj["verified"] is True
    assert disj["axiom_count"] == 3
    ab3 = res["ablations"]["sample_marginals_without_genome_linkage"]
    assert ab3["linked_within_frechet_bounds_all_samples"] is True
    assert ab3["nonzero_ambiguity_samples"] == 549


def test_conformal_intervals_and_mags_ontology_api():
    panel = pd.read_csv(DATA / "sample_evidence.tsv", sep="\t")
    bt = backtest(panel)
    ridge_preds = [r for r in bt["predictions"] if r["model"] == "geography_season_ridge"]
    assert sum(r["covered_90"] for r in ridge_preds) == 237
    assert len(ridge_preds) == 278
    assert all(r["lower_90_log1p_ge"] <= r["predicted_log1p_ge"] <= r["upper_90_log1p_ge"] for r in bt["predictions"])
    changed = panel.copy()
    changed.loc[changed.trip_num.eq(5), "total_h2_marker_ge"] *= 100
    after = backtest(changed)
    before_q = [r["q90_log1p_ge"] for r in bt["predictions"] if r["campaign"] == 5]
    after_q = [r["q90_log1p_ge"] for r in after["predictions"] if r["campaign"] == 5]
    assert before_q == after_q

    client = TestClient(app)
    mags_resp = client.get("/v3/mags?candidates_only=true&limit=50").json()
    assert mags_resp["h2_cbb_candidate_mags"] == 136
    assert len(mags_resp["mags"]) == 50
    assert all(m["h2_cbb_all_steps"] for m in mags_resp["mags"])
    assert "K01601" in mags_resp["mags"][0]["cbb_kos"]
    ont_resp = client.get("/v3/ontology").json()
    assert ont_resp["consistent"] is True


def test_v05_external_source_tampering_rejected(tmp_path):
    import sys
    repo_root = Path(__file__).resolve().parents[1]
    sys.path.insert(0, str(repo_root / "reproducibility/v05"))
    from external_transfer import verified

    ext_root = repo_root / "data/external"
    (tmp_path / "pmc-supplement-manifest.json").write_bytes((ext_root / "pmc-supplement-manifest.json").read_bytes())
    (tmp_path / "buzzard2022").mkdir()
    src = ext_root / "buzzard2022/peerj-10-12966-s007.csv"
    dst = tmp_path / "buzzard2022/peerj-10-12966-s007.csv"
    dst.write_bytes(src.read_bytes())
    assert verified("buzzard2022", "peerj-10-12966-s007.csv", data_dir=tmp_path) == dst
    # Same-length byte tampering must be rejected by digest verification.
    body = bytearray(src.read_bytes())
    body[-5] ^= 0x01
    dst.write_bytes(bytes(body))
    with pytest.raises(ValueError, match="digest mismatch"):
        verified("buzzard2022", "peerj-10-12966-s007.csv", data_dir=tmp_path)


def test_v05_grouped_holdout_boundaries_and_no_test_leakage():
    import sys
    repo_root = Path(__file__).resolve().parents[1]
    sys.path.insert(0, str(repo_root / "reproducibility/v05"))
    from external_transfer import evaluate_sonoran_fold, evaluate_temperate_folds

    obs = pd.read_csv(repo_root / "reviews/v05/sonoran_field_observations.tsv", sep="\t")
    held = "2018-08-02"
    train = obs[obs.campaign.ne(held)].copy()
    test = obs[obs.campaign.eq(held)].copy()
    # Overlapping campaign rows between train and test must be rejected.
    with pytest.raises(ValueError, match="disjoint"):
        evaluate_sonoran_fold(obs, test)
    preds_before, fit_before = evaluate_sonoran_fold(train, test)
    perturbed_test = test.copy()
    perturbed_test["y"] = perturbed_test["y"] * -50.0 + 100.0
    preds_after, fit_after = evaluate_sonoran_fold(train, perturbed_test)
    assert fit_before == fit_after
    for model in preds_before:
        assert np.array_equal(preds_before[model], preds_after[model])

    rates = pd.read_csv(repo_root / "reviews/v03/external_bottle_rates.tsv", sep="\t")
    assert (rates.bottle_k_min < 0).any()
    temp_before = evaluate_temperate_folds(rates)
    perturbed_rates = rates.copy()
    perturbed_rates.loc[perturbed_rates.soil.eq("Meadow"), "bottle_k_min"] += 10.0
    temp_after = evaluate_temperate_folds(perturbed_rates)
    meadow_before = temp_before.loc[temp_before.heldout_soil.eq("Meadow"), "predicted_k_min"].to_numpy()
    meadow_after = temp_after.loc[temp_after.heldout_soil.eq("Meadow"), "predicted_k_min"].to_numpy()
    assert np.array_equal(meadow_before, meadow_after)


def test_v05_native_unit_conversions_and_unquantified_endpoints():
    repo_root = Path(__file__).resolve().parents[1]
    # Incompatible external/internal units must not be pooled.
    areal = Observation(1.5, "nmol H2 m-2 s-1", "collar", "buzzard2022", "areal_flux", "field")
    bottle_mass = Observation(0.8, "nmol h-1 g-1", "soil_dry_g", "tribbia2026", "bottle_rate", "wet")
    bottle_k = Observation(0.02, "min-1", "headspace", "reji2025", "first_order_rate", "forest")
    with pytest.raises(ValueError):
        additive_total([areal, bottle_mass], disjoint=True)
    with pytest.raises(ValueError):
        additive_total([bottle_mass, bottle_k], disjoint=True)

    # Namib dry-soil gravimetric water addition and unquantified dry endpoints.
    namib_pred = pd.read_csv(repo_root / "reviews/v05/namib_conditional_predictions.tsv", sep="\t")
    dry_soil_g = 2.0 / (1.0 + namib_pred.dry_gwc)
    water_g = 2.0 - dry_soil_g + 0.4
    assert np.allclose(namib_pred.wet_gwc, water_g / dry_soil_g, rtol=1e-12, atol=1e-12)
    namib_contrasts = pd.read_csv(repo_root / "reviews/v05/namib_hydration_contrasts.tsv", sep="\t")
    unquantified = namib_contrasts.loc[namib_contrasts.n_dry_rates.eq(0)]
    assert len(unquantified) == 8
    assert unquantified[["dry_mean", "fold", "fold_lcl", "fold_ucl"]].isna().all().all()
