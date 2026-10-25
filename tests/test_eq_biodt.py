"""Comprehensive verification suite for the Empty Quarter Biodiversity Digital Twin (EQ-BioDT)."""

from __future__ import annotations

from pathlib import Path

import pytest
from fastapi.testclient import TestClient

from eq_biodt.api.main import app
from eq_biodt.data import (
    CanonicalDataRegistry,
    OntologyRegistry,
    TwinStateCubeBuilder,
)
from eq_biodt.interventions import (
    CloudSeedingIntervention,
    DesalinationBrineIntervention,
    MegafaunaGrazingIntervention,
    ParetoPortfolioOptimizer,
    RefugiumThermalIntervention,
    RevegetationIntervention,
    SeismicHydrocarbonIntervention,
)
from eq_biodt.models import (
    CarbonBudgetEngine,
    CausalOntologicalTwinEngine,
    CommunityMetabolicInputs,
    MetabolicCommunityEngine,
    MoisturePulseDiffusionModel,
    SpatialJSDMViromeModel,
)
from eq_biodt.standards import (
    B3OccurrenceCubeBuilder,
    DwCEMoFExporter,
    FDOProvenanceRecorder,
)
from eq_biodt.validation import TwinValidationSuite


def test_canonical_data_and_ontology_registry() -> None:
    reg = CanonicalDataRegistry()
    report = reg.verify_all_sources()
    assert len(report) == 23

    ont = reg.ontology
    # Verify programmatic ChEBI & RAK quality/value IRIs
    assert ont.get_chebi_iri("Si") == "http://purl.obolibrary.org/obo/CHEBI_27573"
    assert ont.get_chebi_iri("V") == "http://purl.obolibrary.org/obo/CHEBI_27698"
    q_si, v_si = ont.get_rak_analyte_iris("Si")
    assert q_si == "https://rubalkhali.science/kb/RAK_0000100"
    assert v_si == "https://rubalkhali.science/kb/RAK_0000500"
    q_ti, v_ti = ont.get_rak_analyte_iris("Ti")
    assert q_ti == "https://rubalkhali.science/kb/RAK_0000106"
    assert v_ti == "https://rubalkhali.science/kb/RAK_0000506"
    q_al, v_al = ont.get_rak_analyte_iris("Al")
    assert q_al == "https://rubalkhali.science/kb/RAK_0000107"
    assert v_al == "https://rubalkhali.science/kb/RAK_0000507"
    q_le, v_le = ont.get_rak_analyte_iris("LE")
    assert q_le == "https://rubalkhali.science/kb/RAK_0000032"
    assert v_le == "https://rubalkhali.science/kb/RAK_0000033"

    # Verify ENVO terms from biome_codes.yml
    assert ont.get_envo_iri("desert biome") == "http://purl.obolibrary.org/obo/ENVO_01000179"
    assert ont.get_envo_iri("oilspill") == "http://purl.obolibrary.org/obo/ENVO_00002061"
    assert ont.get_envo_iri("fecal material") == "http://purl.obolibrary.org/obo/ENVO_00002003"

    # Verify canonical auxiliary sites and Trip-1 alias resolution (61..64 -> RAK_1000067..1000070)
    assert ont.site_iris["Site oill spill"] == "https://rubalkhali.science/kb/RAK_1000061"
    assert ont.site_iris["Site Camel1"] == "https://rubalkhali.science/kb/RAK_1000064"
    assert ont.get_site_iri("61") == "https://rubalkhali.science/kb/RAK_1000067"
    assert ont.get_site_iri("64") == "https://rubalkhali.science/kb/RAK_1000070"

    # Verify programmatic NCBITaxon domain and phylum IRIs from taxonomy_mapping_decisions.tsv
    assert ont.get_taxon_iri("Bacteria") == "http://purl.obolibrary.org/obo/NCBITaxon_2"
    assert ont.get_taxon_iri("Archaea") == "http://purl.obolibrary.org/obo/NCBITaxon_2157"
    assert ont.get_taxon_iri("Actinomycetota") == "http://purl.obolibrary.org/obo/NCBITaxon_201174"
    assert ont.get_taxon_iri("Thermoplasmatota") == "http://purl.obolibrary.org/obo/NCBITaxon_2283796"
    assert ont.get_taxon_iri("Halobacteriota") == "http://purl.obolibrary.org/obo/NCBITaxon_3100166"

    # Verify aquatic water chemistry loader from Table S2
    aq = reg.load_aquatic_water_chemistry()
    assert len(aq) == 27
    assert {"59-MR-WR1", "59-MR-WR2", "60-MR-WR1"} <= set(aq.columns)


def test_state_cubes_and_eu_biodt_standards(tmp_path: Path) -> None:
    builder = TwinStateCubeBuilder()
    sites_df = builder.build_sites_cube()
    comp_df = builder.build_compartment_state_cube()
    visit_df = builder.build_visit_compartment_cube()
    grid_df = builder.build_era5_h3_grid()

    assert len(sites_df) == 70  # 60 primary + 10 auxiliary sites
    assert len(comp_df) == 180  # 60 primary sites x 3 compartments
    assert len(visit_df) > 400  # Longitudinal visit-compartment cube across Trips 1..5
    assert len(grid_df) > 200
    assert "in_empirical_support" in grid_df.columns
    assert "abstention_reason" in grid_df.columns
    assert "h3_5" in sites_df.columns and "h3_6" in sites_df.columns
    # Verify PlantRoot normalization recovered 59/60 RootAdjacent WGS observations
    root_df = comp_df[comp_df["compartment_norm"] == "RootAdjacent"]
    assert int(root_df["wgs_observed"].sum()) == 59
    assert int(comp_df["wgs_observed"].sum()) == 155
    # Verify 1,227 primary-site QC-passed 16S profiles (372 Surface, 381 ShallowSubsurface, 474 RootAdjacent)
    assert int(comp_df["amplicon_16s_profile_count"].sum()) == 1227
    assert int(comp_df.loc[comp_df["compartment_norm"] == "Surface", "amplicon_16s_profile_count"].sum()) == 372
    assert int(comp_df.loc[comp_df["compartment_norm"] == "ShallowSubsurface", "amplicon_16s_profile_count"].sum()) == 381
    assert int(comp_df.loc[comp_df["compartment_norm"] == "RootAdjacent", "amplicon_16s_profile_count"].sum()) == 474
    assert "prescribed_portfolio" in grid_df.columns
    assert int((grid_df["prescribed_portfolio"] != "ABSTAIN").sum()) == 128
    # Verify unimputed *_raw columns and empirical phylum columns exist
    assert "hydB_1l_ge_raw" in comp_df.columns
    assert "actinomycetota_rel_abund" in comp_df.columns
    assert "actino_16s_rel_abund" in comp_df.columns
    assert "unmapped_read_pct" in comp_df.columns
    # Verify Sites 62 and 64 have integer site_num in sites_df
    assert 62 in set(sites_df["site_num"].dropna().astype(int))
    assert 64 in set(sites_df["site_num"].dropna().astype(int))

    # Test DwC-A + eMoF export
    dwc = DwCEMoFExporter(builder)
    dwc_files = dwc.export_archive(tmp_path / "dwca")
    assert dwc_files["event"].exists()
    assert dwc_files["occurrence"].exists()
    assert dwc_files["emof"].exists()

    # Test B-Cubed occurrence cube export
    b3 = B3OccurrenceCubeBuilder(builder)
    b3_files = b3.export_b3_datapackage(tmp_path / "b3")
    assert b3_files["datapackage_json"].exists()
    assert b3_files["cube_parquet"].exists()


def test_two_layer_causal_ontological_engine() -> None:
    engine = CausalOntologicalTwinEngine()

    # Verify formal OWL 2 DL specification compiled into OntologicalPolytope
    assert len(engine.polytope.partitions) == 5
    assert len(engine.polytope.subsumptions) == 13
    # Verify fail-closed behavior on empty dict or missing coordinates
    assert engine.polytope.compute_ocs([{}]) == 0.0
    assert any(v.startswith("EmptyState:") for v in engine.polytope.check_violations({}))
    assert any(v.startswith("MissingCoordinate:") for v in engine.polytope.check_violations({"hydB_1l_ge": 0.5}))

    # Intentionally inconsistent raw prediction (subsumption violation: hydB_1l > total_h2ase)
    violated = {
        "actinomycetota_rel_abund": 0.60,
        "thermoplasmatota_rel_abund": 0.30,
        "halobacteriota_rel_abund": 0.25,
        "other_prokaryotes_rel_abund": 0.10,
        "total_prokaryotic_rel_abund": 0.80,  # sum of children (1.25) > parent (0.80)
        "hydB_1l_ge": 0.90,
        "hydB_1h_ge": 0.20,
        "hydB_2a_ge": 0.05,
        "total_h2ase_ge": 0.40,               # violated!
        "coxL_formI_ge": 1.10,
        "total_trace_gas_scavenger_ge": 0.90, # violated!
        "h2_oxidation_power_uw_m3": 20.0,
        "co_oxidation_power_uw_m3": 2.5,
        "total_trace_gas_power_uw_m3": 10.0,  # violated!
        "viable_dna_fraction": 0.60,
        "relic_dna_fraction": 0.60,
        "total_dna_fraction": 1.0,
    }
    assert engine.polytope.compute_ocs([violated]) == 0.0

    # Project through constitutive layer -> OCS must equal 1.0
    repaired = engine.evaluate_constitutive_layer(violated)
    assert repaired["ocs_score"] == 1.0
    assert len(engine.polytope.check_violations(repaired)) == 0

    # Verify NaN / inf / empty rejection in OntologicalPolytope and evaluate_constitutive_layer
    nan_state = dict(violated)
    nan_state["hydB_1l_ge"] = float("nan")
    assert any(v.startswith("NonFinite:") for v in engine.polytope.check_violations(nan_state))
    assert engine.polytope.compute_ocs([nan_state]) == 0.0
    with pytest.raises(ValueError, match="non-finite|NaN"):
        engine.polytope.project(nan_state)
    with pytest.raises(ValueError, match="empty"):
        engine.polytope.project({})
    with pytest.raises(ValueError, match="empty"):
        engine.evaluate_constitutive_layer({})
    with pytest.raises(ValueError, match="Non-finite"):
        engine.evaluate_constitutive_layer({"hydB_1l_ge": float("nan")})

    # Test Proposition 2 Agreement-Weighted Class Transport with DerSimonian-Laird tau2 & abundance weights
    transport = engine.agreement_weighted_class_transport(
        {
            "group_1l_NiFe": (-0.79, 0.08),
            "group_1h_NiFe": (-0.42, 0.12),
            "group_2a_NiFe": (-0.31, 0.15),
        },
        baseline_abundances={"group_1l_NiFe": 0.65, "group_1h_NiFe": 0.20, "group_2a_NiFe": 0.15},
    )
    assert transport["sign_concordance"] == 1.0
    assert transport["pooled_effect"] < 0.0
    assert transport["dersimonian_laird_tau2"] >= 0.0
    assert transport["abundance_weighted_parent_log2fc"] < 0.0


def test_moisture_pulse_and_highs_fba_models() -> None:
    pulse = MoisturePulseDiffusionModel()
    traj = pulse.simulate_pulse_trajectory(rain_pulse_mm=8.0, horizon_days=14)
    # Distributed-lag peak should occur within days 2-3 (1-4 d window)
    assert 1 <= traj["peak_day"] <= 4
    assert traj["delta_shannon_peak"] > 0.0

    # Extreme saturation (35 mm) throttles H2 oxidation power compared to moderate pulse (6 mm)
    mod = pulse.simulate_pulse_trajectory(rain_pulse_mm=6.0)
    sat = pulse.simulate_pulse_trajectory(rain_pulse_mm=35.0)
    assert sat["min_h2_power_uw_m3"] < mod["min_h2_power_uw_m3"]

    # Mechanical compaction (phi: 0.40 -> 0.22) must monotonically reduce diffusivity and H2 energy
    uncompacted = pulse.simulate_pulse_trajectory(rain_pulse_mm=0.0, total_porosity=0.40)
    compacted = pulse.simulate_pulse_trajectory(rain_pulse_mm=0.0, total_porosity=0.22)
    assert compacted["min_h2_power_uw_m3"] < uncompacted["min_h2_power_uw_m3"]
    assert compacted["integrated_14d_h2_energy_j_m3"] < uncompacted["integrated_14d_h2_energy_j_m3"]

    # HiGHS FBA: Site 60 brackish Ca-SO4 vs Site 59 hypersaline Na-Cl
    fba = MetabolicCommunityEngine()
    s60 = fba.solve_community_fba(
        CommunityMetabolicInputs(
            salinity_ms_cm=10.65,
            sulfate_so4_bound=0.45,
            rbcl_hyab_plasmid_inoc=0.8,
        )
    )
    s59 = fba.solve_community_fba(
        CommunityMetabolicInputs(
            h2_uptake_bound=0.14,
            co_uptake_bound=0.03,
            salinity_ms_cm=210.0,
            sulfate_so4_bound=0.05,
            rbcl_hyab_plasmid_inoc=0.0,
        )
    )
    assert s60["solver_status"] == "optimal"
    assert s60["hard_maintenance_feasible"] is True
    assert s59["hard_maintenance_feasible"] is False
    assert s59["solver_status"] == "infeasible_maintenance_deficit"
    assert s59["native_highs_status"] == 2
    assert s59["total_biomass_flux"] == 0.0
    assert s60["autotrophic_carbon_fixation_flux"] > s59["autotrophic_carbon_fixation_flux"]
    assert s60["dark_proteome_active_fraction"] > s59["dark_proteome_active_fraction"]

    # Dormancy fallback with positive DOC bound must keep heterotrophic growth at 0 and balance catabolic ATP
    s59_doc = fba.solve_community_fba(
        CommunityMetabolicInputs(
            h2_uptake_bound=0.14,
            co_uptake_bound=0.03,
            root_exudate_doc_bound=0.25,
            exogenous_nh4_no3_bound=0.10,
            salinity_ms_cm=210.0,
            sulfate_so4_bound=0.05,
        )
    )
    assert s59_doc["hard_maintenance_feasible"] is False
    assert s59_doc["heterotrophic_doc_flux"] == 0.0
    assert s59_doc["total_biomass_flux"] == 0.0
    expected_cat_atp = 1.4 * 0.14 + 0.9 * 0.03 + 0.5 * 0.10 + 0.8 * min(0.05 * 2.718281828**(-(210.0 - 25.0)/85.0), 0.50 * 0.14 + 0.35 * 0.03)
    assert abs(s59_doc["atpm_maintenance_drain"] - expected_cat_atp) < 1e-4

    # Sulfate without electron donors (H2=0, CO=0, DOC=0, NH4=0) must produce zero growth
    so4_no_donor = fba.solve_community_fba(
        CommunityMetabolicInputs(
            h2_uptake_bound=0.0,
            co_uptake_bound=0.0,
            root_exudate_doc_bound=0.0,
            exogenous_nh4_no3_bound=0.0,
            sulfate_so4_bound=0.50,
            salinity_ms_cm=1.0,
            rbcl_hyab_plasmid_inoc=0.0,
        )
    )
    assert so4_no_donor["total_biomass_flux"] <= 1e-7
    assert so4_no_donor["sulfur_redox_flux"] <= 1e-7

    # Spatial JSDM Virome: null intervention must preserve baseline viral and plasmid richness
    jsdm = SpatialJSDMViromeModel()
    res_null = jsdm.predict_virome_mobilome_state(
        shannon_16s=5.8,
        compartment="Surface",
        map_mm=18.5,
        xrf_evaporite_index=-4.2,
        vegetation_cover_frac=0.0,
        exogenous_manure_load=0.0,
        plasmid_bioaug_dose=0.0,
        baseline_viral_richness=18.0,
        baseline_plasmid_richness=31.0,
        baseline_shannon_16s=5.8,
    )
    assert abs(res_null["viral_richness_votu"] - 18.0) < 1e-9
    assert abs(res_null["plasmid_richness_potu"] - 31.0) < 1e-9


def test_all_six_interventions_and_pareto_optimizer() -> None:
    builder = TwinStateCubeBuilder()
    comp_df = builder.build_compartment_state_cube()
    site35 = comp_df[(comp_df["site_num"] == 35) & (comp_df["compartment_norm"] == "Surface")].iloc[0].to_dict()

    r1 = CloudSeedingIntervention().simulate(site35, rain_pulse_mm=7.0)
    assert r1["summary"]["ocs_score"] == 1.0

    r2 = RevegetationIntervention().simulate(
        site35, shrub_cover_fraction=0.25, rak_pr_pgpr_dose=0.8, site60_rbcl_hyab_plasmid_dose=0.6
    )
    assert r2["summary"]["ocs_score"] == 1.0
    assert r2["summary"]["biomass_dna_multiplier"] > 2.0
    # Rhizosphere paradox: root-zone hydB_1l_ge is lower than baseline hydB_1l_ge
    assert r2["summary"]["root_zone_hydB_1l_ge"] < r2["summary"]["baseline_hydB_1l_ge"]

    r3 = DesalinationBrineIntervention().simulate(
        site35, discharge_salinity_ms_cm=10.65, ca_so4_to_nacl_ratio=1.5
    )
    assert r3["summary"]["ocs_score"] == 1.0
    assert r3["summary"]["sulfur_geochemical_regime"] == "CA_SO4_SULFUR_FLUX_OASIS"

    r4 = RefugiumThermalIntervention().simulate(
        site35, subsurface_refugium_intact_fraction=0.95, surface_temp_c=60.0, exogenous_choline_betaine_mm=0.8
    )
    assert r4["summary"]["ocs_score"] == 1.0
    assert r4["summary"]["bacillus_60c_survival_index"] > 0.5

    r5 = MegafaunaGrazingIntervention().simulate(
        site35, oryx_head_per_100km2=4.0, gazelle_head_per_100km2=6.0, camel_head_per_100km2=1.0
    )
    assert r5["summary"]["ocs_score"] == 1.0
    assert r5["summary"]["resistome_biosafety_pass"] is True

    r6 = SeismicHydrocarbonIntervention().simulate(
        site35,
        seismic_line_density_km_km2=5.0,
        low_impact_nodal_mode=True,
        hydrocarbon_tph_g_kg=12.0,
        alkb_consortium_bioaug_dose=0.9,
    )
    assert r6["summary"]["ocs_score"] == 1.0
    assert r6["summary"]["tph_biodegradation_pct"] > 50.0
    assert r6["summary"]["xrf_vanadium_chebi_iri"] == "http://purl.obolibrary.org/obo/CHEBI_27698"

    opt = ParetoPortfolioOptimizer().optimize_pareto_front(site35)
    assert len(opt["pareto_frontier"]) >= 2
    assert opt["recommended_policy"]["ocs_score"] == 1.0
    assert opt["recommended_policy"]["probabilistic_robustness"]["pareto_dominance_prob"] > 0.50

    # Test 4x2x2 common-baseline factorial brine grid
    f_grid = DesalinationBrineIntervention().evaluate_factorial_brine_grid(site35)
    assert len(f_grid) == 16
    assert any(r["hard_maintenance_feasible"] is True for r in f_grid)
    assert any(r["hard_maintenance_feasible"] is False for r in f_grid)


def test_fastapi_endpoints() -> None:
    client = TestClient(app)

    r_health = client.get("/api/v1/health")
    assert r_health.status_code == 200
    assert r_health.json()["sites_count"] == 70

    # Primary site & auxiliary site (Site oill spill RAK_1000061)
    r_s60 = client.get("/api/v1/sites/60")
    assert r_s60.status_code == 200
    r_oil = client.get("/api/v1/sites/RAK_1000061")
    assert r_oil.status_code == 200
    assert r_oil.json()["site_label"] == "Site oill spill"

    # Intervention endpoint
    r_int = client.post(
        "/api/v1/intervene/megafauna_grazing",
        json={
            "site_id": "RAK_1000064",
            "compartment": "Surface",
            "parameters": {"oryx_head_per_100km2": 4.5, "camel_head_per_100km2": 1.0},
        },
    )
    assert r_int.status_code == 200
    assert r_int.json()["summary"]["ocs_score"] == 1.0

    # Pareto optimization & carbon-budget endpoints
    r_opt = client.post("/api/v1/optimize", json={"site_id": "35", "compartment": "Surface"})
    assert r_opt.status_code == 200
    assert "recommended_policy" in r_opt.json()
    assert "carbon_budget" in r_opt.json()["recommended_policy"]

    r_cb = client.get("/api/v1/carbon-budget")
    assert r_cb.status_code == 200
    assert r_cb.json()["baseline_observed_state"]["supported_corridor_area_km2"] == 92160.0
    assert len(r_cb.json()["candidate_portfolio_budgets"]) == 6

    # RO-Crate & Web Console
    r_crate = client.get("/api/v1/standards/ro-crate")
    assert r_crate.status_code == 200
    r_ui = client.get("/twin")
    assert r_ui.status_code == 200
    assert "EQ-BioDT" in r_ui.text


def test_five_tier_quantitative_validation_metrics(tmp_path: Path) -> None:
    suite = TwinValidationSuite()
    rep = suite.run_all_validations(output_dir=tmp_path / "val")

    # Tier 1: Spatial Block CV (zero site leakage, n=550 observed guild WGS, n=549 observed mobilome WGS, n=180 site-compartments)
    t1 = rep["tier1_spatial_block_cv"]["targets"]
    assert t1["log2_hydB_1l"]["n"] == 550
    assert t1["log2_hydB_1l"]["spearman_rho"] > 0.50
    assert t1["log2_hydB_1l"]["r2"] > 0.24
    assert len(t1["log2_hydB_1l"]["spearman_rho_cluster_ci95"]) == 2
    assert t1["log1p_viral_richness"]["n"] == 549
    assert t1["shannon_16s_site_compartment"]["n"] == 180
    assert t1["shannon_16s_site_compartment"]["spearman_rho"] > 0.50
    base1 = rep["tier1_spatial_block_cv"]["baseline_comparisons"]
    assert t1["log2_hydB_1l"]["r2"] > base1["log2_hydB_1l_compartment_only"]["r2"]
    assert t1["log1p_viral_richness"]["r2"] > base1["log1p_viral_richness_sequencing_depth_only"]["r2"]

    # Tier 2: Temporal Leave-One-Campaign-Out CV, Space-for-Time Transport, and Trip 5 Prospective Forecast
    t2 = rep["tier2_temporal_loco_cv"]["pooled_metrics"]
    assert t2["log2_hydB_1l"]["n"] == 550
    assert t2["log2_hydB_1l"]["spearman_rho"] > 0.50
    t2_prosp = rep["tier2_temporal_loco_cv"]["trip5_prospective_forecast_benchmark"]
    assert t2_prosp["n_test_visit_compartments"] >= 85
    assert (
        t2_prosp["level_forecasts_log2_hydB_1l"]["longitudinal_twin"]["r2"]
        > t2_prosp["level_forecasts_log2_hydB_1l"]["site_compartment_persistence"]["r2"]
    )
    assert t2_prosp["temporal_anomaly_forecasts_delta_trip5_minus_history"]["longitudinal_twin"]["pearson_r"] > 0.35

    # Tier 3: Causal-Ontological Polytope Projection & MinT Reconciliation (OCS 0.0 -> 1.0, CVM -> 0.0, Theorem 2 MSE reduction)
    t3 = rep["tier3_causal_ontological_projection"]
    assert t3["n_samples_evaluated"] == 550
    assert t3["ocs_unconstrained"] < 0.05
    assert t3["ocs_projected_two_layer"] == 1.0
    assert t3["ocs_nonnegative_mint"] == 1.0
    assert t3["mean_cvm_projected_two_layer"] == 0.0
    assert t3["rmse_projected_two_layer"] < t3["rmse_unconstrained"]
    assert t3["leaf_rmse_nonnegative_mint"] < t3["leaf_rmse_unconstrained"]
    assert t3["rmse_projected_two_layer"] < t3["reconciliation_baselines"]["bottom_up_nonneg_leaf_sum_rmse"]
    assert t3["theorem2_samplewise_error_reduction_rate"] == 1.0

    # Tier 4: Paired Natural-Experiment Counterfactual Calibration & Lineage-Turnover Decomposition
    t4 = rep["tier4_natural_experiment_calibration"]
    assert t4["root_vs_bulk_rhizosphere_paradox"]["within_empirical_95ci"] is True
    assert t4["surface_vs_subsurface_refugium"]["within_empirical_95ci"] is True
    assert t4["aridity_gradient_hydB_1l_vs_MAP"]["sign_and_significance_match"] is True
    assert t4["lineage_turnover_vs_copy_number_decomposition"]["root_vs_bulk_decomposition"]["frac_delta_from_turnover"] > 0.95

    # Tier 5: Synthetic Metagenome Benchmark, 72-MAG GSMM Pilot Audit & HiGHS FBA Thermodynamic Invariant Audit
    t5 = rep["tier5_synthetic_and_metabolic_benchmark"]
    assert t5["synthetic_metagenome_ge_benchmark"]["mean_pearson_r_across_guilds"] > 0.95
    assert t5["metabolic_pilot_72mag_10584run_audit"]["pilot_panel_mags"] == 72
    assert t5["metabolic_pilot_72mag_10584run_audit"]["total_fba_scenarios_run"] == 10584
    assert t5["metabolic_pilot_72mag_10584run_audit"]["uncurated_native_h2_gpr_wired_count"] == 1
    assert t5["highs_fba_thermodynamic_audit"]["closed_exchange_zero_growth_pass"] is True
    assert t5["highs_fba_thermodynamic_audit"]["sulfate_electron_donor_closure_pass"] is True
    assert t5["highs_fba_thermodynamic_audit"]["h2_uptake_monotonicity_25_steps_pass"] is True
    assert t5["highs_fba_thermodynamic_audit"]["salinity_atpm_monotonicity_25_steps_pass"] is True

    # Carbon Budget & Cross-Desert Validation
    tc = rep["carbon_budget_cross_desert_validation"]
    assert tc["all_cross_desert_benchmarks_pass"] is True
    assert all(tc["internal_invariants_audit"].values())


def test_carbon_budget_and_cross_desert_validation() -> None:
    suite = TwinValidationSuite()
    tc = suite.run_carbon_budget_cross_desert_validation()

    # 1. Empirical CBB / RuBisCO community shares across 155 observed WGS site-compartments
    assert 0.17 <= tc["overall_wgs_cbb_community_share_mean"] <= 0.28
    assert 0.17 <= tc["mag_hydB_1l_cbb_cocarriage_fraction"] <= 0.35
    assert tc["empirical_cbb_rubisco_by_compartment"]["ShallowSubsurface"]["share_cbb_mapped_mean"] > 0.24
    assert tc["empirical_cbb_rubisco_by_compartment"]["Surface"]["share_cbb_mapped_mean"] > 0.20

    # 2. Internal thermodynamic, stoichiometric, and physical invariants
    inv = tc["internal_invariants_audit"]
    assert inv["maintenance_infeasible_zero_rubisco_pass"] is True
    assert inv["zero_trace_gas_zero_rubisco_pass"] is True
    assert inv["plasmid_cbb_monotonicity_pass"] is True
    assert inv["shrub_cover_npp_and_sic_monotonicity_pass"] is True
    assert inv["compaction_diffusivity_h2_sink_monotonicity_pass"] is True
    assert inv["domain_of_support_abstention_preserved_pass"] is True

    # 3. Cross-desert validation against all 7 published desert benchmarks
    assert len(tc["cross_desert_benchmark_validation"]) == 7
    assert tc["all_cross_desert_benchmarks_pass"] is True

    # 4. Upscaled corridor (92,160 km^2) and basin (650,000 km^2) portfolio ordering
    ports = {
        p["policy_name"]: p
        for p in tc["corridor_and_basin_budgets"]["candidate_portfolios"]
    }
    p1 = ports["Baseline Hyperarid Do-Nothing"]
    p2 = ports["Integrated Eco-Restoration & Oryx Rewilding (Balanced Optimum)"]
    p4 = ports["High-Canopy Oasis Revegetation + Site 60 Ca-SO4 Irrigation"]
    p6 = ports["Unmitigated Heavy Seeding + Brine Disposal + Camel Overgrazing + Vibroseis"]

    assert (
        p4["total_direct_co2_bound"]["corridor_kt_yr_point"]
        > p2["total_direct_co2_bound"]["corridor_kt_yr_point"]
        > p1["total_direct_co2_bound"]["corridor_kt_yr_point"]
        > p6["total_direct_co2_bound"]["corridor_kt_yr_point"]
    )
    assert p6["direct_microbial_rubisco_co2"]["corridor_kt_yr_point"] == 0.0
    assert (
        p2["indirect_h2_sink_co2eq"]["corridor_kt_yr_point"]
        > p6["indirect_h2_sink_co2eq"]["corridor_kt_yr_point"]
    )
