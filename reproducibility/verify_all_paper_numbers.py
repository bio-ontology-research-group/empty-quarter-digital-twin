"""Automated numerical audit verifying every number reported in `paper/main.tex` and
`paper/supplementary.tex` against the computed outputs of the EQ-BioDT engines and
validation suite.

Usage:
    uv run python reproducibility/verify_all_paper_numbers.py [--standalone]
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any, List, Tuple

import numpy as np
import pandas as pd

from eq_biodt.api.main import get_runtime
from eq_biodt.data.build_twin_state_cube import TwinStateCubeBuilder
from eq_biodt.data.manifest_registry import CanonicalDataRegistry, CanonicalPaths
from eq_biodt.validation.twin_validator import TwinValidationSuite


class PaperNumberVerifier:
    """Collects and verifies exact numerical assertions across all manuscript sections and tables."""

    def __init__(self) -> None:
        self.checks: List[Tuple[str, Any, Any, bool]] = []

    def check_exact(self, label: str, actual: Any, expected: Any) -> None:
        ok = actual == expected
        self.checks.append((label, actual, expected, ok))
        if not ok:
            raise AssertionError(f"[FAIL] {label}: actual={actual!r} != expected={expected!r}")

    def check_float(self, label: str, actual: float, expected: float, tol: float = 1e-4) -> None:
        ok = math.isfinite(float(actual)) and abs(float(actual) - float(expected)) <= tol
        self.checks.append((label, round(float(actual), 6), expected, ok))
        if not ok:
            raise AssertionError(
                f"[FAIL] {label}: actual={actual:.6f} != expected={expected:.6f} (tol={tol})"
            )


def run_verification(repo_root: Path, standalone: bool = False) -> dict[str, Any]:
    v = PaperNumberVerifier()

    if standalone:
        paths = CanonicalPaths(
            software_root=Path("/nonexistent_software_root"),
            papers_root=Path("/nonexistent_papers_root"),
            repo_root=repo_root,
        )
    else:
        paths = CanonicalPaths(repo_root=repo_root)

    reg = CanonicalDataRegistry(paths)
    builder = TwinStateCubeBuilder(reg)
    suite = TwinValidationSuite(builder)

    sites_df = builder.build_sites_cube()
    comp_df = builder.build_compartment_state_cube()
    visit_df = builder.build_visit_compartment_cube()
    grid_df = builder.build_era5_h3_grid()
    wgs_panel = suite.build_sample_level_validation_panel(id_floor=60)
    rep = suite.run_all_validations()

    # =========================================================================
    # 1. Sample Universe & State Cube Denominators (Table 1 & Section 2.1)
    # =========================================================================
    v.check_exact("Table 1: Total sites (primary + auxiliary)", len(sites_df), 70)
    v.check_exact(
        "Table 1: Primary transect sites (Sites 1-60)",
        int(sites_df["is_primary_transect"].sum()),
        60,
    )
    v.check_exact(
        "Table 1: Canonical auxiliary sites (RAK_1000061..1000070)",
        int((~sites_df["is_primary_transect"]).sum()),
        10,
    )
    v.check_exact(
        "Table 1: Confirmed Trip-1 numeric site aliases (61-64)",
        len(reg.ontology.site_aliases),
        4,
    )
    v.check_exact("Table 1: Site x compartment state cube rows", len(comp_df), 180)
    v.check_exact(
        "Table 1: Primary-site QC-passed 16S profiles",
        int(comp_df["amplicon_16s_profile_count"].sum()),
        1227,
    )
    v.check_exact(
        "Table 1: Surface 16S profiles",
        int(
            comp_df.loc[
                comp_df["compartment_norm"] == "Surface", "amplicon_16s_profile_count"
            ].sum()
        ),
        372,
    )
    v.check_exact(
        "Table 1: ShallowSubsurface 16S profiles",
        int(
            comp_df.loc[
                comp_df["compartment_norm"] == "ShallowSubsurface", "amplicon_16s_profile_count"
            ].sum()
        ),
        381,
    )
    v.check_exact(
        "Table 1: RootAdjacent 16S profiles",
        int(
            comp_df.loc[
                comp_df["compartment_norm"] == "RootAdjacent", "amplicon_16s_profile_count"
            ].sum()
        ),
        474,
    )
    v.check_exact(
        "Table 1: Observed WGS site-compartment states",
        int(comp_df["wgs_observed"].sum()),
        155,
    )
    v.check_exact(
        "Table 1: Observed WGS Surface states",
        int(comp_df.loc[comp_df["compartment_norm"] == "Surface", "wgs_observed"].sum()),
        45,
    )
    v.check_exact(
        "Table 1: Observed WGS ShallowSubsurface states",
        int(comp_df.loc[comp_df["compartment_norm"] == "ShallowSubsurface", "wgs_observed"].sum()),
        51,
    )
    v.check_exact(
        "Table 1: Observed WGS RootAdjacent states",
        int(comp_df.loc[comp_df["compartment_norm"] == "RootAdjacent", "wgs_observed"].sum()),
        59,
    )
    v.check_exact("Table 1: WGS soil samples with metadata", len(wgs_panel), 551)
    v.check_exact(
        "Table 1: Primary-site WGS soil samples (Sites 1-60)",
        int(wgs_panel["site"].between(1, 60).sum()),
        547,
    )
    v.check_exact(
        "Table 1: Primary-site finite guild GE WGS samples (Sites 1-60)",
        int((wgs_panel["site"].between(1, 60) & wgs_panel["guild_observed"]).sum()),
        546,
    )
    v.check_exact(
        "Table 1: Finite guild GE WGS samples",
        int(wgs_panel["guild_observed"].sum()),
        550,
    )
    v.check_exact(
        "Table 1: Observed mobilome + guild WGS samples",
        int((wgs_panel["mobilome_observed"] & wgs_panel["guild_observed"]).sum()),
        549,
    )
    trip_counts = wgs_panel.groupby("trip").size().to_dict()
    v.check_exact("Section 2.1: Trip 1 WGS samples", int(trip_counts["trip1"]), 210)
    v.check_exact("Section 2.1: Trip 2 WGS samples", int(trip_counts["trip2"]), 37)
    v.check_exact("Section 2.1: Trip 3 WGS samples", int(trip_counts["trip3"]), 128)
    v.check_exact("Section 2.1: Trip 4 WGS samples", int(trip_counts["trip4"]), 85)
    v.check_exact("Section 2.1: Trip 5 WGS samples", int(trip_counts["trip5"]), 91)
    v.check_exact("Section 2.1: Longitudinal visit-compartment rows", len(visit_df), 414)
    v.check_exact(
        "Table 1: Admitted CaCl2 pH measurements",
        int(comp_df["ph_n"].dropna().sum()),
        712,
    )
    v.check_exact(
        "Table 1: Site-compartments with admitted pH",
        int(comp_df["ph_observed"].sum()),
        180,
    )
    v.check_exact(
        "Table 1: Mapped XRF chemical analytes (including LE)",
        len(reg.ontology.xrf_mappings),
        93,
    )
    v.check_exact("Table 1: Corridor 0.25-deg grid cells", len(grid_df), 248)
    v.check_exact(
        "Table 1: Supported corridor grid cells (<=45 km)",
        int(grid_df["in_empirical_support"].sum()),
        128,
    )
    v.check_exact(
        "Table 1: Abstained corridor grid cells (>45 km)",
        int((~grid_df["in_empirical_support"]).sum()),
        120,
    )
    supp_counts = (
        grid_df[grid_df["in_empirical_support"]]["prescribed_portfolio"].value_counts().to_dict()
    )
    v.check_exact("Section 3.3: Prescribed P3 cells in corridor", int(supp_counts.get("P3", 0)), 122)
    v.check_exact("Section 3.3: Prescribed P4 cells in corridor", int(supp_counts.get("P4", 0)), 6)

    # =========================================================================
    # 2. Tier 1: 5-Fold Spatial Block CV (Table 2 & Section 3.1)
    # =========================================================================
    t1 = rep["tier1_spatial_block_cv"]["targets"]
    b1 = rep["tier1_spatial_block_cv"]["baseline_comparisons"]

    v.check_float("Tier 1 hydB_1l RMSE", t1["log2_hydB_1l"]["rmse"], 1.7689)
    v.check_float("Tier 1 hydB_1l MAE", t1["log2_hydB_1l"]["mae"], 1.4201)
    v.check_float("Tier 1 hydB_1l R2", t1["log2_hydB_1l"]["r2"], 0.2501)
    v.check_float("Tier 1 hydB_1l Pearson r", t1["log2_hydB_1l"]["pearson_r"], 0.5011)
    v.check_float("Tier 1 hydB_1l Spearman rho", t1["log2_hydB_1l"]["spearman_rho"], 0.5428)
    v.check_float(
        "Tier 1 hydB_1l rho CI low",
        t1["log2_hydB_1l"]["spearman_rho_cluster_ci95"][0],
        0.4697,
    )
    v.check_float(
        "Tier 1 hydB_1l rho CI high",
        t1["log2_hydB_1l"]["spearman_rho_cluster_ci95"][1],
        0.6221,
    )
    v.check_float(
        "Tier 1 hydB_1l R2 CI low", t1["log2_hydB_1l"]["r2_cluster_ci95"][0], 0.1483
    )
    v.check_float(
        "Tier 1 hydB_1l R2 CI high", t1["log2_hydB_1l"]["r2_cluster_ci95"][1], 0.3378
    )
    v.check_float(
        "Tier 1 Compartment-only R2",
        b1["log2_hydB_1l_compartment_only"]["r2"],
        0.1386,
    )
    v.check_float(
        "Tier 1 Compartment-only rho",
        b1["log2_hydB_1l_compartment_only"]["spearman_rho"],
        0.3469,
    )
    v.check_float(
        "Tier 1 Climate-only R2", b1["log2_hydB_1l_climate_only"]["r2"], 0.0583
    )
    v.check_float(
        "Tier 1 Climate-only rho",
        b1["log2_hydB_1l_climate_only"]["spearman_rho"],
        0.2431,
    )
    v.check_float("Tier 1 total_h2ase RMSE", t1["log2_total_h2ase_ge"]["rmse"], 1.4472)
    v.check_float("Tier 1 total_h2ase R2", t1["log2_total_h2ase_ge"]["r2"], 0.2188)
    v.check_float(
        "Tier 1 total_h2ase Spearman rho",
        t1["log2_total_h2ase_ge"]["spearman_rho"],
        0.5191,
    )
    v.check_float("Tier 1 viral_richness RMSE", t1["log1p_viral_richness"]["rmse"], 0.2103)
    v.check_float("Tier 1 viral_richness R2", t1["log1p_viral_richness"]["r2"], 0.3784)
    v.check_float(
        "Tier 1 viral_richness Spearman rho",
        t1["log1p_viral_richness"]["spearman_rho"],
        0.3388,
    )
    v.check_float(
        "Tier 1 Depth-only viral R2",
        b1["log1p_viral_richness_sequencing_depth_only"]["r2"],
        0.1410,
    )
    v.check_float(
        "Tier 1 plasmid_richness RMSE", t1["log1p_plasmid_richness"]["rmse"], 0.1827
    )
    v.check_float("Tier 1 plasmid_richness R2", t1["log1p_plasmid_richness"]["r2"], 0.3463)
    v.check_float(
        "Tier 1 plasmid_richness Spearman rho",
        t1["log1p_plasmid_richness"]["spearman_rho"],
        0.2996,
    )
    v.check_float(
        "Tier 1 shannon_16s RMSE",
        t1["shannon_16s_site_compartment"]["rmse"],
        0.8385,
    )
    v.check_float(
        "Tier 1 shannon_16s R2", t1["shannon_16s_site_compartment"]["r2"], 0.2453
    )
    v.check_float(
        "Tier 1 shannon_16s Spearman rho",
        t1["shannon_16s_site_compartment"]["spearman_rho"],
        0.5470,
    )

    # =========================================================================
    # 3. Tier 2: LOCO CV, Space-for-Time, and Trip 5 Forecast (Table 2 & Suppl. S1)
    # =========================================================================
    t2 = rep["tier2_temporal_loco_cv"]
    v.check_float(
        "Tier 2A LOCO hydB_1l RMSE", t2["pooled_metrics"]["log2_hydB_1l"]["rmse"], 1.6842
    )
    v.check_float(
        "Tier 2A LOCO hydB_1l R2", t2["pooled_metrics"]["log2_hydB_1l"]["r2"], 0.3202
    )
    v.check_float(
        "Tier 2A LOCO hydB_1l Pearson r",
        t2["pooled_metrics"]["log2_hydB_1l"]["pearson_r"],
        0.5660,
    )
    v.check_float(
        "Tier 2A LOCO hydB_1l Spearman rho",
        t2["pooled_metrics"]["log2_hydB_1l"]["spearman_rho"],
        0.5786,
    )
    v.check_float(
        "Tier 2A LOCO Trip 1 R2",
        t2["per_campaign_breakdown"]["trip1"]["log2_hydB_1l"]["r2"],
        0.3364,
    )
    v.check_float(
        "Tier 2A LOCO Trip 2 R2",
        t2["per_campaign_breakdown"]["trip2"]["log2_hydB_1l"]["r2"],
        0.2946,
    )
    v.check_float(
        "Tier 2A LOCO Trip 3 R2",
        t2["per_campaign_breakdown"]["trip3"]["log2_hydB_1l"]["r2"],
        0.2566,
    )
    v.check_float(
        "Tier 2A LOCO Trip 4 R2",
        t2["per_campaign_breakdown"]["trip4"]["log2_hydB_1l"]["r2"],
        0.3109,
    )
    v.check_float(
        "Tier 2A LOCO Trip 5 R2",
        t2["per_campaign_breakdown"]["trip5"]["log2_hydB_1l"]["r2"],
        0.3219,
    )

    t5_f = t2["trip5_prospective_forecast_benchmark"]
    v.check_exact("Tier 2C Train visit-compartments", t5_f["n_train_visit_compartments"], 325)
    v.check_exact("Tier 2C Test visit-compartments", t5_f["n_test_visit_compartments"], 89)
    v.check_exact(
        "Tier 2C Test with prior history",
        t5_f["n_test_with_prior_site_compartment_history"],
        85,
    )
    lvl = t5_f["level_forecasts_log2_hydB_1l"]
    v.check_float("Tier 2C Longitudinal twin level RMSE", lvl["longitudinal_twin"]["rmse"], 1.8363)
    v.check_float("Tier 2C Longitudinal twin level R2", lvl["longitudinal_twin"]["r2"], 0.2991)
    v.check_float(
        "Tier 2C Longitudinal twin level Pearson r",
        lvl["longitudinal_twin"]["pearson_r"],
        0.6170,
    )
    v.check_float(
        "Tier 2C Longitudinal twin level Spearman rho",
        lvl["longitudinal_twin"]["spearman_rho"],
        0.6188,
    )
    v.check_float(
        "Tier 2C Persistence level RMSE",
        lvl["site_compartment_persistence"]["rmse"],
        1.9753,
    )
    v.check_float(
        "Tier 2C Persistence level R2", lvl["site_compartment_persistence"]["r2"], 0.1890
    )
    v.check_float(
        "Tier 2C Env ridge level R2",
        lvl["environmental_ridge_spatial_only"]["r2"],
        0.2897,
    )
    v.check_float(
        "Tier 2C Random forest level R2", lvl["random_forest_baseline"]["r2"], 0.3244
    )
    v.check_float(
        "Tier 2C Gradient boosting level R2",
        lvl["gradient_boosting_baseline"]["r2"],
        0.3308,
    )
    anom = t5_f["temporal_anomaly_forecasts_delta_trip5_minus_history"]
    v.check_float("Tier 2C Longitudinal twin Delta R2", anom["longitudinal_twin"]["r2"], 0.0551)
    v.check_float(
        "Tier 2C Longitudinal twin Delta Pearson r",
        anom["longitudinal_twin"]["pearson_r"],
        0.4094,
    )
    v.check_float(
        "Tier 2C Longitudinal twin Delta Spearman rho",
        anom["longitudinal_twin"]["spearman_rho"],
        0.3753,
    )

    # =========================================================================
    # 4. Tier 3: OWL 2 Polytope Projection & MinT Reconciliation (Table 2)
    # =========================================================================
    t3 = rep["tier3_causal_ontological_projection"]
    v.check_exact("Tier 3 n_samples_evaluated", t3["n_samples_evaluated"], 550)
    v.check_float("Tier 3 OCS unconstrained", t3["ocs_unconstrained"], 0.0000)
    v.check_float("Tier 3 OCS projected Euclidean", t3["ocs_projected_two_layer"], 1.0000)
    v.check_float("Tier 3 OCS Nonnegative MinT", t3["ocs_nonnegative_mint"], 1.0000)
    v.check_float("Tier 3 Mean CVM unconstrained", t3["mean_cvm_unconstrained"], 0.2634)
    v.check_float("Tier 3 Mean CVM projected", t3["mean_cvm_projected_two_layer"], 0.0000)
    v.check_float("Tier 3 All-node RMSE unconstrained", t3["rmse_unconstrained"], 0.4708)
    v.check_float("Tier 3 All-node RMSE Euclidean", t3["rmse_projected_two_layer"], 0.4624)
    v.check_float("Tier 3 All-node RMSE Nonnegative MinT", t3["rmse_nonnegative_mint"], 0.4633)
    v.check_float("Tier 3 Leaf RMSE unconstrained", t3["leaf_rmse_unconstrained"], 0.3315)
    v.check_float("Tier 3 Leaf RMSE Euclidean", t3["leaf_rmse_projected_two_layer"], 0.3192)
    v.check_float("Tier 3 Leaf RMSE Nonnegative MinT", t3["leaf_rmse_nonnegative_mint"], 0.3161)
    v.check_float("Tier 3 Aggregate RMSE unconstrained", t3["aggregate_rmse_unconstrained"], 0.6672)
    v.check_float("Tier 3 Aggregate RMSE Euclidean", t3["aggregate_rmse_projected_two_layer"], 0.6617)
    v.check_float("Tier 3 MSE reduction pct", t3["mse_reduction_pct"], 3.51, tol=0.02)
    v.check_float(
        "Tier 3 Theorem 2 pass rate",
        t3["theorem2_samplewise_error_reduction_rate"],
        1.0000,
    )
    v.check_float(
        "Tier 3 Strict samplewise improvement rate",
        t3["strict_samplewise_error_reduction_rate"],
        1.0000,
    )

    # =========================================================================
    # 5. Tier 4 & Tier 5 Benchmarks (Table 2, Suppl. Tables S2-S3)
    # =========================================================================
    t4 = rep["tier4_natural_experiment_calibration"]
    rvb = t4["root_vs_bulk_rhizosphere_paradox"]
    v.check_exact("Tier 4 Root vs bulk paired visits", rvb["n_paired_visits"], 106)
    v.check_exact("Tier 4 Root vs bulk paired sites", rvb["n_paired_sites"], 47)
    v.check_float("Tier 4 Root vs bulk empirical log2FC", rvb["empirical_median_log2fc_hydB_1l"], -2.2747)
    v.check_float("Tier 4 Root vs bulk twin log2FC", rvb["twin_predicted_median_log2fc_hydB_1l"], -2.0969)
    v.check_float(
        "Tier 4 Root vs bulk lineage turnover fraction",
        t4["lineage_turnover_vs_copy_number_decomposition"]["root_vs_bulk_decomposition"][
            "frac_delta_from_turnover"
        ],
        0.9745,
    )
    svb = t4["surface_vs_subsurface_refugium"]
    v.check_exact("Tier 4 Surface vs subsurface visits", svb["n_paired_visits"], 68)
    v.check_exact("Tier 4 Surface vs subsurface sites", svb["n_paired_sites"], 36)
    v.check_float("Tier 4 Surface vs subsurface empirical log2FC", svb["empirical_median_log2fc_hydB_1l"], -1.0488)
    v.check_float("Tier 4 Surface vs subsurface twin log2FC", svb["twin_predicted_median_log2fc_hydB_1l"], -0.7562)
    arid = t4["aridity_gradient_hydB_1l_vs_MAP"]
    v.check_exact("Tier 4 Aridity gradient observed sites", arid["n_sites"], 51)
    v.check_float("Tier 4 Aridity gradient Spearman rho", arid["cube_spearman_rho"], -0.5137)

    t5 = rep["tier5_synthetic_and_metabolic_benchmark"]
    v.check_float(
        "Tier 5 Synthetic mean Pearson r",
        t5["synthetic_metagenome_ge_benchmark"]["mean_pearson_r_across_guilds"],
        0.9649,
    )
    v.check_float(
        "Tier 5 Synthetic mean Spearman rho",
        t5["synthetic_metagenome_ge_benchmark"]["mean_spearman_rho_across_guilds"],
        0.9580,
    )
    v.check_float(
        "Tier 5 Synthetic root vs bulk R2",
        t5["synthetic_metagenome_ge_benchmark"]["root_vs_bulk_synthetic_contrast_hydB_1l"]["r2"],
        0.9888,
    )
    v.check_exact(
        "Tier 5 72-MAG GSMM pilot count",
        t5["metabolic_pilot_72mag_10584run_audit"]["pilot_panel_mags"],
        72,
    )
    v.check_exact(
        "Tier 5 10,584 FBA runs",
        t5["metabolic_pilot_72mag_10584run_audit"]["total_fba_scenarios_run"],
        10584,
    )

    # =========================================================================
    # 6. Six Biophysical Interventions & Pareto Portfolios (Table 3 & Suppl. S4-S5)
    # =========================================================================
    rt = get_runtime()
    row30 = rt.resolve_site_compartment_row("30", "Surface")
    row59 = rt.resolve_site_compartment_row("59", "Surface")
    row60 = rt.resolve_site_compartment_row("60", "Surface")

    v.check_float("Table 3 S30 baseline Shannon", float(row30["shannon_16s"]), 6.1289)
    v.check_float("Table 3 S30 baseline hydB_1l GE", float(row30["hydB_1l_ge"]), 2.1059)
    v.check_float("Table 3 S60 baseline Shannon", float(row60["shannon_16s"]), 4.5249)
    v.check_float("Table 3 S59 baseline Shannon", float(row59["shannon_16s"]), 2.7346)

    c6 = rt.cloud_sim.simulate(row30, rain_pulse_mm=6.0)["summary"]
    c25 = rt.cloud_sim.simulate(row30, rain_pulse_mm=25.0)["summary"]
    v.check_float("Table 3 Cloud +6mm Delta Shannon", c6["delta_shannon_16s"], 0.3473)
    v.check_float("Table 3 Cloud +6mm peak Shannon", c6["counterfactual_peak_shannon_16s"], 6.4762)
    v.check_float("Table 3 Cloud +6mm min H2 power", c6["min_h2_power_uw_m3"], 40.1066)
    v.check_float("Table 3 Cloud +25mm Delta Shannon", c25["delta_shannon_16s"], 0.6972)
    v.check_float("Table 3 Cloud +25mm peak Shannon", c25["counterfactual_peak_shannon_16s"], 6.8261)
    v.check_float("Table 3 Cloud +25mm min H2 power", c25["min_h2_power_uw_m3"], 10.4208)

    rv_un = rt.reveg_sim.simulate(
        row30, shrub_cover_fraction=0.30, rak_pr_pgpr_dose=0.0, site60_rbcl_hyab_plasmid_dose=0.0
    )["summary"]
    rv_co = rt.reveg_sim.simulate(
        row30, shrub_cover_fraction=0.30, rak_pr_pgpr_dose=0.85, site60_rbcl_hyab_plasmid_dose=0.75
    )["summary"]
    v.check_float("Table 3 Reveg uninoc DNA mult", rv_un["biomass_dna_multiplier"], 3.4660)
    v.check_float("Table 3 Reveg uninoc landscape hydB_1l", rv_un["landscape_hydB_1l_ge"], 1.6300)
    v.check_float("Table 3 Reveg uninoc biomass-weighted hydB_1l", rv_un["biomass_weighted_hydB_1l_ge"], 0.8399)
    v.check_float("Table 3 Reveg uninoc abs density mult", rv_un["absolute_hydB_1l_density_multiplier"], 1.3822)
    v.check_float("Table 3 Reveg coinoc DNA mult", rv_co["biomass_dna_multiplier"], 4.2900)
    v.check_float("Table 3 Reveg coinoc landscape hydB_1l", rv_co["landscape_hydB_1l_ge"], 1.6087)
    v.check_float("Table 3 Reveg coinoc biomass-weighted hydB_1l", rv_co["biomass_weighted_hydB_1l_ge"], 0.7189)
    v.check_float("Table 3 Reveg coinoc abs density mult", rv_co["absolute_hydB_1l_density_multiplier"], 1.4644)
    v.check_float("Table 3 Reveg coinoc C fixation", rv_co["autotrophic_carbon_fixation_flux"], 1.35875)

    d60 = rt.desal_sim.simulate(
        row60, discharge_salinity_ms_cm=45.0, ca_so4_to_nacl_ratio=2.4, sphingomonas_ect_kdp_plasmid_dose=0.8
    )["summary"]
    d59 = rt.desal_sim.simulate(
        row59, discharge_salinity_ms_cm=180.0, ca_so4_to_nacl_ratio=0.15, sphingomonas_ect_kdp_plasmid_dose=0.0
    )["summary"]
    v.check_float("Table 3 Desal S60 Shannon", d60["counterfactual_shannon_16s"], 4.5883)
    v.check_float("Table 3 Desal S60 C fixation", d60["autotrophic_carbon_fixation_flux"], 1.4724)
    v.check_float("Table 3 Desal S60 sulfur flux", d60["sulfur_redox_flux"], 0.8291)
    v.check_float("Table 3 Desal S59 Shannon", d59["counterfactual_shannon_16s"], 1.3500)
    v.check_float("Table 3 Desal S59 C fixation", d59["autotrophic_carbon_fixation_flux"], 0.0000)
    v.check_float("Table 3 Desal S59 ATPM deficit", d59["atpm_maintenance_deficit"], 1.5837)

    rf95 = rt.refugium_sim.simulate(
        row30, subsurface_refugium_intact_fraction=0.95, surface_temp_c=58.0, exogenous_choline_betaine_mm=0.5
    )["summary"]
    rf20 = rt.refugium_sim.simulate(
        row30, subsurface_refugium_intact_fraction=0.20, surface_temp_c=58.0, exogenous_choline_betaine_mm=0.0
    )["summary"]
    v.check_float("Table 3 Refugium 95% Shannon", rf95["counterfactual_shannon_16s"], 6.1564)
    v.check_float("Table 3 Refugium 95% Bacillus survival", rf95["bacillus_60c_survival_index"], 0.5759)
    v.check_float("Table 3 Refugium 20% Shannon", rf20["counterfactual_shannon_16s"], 5.4439)
    v.check_float("Table 3 Refugium 20% Bacillus survival", rf20["bacillus_60c_survival_index"], 0.4122)

    mg_wild = rt.megafauna_sim.simulate(
        row30, oryx_head_per_100km2=4.0, gazelle_head_per_100km2=8.0, camel_head_per_100km2=0.0
    )["summary"]
    mg_camel = rt.megafauna_sim.simulate(
        row30,
        oryx_head_per_100km2=0.0,
        gazelle_head_per_100km2=0.0,
        camel_head_per_100km2=25.0,
        rotational_water_point_management=False,
    )["summary"]
    v.check_float("Table 3 Megafauna wild Shannon", mg_wild["counterfactual_shannon_16s"], 6.3159)
    v.check_float("Table 3 Megafauna wild porosity", mg_wild["effective_soil_porosity"], 0.3949)
    v.check_float("Table 3 Megafauna wild ARG risk", mg_wild["resistome_arg_risk_index"], 0.0753)
    v.check_float("Table 3 Megafauna camel Shannon", mg_camel["counterfactual_shannon_16s"], 4.3089)
    v.check_float("Table 3 Megafauna camel porosity", mg_camel["effective_soil_porosity"], 0.2200)
    v.check_float("Table 3 Megafauna camel ARG risk", mg_camel["resistome_arg_risk_index"], 0.6650)

    ss_nod = rt.seismic_sim.simulate(
        row30,
        seismic_line_density_km_km2=6.0,
        cut_line_grading_fraction=0.03,
        low_impact_nodal_mode=True,
        hydrocarbon_tph_g_kg=12.0,
        alkb_consortium_bioaug_dose=0.9,
    )["summary"]
    ss_vib = rt.seismic_sim.simulate(
        row30,
        seismic_line_density_km_km2=12.0,
        cut_line_grading_fraction=0.65,
        low_impact_nodal_mode=False,
        hydrocarbon_tph_g_kg=25.0,
        alkb_consortium_bioaug_dose=0.0,
    )["summary"]
    v.check_float("Table 3 Seismic nodal Shannon", ss_nod["counterfactual_shannon_16s"], 6.1045)
    v.check_float("Table 3 Seismic nodal H2 power", ss_nod["h2_oxidation_power_uw_m3"], 44.9109)
    v.check_float("Table 3 Seismic nodal TPH deg pct", ss_nod["tph_biodegradation_pct"], 76.97)
    v.check_float("Table 3 Seismic vibroseis Shannon", ss_vib["counterfactual_shannon_16s"], 4.3934)
    v.check_float("Table 3 Seismic vibroseis H2 power", ss_vib["h2_oxidation_power_uw_m3"], 7.1005)
    v.check_float("Table 3 Seismic vibroseis TPH deg pct", ss_vib["tph_biodegradation_pct"], 20.15)

    # Pareto portfolios at S30 (Suppl. Table S5)
    opt = rt.pareto_opt.optimize_pareto_front(row30)
    by_name = {ev["policy"]["name"]: ev for ev in opt["all_evaluated_portfolios"]}
    p3 = by_name["Low-Water Chemosynthetic Crust & Refugium Conservation"]
    p2 = by_name["Integrated Eco-Restoration & Oryx Rewilding (Balanced Optimum)"]
    p4 = by_name["High-Canopy Oasis Revegetation + Site 60 Ca-SO4 Irrigation"]
    p1 = by_name["Baseline Hyperarid Do-Nothing"]
    p5 = by_name["Active Energy-Corridor Nodal Seismic + Oil-Spill Bioremediation"]
    p6 = by_name["Unmitigated Heavy Seeding + Brine Disposal + Camel Overgrazing + Vibroseis"]

    v.check_float("Suppl. Table S5 P3 utility", p3["desdeo_utility_score"], 0.9481)
    v.check_float("Suppl. Table S5 P3 dominance prob", p3["probabilistic_robustness"]["pareto_dominance_prob"], 1.000)
    v.check_float("Suppl. Table S5 P3 rank1 prob", p3["probabilistic_robustness"]["rank1_prob"], 0.310)
    v.check_float("Suppl. Table S5 P2 utility", p2["desdeo_utility_score"], 0.9432)
    v.check_float("Suppl. Table S5 P2 rank1 prob", p2["probabilistic_robustness"]["rank1_prob"], 0.245)
    v.check_float("Suppl. Table S5 P4 utility", p4["desdeo_utility_score"], 0.9396)
    v.check_float("Suppl. Table S5 P4 rank1 prob", p4["probabilistic_robustness"]["rank1_prob"], 0.305)
    v.check_float("Suppl. Table S5 P1 utility", p1["desdeo_utility_score"], 0.8367)
    v.check_float("Suppl. Table S5 P5 utility", p5["desdeo_utility_score"], 0.9183)
    v.check_float("Suppl. Table S5 P6 utility", p6["desdeo_utility_score"], -0.0525)

    # =========================================================================
    # 7. Four-Pathway Carbon Budget & Cross-Desert Validation (Table 4 & Suppl. S6)
    # =========================================================================
    tc = rep["carbon_budget_cross_desert_validation"]
    v.check_exact("Section 3.4 All 7 cross-desert benchmarks pass", tc["all_cross_desert_benchmarks_pass"], True)
    v.check_float("Section 3.4 Overall WGS CBB share", tc["overall_wgs_cbb_community_share_mean"], 0.2293)
    v.check_float("Section 3.4 Overall WGS CBB+prkB share", tc["overall_wgs_cbb_prkb_community_share_mean"], 0.1321)
    v.check_float("Section 3.4 Overall WGS 3-gene CBB share", tc["overall_wgs_cbb_3gene_community_share_mean"], 0.1123)
    v.check_float("Section 3.4 Wetted active dark CO2 fixation", tc["dark_co2_fixation_wetted_active_horizon_nmol_gdw_d"], 0.04824)
    v.check_float("Section 3.4 Annual dry column dark CO2 fixation", tc["dark_co2_fixation_annual_dry_column_nmol_gdw_d"], 0.000397)

    cb = tc["corridor_and_basin_budgets"]
    base_cb = cb["baseline_observed"]
    atlas_cb = cb["prescribed_corridor_atlas"]
    ports_cb = {p["policy_name"]: p for p in cb["candidate_portfolios"]}

    v.check_float("Table 4 Baseline direct CO2 corridor kt/yr", base_cb["total_direct_co2_bound"]["corridor_kt_yr_point"], 801.60, tol=0.01)
    v.check_float("Table 4 Baseline persistent CO2 corridor kt/yr", base_cb["net_persistent_co2_sequestered"]["corridor_kt_yr_point"], 140.88, tol=0.01)
    v.check_float("Table 4 Baseline total CO2eq corridor kt/yr", base_cb["total_co2eq_bound"]["corridor_kt_yr_point"], 811.30, tol=0.01)
    v.check_float("Table 4 Baseline total CO2eq basin kt/yr", base_cb["total_co2eq_bound"]["basin_kt_yr_point"], 1496.05, tol=0.01)

    v.check_float("Table 4 Prescribed Atlas direct CO2 corridor kt/yr", atlas_cb["total_direct_co2_bound"]["corridor_kt_yr_point"], 5192.04, tol=0.01)
    v.check_float("Table 4 Prescribed Atlas persistent CO2 corridor kt/yr", atlas_cb["net_persistent_co2_sequestered"]["corridor_kt_yr_point"], 814.52, tol=0.01)
    v.check_float("Table 4 Prescribed Atlas total CO2eq corridor kt/yr", atlas_cb["total_co2eq_bound"]["corridor_kt_yr_point"], 5203.77, tol=0.01)
    v.check_float("Table 4 Prescribed Atlas total CO2eq basin kt/yr", atlas_cb["total_co2eq_bound"]["basin_kt_yr_point"], 5888.53, tol=0.01)

    p1_cb = ports_cb["Baseline Hyperarid Do-Nothing"]
    p2_cb = ports_cb["Integrated Eco-Restoration & Oryx Rewilding (Balanced Optimum)"]
    p3_cb = ports_cb["Low-Water Chemosynthetic Crust & Refugium Conservation"]
    p4_cb = ports_cb["High-Canopy Oasis Revegetation + Site 60 Ca-SO4 Irrigation"]
    p5_cb = ports_cb["Active Energy-Corridor Nodal Seismic + Oil-Spill Bioremediation"]
    p6_cb = ports_cb["Unmitigated Heavy Seeding + Brine Disposal + Camel Overgrazing + Vibroseis"]

    v.check_float("Table 4 P1 total CO2eq corridor kt/yr", p1_cb["total_co2eq_bound"]["corridor_kt_yr_point"], 1118.669, tol=0.01)
    v.check_float("Table 4 P1 persistent CO2 corridor kt/yr", p1_cb["net_persistent_co2_sequestered"]["corridor_kt_yr_point"], 188.140, tol=0.01)
    v.check_float("Table 4 P3 total CO2eq corridor kt/yr", p3_cb["total_co2eq_bound"]["corridor_kt_yr_point"], 4597.953, tol=0.01)
    v.check_float("Table 4 P3 persistent CO2 corridor kt/yr", p3_cb["net_persistent_co2_sequestered"]["corridor_kt_yr_point"], 721.491, tol=0.01)
    v.check_float("Table 4 P2 total CO2eq corridor kt/yr", p2_cb["total_co2eq_bound"]["corridor_kt_yr_point"], 8421.650, tol=0.01)
    v.check_float("Table 4 P2 persistent CO2 corridor kt/yr", p2_cb["net_persistent_co2_sequestered"]["corridor_kt_yr_point"], 1307.091, tol=0.01)
    v.check_float("Table 4 P4 total CO2eq corridor kt/yr", p4_cb["total_co2eq_bound"]["corridor_kt_yr_point"], 14713.164, tol=0.01)
    v.check_float("Table 4 P4 persistent CO2 corridor kt/yr", p4_cb["net_persistent_co2_sequestered"]["corridor_kt_yr_point"], 2273.738, tol=0.01)
    v.check_float("Table 4 P5 total CO2eq corridor kt/yr", p5_cb["total_co2eq_bound"]["corridor_kt_yr_point"], 6009.729, tol=0.01)
    v.check_float("Table 4 P5 persistent CO2 corridor kt/yr", p5_cb["net_persistent_co2_sequestered"]["corridor_kt_yr_point"], 937.107, tol=0.01)
    v.check_float("Table 4 P6 total CO2eq corridor kt/yr", p6_cb["total_co2eq_bound"]["corridor_kt_yr_point"], 131.899, tol=0.01)
    v.check_float("Table 4 P6 persistent CO2 corridor kt/yr", p6_cb["net_persistent_co2_sequestered"]["corridor_kt_yr_point"], 40.955, tol=0.01)
    v.check_float("Table 4 P6 RuBisCO CO2 corridor kt/yr", p6_cb["direct_microbial_rubisco_co2"]["corridor_kt_yr_point"], 0.00, tol=0.001)

    # =========================================================================
    # 8. 12 Focal Bacterial Genera & Arabian Desert WebGL Grid (Section 3.5)
    # =========================================================================
    taxon_json = json.loads((repo_root / "data/processed/taxon_responses.json").read_text())
    grid_summary = json.loads((repo_root / "data/processed/arabian_desert_grid_summary.json").read_text())

    v.check_exact("Section 3.5 Total focal bacterial genera", len(taxon_json["taxa"]), 12)
    v.check_exact(
        "Section 3.5 Hallmark desert genera count",
        sum(1 for t in taxon_json["taxa"] if t["hallmark"]),
        7,
    )
    v.check_exact(
        "Section 3.5 Root-adjacent indicator genera count",
        sum(1 for t in taxon_json["taxa"] if t["plant_associated"]),
        5,
    )
    csp_idx = [t["name"] for t in taxon_json["taxa"]].index("CSP1-2")
    flav_idx = [t["name"] for t in taxon_json["taxa"]].index("Flavisolibacter")
    v.check_float(
        "Section 3.5 CSP1-2 surface 5-block spatial CV R2",
        taxon_json["compartments"]["surface"]["spatial_cv_r2"][csp_idx],
        0.4436,
        tol=1e-3,
    )
    v.check_float(
        "Section 3.5 CSP1-2 deep 5-block spatial CV R2",
        taxon_json["compartments"]["deep"]["spatial_cv_r2"][csp_idx],
        0.4517,
        tol=1e-3,
    )
    v.check_float(
        "Section 3.5 Flavisolibacter deep 5-block spatial CV R2",
        taxon_json["compartments"]["deep"]["spatial_cv_r2"][flav_idx],
        0.4763,
        tol=1e-3,
    )
    v.check_exact(
        "Section 3.5 Medium grid total pixels (0.05 deg, 360x480)",
        grid_summary["resolutions"]["medium"]["total_pixels"],
        172800,
    )
    v.check_exact(
        "Section 3.5 Medium grid Arabian Desert mask pixels",
        grid_summary["resolutions"]["medium"]["in_arabian_desert_mask_pixels"],
        63902,
    )
    v.check_exact(
        "Section 3.5 Medium grid southern transect envelope pixels",
        grid_summary["resolutions"]["medium"]["in_southern_transect_envelope_pixels"],
        2916,
    )
    v.check_exact(
        "Section 3.5 Fine grid total pixels (0.04 deg, 450x600)",
        grid_summary["resolutions"]["fine"]["total_pixels"],
        270000,
    )
    v.check_exact(
        "Section 3.5 Fine grid Arabian Desert mask pixels",
        grid_summary["resolutions"]["fine"]["in_arabian_desert_mask_pixels"],
        99870,
    )
    v.check_exact(
        "Section 3.5 Fine grid southern transect envelope pixels",
        grid_summary["resolutions"]["fine"]["in_southern_transect_envelope_pixels"],
        4572,
    )

    return {
        "total_assertions_checked": len(v.checks),
        "all_passed": all(c[3] for c in v.checks),
        "standalone_mode": standalone,
        "checks": [
            {"label": c[0], "actual": c[1], "expected": c[2], "passed": c[3]} for c in v.checks
        ],
    }


def main() -> None:
    repo_root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--standalone",
        action="store_true",
        help="Force standalone verification using only repository-bundled files.",
    )
    args = parser.parse_args()
    res = run_verification(repo_root, standalone=args.standalone)
    out_file = repo_root / "reproducibility" / "outputs" / "paper_numbers_audit.json"
    out_file.parent.mkdir(parents=True, exist_ok=True)
    out_file.write_text(json.dumps(res, indent=2) + "\n", encoding="utf-8")
    print(
        f"VERIFIED {res['total_assertions_checked']} MANUSCRIPT ASSERTIONS "
        f"(standalone={res['standalone_mode']}): ALL PASSED -> {out_file}"
    )


if __name__ == "__main__":
    main()
