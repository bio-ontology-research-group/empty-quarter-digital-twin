"""Export the verified EQ-BioDT state cube, 248-cell corridor atlas, six intervention curves,
six portfolios, four-pathway carbon budgets, 12 focal bacterial genera, Arabian Desert raster
grid summary, and 5-tier validation report to a compact JSON bundle for the interactive
web application at https://rubalkhali.science/twin (and https://rubalkhali.science/model).

Usage:
    uv run python reproducibility/export_web_twin_bundle.py [--output path/to/eq_biodt_bundle.json]
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from eq_biodt.interventions.cloud_seeding import CloudSeedingIntervention
from eq_biodt.interventions.desalination_brine import DesalinationBrineIntervention
from eq_biodt.interventions.megafauna_grazing import MegafaunaGrazingIntervention
from eq_biodt.interventions.pareto_optimizer import ParetoPortfolioOptimizer
from eq_biodt.interventions.refugium_thermal import RefugiumThermalIntervention
from eq_biodt.interventions.revegetation import RevegetationIntervention
from eq_biodt.interventions.seismic_hydrocarbon import SeismicHydrocarbonIntervention
from eq_biodt.models.carbon_budget_engine import CarbonBudgetEngine

PORTFOLIO_SHORT = {
    "Baseline Hyperarid Do-Nothing": ("P1", "P1: Current-use baseline"),
    "Integrated Eco-Restoration & Oryx Rewilding (Balanced Optimum)": (
        "P2",
        "P2: Integrated eco-restoration & Oryx",
    ),
    "Low-Water Chemosynthetic Crust & Refugium Conservation": (
        "P3",
        "P3: Low-water crust & refugium conservation",
    ),
    "High-Canopy Oasis Revegetation + Site 60 Ca-SO4 Irrigation": (
        "P4",
        "P4: High-canopy Ca-SO4 oasis",
    ),
    "Active Energy-Corridor Nodal Seismic + Oil-Spill Bioremediation": (
        "P5",
        "P5: Nodal seismic & bioremediation",
    ),
    "Unmitigated Heavy Seeding + Brine Disposal + Camel Overgrazing + Vibroseis": (
        "P6",
        "P6: Unmitigated disturbance",
    ),
}


def clean_json(obj: Any) -> Any:
    if isinstance(obj, dict):
        return {str(k): clean_json(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [clean_json(v) for v in obj]
    if isinstance(obj, (np.integer, int)):
        return int(obj)
    if isinstance(obj, (np.floating, float)):
        val = float(obj)
        return None if (math.isnan(val) or math.isinf(val)) else round(val, 6)
    if isinstance(obj, (np.bool_, bool)):
        return bool(obj)
    return obj


def export_web_bundle(repo_root: Path, out_path: Path) -> dict[str, Any]:
    sites_df = pd.read_parquet(repo_root / "data/processed/eq_sites_cube.parquet")
    comp_df = pd.read_parquet(repo_root / "data/processed/eq_compartment_state_cube.parquet")
    grid_df = pd.read_parquet(repo_root / "data/processed/eq_era5_h3_corridor_grid.parquet")
    val_report = json.loads((repo_root / "data/validation/validation_report.json").read_text())

    opt = ParetoPortfolioOptimizer()
    carbon_engine = CarbonBudgetEngine()

    comp_by_site: dict[int, dict[str, Any]] = {}
    for _, row in comp_df.iterrows():
        snum = int(row["site_num"])
        cname = str(row["compartment_norm"])
        h1l = float(row.get("hydB_1l_ge", 0.0))
        h1h = float(row.get("hydB_1h_ge", 0.0))
        h2a = float(row.get("hydB_2a_ge", 0.0))
        coxl = float(row.get("coxL_formI_ge", 0.0))
        tot_h2ase = h1l + h1h + h2a
        h2_pw = 29.8 * h1l + 8.5 * (h1h + h2a)
        co_pw = 1.95 * coxl
        s1l = max(1e-9, float(row.get("share_1l_mapped", 0.0)))
        s1l_cbb = float(row.get("share_1l_cbb_mapped", 0.0))
        s1l_3g = float(row.get("share_1l_cbb_3gene_mapped", 0.0))
        comp_by_site.setdefault(snum, {})[cname] = {
            "wgs_observed": bool(row.get("wgs_observed", False)),
            "ph_observed": bool(row.get("ph_observed", True)),
            "amplicon_16s_observed": bool(row.get("amplicon_16s_observed", True)),
            "wgs_sample_count": int(row.get("wgs_sample_count", 0)),
            "amplicon_16s_profile_count": int(row.get("amplicon_16s_profile_count", 0)),
            "shannon_16s": float(row.get("shannon_16s", 0.0)),
            "hydB_1l_ge": h1l,
            "hydB_1h_ge": h1h,
            "hydB_2a_ge": h2a,
            "coxL_formI_ge": coxl,
            "total_h2ase_ge": tot_h2ase,
            "total_trace_gas_power_uw_m3": h2_pw + co_pw,
            "h2_oxidation_power_uw_m3": h2_pw,
            "co_oxidation_power_uw_m3": co_pw,
            "actinomycetota_16s": float(row.get("actino_16s_rel_abund", 0.0)),
            "thermoplasmatota_16s": float(row.get("thermoplasmatota_16s_rel_abund", 0.0)),
            "halobacteriota_16s": float(row.get("halobacteriota_16s_rel_abund", 0.0)),
            "actinomycetota_mag": float(row.get("actinomycetota_rel_abund", 0.0)),
            "thermoplasmatota_mag": float(row.get("thermoplasmatota_rel_abund", 0.0)),
            "halobacteriota_mag": float(row.get("halobacteriota_rel_abund", 0.0)),
            "cbb_rbcl_share": float(row.get("share_cbb_mapped", 0.0)),
            "cbb_rbcl_prk_share": float(row.get("share_cbb_prk_mapped", 0.0)),
            "cbb_strict_3gene_share": float(row.get("share_cbb_3gene_mapped", 0.0)),
            "cbb_cond_1l_share": (
                float(np.clip(s1l_cbb / s1l, 0.0, 1.0))
                if float(row.get("share_1l_mapped", 0.0)) > 1e-6
                else 0.0
            ),
            "cbb_strict_3gene_cond_1l": (
                float(np.clip(s1l_3g / s1l, 0.0, 1.0))
                if float(row.get("share_1l_mapped", 0.0)) > 1e-6
                else 0.0
            ),
            "ph_cacl2": float(row.get("ph_mean", 7.8)),
            "dna_conc_ng_ul": float(row.get("dna_conc_mean", 0.02)),
            "xrf_evaporite_index": float(row.get("xrf_evaporite_index", 0.0)),
            "viral_votu_richness": float(row.get("viral_richness_mean", 0.0)),
            "plasmid_potu_richness": float(row.get("plasmid_richness_mean", 0.0)),
            "unmapped_read_pct": float(row.get("unmapped_read_pct", 0.0)),
            "ocs": 1.0,
        }

    sites_list = []
    for _, s in sites_df.iterrows():
        snum = int(s["site_num"]) if pd.notna(s.get("site_num")) else None
        sites_list.append(
            {
                "site_id": str(s.get("site_id", f"S{snum:02d}" if snum else "")),
                "site_num": snum,
                "site_label": str(s.get("site_label", f"Site {snum}")),
                "site_iri": str(s.get("site_iri", "")),
                "cohort": "primary" if bool(s.get("is_primary_transect", True)) else "auxiliary",
                "lat": float(s["lat"]),
                "lon": float(s["lon"]),
                "elevation_m": (
                    float(s.get("elevation_m", 0.0)) if pd.notna(s.get("elevation_m")) else 0.0
                ),
                "habitat_regime": str(s.get("habitat_regime", "aeolian_desert_soil")),
                "transect_zone": str(s.get("transect_zone", "central_hyperarid")),
                "map_mm_yr": float(s.get("map_mm", 20.0)),
                "era5_map_mm_yr": float(s.get("map_era5_mm", 16.0)),
                "mat_c": float(s.get("mat_c", 28.0)),
                "h3_5": str(s.get("h3_5", "")),
                "compartments": comp_by_site.get(snum, {}) if snum is not None else {},
            }
        )

    corridor_cells = []
    prescribed_counts = grid_df["prescribed_portfolio"].value_counts().to_dict()
    for idx, grow in grid_df.iterrows():
        p_code = str(grow.get("prescribed_portfolio", "ABSTAIN"))
        corridor_cells.append(
            {
                "cell_idx": int(idx),
                "era5_cell": str(grow.get("era5_cell", f"C{idx:03d}")),
                "h3_5": str(grow.get("h3_5", "")),
                "lat": float(grow["lat"]),
                "lon": float(grow["lon"]),
                "map_mm_yr": float(grow.get("map_mm", 20.0)),
                "mat_c": float(grow.get("mat_c", 28.0)),
                "elevation_m": (
                    float(grow.get("elevation_m", 150.0))
                    if pd.notna(grow.get("elevation_m"))
                    else 150.0
                ),
                "dist_km": float(grow.get("corridor_distance_km", 0.0)),
                "nearest_site_num": int(grow.get("nearest_primary_site", 30)),
                "transect_zone": str(grow.get("transect_zone", "central_hyperarid")),
                "in_support": bool(grow["in_empirical_support"]),
                "portfolio_id": p_code,
                "abstention_reason": str(grow.get("abstention_reason", "")),
            }
        )

    s30_row = (
        comp_df[(comp_df["site_num"] == 30) & (comp_df["compartment_norm"] == "Surface")]
        .iloc[0]
        .to_dict()
    )
    s60_row = (
        comp_df[(comp_df["site_num"] == 60) & (comp_df["compartment_norm"] == "Surface")]
        .iloc[0]
        .to_dict()
    )

    cs_engine = CloudSeedingIntervention()
    cloud_seeding_curve = []
    for rain_mm in [0.5, 2.0, 3.5, 6.0, 8.0, 10.0, 12.0, 15.0, 20.0, 25.0, 30.0]:
        sim = cs_engine.simulate(s30_row, rain_pulse_mm=float(rain_mm), horizon_days=14)
        s = sim["summary"]
        cloud_seeding_curve.append(
            {
                "rain_mm": rain_mm,
                "corridor_water_million_m3": round(rain_mm * 92.16, 1),
                "delta_shannon": round(float(s["delta_shannon_16s"]), 3),
                "shannon": round(float(s["counterfactual_peak_shannon_16s"]), 3),
                "min_h2_power_retention_pct": round(
                    100.0 * float(s["h2_power_retention_ratio"]), 1
                ),
                "min_h2_power_uW_m3": round(float(s["min_h2_power_uw_m3"]), 2),
                "viable_dna_fraction": round(float(s["viable_dna_fraction_peak"]), 3),
                "advisory": str(s["management_advisory"]),
            }
        )

    sim_6mm = cs_engine.simulate(s30_row, rain_pulse_mm=6.0, horizon_days=14)
    sim_15mm = cs_engine.simulate(s30_row, rain_pulse_mm=15.0, horizon_days=14)
    sim_25mm = cs_engine.simulate(s30_row, rain_pulse_mm=25.0, horizon_days=14)
    moisture_pulse_14d = []
    for d in range(len(sim_6mm["trajectory"])):
        t6 = sim_6mm["trajectory"][d]
        t15 = sim_15mm["trajectory"][d]
        t25 = sim_25mm["trajectory"][d]
        moisture_pulse_14d.append(
            {
                "day": int(t6["day"]),
                "vwc_6mm": round(float(t6["vwc_theta"]), 4),
                "vwc_15mm": round(float(t15["vwc_theta"]), 4),
                "vwc_25mm": round(float(t25["vwc_theta"]), 4),
                "rel_diff_6mm": round(float(t6["relative_gas_diffusivity"]), 3),
                "rel_diff_15mm": round(float(t15["relative_gas_diffusivity"]), 3),
                "rel_diff_25mm": round(float(t25["relative_gas_diffusivity"]), 3),
                "h2_power_6mm": round(float(t6["h2_oxidation_power_uw_m3"]), 2),
                "h2_power_15mm": round(float(t15["h2_oxidation_power_uw_m3"]), 2),
                "h2_power_25mm": round(float(t25["h2_oxidation_power_uw_m3"]), 2),
                "shannon_6mm": round(float(t6["shannon_16s"]), 3),
                "shannon_25mm": round(float(t25["shannon_16s"]), 3),
                "viable_dna_6mm": round(float(t6["viable_dna_fraction"]), 3),
                "viable_dna_25mm": round(float(t25["viable_dna_fraction"]), 3),
            }
        )

    rev_engine = RevegetationIntervention()
    revegetation_curve = []
    for cover_pct in [4, 10, 14, 20, 24, 30, 40, 50, 60]:
        f_cov = cover_pct / 100.0
        sim_inoc = rev_engine.simulate(
            s30_row,
            shrub_cover_fraction=f_cov,
            rak_pr_pgpr_dose=0.85,
            site60_rbcl_hyab_plasmid_dose=0.75,
        )["summary"]
        sim_uninoc = rev_engine.simulate(
            s30_row,
            shrub_cover_fraction=f_cov,
            rak_pr_pgpr_dose=0.0,
            site60_rbcl_hyab_plasmid_dose=0.0,
        )["summary"]
        cb_inoc = carbon_engine.compute_state_carbon_budget(
            s30_row,
            shrub_cover_fraction=f_cov,
            rak_pr_pgpr_dose=0.85,
            rbcl_hyab_plasmid_dose=0.75,
        )
        revegetation_curve.append(
            {
                "shrub_cover_pct": cover_pct,
                "area_weighted_1l_ge_inoc": round(float(sim_inoc["landscape_hydB_1l_ge"]), 3),
                "area_weighted_1l_ge_uninoc": round(float(sim_uninoc["landscape_hydB_1l_ge"]), 3),
                "biomass_weighted_1l_ge_inoc": round(
                    float(sim_inoc["biomass_weighted_hydB_1l_ge"]), 3
                ),
                "biomass_weighted_1l_ge_uninoc": round(
                    float(sim_uninoc["biomass_weighted_hydB_1l_ge"]), 3
                ),
                "abs_density_multiplier_inoc": round(
                    float(sim_inoc["absolute_hydB_1l_density_multiplier"]), 3
                ),
                "abs_density_multiplier_uninoc": round(
                    float(sim_uninoc["absolute_hydB_1l_density_multiplier"]), 3
                ),
                "gross_co2eq_t_km2_yr": round(
                    float(cb_inoc["total_co2eq_bound_t_km2_yr"]["point"]), 2
                ),
                "persistent_co2_t_km2_yr": round(
                    float(cb_inoc["net_persistent_co2_sequestered_t_km2_yr"]["point"]), 2
                ),
            }
        )

    brine_engine = DesalinationBrineIntervention()
    brine_grid = brine_engine.evaluate_factorial_brine_grid(s60_row)

    ref_engine = RefugiumThermalIntervention()
    refugium_curve = []
    for intact_frac in [1.0, 0.85, 0.70, 0.50, 0.35, 0.15, 0.0]:
        sim_r = ref_engine.simulate(
            s30_row,
            subsurface_refugium_intact_fraction=intact_frac,
            exogenous_choline_betaine_mm=0.5 if intact_frac >= 0.85 else 0.0,
        )["summary"]
        refugium_curve.append(
            {
                "intact_refugium_pct": round(intact_frac * 100, 0),
                "effective_temp_c": round(float(sim_r["effective_horizon_temp_c"]), 1),
                "bacillus_survival_index": round(float(sim_r["bacillus_60c_survival_index"]), 3),
                "archaeal_c_fixation_flux": round(
                    float(sim_r["archaeal_3hp4hb_c_fixation_flux"]), 3
                ),
                "shannon": round(float(sim_r["counterfactual_shannon_16s"]), 3),
            }
        )

    meg_engine = MegafaunaGrazingIntervention()
    megafauna_scenarios = []
    for label, oryx, gaz, camel, rot in [
        ("Current low-density baseline (P1)", 0.5, 1.0, 2.0, False),
        ("Native Oryx + Gazelle rewilding (P2)", 4.0, 7.0, 1.0, True),
        ("Crust & refugium conservation (P3)", 3.0, 5.0, 0.5, True),
        ("Severe camel overgrazing (P6: 25/100 km²)", 0.0, 0.0, 25.0, False),
    ]:
        sim_m = meg_engine.simulate(
            s30_row,
            oryx_head_per_100km2=oryx,
            gazelle_head_per_100km2=gaz,
            camel_head_per_100km2=camel,
            rotational_water_point_management=rot,
        )["summary"]
        megafauna_scenarios.append(
            {
                "label": label,
                "oryx_100km2": oryx,
                "gazelle_100km2": gaz,
                "camel_100km2": camel,
                "rotational": rot,
                "effective_porosity": round(float(sim_m["effective_soil_porosity"]), 3),
                "rel_gas_diffusivity": round(float(sim_m["relative_gas_diffusivity"]), 3),
                "nitrification_flux": round(float(sim_m["nitrification_flux"]), 3),
                "shannon": round(float(sim_m["counterfactual_shannon_16s"]), 3),
                "arg_burden": round(float(sim_m["resistome_arg_risk_index"]), 3),
            }
        )

    seis_engine = SeismicHydrocarbonIntervention()
    seismic_scenarios = []
    for label, line_d, grading, nodal, tph, alkb in [
        ("Low-impact wireless nodal (P2/P3: no spill)", 1.5, 0.0, True, 0.0, 0.2),
        (
            "Active energy corridor nodal + alkB/ladA bioremediation (P5)",
            8.0,
            0.04,
            True,
            15.0,
            0.95,
        ),
        ("Conventional cut-line grading + unmitigated spill (P6)", 12.0, 0.65, False, 25.0, 0.0),
    ]:
        sim_s = seis_engine.simulate(
            s30_row,
            seismic_line_density_km_km2=line_d,
            cut_line_grading_fraction=grading,
            low_impact_nodal_mode=nodal,
            hydrocarbon_tph_g_kg=tph,
            alkb_consortium_bioaug_dose=alkb,
        )["summary"]
        seismic_scenarios.append(
            {
                "label": label,
                "grading_pct": round(grading * 100, 0),
                "tph_g_kg": tph,
                "alkb_dose": alkb,
                "tph_degradation_pct": round(float(sim_s["tph_biodegradation_pct"]), 1),
                "residual_tph_g_kg": round(float(sim_s["residual_tph_g_kg"]), 2),
                "h2_power_uW_m3": round(float(sim_s["h2_oxidation_power_uw_m3"]), 2),
                "shannon": round(float(sim_s["counterfactual_shannon_16s"]), 3),
            }
        )

    pareto_res = opt.optimize_pareto_front(s30_row)
    portfolios_formatted = []
    for ev in pareto_res["all_evaluated_portfolios"]:
        pname = ev["policy"]["name"]
        code, short_lbl = PORTFOLIO_SHORT.get(pname, ("P?", pname))
        portfolios_formatted.append(
            {
                "code": code,
                "label": short_lbl,
                "full_name": pname,
                "is_pareto": bool(ev["is_pareto_optimal"]),
                "utility": round(float(ev["desdeo_utility_score"]), 4),
                "ocs": float(ev["ocs_score"]),
                "biosafety_pass": bool(ev["resistome_biosafety_pass"]),
                "objectives": ev["objectives"],
                "robustness": ev["probabilistic_robustness"],
                "policy": ev["policy"],
            }
        )
    portfolios_formatted.sort(key=lambda x: x["code"])

    carbon_val = val_report["carbon_budget_cross_desert_validation"]

    bundle = clean_json(
        {
            "version": "eq-biodt-0.2.0",
            "summary": {
                "n_sites_total": 70,
                "n_primary_sites": 60,
                "n_auxiliary_sites": 10,
                "n_compartment_states": 180,
                "n_observed_wgs_states": 155,
                "n_wgs_libraries": 551,
                "n_finite_guild_wgs": 550,
                "n_primary_wgs_cohort": 546,
                "n_16s_qc_primary_profiles": 1227,
                "n_16s_qc_total_profiles": 1237,
                "n_corridor_grid_cells": 248,
                "n_supported_corridor_cells": 128,
                "n_abstained_cells": 120,
                "supported_corridor_area_km2": 92160.0,
                "sand_sea_polygon_area_km2": 500000.0,
                "basin_total_area_km2": 650000.0,
                "prescribed_atlas_counts": {str(k): int(v) for k, v in prescribed_counts.items()},
            },
            "sites": sites_list,
            "corridor_cells": corridor_cells,
            "interventions": {
                "cloud_seeding_curve": cloud_seeding_curve,
                "moisture_pulse_14d": moisture_pulse_14d,
                "revegetation_curve": revegetation_curve,
                "brine_grid": brine_grid,
                "refugium_curve": refugium_curve,
                "megafauna_scenarios": megafauna_scenarios,
                "seismic_scenarios": seismic_scenarios,
            },
            "portfolios": portfolios_formatted,
            "carbon_budget": carbon_val,
            "validation_tiers": {
                "tier1_spatial_cv": val_report["tier1_spatial_block_cv"],
                "tier2_temporal": val_report["tier2_temporal_loco_cv"],
                "tier3_polytope": val_report["tier3_causal_ontological_projection"],
                "tier4_contrasts": val_report["tier4_natural_experiment_calibration"],
                "tier5_invariants": val_report["tier5_synthetic_and_metabolic_benchmark"],
            },
        }
    )

    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(bundle, indent=2) + "\n", encoding="utf-8")
    return bundle


def main() -> None:
    repo_root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        type=Path,
        default=repo_root / "reproducibility" / "outputs" / "eq_biodt_bundle.json",
    )
    args = parser.parse_args()
    export_web_bundle(repo_root, args.output)
    print(f"Wrote {args.output} ({args.output.stat().st_size / 1024:.1f} KB)")


if __name__ == "__main__":
    main()
