"""Regenerate the intervention and Pareto values reported in main.tex (Site 30, 59, 60).

Run from the repository root: `uv run python paper/verify_manuscript_numbers.py`.
"""

import json

import pandas as pd

from eq_biodt.api.main import get_runtime
from eq_biodt.data.manifest_registry import CanonicalDataRegistry


def show(label, res):
    print(label, json.dumps(res["summary"], default=str))


def main():
    rt = get_runtime()
    row30 = rt.resolve_site_compartment_row("30", "Surface")
    keys = [
        "shannon_16s",
        "hydB_1l_ge",
        "plasmid_richness_mean",
        "viral_richness_mean",
        "thermoplasmatota_rel_abund",
        "xrf_evaporite_index",
        "map_mm",
    ]
    print("S30 baseline", {k: row30.get(k) for k in keys})

    show("cloud_6mm", rt.cloud_sim.simulate(row30, rain_pulse_mm=6.0))
    show("cloud_25mm", rt.cloud_sim.simulate(row30, rain_pulse_mm=25.0))

    show(
        "reveg_uninoculated",
        rt.reveg_sim.simulate(
            row30, shrub_cover_fraction=0.30, rak_pr_pgpr_dose=0.0, site60_rbcl_hyab_plasmid_dose=0.0
        ),
    )
    show(
        "reveg_coinoculated",
        rt.reveg_sim.simulate(
            row30, shrub_cover_fraction=0.30, rak_pr_pgpr_dose=0.85, site60_rbcl_hyab_plasmid_dose=0.75
        ),
    )

    row60 = rt.resolve_site_compartment_row("60", "Surface")
    row59 = rt.resolve_site_compartment_row("59", "Surface")
    print("baseline shannon S60", row60.get("shannon_16s"), "S59", row59.get("shannon_16s"))
    show(
        "desal_S60",
        rt.desal_sim.simulate(
            row60, discharge_salinity_ms_cm=45.0, ca_so4_to_nacl_ratio=2.4, sphingomonas_ect_kdp_plasmid_dose=0.8
        ),
    )
    show(
        "desal_S59",
        rt.desal_sim.simulate(
            row59, discharge_salinity_ms_cm=180.0, ca_so4_to_nacl_ratio=0.15, sphingomonas_ect_kdp_plasmid_dose=0.0
        ),
    )

    show(
        "refugium_95",
        rt.refugium_sim.simulate(
            row30, subsurface_refugium_intact_fraction=0.95, surface_temp_c=58.0, exogenous_choline_betaine_mm=0.5
        ),
    )
    show("refugium_20_default_betaine", rt.refugium_sim.simulate(row30, subsurface_refugium_intact_fraction=0.20, surface_temp_c=58.0))
    show(
        "refugium_20_no_betaine",
        rt.refugium_sim.simulate(
            row30, subsurface_refugium_intact_fraction=0.20, surface_temp_c=58.0, exogenous_choline_betaine_mm=0.0
        ),
    )

    show(
        "mega_wild_no_camel",
        rt.megafauna_sim.simulate(row30, oryx_head_per_100km2=4.0, gazelle_head_per_100km2=8.0, camel_head_per_100km2=0.0),
    )
    show("mega_wild_default_camel", rt.megafauna_sim.simulate(row30, oryx_head_per_100km2=4.0, gazelle_head_per_100km2=8.0))
    show(
        "mega_camel",
        rt.megafauna_sim.simulate(
            row30,
            oryx_head_per_100km2=0.0,
            gazelle_head_per_100km2=0.0,
            camel_head_per_100km2=25.0,
            rotational_water_point_management=False,
        ),
    )

    show(
        "seismic_nodal",
        rt.seismic_sim.simulate(
            row30,
            seismic_line_density_km_km2=6.0,
            cut_line_grading_fraction=0.03,
            low_impact_nodal_mode=True,
            hydrocarbon_tph_g_kg=12.0,
            alkb_consortium_bioaug_dose=0.9,
        ),
    )
    show(
        "seismic_vibroseis",
        rt.seismic_sim.simulate(
            row30,
            seismic_line_density_km_km2=12.0,
            cut_line_grading_fraction=0.65,
            low_impact_nodal_mode=False,
            hydrocarbon_tph_g_kg=25.0,
            alkb_consortium_bioaug_dose=0.0,
        ),
    )

    opt = rt.pareto_opt.optimize_pareto_front(row30)
    for ev in opt["all_evaluated_portfolios"]:
        print(
            "PARETO",
            ev["policy"]["name"],
            ev["is_pareto_optimal"],
            ev["desdeo_utility_score"],
            json.dumps(ev["objectives"]),
            ev["ocs_score"],
        )

    reg = CanonicalDataRegistry()
    print("xrf analytes mapped", len(reg.ontology.xrf_mappings))
    cube = pd.read_parquet("data/processed/eq_compartment_state_cube.parquet")
    print(
        "cube rows",
        len(cube),
        "ph_n notna",
        int(cube["ph_n"].notna().sum()),
        "sum ph_n",
        cube["ph_n"].sum(),
        "shannon notna",
        int(cube["shannon_16s"].notna().sum()),
    )
    print("compartments", cube["compartment_norm"].value_counts().to_dict())
    print("cube cols", list(cube.columns))
    sites = pd.read_parquet("data/processed/eq_sites_cube.parquet")
    print("sites", len(sites), "MAP", sites["map_mm"].min(), sites["map_mm"].max(), "MAT", sites["mat_c"].min(), sites["mat_c"].max())
    print("grid cells", len(pd.read_parquet("data/processed/eq_era5_h3_corridor_grid.parquet")))
    wgs = reg.load_wgs_sample_metadata()
    print("wgs", len(wgs), wgs.columns.tolist())
    if "trip" in wgs.columns:
        print("wgs per trip", wgs.groupby("trip").size().to_dict())
    led = reg.load_sample_ledger()
    print("ledger", len(led), led["campaign_role"].value_counts().to_dict())


if __name__ == "__main__":
    main()
