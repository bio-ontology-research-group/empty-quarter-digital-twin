"""Master end-to-end reproducibility pipeline for the Empty Quarter Biodiversity Digital Twin
(`EQ-BioDT`) manuscript.

Regenerates every validation metric, scenario trajectory, multi-objective Pareto portfolio,
dimensional carbon budget, 12-genus taxon response model, LaTeX table fragment, vector PDF
figure, and web twin bundle from the repository's self-contained data cubes, followed by an
automated 132-assertion numerical audit against the manuscript text and tables.

Usage:
    uv run python reproducibility/reproduce_all.py [--standalone]
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any, Dict

import numpy as np
import pandas as pd
from sklearn.linear_model import Ridge
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT / "src") not in sys.path:
    sys.path.insert(0, str(REPO_ROOT / "src"))
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from eq_biodt.api.main import get_runtime
from eq_biodt.data.build_twin_state_cube import TwinStateCubeBuilder
from eq_biodt.data.manifest_registry import CanonicalDataRegistry, CanonicalPaths
from eq_biodt.validation.twin_validator import TwinValidationSuite
from reproducibility.export_web_twin_bundle import export_web_bundle
from reproducibility.verify_all_paper_numbers import run_verification
from scripts.generate_figures import generate_fig1, generate_fig2, generate_fig3

PORTFOLIO_LABELS = {
    "Baseline Hyperarid Do-Nothing": ("P1", "Do-nothing baseline"),
    "Integrated Eco-Restoration & Oryx Rewilding (Balanced Optimum)": (
        "P2",
        "Integrated eco-restoration and oryx rewilding",
    ),
    "Low-Water Chemosynthetic Crust & Refugium Conservation": (
        "P3",
        "Low-water chemosynthetic crust and refugium conservation",
    ),
    "High-Canopy Oasis Revegetation + Site 60 Ca-SO4 Irrigation": (
        "P4",
        "High-canopy oasis revegetation with Ca-SO$_4$ irrigation",
    ),
    "Active Energy-Corridor Nodal Seismic + Oil-Spill Bioremediation": (
        "P5",
        "Energy corridor with nodal seismic and oil-spill bioremediation",
    ),
    "Unmitigated Heavy Seeding + Brine Disposal + Camel Overgrazing + Vibroseis": (
        "P6",
        "Unmitigated seeding, brine, camel overgrazing and Vibroseis",
    ),
}


def refit_12_genus_ridge_models(repo_root: Path) -> Dict[str, Any]:
    """Refit the 12 focal bacterial genera (7 Hallmark + 5 root-adjacent) site-mean CLR ridge
    regression models, 200-bootstrap distributions, and 5-block spatial CV R^2 from the bundled
    180-row site x compartment panel.
    """
    panel_path = repo_root / "reproducibility" / "data" / "taxon_12genera_site_clr_panel.tsv"
    frames = pd.read_csv(panel_path, sep="\t")
    ref_json = json.loads((repo_root / "data" / "processed" / "taxon_responses.json").read_text())
    taxa = [t["name"] for t in ref_json["taxa"]]
    plant = [t["name"] for t in ref_json["taxa"] if t["plant_associated"]]
    comp_map = {"surface": "Surface", "deep": "Deep", "rhizosphere": "Rhizosphere"}

    def _fit(x: np.ndarray, y: np.ndarray):
        model = make_pipeline(StandardScaler(), Ridge(alpha=10.0)).fit(x, y)
        scaler, ridge = model.steps[0][1], model.steps[1][1]
        beta = ridge.coef_ / scaler.scale_
        return ridge.intercept_ - beta @ scaler.mean_, beta

    compartments_out: Dict[str, Any] = {}
    for comp_key, comp_label in comp_map.items():
        frame = frames[frames["Type"] == comp_label].sort_values("longitude")
        x = frame[["MAT", "rain", "S"]].to_numpy(dtype=float)
        y = frame[taxa].to_numpy(dtype=float)
        intercept, beta = _fit(x, y)
        cv = np.empty_like(y)
        for test in np.array_split(np.arange(len(x)), 5):
            train = np.setdiff1d(np.arange(len(x)), test)
            a, b = _fit(x[train], y[train])
            cv[test] = x[test] @ b.T + a
        r2 = 1.0 - ((y - cv) ** 2).sum(axis=0) / ((y - y.mean(axis=0)) ** 2).sum(axis=0)
        rng = np.random.default_rng(20261004)
        draws = []
        for _ in range(200):
            sample = rng.integers(0, len(x), len(x))
            _, b = _fit(x[sample], y[sample])
            draws.append(b.tolist())
        plant_idx = [taxa.index(t) for t in plant]
        py, pcv = y[:, plant_idx].mean(axis=1), cv[:, plant_idx].mean(axis=1)
        plant_r2 = 1.0 - ((py - pcv) ** 2).sum() / ((py - py.mean()) ** 2).sum()
        compartments_out[comp_key] = {
            "n_sites": int(len(frame)),
            "intercept": intercept.tolist(),
            "beta": beta.tolist(),
            "beta_bootstrap": draws,
            "spatial_cv_r2": r2.tolist(),
            "plant_spatial_cv_r2": float(plant_r2),
            "training_min": x.min(axis=0).tolist(),
            "training_max": x.max(axis=0).tolist(),
        }
        # Verify exact match with reference taxon_responses.json
        assert np.allclose(
            r2, np.asarray(ref_json["compartments"][comp_key]["spatial_cv_r2"]), atol=1e-6
        )

    return {
        "taxa": ref_json["taxa"],
        "compartments": compartments_out,
    }


def write_table_s5_tex(repo_root: Path, opt_site30: Dict[str, Any]) -> Path:
    ranked = opt_site30["all_evaluated_portfolios"]
    rank_of = {ev["policy"]["name"]: i + 1 for i, ev in enumerate(ranked)}
    rows = []
    for ev in sorted(ranked, key=lambda e: PORTFOLIO_LABELS[e["policy"]["name"]][0]):
        name = ev["policy"]["name"]
        code, label = PORTFOLIO_LABELS[name]
        o = ev["objectives"]
        rob = ev["probabilistic_robustness"]
        rows.append(
            f"{code} & {label} & {'yes' if ev['is_pareto_optimal'] else 'no'} & "
            f"{o['J1_biodiversity_refugium_max']:.4f} & "
            f"{o['J2_chemosynthesis_cn_flux_max']:.4f} & "
            f"{o['J3_salinity_compaction_risk_min']:.4f} & "
            f"{o['J4_water_energy_arg_risk_min']:.4f} & "
            f"{ev['desdeo_utility_score']:.4f} & "
            f"{rob['pareto_dominance_prob']:.3f} & "
            f"{rob['rank1_prob']:.3f} & "
            f"{rank_of[name]} & {ev['ocs_score']:.1f} \\\\"
        )
    out = repo_root / "paper" / "tables" / "table_s3_scores.tex"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text("\n".join(rows) + "\n", encoding="utf-8")
    return out


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--standalone",
        action="store_true",
        help="Run using only repository-bundled files (no external workspace paths).",
    )
    args = parser.parse_args()

    out_dir = REPO_ROOT / "reproducibility" / "outputs"
    out_dir.mkdir(parents=True, exist_ok=True)

    if args.standalone:
        paths = CanonicalPaths(
            software_root=Path("/nonexistent_software_root"),
            papers_root=Path("/nonexistent_papers_root"),
            repo_root=REPO_ROOT,
        )
    else:
        paths = CanonicalPaths(repo_root=REPO_ROOT)

    print("[1/6] Loading state cubes & running Five-Tier Quantitative Validation Suite...")
    reg = CanonicalDataRegistry(paths)
    builder = TwinStateCubeBuilder(reg)
    suite = TwinValidationSuite(builder)
    val_report = suite.run_all_validations(output_dir=REPO_ROOT / "data" / "validation")
    print("      -> Wrote data/validation/validation_report.json & validation_metrics_summary.tsv")

    print("[2/6] Evaluating all 6 biophysical interventions & Pareto portfolios at S30/S35/S59/S60...")
    rt = get_runtime()
    row30 = rt.resolve_site_compartment_row("30", "Surface")
    row35 = rt.resolve_site_compartment_row("35", "Surface")
    row59 = rt.resolve_site_compartment_row("59", "Surface")
    row60 = rt.resolve_site_compartment_row("60", "Surface")

    opt30 = rt.pareto_opt.optimize_pareto_front(row30)
    opt35 = rt.pareto_opt.optimize_pareto_front(row35)
    comp_df = builder.build_compartment_state_cube()
    grid_df = builder.build_era5_h3_grid()
    atlas_summary = rt.pareto_opt.evaluate_corridor_atlas(comp_df, grid_df)
    factorial_brine = rt.desal_sim.evaluate_factorial_brine_grid(row60)

    scenarios_summary = {
        "pareto_site30": opt30,
        "pareto_site35": opt35,
        "corridor_atlas_summary": atlas_summary,
        "factorial_brine_site60_baseline": factorial_brine,
        "single_interventions": {
            "cloud_seeding_6mm_s30": rt.cloud_sim.simulate(row30, rain_pulse_mm=6.0)["summary"],
            "cloud_seeding_25mm_s30": rt.cloud_sim.simulate(row30, rain_pulse_mm=25.0)["summary"],
            "revegetation_uninoc_s30": rt.reveg_sim.simulate(
                row30,
                shrub_cover_fraction=0.30,
                rak_pr_pgpr_dose=0.0,
                site60_rbcl_hyab_plasmid_dose=0.0,
            )["summary"],
            "revegetation_coinoc_s30": rt.reveg_sim.simulate(
                row30,
                shrub_cover_fraction=0.30,
                rak_pr_pgpr_dose=0.85,
                site60_rbcl_hyab_plasmid_dose=0.75,
            )["summary"],
            "desalination_ca_so4_s60": rt.desal_sim.simulate(
                row60,
                discharge_salinity_ms_cm=45.0,
                ca_so4_to_nacl_ratio=2.4,
                sphingomonas_ect_kdp_plasmid_dose=0.8,
            )["summary"],
            "desalination_nacl_s59": rt.desal_sim.simulate(
                row59,
                discharge_salinity_ms_cm=180.0,
                ca_so4_to_nacl_ratio=0.15,
                sphingomonas_ect_kdp_plasmid_dose=0.0,
            )["summary"],
            "refugium_intact_95_s30": rt.refugium_sim.simulate(
                row30,
                subsurface_refugium_intact_fraction=0.95,
                surface_temp_c=58.0,
                exogenous_choline_betaine_mm=0.5,
            )["summary"],
            "refugium_graded_20_s30": rt.refugium_sim.simulate(
                row30,
                subsurface_refugium_intact_fraction=0.20,
                surface_temp_c=58.0,
                exogenous_choline_betaine_mm=0.0,
            )["summary"],
            "megafauna_rewilding_s30": rt.megafauna_sim.simulate(
                row30,
                oryx_head_per_100km2=4.0,
                gazelle_head_per_100km2=8.0,
                camel_head_per_100km2=0.0,
            )["summary"],
            "megafauna_camel_overgrazing_s30": rt.megafauna_sim.simulate(
                row30,
                oryx_head_per_100km2=0.0,
                gazelle_head_per_100km2=0.0,
                camel_head_per_100km2=25.0,
                rotational_water_point_management=False,
            )["summary"],
            "seismic_nodal_bioaug_s30": rt.seismic_sim.simulate(
                row30,
                seismic_line_density_km_km2=6.0,
                cut_line_grading_fraction=0.03,
                low_impact_nodal_mode=True,
                hydrocarbon_tph_g_kg=12.0,
                alkb_consortium_bioaug_dose=0.9,
            )["summary"],
            "seismic_vibroseis_unmitigated_s30": rt.seismic_sim.simulate(
                row30,
                seismic_line_density_km_km2=12.0,
                cut_line_grading_fraction=0.65,
                low_impact_nodal_mode=False,
                hydrocarbon_tph_g_kg=25.0,
                alkb_consortium_bioaug_dose=0.0,
            )["summary"],
        },
    }
    (REPO_ROOT / "data" / "validation" / "all_scenarios_summary.json").write_text(
        json.dumps(scenarios_summary, indent=2) + "\n", encoding="utf-8"
    )
    table_s5_path = write_table_s5_tex(REPO_ROOT, opt30)
    print(f"      -> Wrote {table_s5_path.relative_to(REPO_ROOT)}")

    print("[3/6] Refitting 12-genus (7 Hallmark + 5 root-adjacent) CLR ridge models & 5-block CV...")
    taxon_fit = refit_12_genus_ridge_models(REPO_ROOT)
    (out_dir / "taxon_12genera_refit_summary.json").write_text(
        json.dumps(taxon_fit, indent=2) + "\n", encoding="utf-8"
    )
    print("      -> Verified exact match with data/processed/taxon_responses.json")

    print("[4/6] Exporting interactive WebGL Digital Twin bundle...")
    bundle_path = out_dir / "eq_biodt_bundle.json"
    export_web_bundle(REPO_ROOT, bundle_path)
    print(f"      -> Wrote {bundle_path.relative_to(REPO_ROOT)}")

    print("[5/6] Regenerating publication vector PDF figures (Figs. 1-3)...")
    figures_dir = REPO_ROOT / "paper" / "figures"
    generate_fig1(REPO_ROOT, figures_dir)
    generate_fig2(REPO_ROOT, figures_dir)
    generate_fig3(REPO_ROOT, figures_dir)
    print(f"      -> Wrote Fig 1, Fig 2, Fig 3 to {figures_dir.relative_to(REPO_ROOT)}")

    print("[6/6] Running automated numerical assertion audit against manuscript claims...")
    audit = run_verification(REPO_ROOT, standalone=args.standalone)
    audit_path = out_dir / "paper_numbers_audit.json"
    audit_path.write_text(json.dumps(audit, indent=2) + "\n", encoding="utf-8")
    pd.DataFrame(audit["checks"]).to_csv(
        out_dir / "reproducibility_summary.tsv", sep="\t", index=False
    )
    print(
        f"SUCCESS: Reproduced and verified {audit['total_assertions_checked']} manuscript numbers "
        f"(all_passed={audit['all_passed']})."
    )


if __name__ == "__main__":
    main()
