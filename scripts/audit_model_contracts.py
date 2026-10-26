"""Reproduce scientific-contract checks independently of the numerical fixtures.

Run: uv run python scripts/audit_model_contracts.py
Failures are reported as findings, without changing model parameters or outputs.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

import numpy as np
import pandas as pd

from eq_biodt.models.causal_ontological_twin import CausalOntologicalTwinEngine, OntologicalPolytope
from eq_biodt.models.carbon_budget_engine import CarbonBudgetEngine
from eq_biodt.models.metabolic_community_engine import CommunityMetabolicInputs, MetabolicCommunityEngine
from eq_biodt.models.spatial_jsdm_virome import SpatialJSDMViromeModel
from eq_biodt.interventions.pareto_optimizer import ParetoPortfolioOptimizer


def audit():
    root = Path(__file__).resolve().parents[1]
    findings = []

    def record(name, passed, evidence):
        findings.append(dict(check=name, status="pass" if passed else "fail", evidence=evidence))

    engine = CausalOntologicalTwinEngine()
    poly = engine.polytope
    record("OWL module dimensions", (len(poly.variable_names), len(poly.subsumptions), len(poly.partitions)) == (17, 13, 5),
           dict(variables=len(poly.variable_names), subsumptions=len(poly.subsumptions), partitions=len(poly.partitions)))
    for label, call in [
        ("Projection rejects missing coordinates", lambda: poly.project({"hydB_1l_ge": 1.0})),
        ("Constitutive evaluation rejects missing leaf coordinates", lambda: engine.evaluate_constitutive_layer({"hydB_1l_ge": 1.0})),
        ("Projection rejects indefinite metric", lambda: poly.project(dict.fromkeys(poly.variable_names, 0.0), -np.eye(17))),
    ]:
        try:
            value = call()
            record(label, False, {"accepted": True, "output_coordinates": len(value)})
        except ValueError as exc:
            record(label, True, str(exc))

    with patch("eq_biodt.models.metabolic_community_engine.linprog", return_value=SimpleNamespace(success=False, status=4, x=None)):
        result = MetabolicCommunityEngine().solve_community_fba(CommunityMetabolicInputs())
        record("Numerical LP failure is distinct from biological dormancy", result["metabolic_viability_state"] != "MAINTENANCE_DEFICIT_DORMANCY", result)

    result = MetabolicCommunityEngine().solve_community_fba(CommunityMetabolicInputs(
        h2_uptake_bound=0, co_uptake_bound=0, root_exudate_doc_bound=.1,
        exogenous_nh4_no3_bound=0, sulfate_so4_bound=0, salinity_ms_cm=216.2))
    record("Dormancy fallback no longer creates ATP from unconsumed DOC", result["atpm_maintenance_drain"] == 0, result)

    result = SpatialJSDMViromeModel().predict_virome_mobilome_state(4.0, "Surface", 25, 0,
        baseline_viral_richness=0, baseline_plasmid_richness=0)
    record("Zero-intervention virome identity includes zero baselines", result["viral_richness_votu"] == 0 and result["plasmid_richness_potu"] == 0, result)

    frame = pd.read_parquet(root / "data/processed/eq_compartment_state_cube.parquet")
    grid = pd.read_parquet(root / "data/processed/eq_era5_h3_corridor_grid.parquet")
    row = frame[(frame.site_num == 30) & (frame.compartment_norm == "Surface")].iloc[0].to_dict()
    carbon = CarbonBudgetEngine()
    policy = dict(ParetoPortfolioOptimizer.CANDIDATE_PORTFOLIOS[1])
    original = carbon._evaluate_policy_on_cell(row, row["map_mm"], policy)
    policy["rain_pulse_mm"] = 100
    changed = carbon._evaluate_policy_on_cell(row, row["map_mm"], policy)
    record("Regional carbon responds to portfolio rainfall", original != changed,
           {"rainfall_mm": [6.5, 100], "outputs_identical": original == changed})
    zero = dict(row, share_1l_mapped=.2, share_1l_cbb_mapped=0., share_1l_cbb_3gene_mapped=0.)
    result = carbon.compute_state_carbon_budget(zero)
    record("Zero rbcL co-carriage preserves zero microbial fixation", result["direct_microbial_rubisco_t_co2_km2_yr"]["point"] == 0, result["direct_microbial_rubisco_t_co2_km2_yr"])

    budgets = ParetoPortfolioOptimizer().evaluate_basin_carbon_budgets(frame, grid)
    p2 = budgets["candidate_portfolio_budgets"][1]
    record("Carbon pathway arithmetic", abs(sum(p2[k]["corridor_kt_yr_point"] for k in ["vegetation_npp_co2", "pedogenic_carbonate_co2", "direct_microbial_rubisco_co2", "indirect_h2_sink_co2eq"]) - p2["total_co2eq_bound"]["corridor_kt_yr_point"]) < .01,
           {"P2_corridor": p2, "vegetation_fraction_of_total": p2["vegetation_npp_co2"]["corridor_kt_yr_point"] / p2["total_co2eq_bound"]["corridor_kt_yr_point"]})
    surface = frame[frame.compartment_norm == "Surface"].set_index("site_num")
    assigned = grid[grid.in_empirical_support].nearest_primary_site.map(surface.wgs_observed)
    record("All supported grid baselines have observed surface WGS", bool(assigned.all()),
           {"supported_cells": len(assigned), "cells_assigned_imputed_surface_WGS": int((~assigned).sum())})

    snapshot = {}
    for path in sorted((root / "src").rglob("*.py")):
        snapshot[str(path.relative_to(root))] = hashlib.sha256(path.read_bytes()).hexdigest()
    return {"findings": findings, "source_sha256": snapshot}


if __name__ == "__main__":
    report = audit()
    path = Path(__file__).resolve().parents[1] / "reviews/validation-2026-10-26/model-contracts.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    for result in report["findings"]:
        print(f"{result['status'].upper()}: {result['check']}")
    print(path)
