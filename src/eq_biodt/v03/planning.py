"""Mechanism-discrimination experiments and scenario sensitivity, without probabilities."""
from __future__ import annotations

from dataclasses import asdict, replace
from itertools import combinations, product

import numpy as np

from .process import Parameters, State, gas_exchange, simulate


def parameter_ensemble(base: Parameters) -> list[Parameters]:
    """Explicit sensitivity grid; its members are not posterior samples."""
    return [replace(base, retention_b=b, dry_threshold_mpa=dry,
                    activation_shape=shape, k_h2_s=base.k_h2_s*rate,
                    k_co_s=base.k_co_s*rate)
            for b, dry, shape, rate in product((3., 5.), (-70., -100.), (.5, 2.), (.1, 1., 10.))]


def compare(initial: State, scenarios: list[dict], members: list[Parameters]) -> dict:
    if not scenarios or not members:
        raise ValueError("Scenarios and explicit parameter members required")
    names = [s["name"] for s in scenarios]
    if len(set(names)) != len(names):
        raise ValueError("Scenario names must be unique")
    if len({len(s["forcing"]) for s in scenarios}) != 1:
        raise ValueError("Compare scenarios over the same time interval")
    outcomes = []
    fronts = []
    for i,p in enumerate(members):
        rows = []
        for s in scenarios:
            if set(s)-{"name", "forcing", "disturbance_m2", "porosity"}:
                raise ValueError("Unsupported scenario field")
            cost = s["disturbance_m2"]
            if not np.isfinite(cost) or cost < 0:
                raise ValueError("Explicit nonnegative disturbed area required")
            pp = replace(p, porosity=s["porosity"]) if "porosity" in s else p
            result = simulate(initial, s["forcing"], pp)["summary"]
            rows.append({"member": i, "scenario": s["name"], **result, "disturbance_m2": cost})
        valid = all(r["h2_uptake_mol_m2"] is not None for r in rows)
        front = []
        if valid:
            # All four objectives are explicit: uptake up; water, salt, disturbance down.
            values = np.array([[-r["h2_uptake_mol_m2"], r["water_input_mm"],
                                r["final_salt_mol_m2"], r["disturbance_m2"]] for r in rows])
            for j, row in enumerate(rows):
                dominated = any(np.all(values[k] <= values[j]) and np.any(values[k] < values[j])
                                for k in range(len(rows)) if k != j)
                if not dominated:
                    front.append(row["scenario"])
        fronts.append({"member": i, "pareto": front if valid else None})
        outcomes.extend(rows)
    reversals = []
    for a,b in combinations(names, 2):
        differences = []
        for i in range(len(members)):
            aa = next(r for r in outcomes if r["member"] == i and r["scenario"] == a)
            bb = next(r for r in outcomes if r["member"] == i and r["scenario"] == b)
            if aa["h2_uptake_mol_m2"] is not None and bb["h2_uptake_mol_m2"] is not None:
                differences.append(aa["h2_uptake_mol_m2"]-bb["h2_uptake_mol_m2"])
        reversals.append({"a": a, "b": b, "uptake_difference_range_mol_m2":
                          [min(differences), max(differences)] if differences else None,
                          "rank_reversal": bool(differences and min(differences) < 0 < max(differences))})
    complete = all(f["pareto"] is not None for f in fronts)
    persistent = sorted(set.intersection(*(set(f["pareto"]) for f in fronts))) if complete else None
    return {"status": "conditional_sensitivity_comparison", "parameters": [asdict(p) for p in members],
            "objectives": {"h2_uptake": "maximise", "water_input": "minimise",
                           "residual_salt": "minimise", "disturbed_area": "minimise"},
            "interpretation": "Equal enumeration is sensitivity, not probability. No biodiversity desirability score.",
            "outcomes": outcomes, "pareto_by_member": fronts,
            "pareto_in_every_member": persistent, "pairwise_uptake_comparisons": reversals}


def experiment(p: Parameters, *, replicates: int = 3, seed: int = 42,
               measurement_noise_mol_m2_s: float = 1e-9) -> dict:
    if replicates < 2 or not np.isfinite(measurement_noise_mol_m2_s) or measurement_noise_mol_m2_s <= 0:
        raise ValueError("At least two blocks and positive assumed measurement noise required")
    cells = []
    for psi, porosity in product((-150., -85., -30., -3., -.3, -.03), (.3, .4, .5)):
        pp = replace(p, porosity=porosity)
        theta = porosity*(psi/pp.psi_entry_mpa)**(-1/pp.retention_b)
        active = gas_exchange(theta, 25., pp, measured_matric_mpa=psi)["h2_mol_m2_s"]
        alternative = gas_exchange(theta, 25., replace(pp, mechanism="diffusion_only"),
                                   measured_matric_mpa=psi)["h2_mol_m2_s"]
        cells.append({"target_water_potential_mpa": psi, "target_porosity": porosity,
                      "illustrative_theta": theta, "activation_model_h2": active,
                      "diffusion_only_model_h2": alternative,
                      "disagreement_to_assumed_noise": abs(active-alternative)/measurement_noise_mol_m2_s,
                      "role": "held_out_condition" if psi == -3. else "fit_condition"})
    rng = np.random.default_rng(seed)
    schedule = []
    for block in range(replicates):
        for order,i in enumerate(rng.permutation(len(cells))):
            schedule.append({"block": block+1, "order": order+1, "cell": int(i), **cells[i]})
    return {"design": "randomised complete block moisture by measured porosity factorial",
            "parameters": asdict(p), "seed": seed, "replicates": replicates,
            "assumed_measurement_noise_mol_m2_s": measurement_noise_mol_m2_s,
            "cells": cells, "schedule": schedule,
            "measurement_protocol": [
                "Use independent soil units per treatment within spatial blocks; preserve an intact-soil control.",
                "Measure water potential, retention curves, bulk density/porosity and temperature; theta is illustrative.",
                "Measure H2, CO and O2 flux with chamber blanks and killed controls; report limits of detection.",
                "Collect WGS plus absolute biomass/cell estimates and RNA/protein or isotope activity measurements.",
                "Measure at baseline and 1, 6, 24, 72 hours and 7 days after wetting; repeated reads are not replicates.",
                "Fit kinetics using fit conditions only; keep every -3 MPa unit withheld for predictive evaluation.",
                "Record salinity and assay methods; high ionic strength requires nonideal speciation.",
                "Final replication requires pilot variance; this design does not claim statistical power."],
            "selection_criterion": "Full factorial with reported model disagreement; no unmeasured information-gain claim."}
