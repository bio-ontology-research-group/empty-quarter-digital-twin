"""Intervention 1: Cloud Seeding & Artificial Precipitation Pulse Simulator.

Simulates artificial rainfall pulses (Delta P in [0, 35] mm) across Empty Quarter sites,
coupling:
- 1-4 day distributed-lag richness resuscitation (`empty-quarter-rainfall`)
- PMAxx viable vs extracellular relic DNA turnover (`empty-quarter-pmaxx`)
- Millington-Quirk H2/CO diffusion throttling (`empty-quarter-lifeonair`)
- Evaporite sabkha salt-crust risk (`empty-quarter-map`)
- Two-layer causal-ontological projection (OCS == 1.0) and FDO provenance.
"""

from __future__ import annotations

from typing import Any, Dict, Optional

import numpy as np

from eq_biodt.models.causal_ontological_twin import CausalOntologicalTwinEngine
from eq_biodt.models.moisture_pulse_diffusion import MoisturePulseDiffusionModel
from eq_biodt.standards.fdo_provenance import FDOProvenanceRecorder


class CloudSeedingIntervention:
    """Counterfactual simulator for cloud-seeding rainfall pulses (`do(RainPulse_mm = p)`)."""

    def __init__(
        self,
        pulse_model: Optional[MoisturePulseDiffusionModel] = None,
        causal_engine: Optional[CausalOntologicalTwinEngine] = None,
    ) -> None:
        self.pulse_model = pulse_model or MoisturePulseDiffusionModel()
        self.causal_engine = causal_engine or CausalOntologicalTwinEngine()

    def simulate(
        self,
        site_row: Dict[str, Any],
        rain_pulse_mm: float = 8.0,
        horizon_days: int = 14,
    ) -> Dict[str, Any]:
        """Run a cloud-seeding counterfactual at a given site/compartment row."""
        base_shannon = float(site_row.get("shannon_16s", 4.85))
        base_h1l = float(site_row.get("hydB_1l_ge", 0.61))
        base_h1h = float(site_row.get("hydB_1h_ge", 0.12))
        base_h2a = float(site_row.get("hydB_2a_ge", 0.03))
        base_cox = float(site_row.get("coxL_formI_ge", 1.15))
        evap_idx = float(site_row.get("xrf_evaporite_index", 0.0))
        site_iri = str(site_row.get("site_iri", "https://rubalkhali.science/kb/RAK_1000001"))

        pulse_res = self.pulse_model.simulate_pulse_trajectory(
            rain_pulse_mm=rain_pulse_mm,
            baseline_shannon=base_shannon,
            baseline_hydB_1l_ge=base_h1l,
            baseline_coxL_ge=base_cox,
            xrf_evaporite_index=evap_idx,
            horizon_days=horizon_days,
        )

        # Extract peak-window state (day 2-3) and project onto ontological polytope
        peak_day_idx = min(len(pulse_res["trajectory"]) - 1, max(1, pulse_res["peak_day"]))
        peak_pt = pulse_res["trajectory"][peak_day_idx]

        raw_leaf_state = {
            "actinomycetota_rel_abund": float(site_row.get("actinomycetota_rel_abund", 0.48))
            * (1.0 - 0.05 * np.log1p(rain_pulse_mm / 5.0)),
            "thermoplasmatota_rel_abund": float(site_row.get("thermoplasmatota_rel_abund", 0.04))
            * (1.0 + 0.08 * np.log1p(rain_pulse_mm / 5.0)),
            "halobacteriota_rel_abund": float(site_row.get("halobacteriota_rel_abund", 0.02))
            * (1.0 + 0.12 * float(pulse_res["salt_crust_penalty"])),
            "hydB_1l_ge": float(peak_pt["hydB_1l_ge"]),
            "hydB_1h_ge": base_h1h * (1.0 + 0.06 * np.log1p(rain_pulse_mm / 5.0)),
            "hydB_2a_ge": base_h2a * (1.0 + 0.08 * np.log1p(rain_pulse_mm / 5.0)),
            "coxL_formI_ge": float(peak_pt["coxL_formI_ge"]),
            "h2_oxidation_power_uw_m3": float(peak_pt["h2_oxidation_power_uw_m3"]),
            "co_oxidation_power_uw_m3": float(peak_pt["co_oxidation_power_uw_m3"]),
            "viable_dna_fraction": float(peak_pt["viable_dna_fraction"]),
            "shannon_16s": float(peak_pt["shannon_16s"]),
        }
        projected_state = self.causal_engine.evaluate_constitutive_layer(raw_leaf_state)

        # Version-ambiguity interval across convective downpour (runoff loss, eta=0.55)
        # vs stratiform drizzle (high infiltration, eta=0.92) propagated through the pulse model,
        # plus linear Proposition 1 bounds for comparison.
        pulse_conv = self.pulse_model.simulate_pulse_trajectory(
            rain_pulse_mm=rain_pulse_mm,
            baseline_shannon=base_shannon,
            baseline_hydB_1l_ge=base_h1l,
            baseline_coxL_ge=base_cox,
            xrf_evaporite_index=evap_idx,
            infiltration_efficiency=0.55,
            horizon_days=horizon_days,
        )
        pulse_strat = self.pulse_model.simulate_pulse_trajectory(
            rain_pulse_mm=rain_pulse_mm,
            baseline_shannon=base_shannon,
            baseline_hydB_1l_ge=base_h1l,
            baseline_coxL_ge=base_cox,
            xrf_evaporite_index=evap_idx,
            infiltration_efficiency=0.92,
            horizon_days=horizon_days,
        )
        nonlinear_shannon_bounds = (
            min(float(pulse_conv["delta_shannon_peak"]), float(pulse_strat["delta_shannon_peak"])),
            max(float(pulse_conv["delta_shannon_peak"]), float(pulse_strat["delta_shannon_peak"])),
        )
        nonlinear_h2_bounds = (
            min(float(pulse_conv["min_h2_power_uw_m3"]), float(pulse_strat["min_h2_power_uw_m3"])),
            max(float(pulse_conv["min_h2_power_uw_m3"]), float(pulse_strat["min_h2_power_uw_m3"])),
        )
        amb_interval = self.causal_engine.compute_version_ambiguity_interval(
            target_parent_total=rain_pulse_mm,
            leaf_sensitivities={
                "stratiform_drizzle_infiltration": 0.032,
                "convective_downpour_runoff": 0.011,
            },
            min_leaf_shares={"stratiform_drizzle_infiltration": 0.25, "convective_downpour_runoff": 0.15},
        )

        # Optimal window guidance derived directly from realized H2 retention >= 65% and ΔH' >= 0.15
        base_h2_power = float(pulse_res["trajectory"][0]["h2_oxidation_power_uw_m3"])
        h2_retention_ratio = float(pulse_res["min_h2_power_uw_m3"]) / max(1e-4, base_h2_power)
        delta_shannon = float(projected_state["shannon_16s"] - base_shannon)
        management_advisory = (
            "OPTIMAL_PULSE: Resuscitates dormant diversity while preserving >=65% atmospheric H2 oxidation."
            if (delta_shannon >= 0.15 and h2_retention_ratio >= 0.65 and pulse_res["salt_crust_penalty"] < 0.02)
            else (
                "WARNING_SALT_CRUST_OR_ANOXIA: Pulse throttles Millington-Quirk H2 diffusion (<65% retention) or mobilizes sabkha evaporite salts."
                if (h2_retention_ratio < 0.65 or pulse_res["salt_crust_penalty"] >= 0.02)
                else "SUB_THRESHOLD_PULSE: Minimal microbial resuscitation beyond baseline."
            )
        )

        summary = {
            "intervention": "cloud_seeding",
            "site_iri": site_iri,
            "rain_pulse_mm": float(rain_pulse_mm),
            "effective_infiltrated_mm": pulse_res["effective_infiltrated_mm"],
            "peak_resuscitation_day": pulse_res["peak_day"],
            "baseline_shannon_16s": round(base_shannon, 4),
            "counterfactual_peak_shannon_16s": round(projected_state["shannon_16s"], 4),
            "delta_shannon_16s": round(delta_shannon, 4),
            "baseline_h2_power_uw_m3": round(base_h2_power, 4),
            "min_h2_power_uw_m3": round(float(pulse_res["min_h2_power_uw_m3"]), 4),
            "integrated_14d_h2_energy_j_m3": pulse_res["integrated_14d_h2_energy_j_m3"],
            "h2_power_retention_ratio": round(h2_retention_ratio, 4),
            "salt_crust_penalty": pulse_res["salt_crust_penalty"],
            "viable_dna_fraction_peak": round(projected_state["viable_dna_fraction"], 4),
            "ocs_score": projected_state["ocs_score"],
            "version_ambiguity_shannon_bounds": [round(amb_interval[0], 4), round(amb_interval[1], 4)],
            "nonlinear_trajectory_shannon_bounds": [
                round(nonlinear_shannon_bounds[0], 4),
                round(nonlinear_shannon_bounds[1], 4),
            ],
            "nonlinear_trajectory_min_h2_power_bounds": [
                round(nonlinear_h2_bounds[0], 4),
                round(nonlinear_h2_bounds[1], 4),
            ],
            "management_advisory": management_advisory,
        }

        fdo_crate = FDOProvenanceRecorder.build_run_crate(
            intervention_name="cloud_seeding",
            site_iri=site_iri,
            parameters={"rain_pulse_mm": rain_pulse_mm, "horizon_days": horizon_days},
            outputs=summary,
            ocs_score=projected_state["ocs_score"],
            version_ambiguity_interval=amb_interval,
        )

        return {
            "summary": summary,
            "projected_state": projected_state,
            "trajectory": pulse_res["trajectory"],
            "fdo_provenance": fdo_crate["fdo_kernel_metadata"],
        }
