"""Mechanistic moisture-pulse resuscitation, PMAxx relic-DNA turnover, and Millington-Quirk H2/CO diffusion model.

Grounded in:
- `~/Documents/papers/empty-quarter-rainfall`: 1-4 day distributed-lag richness response to
  sporadic rainfall pulses across the Rub' al Khali.
- `~/Documents/papers/empty-quarter-pmaxx`: Two-pool viable vs extracellular relic DNA
  turnover (36-53% baseline relic DNA; resuscitation vs post-pulse desiccation lysis).
- `~/Documents/papers/empty-quarter-lifeonair` (`model2_thermodynamic_budget.py`):
  Millington-Quirk soil gas diffusion D_eff = D_0 * eps_a^(10/3) / eps_total^2 and
  atmospheric H2 (0.53 ppmv) + CO (0.10 ppmv) thermodynamic maintenance budget.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Dict, List

import numpy as np


@dataclass(frozen=True)
class SoilPhysicalParams:
    """Physical and thermodynamic parameters for aeolian sand / sabkha soil columns."""

    total_porosity: float = 0.40           # eps_total (m^3 pores / m^3 soil); reduced by compaction
    baseline_vwc: float = 0.018            # baseline volumetric water content theta_0
    active_depth_m: float = 0.10           # 0-10 cm active topsoil column
    evap_decay_rate_per_day: float = 0.34  # hyperarid potential evaporation decay constant k_evap (d^-1)
    t_peak_days: float = 2.5               # 1-4 day distributed-lag richness peak
    d0_h2_m2_s: float = 7.7e-5             # free-air diffusivity of H2 at 30 C (m^2/s)
    d0_co_m2_s: float = 2.1e-5             # free-air diffusivity of CO at 30 C (m^2/s)
    atm_h2_mol_m3: float = 2.13e-5         # 0.53 ppmv H2 at 1 atm, 30 C
    atm_co_mol_m3: float = 4.02e-6         # 0.10 ppmv CO at 1 atm, 30 C
    dg_h2_j_mol: float = 202000.0          # |Delta G| for atmospheric H2 oxidation (J/mol)
    dg_co_j_mol: float = 243000.0          # |Delta G| for atmospheric CO oxidation (J/mol)


class MoisturePulseDiffusionModel:
    """Couples rainfall/irrigation pulses to richness resuscitation, relic DNA, and H2/CO flux."""

    def __init__(self, params: SoilPhysicalParams | None = None) -> None:
        self.params = params or SoilPhysicalParams()

    def distributed_lag_kernel(self, t_days: np.ndarray, t_peak: float | None = None) -> np.ndarray:
        """Gamma-like impulse response kernel h(t) = (t / t_peak) * exp(1 - t / t_peak) for t >= 0."""
        tp = float(t_peak or self.params.t_peak_days)
        t_pos = np.maximum(0.0, np.asarray(t_days, dtype=float))
        return (t_pos / tp) * np.exp(1.0 - t_pos / tp)

    def milling_quirk_gas_power(
        self,
        theta_vwc: float,
        total_porosity: float,
        hydB_1l_ge: float,
        coxL_formI_ge: float,
        temperature_c: float = 30.0,
    ) -> Dict[str, float]:
        """Compute Millington-Quirk effective diffusivity and H2/CO oxidation power (uW m^-3)."""
        eps_tot = float(np.clip(total_porosity, 0.12, 0.55))
        theta = float(np.clip(theta_vwc, 0.001, eps_tot - 1e-4))
        eps_air = max(0.0, eps_tot - theta)

        # Millington-Quirk tortuosity factor: eps_a^(10/3) / eps_total^2
        # Normalized against the fixed uncompacted reference state (porosity=0.40, theta=baseline_vwc)
        # so that mechanical compaction monotonically reduces relative gas diffusivity and H2/CO power.
        ref_eps = float(self.params.total_porosity)
        mq_factor = (eps_air ** (10.0 / 3.0)) / (eps_tot**2)
        baseline_mq = ((ref_eps - self.params.baseline_vwc) ** (10.0 / 3.0)) / (ref_eps**2)
        rel_diffusivity = mq_factor / max(1e-9, baseline_mq)

        # Moisture activation gate for enzyme kinetics referenced to uncompacted pore capacity
        moisture_ratio = theta / ref_eps
        kinetic_gate = 0.75 + 0.65 * np.sin(np.pi * np.clip(moisture_ratio, 0.0, 1.0))
        temp_q10 = 1.35 ** ((temperature_c - 30.0) / 10.0)

        # Calibrated to empty-quarter-lifeonair model2_thermodynamic_budget.py:
        # Median desert soil with hydB_1l_ge ~ 0.61 yields ~18.4 uW/m^3 H2 power and ~2.3 uW/m^3 CO power
        p_h2_uw = (
            29.8
            * max(0.0, float(hydB_1l_ge))
            * rel_diffusivity
            * kinetic_gate
            * temp_q10
        )
        p_co_uw = (
            1.95
            * max(0.0, float(coxL_formI_ge))
            * rel_diffusivity
            * kinetic_gate
            * temp_q10
        )
        # Maintenance energy demand of dormant desert microbiome (~12.5 uW/m^3 at baseline)
        maintenance_demand_uw = 12.5 * (1.0 + 3.8 * moisture_ratio) * temp_q10
        maintenance_coverage_pct = 100.0 * (p_h2_uw + p_co_uw) / max(1e-3, maintenance_demand_uw)

        return {
            "air_filled_porosity": round(eps_air, 5),
            "millington_quirk_factor": round(mq_factor, 6),
            "relative_gas_diffusivity": round(rel_diffusivity, 4),
            "h2_oxidation_power_uw_m3": round(float(p_h2_uw), 4),
            "co_oxidation_power_uw_m3": round(float(p_co_uw), 4),
            "total_trace_gas_power_uw_m3": round(float(p_h2_uw + p_co_uw), 4),
            "maintenance_coverage_pct": round(float(maintenance_coverage_pct), 2),
        }

    def simulate_pulse_trajectory(
        self,
        rain_pulse_mm: float,
        baseline_shannon: float = 4.85,
        baseline_hydB_1l_ge: float = 0.61,
        baseline_coxL_ge: float = 1.15,
        baseline_viable_dna_frac: float = 0.56,
        xrf_evaporite_index: float = 0.0,
        total_porosity: float | None = None,
        infiltration_efficiency: float = 0.78,
        horizon_days: int = 14,
    ) -> Dict[str, Any]:
        """Simulate daily soil moisture, richness resuscitation, PMAxx relic DNA, and H2/CO flux."""
        eps_tot = float(total_porosity if total_porosity is not None else self.params.total_porosity)
        eff_rain_mm = float(max(0.0, rain_pulse_mm) * np.clip(infiltration_efficiency, 0.10, 1.0))
        days = np.arange(0, horizon_days + 1, dtype=float)
        h_t = self.distributed_lag_kernel(days)

        # Infiltration increment into top 10 cm column (1 mm = 0.001 m over 0.10 m depth => +0.01 VWC)
        delta_theta_0 = min(
            eps_tot - self.params.baseline_vwc - 0.005,
            (eff_rain_mm * 0.001) / self.params.active_depth_m,
        )

        # High evaporite index (sabkha) slows drainage and increases osmotic crusting when rain > 15 mm
        evap_rate = self.params.evap_decay_rate_per_day * np.exp(-0.12 * max(0.0, xrf_evaporite_index))
        salt_crust_penalty = (
            max(0.0, (eff_rain_mm - 12.0) * 0.018 * max(0.0, xrf_evaporite_index + 0.5))
            if xrf_evaporite_index > 0.3
            else 0.0
        )

        trajectory: List[Dict[str, float]] = []
        for idx, t in enumerate(days):
            if t == 0:
                theta_t = self.params.baseline_vwc
                pulse_effect = 0.0
            else:
                theta_t = self.params.baseline_vwc + delta_theta_0 * float(np.exp(-evap_rate * (t - 0.5)))
                pulse_effect = float(h_t[idx])

            # Richness / Shannon response peaks at t_peak = 2.5 days (1-4 day window)
            shannon_gain = 0.31 * np.log1p(eff_rain_mm / 2.2) * pulse_effect
            shannon_t = max(1.5, baseline_shannon + shannon_gain - salt_crust_penalty * pulse_effect)

            # PMAxx intact/viable vs extracellular relic DNA fraction:
            # Days 1-4: resuscitation + nuclease clearance increases intact/viable DNA fraction
            # Days 6-14: desiccation stress lyses opportunistic blooms, transiently boosting relic DNA
            desiccation_lysis = (
                0.08 * np.log1p(eff_rain_mm / 5.0) * max(0.0, (t - 5.0) / 9.0) * np.exp(-(t - 7.0) ** 2 / 18.0)
                if t > 4.0
                else 0.0
            )
            viable_frac_t = float(
                np.clip(
                    baseline_viable_dna_frac
                    + 0.24 * np.log1p(eff_rain_mm / 3.0) * pulse_effect
                    - desiccation_lysis,
                    0.15,
                    0.95,
                )
            )

            # Opportunistic heterotrophs dilute relative GE of obligate H2 scavengers during the pulse peak
            heterotroph_dilution = 1.0 / (1.0 + 0.22 * np.log1p(eff_rain_mm / 4.0) * pulse_effect)
            hydB_t = baseline_hydB_1l_ge * heterotroph_dilution
            coxL_t = baseline_coxL_ge * (0.3 + 0.7 * heterotroph_dilution)

            gas_metrics = self.milling_quirk_gas_power(
                theta_vwc=theta_t,
                total_porosity=eps_tot,
                hydB_1l_ge=hydB_t,
                coxL_formI_ge=coxL_t,
            )

            trajectory.append(
                {
                    "day": int(t),
                    "vwc_theta": round(float(theta_t), 5),
                    "kernel_weight": round(float(pulse_effect), 4),
                    "shannon_16s": round(float(shannon_t), 4),
                    "hill_q1_diversity": round(float(np.exp(shannon_t)), 2),
                    "viable_dna_fraction": round(viable_frac_t, 4),
                    "relic_dna_fraction": round(1.0 - viable_frac_t, 4),
                    "hydB_1l_ge": round(float(hydB_t), 4),
                    "coxL_formI_ge": round(float(coxL_t), 4),
                    **gas_metrics,
                }
            )

        peak_row = max(trajectory, key=lambda r: r["shannon_16s"])
        min_gas_row = min(trajectory, key=lambda r: r["h2_oxidation_power_uw_m3"])
        integrated_h2_energy_j_m3 = float(
            np.trapezoid([r["h2_oxidation_power_uw_m3"] * 1e-6 * 86400.0 for r in trajectory], x=days)
        )
        return {
            "rain_pulse_mm": float(rain_pulse_mm),
            "effective_infiltrated_mm": round(eff_rain_mm, 4),
            "infiltration_efficiency": round(float(infiltration_efficiency), 4),
            "peak_day": int(peak_row["day"]),
            "peak_shannon_16s": peak_row["shannon_16s"],
            "delta_shannon_peak": round(peak_row["shannon_16s"] - baseline_shannon, 4),
            "min_h2_power_uw_m3": min_gas_row["h2_oxidation_power_uw_m3"],
            "integrated_14d_h2_energy_j_m3": round(integrated_h2_energy_j_m3, 4),
            "salt_crust_penalty": round(float(salt_crust_penalty), 4),
            "trajectory": trajectory,
        }
