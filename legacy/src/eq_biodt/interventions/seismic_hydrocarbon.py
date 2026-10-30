"""Intervention 6: Seismic Exploration (Vibroseis Compaction vs Nodal Surveys) & Hydrocarbon Bioremediation.

Grounded in:
- Canonical sites `Site oill spill` (`https://rubalkhali.science/kb/RAK_1000061`, 19.40536 N, 45.72708 E),
  `Site road (location 1)` (`RAK_1000068`), `Site road (location 2)` (`RAK_1000069`), verified
  ENVO term `ENVO_00002061` (`oil contaminated soil`), and verified XRF Vanadium marker
  `CHEBI_27698` (`V` in `xrf_chemical_mapping.yml`).
- `~/Documents/papers/empty-quarter-lifeonair`: Vibroseis axle compaction reduces total
  soil porosity eps_total, nonlinearly throttling Millington-Quirk H2/CO diffusivity
  D_eff = D_0 * eps_a^(10/3) / eps_total^2.
- `~/Documents/papers/empty-quarter-depth`: Bulldozer seismic cut-line grading strips the
  0-5 cm surface layer and exposes the 5-10 cm Thermoplasmatota/Nitrospira refugium.
- `~/Documents/papers/empty-quarter-enrichments`: Indigenous alkane-degrading (alkB/ladA)
  Actinomycetota/Pseudomonadota enrichment consortia for crude-oil bioremediation.
"""

from __future__ import annotations

from typing import Any, Dict, Optional

import numpy as np

from eq_biodt.models.causal_ontological_twin import CausalOntologicalTwinEngine
from eq_biodt.models.moisture_pulse_diffusion import MoisturePulseDiffusionModel
from eq_biodt.standards.fdo_provenance import FDOProvenanceRecorder


class SeismicHydrocarbonIntervention:
    """Counterfactual simulator for 2D/3D seismic exploration impacts and oil-spill bioremediation."""

    def __init__(
        self,
        causal_engine: Optional[CausalOntologicalTwinEngine] = None,
        diffusion_model: Optional[MoisturePulseDiffusionModel] = None,
    ) -> None:
        self.causal_engine = causal_engine or CausalOntologicalTwinEngine()
        self.diffusion_model = diffusion_model or MoisturePulseDiffusionModel()

    def simulate(
        self,
        site_row: Dict[str, Any],
        seismic_line_density_km_km2: float = 4.0,
        vibroseis_axle_passes: int = 6,
        cut_line_grading_fraction: float = 0.25,
        low_impact_nodal_mode: bool = False,
        hydrocarbon_tph_g_kg: float = 0.0,
        alkb_consortium_bioaug_dose: float = 0.0,
        remediation_horizon_days: int = 90,
    ) -> Dict[str, Any]:
        """Simulate seismic survey disturbance and/or hydrocarbon spill bioremediation."""
        line_dens = float(np.clip(seismic_line_density_km_km2, 0.0, 25.0))
        passes = int(np.clip(vibroseis_axle_passes, 0, 50))
        grading = float(np.clip(cut_line_grading_fraction, 0.0, 1.0))
        tph_initial = float(np.clip(hydrocarbon_tph_g_kg, 0.0, 120.0))
        bioaug = float(np.clip(alkb_consortium_bioaug_dose, 0.0, 1.0))
        days = int(np.clip(remediation_horizon_days, 1, 730))

        base_shannon = float(site_row.get("shannon_16s", 4.85))
        base_h1l = float(site_row.get("hydB_1l_ge", 0.61))
        base_h1h = float(site_row.get("hydB_1h_ge", 0.12))
        base_h2a = float(site_row.get("hydB_2a_ge", 0.03))
        base_cox = float(site_row.get("coxL_formI_ge", 1.15))
        base_v_xrf = float(site_row.get("V", 0.008))
        site_iri = str(site_row.get("site_iri", "https://rubalkhali.science/kb/RAK_1000061"))

        # 1. Seismic compaction of soil porosity eps_total
        # Low-impact nodal/drone seismic reduces axle compaction & bulldozer grading by 82%
        impact_mitigation = 0.18 if low_impact_nodal_mode else 1.0
        compaction_intensity = impact_mitigation * (
            0.014 * np.log1p(passes) * (line_dens / (3.0 + line_dens)) + 0.11 * grading
        )
        effective_porosity = float(np.clip(0.40 - compaction_intensity, 0.19, 0.40))

        # 5-10 cm subsurface refugium integrity (destroyed along bulldozer cut-lines)
        refugium_integrity = float(np.clip(1.0 - 0.85 * grading * impact_mitigation, 0.10, 1.0))

        # 2. Hydrocarbon spill & XRF Vanadium Switch (CHEBI_27698)
        # First-order alkB/ladA biodegradation rate k_deg (d^-1):
        # Natural attenuation in hyperarid sand is slow (k ~ 0.0025 d^-1); indigenous
        # Actinomycetota/Pseudomonadota enrichment consortium accelerates k up to 0.018 d^-1.
        k_deg = 0.0025 + 0.0155 * bioaug * (effective_porosity / 0.40)
        residual_tph = tph_initial * float(np.exp(-k_deg * days))
        tph_removal_pct = (
            100.0 * (tph_initial - residual_tph) / tph_initial if tph_initial > 0.0 else 100.0
        )
        counterfactual_vanadium_pct = base_v_xrf + 0.0042 * tph_initial

        # Hydrocarbon pore-clogging further reduces effective air-filled porosity
        pore_oil_clog = min(0.08, 0.0025 * residual_tph)
        gas_res = self.diffusion_model.milling_quirk_gas_power(
            theta_vwc=0.018 + pore_oil_clog,
            total_porosity=effective_porosity,
            hydB_1l_ge=base_h1l * np.exp(-0.045 * residual_tph) * (0.5 + 0.5 * refugium_integrity),
            coxL_formI_ge=base_cox * np.exp(-0.030 * residual_tph),
        )

        # Taxonomic & guild shifts:
        # - Thermoplasmatota drops sharply when 5-10 cm refugium is graded or contaminated by TPH
        therm_cf = (
            float(site_row.get("thermoplasmatota_rel_abund", 0.045))
            * refugium_integrity
            * np.exp(-0.06 * residual_tph)
        )
        # - Alkane-degrading Actinomycetota (Rhodococcus/Dietzia) expand during active bioremediation
        act_cf = float(
            np.clip(
                float(site_row.get("actinomycetota_rel_abund", 0.48))
                * (1.0 + 0.15 * bioaug * min(1.0, tph_initial / 10.0) - 0.12 * (1.0 - refugium_integrity)),
                0.15,
                0.75,
            )
        )

        # Diversity response
        seismic_shannon_loss = 0.85 * (1.0 - refugium_integrity) + 1.8 * (0.40 - effective_porosity)
        tph_shannon_loss = 0.055 * residual_tph
        remediation_bonus = 0.22 * bioaug * (1.0 - np.exp(-0.1 * tph_initial))
        shannon_cf = max(1.8, base_shannon - seismic_shannon_loss - tph_shannon_loss + remediation_bonus)

        raw_leaf_state = {
            "actinomycetota_rel_abund": act_cf,
            "thermoplasmatota_rel_abund": therm_cf,
            "halobacteriota_rel_abund": float(site_row.get("halobacteriota_rel_abund", 0.02)),
            "hydB_1l_ge": base_h1l * gas_res["relative_gas_diffusivity"] * np.exp(-0.045 * residual_tph),
            "hydB_1h_ge": base_h1h * (1.0 + 0.25 * bioaug),
            "hydB_2a_ge": base_h2a * (1.0 + 0.15 * bioaug),
            "coxL_formI_ge": base_cox * gas_res["relative_gas_diffusivity"] * np.exp(-0.030 * residual_tph),
            "h2_oxidation_power_uw_m3": gas_res["h2_oxidation_power_uw_m3"],
            "co_oxidation_power_uw_m3": gas_res["co_oxidation_power_uw_m3"],
            "viable_dna_fraction": float(np.clip(0.56 * refugium_integrity + 0.18 * bioaug, 0.15, 0.90)),
            "shannon_16s": shannon_cf,
        }
        projected_state = self.causal_engine.evaluate_constitutive_layer(raw_leaf_state)

        amb_interval = self.causal_engine.compute_version_ambiguity_interval(
            target_parent_total=max(0.0, line_dens),
            leaf_sensitivities={
                "low_impact_nodal_seismic_with_alkB_bioaug": -0.008,
                "conventional_vibroseis_bulldozer_cutlines": -0.064,
            },
        )

        advisory = (
            "LOW_IMPACT_SEISMIC_AND_REMEDIATED: Wireless nodal survey preserves Millington-Quirk H2 diffusivity and 5-10 cm Thermoplasmatota refugium."
            if (effective_porosity >= 0.35 and refugium_integrity >= 0.80 and residual_tph <= 2.5)
            else "WARNING_POROSITY_COMPACTION_OR_HYDROCARBON_TOXICITY: Switch to nodal/drone receivers, avoid bulldozer cut-line grading, and inoculate alkB/ladA enrichment consortium."
        )

        summary = {
            "intervention": "seismic_hydrocarbon",
            "site_iri": site_iri,
            "seismic_line_density_km_km2": round(line_dens, 2),
            "vibroseis_axle_passes": passes,
            "cut_line_grading_fraction": round(grading, 3),
            "low_impact_nodal_mode": low_impact_nodal_mode,
            "effective_soil_porosity": round(effective_porosity, 4),
            "subsurface_refugium_integrity": round(refugium_integrity, 4),
            "relative_gas_diffusivity": gas_res["relative_gas_diffusivity"],
            "initial_tph_g_kg": round(tph_initial, 3),
            "residual_tph_g_kg": round(residual_tph, 3),
            "tph_biodegradation_pct": round(tph_removal_pct, 2),
            "xrf_vanadium_marker_pct": round(counterfactual_vanadium_pct, 5),
            "xrf_vanadium_chebi_iri": "http://purl.obolibrary.org/obo/CHEBI_27698",
            "counterfactual_shannon_16s": round(projected_state["shannon_16s"], 4),
            "h2_oxidation_power_uw_m3": round(projected_state["h2_oxidation_power_uw_m3"], 4),
            "ocs_score": projected_state["ocs_score"],
            "management_advisory": advisory,
        }

        fdo_crate = FDOProvenanceRecorder.build_run_crate(
            intervention_name="seismic_hydrocarbon",
            site_iri=site_iri,
            parameters={
                "seismic_line_density_km_km2": line_dens,
                "vibroseis_axle_passes": passes,
                "cut_line_grading_fraction": grading,
                "low_impact_nodal_mode": low_impact_nodal_mode,
                "hydrocarbon_tph_g_kg": tph_initial,
                "alkb_consortium_bioaug_dose": bioaug,
                "remediation_horizon_days": days,
            },
            outputs=summary,
            ocs_score=projected_state["ocs_score"],
            version_ambiguity_interval=amb_interval,
        )

        return {
            "summary": summary,
            "projected_state": projected_state,
            "fdo_provenance": fdo_crate["fdo_kernel_metadata"],
        }
