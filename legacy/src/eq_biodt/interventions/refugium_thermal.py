"""Intervention 4: Subsurface 5-10 cm Refugium Protection & 60 C Thermal/Betaine Shielding.

Grounded in:
- `~/Documents/papers/empty-quarter-depth`: The first 5-10 cm subsurface horizon acts as
  the primary thermal/UV refugium, enriching Thermoplasmatota by 2.26x, Nitrospira, and
  archaeal 3-HP/4-HB carbon fixation relative to 0 cm surface sand.
- `~/Documents/papers/empty-quarter-models-thermophiles`: At 60 C surface extremes,
  thermophilic Bacillus survival requires exogenous choline/glycine-betaine uptake
  (OpuA/GbsAB pathway) to offset thermal protein-repair ATP maintenance.
"""

from __future__ import annotations

from typing import Any, Dict, Optional

import numpy as np

from eq_biodt.models.causal_ontological_twin import CausalOntologicalTwinEngine
from eq_biodt.models.metabolic_community_engine import (
    CommunityMetabolicInputs,
    MetabolicCommunityEngine,
)
from eq_biodt.standards.fdo_provenance import FDOProvenanceRecorder


class RefugiumThermalIntervention:
    """Counterfactual simulator for 5-10 cm subsurface refugium protection & betaine thermal priming."""

    def __init__(
        self,
        causal_engine: Optional[CausalOntologicalTwinEngine] = None,
        metabolic_engine: Optional[MetabolicCommunityEngine] = None,
    ) -> None:
        self.causal_engine = causal_engine or CausalOntologicalTwinEngine()
        self.metabolic_engine = metabolic_engine or MetabolicCommunityEngine()

    def simulate(
        self,
        site_row: Dict[str, Any],
        subsurface_refugium_intact_fraction: float = 0.95,
        surface_temp_c: float = 58.0,
        shade_shield_fraction: float = 0.4,
        exogenous_choline_betaine_mm: float = 0.5,
    ) -> Dict[str, Any]:
        """Simulate subsurface 5-10 cm preservation, thermal shading, and betaine priming."""
        intact = float(np.clip(subsurface_refugium_intact_fraction, 0.0, 1.0))
        t_surf = float(np.clip(surface_temp_c, 20.0, 68.0))
        shade = float(np.clip(shade_shield_fraction, 0.0, 1.0))
        betaine = float(np.clip(exogenous_choline_betaine_mm, 0.0, 5.0))

        # Effective soil horizon temperature damped by 5-10 cm depth refugium (-11 C) and shade (-7 C)
        effective_temp_c = t_surf - 11.0 * intact - 7.0 * shade
        thermal_over_50 = max(0.0, effective_temp_c - 48.0)

        base_shannon = float(site_row.get("shannon_16s", 5.05))
        base_therm = float(site_row.get("thermoplasmatota_rel_abund", 0.072))
        base_h1l = float(site_row.get("hydB_1l_ge", 0.61))
        base_h1h = float(site_row.get("hydB_1h_ge", 0.12))
        base_h2a = float(site_row.get("hydB_2a_ge", 0.03))
        base_cox = float(site_row.get("coxL_formI_ge", 1.15))
        site_iri = str(site_row.get("site_iri", "https://rubalkhali.science/kb/RAK_1000001"))

        fba_res = self.metabolic_engine.solve_community_fba(
            CommunityMetabolicInputs(
                h2_uptake_bound=1.8 * (0.5 + 0.5 * intact),
                co_uptake_bound=0.35 * (0.5 + 0.5 * intact),
                exogenous_nh4_no3_bound=0.09 * intact,
                temperature_c=effective_temp_c,
                choline_betaine_bound=betaine,
            )
        )

        # Thermoplasmatota requires intact 5-10 cm subsurface refugium (2.26x enrichment vs exposed surface)
        therm_cf = base_therm * (0.44 + 0.56 * intact) * np.exp(-0.06 * thermal_over_50)
        # Thermophilic Bacillus / Firmicutes survival at 60 C rescued by GbsAB/OpuA choline/betaine
        bacillus_60c_viability = float(
            np.clip(
                np.exp(-0.09 * max(0.0, t_surf - 45.0) / (1.0 + 1.6 * betaine + 0.8 * shade)),
                0.08,
                0.99,
            )
        )

        shannon_cf = (
            base_shannon
            - 0.65 * (1.0 - intact)
            - 0.045 * thermal_over_50 / (1.0 + 1.2 * betaine)
            + 0.15 * shade
        )

        thermal_h2_atten = float(np.exp(-0.022 * thermal_over_50 / (1.0 + 1.2 * betaine)))
        h1l_cf = base_h1l * (0.65 + 0.35 * intact) * thermal_h2_atten
        h1h_cf = base_h1h * (0.75 + 0.25 * intact) * thermal_h2_atten
        cox_cf = base_cox * (0.70 + 0.30 * intact) * thermal_h2_atten

        raw_leaf_state = {
            "actinomycetota_rel_abund": float(site_row.get("actinomycetota_rel_abund", 0.48))
            * (0.85 + 0.15 * intact),
            "thermoplasmatota_rel_abund": therm_cf,
            "halobacteriota_rel_abund": float(site_row.get("halobacteriota_rel_abund", 0.02)),
            "hydB_1l_ge": h1l_cf,
            "hydB_1h_ge": h1h_cf,
            "hydB_2a_ge": base_h2a,
            "coxL_formI_ge": cox_cf,
            "h2_oxidation_power_uw_m3": 29.8 * h1l_cf + 8.5 * (h1h_cf + base_h2a),
            "co_oxidation_power_uw_m3": 1.95 * cox_cf,
            "viable_dna_fraction": float(np.clip(0.45 + 0.35 * bacillus_60c_viability * intact, 0.18, 0.92)),
            "shannon_16s": max(2.2, shannon_cf),
        }
        projected_state = self.causal_engine.evaluate_constitutive_layer(raw_leaf_state)

        amb_interval = self.causal_engine.compute_version_ambiguity_interval(
            target_parent_total=betaine,
            leaf_sensitivities={
                "glycine_betaine_direct_OpuA": 0.24,
                "choline_precursor_GbsAB_oxidation": 0.17,
            },
        )

        advisory = (
            "REFUGIUM_PROTECTED: Intact 5-10 cm horizon + betaine priming sustains Thermoplasmatota 3-HP/4-HB carbon fixation and 60 C Bacillus survival."
            if intact >= 0.8 and bacillus_60c_viability >= 0.65
            else "WARNING_REFUGIUM_EXPOSURE: Topsoil scraping (<5-10 cm) or unmitigated 60 C thermal stress depletes subsurface Thermoplasmatota and Nitrospira."
        )

        summary = {
            "intervention": "refugium_thermal",
            "site_iri": site_iri,
            "subsurface_refugium_intact_fraction": round(intact, 3),
            "surface_temp_c": round(t_surf, 2),
            "effective_horizon_temp_c": round(effective_temp_c, 2),
            "exogenous_choline_betaine_mm": round(betaine, 3),
            "thermoplasmatota_rel_abund": round(projected_state["thermoplasmatota_rel_abund"], 4),
            "bacillus_60c_survival_index": round(bacillus_60c_viability, 4),
            "archaeal_3hp4hb_c_fixation_flux": fba_res["archaeal_3hp4hb_flux"],
            "betaine_osmoprotectant_flux": fba_res["betaine_osmoprotectant_flux"],
            "counterfactual_shannon_16s": round(projected_state["shannon_16s"], 4),
            "ocs_score": projected_state["ocs_score"],
            "management_advisory": advisory,
        }

        fdo_crate = FDOProvenanceRecorder.build_run_crate(
            intervention_name="refugium_thermal",
            site_iri=site_iri,
            parameters={
                "subsurface_refugium_intact_fraction": intact,
                "surface_temp_c": t_surf,
                "shade_shield_fraction": shade,
                "exogenous_choline_betaine_mm": betaine,
            },
            outputs=summary,
            ocs_score=projected_state["ocs_score"],
            version_ambiguity_interval=amb_interval,
        )

        return {
            "summary": summary,
            "projected_state": projected_state,
            "fba_fluxes": fba_res,
            "fdo_provenance": fdo_crate["fdo_kernel_metadata"],
        }
