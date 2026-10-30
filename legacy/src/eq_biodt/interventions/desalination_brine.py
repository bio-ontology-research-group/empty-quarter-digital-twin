"""Intervention 3: Desalination Brine Management vs Brackish Oasis Irrigation Simulator.

Grounded in the Umm as-Samim paired natural experiment and Table S2 measured water chemistry:
- `~/Documents/papers/empty-quarter-tiannyu` & `empty-quarter-saline`:
  - Site 60 (`RAK_1000060`, `60-MR-WR1`): `EC = 10.65 mS/cm`, `pH = 7.68`, `Ca2+ = 1,066.4 mg/L`,
    `SO4^2- = 2,229.0 mg/L`, `Na+ = 1,548.9 mg/L`, `Cl- = 2,221.4 mg/L` ("Ca-SO4 sulfur flux"
    brackish artesian pool supporting the 541 kb Allorhizobium rbcL-hyaB plasmid, diazotrophic
    Hydrogenophaga/Rubrivivax, and high viral/functional diversity).
  - Site 59 (`RAK_1000059`, `59-MR-WR1` & `59-MR-WR2`): `EC = 216.2-231.4 mS/cm`, `pH = 7.92-8.01`,
    `Na+ = 71,044.9-79,896.2 mg/L`, `Cl- = 103,722-121,069 mg/L`, `Ca2+ = 1,241.7-1,371.9 mg/L`,
    `SO4^2- = 11,074.6-11,923.9 mg/L` ("Na-Cl hypersaline sulfur sink" regime causing hard
    maintenance ATP deficit and collapse into Halobacteriota monoculture).
  - Sphingomonas C-27.9 339 kb ectABC-kdpFABC-proP osmoregulation megaplasmid.
"""

from __future__ import annotations

from typing import Any, Dict, List, Optional

import numpy as np

from eq_biodt.models.causal_ontological_twin import CausalOntologicalTwinEngine
from eq_biodt.models.metabolic_community_engine import (
    CommunityMetabolicInputs,
    MetabolicCommunityEngine,
)
from eq_biodt.models.moisture_pulse_diffusion import MoisturePulseDiffusionModel
from eq_biodt.standards.fdo_provenance import FDOProvenanceRecorder


class DesalinationBrineIntervention:
    """Counterfactual simulator for desalination permeate irrigation vs brine discharge."""

    def __init__(
        self,
        causal_engine: Optional[CausalOntologicalTwinEngine] = None,
        metabolic_engine: Optional[MetabolicCommunityEngine] = None,
        gas_model: Optional[MoisturePulseDiffusionModel] = None,
    ) -> None:
        self.causal_engine = causal_engine or CausalOntologicalTwinEngine()
        self.metabolic_engine = metabolic_engine or MetabolicCommunityEngine()
        self.gas_model = gas_model or MoisturePulseDiffusionModel()

    def simulate(
        self,
        site_row: Dict[str, Any],
        discharge_salinity_ms_cm: float = 10.65,
        ca_so4_to_nacl_ratio: float = 1.4,
        sphingomonas_ect_kdp_plasmid_dose: float = 0.5,
        allorhizobium_rbcl_hyab_dose: float = 0.5,
    ) -> Dict[str, Any]:
        """Simulate brine/permeate discharge with specified EC (mS/cm) and Ca-SO4:Na-Cl ratio."""
        ec = float(np.clip(discharge_salinity_ms_cm, 0.5, 250.0))
        so4_ratio = float(np.clip(ca_so4_to_nacl_ratio, 0.05, 3.0))
        ect_dose = float(np.clip(sphingomonas_ect_kdp_plasmid_dose, 0.0, 1.0))
        rbcl_dose = float(np.clip(allorhizobium_rbcl_hyab_dose, 0.0, 1.0))

        base_shannon = float(site_row.get("shannon_16s", 4.85))
        base_h1l = float(site_row.get("hydB_1l_ge", 0.61))
        base_h1h = float(site_row.get("hydB_1h_ge", 0.12))
        base_h2a = float(site_row.get("hydB_2a_ge", 0.03))
        base_cox = float(site_row.get("coxL_formI_ge", 1.15))
        site_iri = str(site_row.get("site_iri", "https://rubalkhali.science/kb/RAK_1000060"))

        # Divalent Ca-SO4 buffering lowers effective monovalent Na-Cl osmotic activity:
        effective_nacl_ec = ec / (1.0 + 0.55 * so4_ratio)
        is_sulfur_flux = (effective_nacl_ec <= 25.0) and (so4_ratio >= 0.6)
        sulfur_regime = (
            "CA_SO4_SULFUR_FLUX_OASIS" if is_sulfur_flux else "NACL_HYPERSALINE_SULFUR_SINK"
        )

        # Solve HiGHS FBA with hard maintenance feasibility for the aquatic/soil boundary community
        fba_res = self.metabolic_engine.solve_community_fba(
            CommunityMetabolicInputs(
                h2_uptake_bound=1.6 * np.exp(-max(0.0, effective_nacl_ec - 12.0) / 65.0),
                co_uptake_bound=0.30 * np.exp(-max(0.0, effective_nacl_ec - 12.0) / 75.0),
                root_exudate_doc_bound=0.04 if is_sulfur_flux else 0.01,
                exogenous_nh4_no3_bound=0.10,
                sulfate_so4_bound=0.45 * min(2.0, so4_ratio),
                salinity_ms_cm=effective_nacl_ec,
                rbcl_hyab_plasmid_inoc=rbcl_dose if effective_nacl_ec < 35.0 else 0.05 * rbcl_dose,
                ect_kdp_plasmid_inoc=ect_dose,
            )
        )

        # Halobacteriota vs Actinomycetota/Thermoplasmatota response to effective Na-Cl EC and ectABC rescue
        halo_bloom = float(
            np.clip(
                0.02
                + 0.76
                / (1.0 + np.exp(-(effective_nacl_ec - 55.0) / 20.0))
                * (1.0 - 0.25 * ect_dose),
                0.01,
                0.85,
            )
        )
        act_abund = float(
            np.clip(
                float(site_row.get("actinomycetota_rel_abund", 0.45))
                * np.exp(-max(0.0, effective_nacl_ec - 15.0) / (60.0 + 35.0 * ect_dose)),
                0.03,
                0.65,
            )
        )
        therm_abund = float(
            np.clip(
                float(site_row.get("thermoplasmatota_rel_abund", 0.04))
                * np.exp(-max(0.0, effective_nacl_ec - 12.0) / 38.0),
                0.002,
                0.20,
            )
        )

        # Shannon diversity: peaks in Ca-SO4 buffered brackish regime, declines under unbuffered Na-Cl
        brackish_oasis_bonus = (
            0.42 * np.exp(-((effective_nacl_ec - 10.65) ** 2) / 95.0) * min(1.5, so4_ratio)
        )
        brine_collapse = (
            2.15
            / (1.0 + np.exp(-(effective_nacl_ec - 60.0) / 22.0))
            * (1.0 - 0.32 * ect_dose)
        )
        shannon_cf = max(1.35, base_shannon + brackish_oasis_bonus - brine_collapse)

        h1l_cf = base_h1l * np.exp(-max(0.0, effective_nacl_ec - 10.0) / 42.0)
        h1h_cf = base_h1h * (
            1.0 + 0.85 * rbcl_dose * np.exp(-max(0.0, effective_nacl_ec - 18.0) / 35.0)
        )
        h2a_cf = base_h2a * (
            1.0 + 1.10 * rbcl_dose * np.exp(-max(0.0, effective_nacl_ec - 18.0) / 35.0)
        )
        cox_cf = base_cox * np.exp(-max(0.0, effective_nacl_ec - 12.0) / 50.0)

        # Unified Millington-Quirk gas diffusion power
        gas_metrics = self.gas_model.milling_quirk_gas_power(
            theta_vwc=0.045 if ec < 50.0 else 0.085,
            total_porosity=0.40 if ec < 80.0 else 0.32,
            hydB_1l_ge=h1l_cf,
            coxL_formI_ge=cox_cf,
        )

        raw_leaf_state = {
            "actinomycetota_rel_abund": act_abund,
            "thermoplasmatota_rel_abund": therm_abund,
            "halobacteriota_rel_abund": halo_bloom,
            "hydB_1l_ge": h1l_cf,
            "hydB_1h_ge": h1h_cf,
            "hydB_2a_ge": h2a_cf,
            "coxL_formI_ge": cox_cf,
            "h2_oxidation_power_uw_m3": gas_metrics["h2_oxidation_power_uw_m3"],
            "co_oxidation_power_uw_m3": gas_metrics["co_oxidation_power_uw_m3"],
            "viable_dna_fraction": float(
                np.clip(0.62 - 0.25 * (ec / 220.0) + 0.12 * ect_dose, 0.20, 0.90)
            ),
            "shannon_16s": shannon_cf,
        }
        projected_state = self.causal_engine.evaluate_constitutive_layer(raw_leaf_state)

        amb_interval = self.causal_engine.compute_version_ambiguity_interval(
            target_parent_total=ec,
            leaf_sensitivities={
                "CaSO4_gypsum_buffered_discharge": -0.003,
                "NaCl_unbuffered_ro_reject_brine": -0.014,
            },
        )

        advisory = (
            "SITE60_SULFUR_FLUX_OASIS: Brackish Ca-SO4 regime (~10.65 mS/cm) sustains Allorhizobium CBB carbon fixation, N2 fixation, and dark-proteome coding diversity."
            if is_sulfur_flux
            else "WARNING_SITE59_SULFUR_SINK_COLLAPSE: Unbuffered Na-Cl brine (>200 mS/cm) drives hard ATP maintenance deficit and Halobacteriota monoculture; dilute below 15 mS/cm and buffer with Ca-SO4."
        )

        summary = {
            "intervention": "desalination_brine",
            "site_iri": site_iri,
            "discharge_salinity_ms_cm": round(ec, 2),
            "effective_nacl_ec_ms_cm": round(effective_nacl_ec, 2),
            "ca_so4_to_nacl_ratio": round(so4_ratio, 3),
            "sulfur_geochemical_regime": sulfur_regime,
            "sphingomonas_ect_kdp_plasmid_dose": round(ect_dose, 3),
            "allorhizobium_rbcl_hyab_dose": round(rbcl_dose, 3),
            "counterfactual_shannon_16s": round(projected_state["shannon_16s"], 4),
            "halobacteriota_rel_abund": round(projected_state["halobacteriota_rel_abund"], 4),
            "hard_maintenance_feasible": fba_res["hard_maintenance_feasible"],
            "metabolic_viability_state": fba_res["metabolic_viability_state"],
            "autotrophic_carbon_fixation_flux": fba_res["autotrophic_carbon_fixation_flux"],
            "sulfur_redox_flux": fba_res["sulfur_redox_flux"],
            "atpm_maintenance_deficit": fba_res["atpm_maintenance_deficit"],
            "dark_proteome_coding_fraction": fba_res["dark_proteome_coding_fraction"],
            "dark_proteome_active_fraction": fba_res["dark_proteome_active_fraction"],
            "ocs_score": projected_state["ocs_score"],
            "management_advisory": advisory,
        }

        fdo_crate = FDOProvenanceRecorder.build_run_crate(
            intervention_name="desalination_brine",
            site_iri=site_iri,
            parameters={
                "discharge_salinity_ms_cm": ec,
                "ca_so4_to_nacl_ratio": so4_ratio,
                "sphingomonas_ect_kdp_plasmid_dose": ect_dose,
                "allorhizobium_rbcl_hyab_dose": rbcl_dose,
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

    def evaluate_factorial_brine_grid(
        self, site_row: Dict[str, Any]
    ) -> List[Dict[str, Any]]:
        """Evaluate a common-baseline 4 x 2 x 2 factorial grid separating EC, Ca-SO4:Na-Cl ratio, and plasmid dose."""
        ec_levels = (10.65, 35.0, 85.0, 216.2)
        ratio_levels = (0.10, 1.40)
        plasmid_levels = (0.0, 0.60)
        results: List[Dict[str, Any]] = []
        for ec in ec_levels:
            for ratio in ratio_levels:
                for p_dose in plasmid_levels:
                    res = self.simulate(
                        site_row,
                        discharge_salinity_ms_cm=ec,
                        ca_so4_to_nacl_ratio=ratio,
                        sphingomonas_ect_kdp_plasmid_dose=p_dose,
                        allorhizobium_rbcl_hyab_dose=p_dose,
                    )["summary"]
                    results.append(res)
        return results
