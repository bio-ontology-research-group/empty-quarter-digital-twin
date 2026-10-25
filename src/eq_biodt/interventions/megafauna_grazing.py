"""Intervention 5: Native Megafauna Reintroduction (Arabian Oryx & Sand Gazelle) vs Camel Grazing.

Grounded in:
- Canonical sites `Site Camel1` (`https://rubalkhali.science/kb/RAK_1000064`, 19.65401 N, 47.68826 E),
  `31urine` (20.294623 N, 49.891405 E), and verified ENVO term `ENVO_00002003` (`fecal material`).
- `~/Documents/papers/empty-quarter-models`: Relieving the 1.0% soil `nifH` nitrogen-fixation
  bottleneck via localized ungulate dung/urine urea/NH4+ and phosphorus pulses to fuel
  Nitrososphaerota/Thermoplasmatota ammonia oxidation, Nitrospira nitrification, and
  archaeal 3-HP/4-HB carbon fixation.
- `~/Documents/papers/empty-quarter-lifeonair`: Hoof compaction of topsoil porosity eps_total
  under high camel stocking throttles Millington-Quirk H2/CO diffusion.
- `~/Documents/papers/empty-quarter-mobilome`: Resistome/ARG safety (`06_resistome_safety`)
  comparing wild native ungulates vs high-density livestock manure.
"""

from __future__ import annotations

from typing import Any, Dict, Optional

import numpy as np

from eq_biodt.models.causal_ontological_twin import CausalOntologicalTwinEngine
from eq_biodt.models.metabolic_community_engine import (
    CommunityMetabolicInputs,
    MetabolicCommunityEngine,
)
from eq_biodt.models.moisture_pulse_diffusion import MoisturePulseDiffusionModel
from eq_biodt.models.spatial_jsdm_virome import SpatialJSDMViromeModel
from eq_biodt.standards.fdo_provenance import FDOProvenanceRecorder


class MegafaunaGrazingIntervention:
    """Counterfactual simulator for Arabian Oryx / Sand Gazelle reintroduction vs Camel grazing."""

    SPECIES_TAXON_IDS = {
        "Oryx leucoryx": "http://purl.obolibrary.org/obo/NCBITaxon_39417",
        "Gazella marica": "http://purl.obolibrary.org/obo/NCBITaxon_1095692",
        "Camelus dromedarius": "http://purl.obolibrary.org/obo/NCBITaxon_9838",
    }

    def __init__(
        self,
        causal_engine: Optional[CausalOntologicalTwinEngine] = None,
        metabolic_engine: Optional[MetabolicCommunityEngine] = None,
        diffusion_model: Optional[MoisturePulseDiffusionModel] = None,
        virome_model: Optional[SpatialJSDMViromeModel] = None,
    ) -> None:
        self.causal_engine = causal_engine or CausalOntologicalTwinEngine()
        self.metabolic_engine = metabolic_engine or MetabolicCommunityEngine()
        self.diffusion_model = diffusion_model or MoisturePulseDiffusionModel()
        self.virome_model = virome_model or SpatialJSDMViromeModel()

    def simulate(
        self,
        site_row: Dict[str, Any],
        oryx_head_per_100km2: float = 3.0,
        gazelle_head_per_100km2: float = 5.0,
        camel_head_per_100km2: float = 1.5,
        rotational_water_point_management: bool = True,
    ) -> Dict[str, Any]:
        """Simulate wildlife reintroduction (Oryx leucoryx, Gazella marica) and camel grazing."""
        oryx = float(np.clip(oryx_head_per_100km2, 0.0, 40.0))
        gazelle = float(np.clip(gazelle_head_per_100km2, 0.0, 60.0))
        camel = float(np.clip(camel_head_per_100km2, 0.0, 80.0))
        rot_factor = 0.55 if rotational_water_point_management else 1.0

        base_shannon = float(site_row.get("shannon_16s", 4.85))
        base_h1l = float(site_row.get("hydB_1l_ge", 0.61))
        base_h1h = float(site_row.get("hydB_1h_ge", 0.12))
        base_h2a = float(site_row.get("hydB_2a_ge", 0.03))
        base_cox = float(site_row.get("coxL_formI_ge", 1.15))
        map_mm = float(site_row.get("map_mm", 28.0))
        evap_idx = float(site_row.get("xrf_evaporite_index", 0.0))
        site_iri = str(site_row.get("site_iri", "https://rubalkhali.science/kb/RAK_1000064"))

        # Standardized Metabolic Equivalents (1 Camel ~ 3.5 Oryx ~ 10 Sand Gazelles)
        wild_ungulate_units = oryx + 0.30 * gazelle
        camel_effective_units = 2.8 * camel * rot_factor

        # 1. Endozoochorous seed dispersal vs overbrowsing canopy loss
        # Wild ungulates at moderate density (1-8 units) disperse Haloxylon/Calligonum seeds;
        # heavy camel browsing strips perennial shrubs.
        seed_dispersal_gain = 0.14 * (wild_ungulate_units / (1.0 + 0.12 * wild_ungulate_units)) * np.exp(
            -0.04 * max(0.0, wild_ungulate_units - 8.0)
        )
        overbrowsing_loss = 0.22 * (camel_effective_units / (8.0 + camel_effective_units))
        net_shrub_cover_frac = float(np.clip(0.08 + seed_dispersal_gain - overbrowsing_loss, 0.01, 0.45))

        # 2. Hoof compaction of 0-10 cm aeolian sand crust (reduces eps_total from 0.40)
        compaction_delta = 0.0035 * camel_effective_units + 0.0008 * wild_ungulate_units
        effective_porosity = float(np.clip(0.40 - compaction_delta, 0.22, 0.40))

        # 3. Dung & urine N/P + choline/betaine pulse (bypasses the 1.0% soil nifH bottleneck)
        bioavailable_n_input = 0.08 + 0.022 * wild_ungulate_units + 0.028 * min(12.0, camel_effective_units)
        dung_choline_betaine = 0.05 * np.log1p(wild_ungulate_units + 0.6 * camel_effective_units)

        gas_res = self.diffusion_model.milling_quirk_gas_power(
            theta_vwc=0.020,
            total_porosity=effective_porosity,
            hydB_1l_ge=base_h1l * (1.0 - 0.35 * net_shrub_cover_frac),
            coxL_formI_ge=base_cox,
        )

        fba_res = self.metabolic_engine.solve_community_fba(
            CommunityMetabolicInputs(
                h2_uptake_bound=1.8 * gas_res["relative_gas_diffusivity"],
                co_uptake_bound=0.35 * gas_res["relative_gas_diffusivity"],
                root_exudate_doc_bound=0.95 * net_shrub_cover_frac + 0.04 * np.log1p(wild_ungulate_units),
                exogenous_nh4_no3_bound=bioavailable_n_input,
                choline_betaine_bound=dung_choline_betaine,
            )
        )

        # Nitrifying & 3-HP/4-HB archaea (Thermoplasmatota / Nitrososphaerota) expand with urea/NH4+ pulses
        # unless soil is heavily compacted by overstocking
        n_stimulation = 1.0 + 0.38 * np.log1p(wild_ungulate_units) - 0.018 * max(0.0, camel_effective_units - 8.0)
        therm_cf = float(
            np.clip(float(site_row.get("thermoplasmatota_rel_abund", 0.045)) * n_stimulation, 0.005, 0.22)
        )

        # Shannon diversity: hump-shaped intermediate disturbance & nutrient enrichment curve
        wild_bonus = 0.26 * (wild_ungulate_units / (2.5 + wild_ungulate_units)) * np.exp(
            -0.025 * max(0.0, wild_ungulate_units - 10.0)
        )
        overgrazing_penalty = 0.028 * max(0.0, camel_effective_units - 5.0)
        shannon_cf = max(2.4, base_shannon + wild_bonus - overgrazing_penalty)

        # Resistome / ARG load (`06_resistome_safety`): wild native Oryx/Gazelle have negligible
        # veterinary antibiotic exposure, whereas high-density camel corrals increase fecal ARG load
        manure_arg_driver = float(np.clip(0.015 * wild_ungulate_units + 0.065 * camel_effective_units, 0.0, 1.5))
        virome_res = self.virome_model.predict_virome_mobilome_state(
            shannon_16s=shannon_cf,
            compartment=str(site_row.get("compartment_norm", "Surface")),
            map_mm=map_mm,
            xrf_evaporite_index=evap_idx,
            vegetation_cover_frac=net_shrub_cover_frac,
            exogenous_manure_load=manure_arg_driver,
            baseline_viral_richness=float(site_row.get("viral_richness_mean", 18.0)),
            baseline_plasmid_richness=float(site_row.get("plasmid_richness_mean", 31.0)),
            baseline_shannon_16s=base_shannon,
        )

        raw_leaf_state = {
            "actinomycetota_rel_abund": float(site_row.get("actinomycetota_rel_abund", 0.48))
            * (1.0 - 0.06 * net_shrub_cover_frac),
            "thermoplasmatota_rel_abund": therm_cf,
            "halobacteriota_rel_abund": float(site_row.get("halobacteriota_rel_abund", 0.02))
            * (1.0 + 0.02 * camel_effective_units),
            "hydB_1l_ge": base_h1l * gas_res["relative_gas_diffusivity"] * (1.0 - 0.25 * net_shrub_cover_frac),
            "hydB_1h_ge": base_h1h * (1.0 + 0.18 * net_shrub_cover_frac),
            "hydB_2a_ge": base_h2a * (1.0 + 0.22 * net_shrub_cover_frac),
            "coxL_formI_ge": base_cox * gas_res["relative_gas_diffusivity"],
            "h2_oxidation_power_uw_m3": gas_res["h2_oxidation_power_uw_m3"],
            "co_oxidation_power_uw_m3": gas_res["co_oxidation_power_uw_m3"],
            "viable_dna_fraction": float(np.clip(0.56 + 0.14 * net_shrub_cover_frac, 0.25, 0.90)),
            "shannon_16s": shannon_cf,
        }
        projected_state = self.causal_engine.evaluate_constitutive_layer(raw_leaf_state)

        amb_interval = self.causal_engine.compute_version_ambiguity_interval(
            target_parent_total=oryx + gazelle + camel,
            leaf_sensitivities={
                "Oryx_leucoryx_nomadic_reintroduction": 0.035,
                "Gazella_marica_seed_disperser": 0.028,
                "Camelus_dromedarius_stationary_grazing": -0.024,
            },
        )

        advisory = (
            "OPTIMAL_ORYX_GAZELLE_REWILDING: Native ungulate N/P pulses bypass the 1% soil nifH bottleneck, boost nitrification & seed dispersal, and pass resistome safety."
            if (oryx + gazelle >= 2.0 and camel_effective_units <= 6.0 and virome_res["resistome_biosafety_pass"])
            else "WARNING_OVERGRAZING_OR_ARG_LOAD: High livestock stocking compacts topsoil porosity, suppresses H2 diffusion, and elevates fecal resistome risk."
        )

        summary = {
            "intervention": "megafauna_grazing",
            "site_iri": site_iri,
            "oryx_head_per_100km2": round(oryx, 2),
            "gazelle_head_per_100km2": round(gazelle, 2),
            "camel_head_per_100km2": round(camel, 2),
            "rotational_water_point_management": rotational_water_point_management,
            "net_shrub_cover_fraction": round(net_shrub_cover_frac, 4),
            "effective_soil_porosity": round(effective_porosity, 4),
            "relative_gas_diffusivity": gas_res["relative_gas_diffusivity"],
            "nitrification_flux": fba_res["nitrification_flux"],
            "archaeal_3hp4hb_c_fixation_flux": fba_res["archaeal_3hp4hb_flux"],
            "counterfactual_shannon_16s": round(projected_state["shannon_16s"], 4),
            "resistome_arg_risk_index": virome_res["resistome_arg_risk_index"],
            "resistome_biosafety_pass": virome_res["resistome_biosafety_pass"],
            "ocs_score": projected_state["ocs_score"],
            "management_advisory": advisory,
        }

        fdo_crate = FDOProvenanceRecorder.build_run_crate(
            intervention_name="megafauna_grazing",
            site_iri=site_iri,
            parameters={
                "oryx_head_per_100km2": oryx,
                "gazelle_head_per_100km2": gazelle,
                "camel_head_per_100km2": camel,
                "rotational_water_point_management": rotational_water_point_management,
            },
            outputs=summary,
            ocs_score=projected_state["ocs_score"],
            version_ambiguity_interval=amb_interval,
        )

        return {
            "summary": summary,
            "projected_state": projected_state,
            "fba_fluxes": fba_res,
            "virome_mobilome": virome_res,
            "fdo_provenance": fdo_crate["fdo_kernel_metadata"],
        }
