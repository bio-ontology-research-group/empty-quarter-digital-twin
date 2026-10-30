"""Intervention 2: Native Shrub Revegetation & RAK-PR / Site 60 Bio-Inoculation Simulator.

Grounded in:
- `~/Documents/papers/empty-quarter-isolate`: Multi-functional desert PGPR isolate
  Peribacillus frigoritolerans RAK-PR (accS, phoD, IAA, siderophore, proline/trehalose).
- `~/Documents/papers/empty-quarter-saline`: Site 60 Allorhizobium 541 kb rbcL-hyaB
  autotrophic plasmid + diazotrophic Hydrogenophaga/Rubrivivax consortium.
- `~/Documents/papers/empty-quarter-lifeonair`: The Rhizosphere Paradox — 10x microbial
  biomass expansion in root zones vs 79% reduction in group 1l [NiFe] H2 scavengers,
  requiring mosaic bare-crust preservation.
- `~/Documents/papers/empty-quarter-viruses`: Plant-associated Caudoviricetes expansion.
"""

from __future__ import annotations

from typing import Any, Dict, Optional

import numpy as np

from eq_biodt.models.carbon_budget_engine import CarbonBudgetEngine
from eq_biodt.models.causal_ontological_twin import CausalOntologicalTwinEngine
from eq_biodt.models.metabolic_community_engine import (
    CommunityMetabolicInputs,
    MetabolicCommunityEngine,
)
from eq_biodt.models.spatial_jsdm_virome import SpatialJSDMViromeModel
from eq_biodt.standards.fdo_provenance import FDOProvenanceRecorder


class RevegetationIntervention:
    """Counterfactual simulator for native shrub planting and PGPR/plasmid inoculation."""

    NATIVE_SHRUB_SPECIES = (
        "Haloxylon persicum",
        "Calligonum crinitum",
        "Zygophyllum mandavillei",
        "Cornulaca arabica",
    )

    def __init__(
        self,
        causal_engine: Optional[CausalOntologicalTwinEngine] = None,
        metabolic_engine: Optional[MetabolicCommunityEngine] = None,
        virome_model: Optional[SpatialJSDMViromeModel] = None,
        carbon_engine: Optional[CarbonBudgetEngine] = None,
    ) -> None:
        self.causal_engine = causal_engine or CausalOntologicalTwinEngine()
        self.metabolic_engine = metabolic_engine or MetabolicCommunityEngine()
        self.virome_model = virome_model or SpatialJSDMViromeModel()
        self.carbon_engine = carbon_engine or CarbonBudgetEngine()

    def simulate(
        self,
        site_row: Dict[str, Any],
        shrub_cover_fraction: float = 0.25,
        rak_pr_pgpr_dose: float = 0.8,
        site60_rbcl_hyab_plasmid_dose: float = 0.6,
        shrub_species: str = "Haloxylon persicum",
    ) -> Dict[str, Any]:
        """Run a revegetation + bio-inoculation counterfactual at a site/compartment."""
        cover = float(np.clip(shrub_cover_fraction, 0.0, 0.85))
        pgpr = float(np.clip(rak_pr_pgpr_dose, 0.0, 1.0))
        plasmid = float(np.clip(site60_rbcl_hyab_plasmid_dose, 0.0, 1.0))

        base_shannon = float(site_row.get("shannon_16s", 4.85))
        base_h1l = float(site_row.get("hydB_1l_ge", 0.61))
        base_h1h = float(site_row.get("hydB_1h_ge", 0.12))
        base_h2a = float(site_row.get("hydB_2a_ge", 0.03))
        base_cox = float(site_row.get("coxL_formI_ge", 1.15))
        base_dna = float(site_row.get("dna_conc_mean", 0.55))
        map_mm = float(site_row.get("map_mm", 28.0))
        evap_idx = float(site_row.get("xrf_evaporite_index", 0.0))
        site_iri = str(site_row.get("site_iri", "https://rubalkhali.science/kb/RAK_1000001"))

        # Rhizosphere paradox: inside root zones, biomass expands ~10.3x (log2 ratio = 3.37) while
        # group 1l [NiFe] H2 scavengers decline by 79.3% (log2FC = -2.2747) due to root-exudate competition.
        # Mosaic planting weights root-zone (cover) and bare inter-shrub crust (1 - cover).
        root_zone_h1l = base_h1l * (1.0 - 0.793 * (0.95 + 0.05 * pgpr))
        bare_zone_h1l = base_h1l
        landscape_h1l = cover * root_zone_h1l + (1.0 - cover) * bare_zone_h1l

        # Site 60 rbcL-hyaB plasmid rescues aerobic/microaerophilic H2 uptake (group 1h/2a) & CBB fixation
        landscape_h1h = base_h1h * (1.0 + 0.55 * cover * (0.4 + 0.6 * plasmid))
        landscape_h2a = base_h2a * (1.0 + 1.20 * cover * plasmid)
        landscape_cox = base_cox * (1.0 - 0.35 * cover + 0.25 * cover * plasmid)

        # Biomass multiplier (1x in bare soil, ~10.34x in wild root zone, up to 12.5x with PGPR+plasmid)
        root_biomass_mult = 1.0 + 9.34 * (0.88 + 0.24 * pgpr + 0.12 * plasmid)
        landscape_biomass_mult = cover * root_biomass_mult + (1.0 - cover) * 1.0
        counterfactual_dna_conc = base_dna * landscape_biomass_mult
        biomass_weighted_h1l = (
            cover * root_biomass_mult * root_zone_h1l + (1.0 - cover) * 1.0 * bare_zone_h1l
        ) / max(1e-9, landscape_biomass_mult)
        absolute_h1l_density_mult = (
            cover * root_biomass_mult * (root_zone_h1l / max(1e-6, base_h1l)) + (1.0 - cover) * 1.0
        )

        # Solve HiGHS FBA for the vegetated + inoculated community
        fba_res = self.metabolic_engine.solve_community_fba(
            CommunityMetabolicInputs(
                h2_uptake_bound=1.8 * (landscape_h1l / max(0.1, base_h1l)),
                co_uptake_bound=0.35 * (landscape_cox / max(0.1, base_cox)),
                root_exudate_doc_bound=1.45 * cover * (1.0 + 0.35 * pgpr),
                exogenous_nh4_no3_bound=0.08 + 0.12 * cover * pgpr,
                sulfate_so4_bound=0.25,
                salinity_ms_cm=max(0.8, 1.8 + 1.5 * max(0.0, evap_idx)),
                rbcl_hyab_plasmid_inoc=plasmid * cover,
            )
        )

        # Predict Shannon diversity & virome/mobilome co-turnover anchored to site baseline
        shannon_cf = base_shannon + 0.38 * cover * (0.6 + 0.4 * pgpr) - 0.25 * max(0.0, cover - 0.55)
        virome_res = self.virome_model.predict_virome_mobilome_state(
            shannon_16s=shannon_cf,
            compartment="RootAdjacent" if cover > 0.15 else str(site_row.get("compartment_norm", "Surface")),
            map_mm=map_mm,
            xrf_evaporite_index=evap_idx,
            vegetation_cover_frac=cover,
            plasmid_bioaug_dose=0.5 * (pgpr + plasmid),
            baseline_viral_richness=float(site_row.get("viral_richness_mean", 18.0)),
            baseline_plasmid_richness=float(site_row.get("plasmid_richness_mean", 31.0)),
            baseline_shannon_16s=base_shannon,
        )

        raw_leaf_state = {
            "actinomycetota_rel_abund": float(site_row.get("actinomycetota_rel_abund", 0.48))
            * (1.0 - 0.12 * cover),
            "thermoplasmatota_rel_abund": float(site_row.get("thermoplasmatota_rel_abund", 0.04))
            * (1.0 + 0.25 * cover),
            "halobacteriota_rel_abund": float(site_row.get("halobacteriota_rel_abund", 0.02))
            * (1.0 - 0.15 * cover * pgpr),
            "hydB_1l_ge": landscape_h1l,
            "hydB_1h_ge": landscape_h1h,
            "hydB_2a_ge": landscape_h2a,
            "coxL_formI_ge": landscape_cox,
            "h2_oxidation_power_uw_m3": 29.8 * landscape_h1l + 8.5 * (landscape_h1h + landscape_h2a),
            "co_oxidation_power_uw_m3": 1.95 * landscape_cox,
            "viable_dna_fraction": float(np.clip(0.56 + 0.28 * cover * (0.6 + 0.4 * pgpr), 0.2, 0.95)),
            "shannon_16s": shannon_cf,
        }
        projected_state = self.causal_engine.evaluate_constitutive_layer(raw_leaf_state)

        amb_interval = self.causal_engine.compute_version_ambiguity_interval(
            target_parent_total=cover,
            leaf_sensitivities={
                "Haloxylon_persicum_deep_taproot": 0.42,
                "Calligonum_crinitum_dune_stabilizer": 0.36,
                "Zygophyllum_mandavillei_halophyte": 0.21,
            },
        )

        advisory = (
            "OPTIMAL_MOSAIC_REVEGETATION: Balances 10x root-zone biomass & viral expansion with inter-shrub atmospheric H2 scavenging."
            if 0.15 <= cover <= 0.45
            else (
                "WARNING_CRUST_LOSS: Shrub cover >45% excessively displaces obligate group 1l H2 scavengers in bare interdune crusts."
                if cover > 0.45
                else "LOW_CANOPY_DENSITY: Modest rhizosphere carbon fixation gain."
            )
        )

        cf_site_row = dict(site_row)
        cf_site_row["hydB_1l_ge"] = projected_state["hydB_1l_ge"]
        cf_site_row["hydB_1h_ge"] = projected_state["hydB_1h_ge"]
        cf_site_row["hydB_2a_ge"] = projected_state["hydB_2a_ge"]
        cf_site_row["coxL_formI_ge"] = projected_state["coxL_formI_ge"]
        carbon_budget = self.carbon_engine.compute_state_carbon_budget(
            cf_site_row,
            shrub_cover_fraction=cover,
            rak_pr_pgpr_dose=pgpr,
            rbcl_hyab_plasmid_dose=plasmid,
            hard_maintenance_feasible=bool(fba_res.get("hard_maintenance_feasible", True)),
            cbb_flux_multiplier=1.0 + 0.65 * plasmid * cover,
            biomass_density_multiplier=landscape_biomass_mult,
        )

        summary = {
            "intervention": "revegetation",
            "site_iri": site_iri,
            "shrub_species": shrub_species,
            "shrub_cover_fraction": round(cover, 3),
            "rak_pr_pgpr_dose": round(pgpr, 3),
            "site60_rbcl_hyab_plasmid_dose": round(plasmid, 3),
            "biomass_dna_multiplier": round(landscape_biomass_mult, 3),
            "counterfactual_dna_conc_ng_ul": round(counterfactual_dna_conc, 4),
            "baseline_hydB_1l_ge": round(base_h1l, 4),
            "root_zone_hydB_1l_ge": round(root_zone_h1l, 4),
            "landscape_hydB_1l_ge": round(projected_state["hydB_1l_ge"], 4),
            "biomass_weighted_hydB_1l_ge": round(biomass_weighted_h1l, 4),
            "absolute_hydB_1l_density_multiplier": round(absolute_h1l_density_mult, 4),
            "total_h2ase_ge": round(projected_state["total_h2ase_ge"], 4),
            "counterfactual_shannon_16s": round(projected_state["shannon_16s"], 4),
            "autotrophic_carbon_fixation_flux": fba_res["autotrophic_carbon_fixation_flux"],
            "nifH_n2_fixation_flux": fba_res["nifH_n2_fixation_flux"],
            "total_direct_co2_bound_t_km2_yr": carbon_budget["total_direct_co2_bound_t_km2_yr"]["point"],
            "total_co2eq_bound_t_km2_yr": carbon_budget["total_co2eq_bound_t_km2_yr"]["point"],
            "viral_richness_votu": virome_res["viral_richness_votu"],
            "plasmid_richness_potu": virome_res["plasmid_richness_potu"],
            "resistome_biosafety_pass": virome_res["resistome_biosafety_pass"],
            "ocs_score": projected_state["ocs_score"],
            "management_advisory": advisory,
        }

        fdo_crate = FDOProvenanceRecorder.build_run_crate(
            intervention_name="revegetation",
            site_iri=site_iri,
            parameters={
                "shrub_species": shrub_species,
                "shrub_cover_fraction": cover,
                "rak_pr_pgpr_dose": pgpr,
                "site60_rbcl_hyab_plasmid_dose": plasmid,
            },
            outputs=summary,
            ocs_score=projected_state["ocs_score"],
            version_ambiguity_interval=amb_interval,
        )

        return {
            "summary": summary,
            "projected_state": projected_state,
            "fba_fluxes": fba_res,
            "carbon_budget": carbon_budget,
            "virome_mobilome": virome_res,
            "fdo_provenance": fdo_crate["fdo_kernel_metadata"],
        }
