"""Dimensional Carbon & CO2-Equivalent Budget Engine for EQ-BioDT.

Quantifies current and counterfactual CO2 binding and indirect greenhouse-gas mitigation
across the Rub' al-Khali (650,000 km^2 total basin; 248 corridor grid cells = 178,560 km^2,
of which 128 cells = 92,160 km^2 lie inside the empirical domain of support) through
four mechanistically distinct pathways:

1. Direct Microbial Dark Carbon Fixation via RuBisCO (CBB cycle, K01601 rbcL + K00855 prkB,
   and archaeal 3-HP/4-HB cycle):
   - Coupled to state-specific `hydB_1l_ge`, `hydB_1h_ge`, `hydB_2a_ge`, `coxL_formI_ge`,
     empirical CBB co-carriage (`share_1l_cbb_mapped`, `share_cbb_mapped`), Millington-Quirk
     gas diffusivity `D_rel`, moisture gate `g(theta)`, and HiGHS FBA maintenance feasibility.
   - Grounded in `lifeonair_4_energetics.py`, Leung et al. (2025) global 25-desert synthesis
     (0.01-1.0 nmol 14CO2 gdw^-1 d^-1; 17-21% CBB-encoding cells), Bay et al. (2021 ISME J /
     Nat. Microbiol.), and Tribbia et al. (2026 AEM, Namib Desert).

2. Indirect Soil H2 Scavenging Greenhouse-Gas Sink (CO2-equivalent via GWP100 = 11 +/- 5):
   - High-affinity [NiFe]-hydrogenases oxidize atmospheric H2 (0.53 ppmv), preserving
     tropospheric OH radicals that oxidize CH4 and reducing stratospheric H2O formation
     (Ocko & Hamburg 2022 ACP; Warwick et al. 2023 ACP; Paulot et al. 2024 ACP).

3. Photosynthetic Shrub & Rhizosphere Soil Organic Carbon Sequestration (Haloxylon persicum,
   Calligonum crinitum, Zygophyllum mandavillei):
   - Scaled by fractional canopy cover `f_shrub` (inferred from 60-day MODIS NDVI at baseline
     or prescribed by revegetation scenarios), mean annual precipitation `MAP`, and
     RAK-PR PGPR + Site 60 `rbcL-hyaB` plasmid rhizosphere inoculation.
   - Calibrated against hyperarid/arid Haloxylon chronosequence and flux benchmarks
     (Su et al. 2007; Wang et al. 2021; Wohlfahrt et al. 2008).

4. Abiotic & Biogenic Pedogenic / Sabkha Carbonate Binding (CaCO3 / CaMg(CO3)2):
   - Coupled to XRF Ca + Mg content, alkaline pH (>7.2), root/microbial respiration HCO3^-
     supply, and divalent Ca-SO4 vs monovalent Na-Cl brine irrigation regime.
   - Calibrated against Schlesinger (1985, Mojave Desert caliche: 0.44-1.54 t CO2 km^-2 yr^-1)
     and Ma et al. (2014, Tarim/Taklamakan alkaline desert basin: up to 26 t CO2 km^-2 yr^-1).
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Dict, List, Optional

import numpy as np
import pandas as pd

from eq_biodt.models.metabolic_community_engine import (
    CommunityMetabolicInputs,
    MetabolicCommunityEngine,
)


@dataclass(frozen=True)
class CarbonBudgetParams:
    """Physical, stoichiometric, and regional scaling constants for the EQ carbon budget."""

    # Spatial denominators
    rub_al_khali_total_area_km2: float = 650000.0
    rub_al_khali_sand_sea_polygon_km2: float = 500000.0
    corridor_cell_area_km2: float = 720.0  # 0.25 deg x 0.25 deg cell at ~20.2 N (~27.8 x 25.9 km)
    active_soil_depth_m: float = 0.10
    top_biocrust_assay_depth_m: float = 0.02
    dry_bulk_density_g_m3: float = 1.6e6   # 1.6 g cm^-3 -> 1.6e5 g dry soil per m^2 in 0-10 cm column
    hyperarid_cells_per_gdw: float = 1.0e6 # lifeonair_4_energetics.py lower-envelope hyperarid cell density
    guild_cell_fraction: float = 0.184     # lifeonair_4_energetics.py measured MAG-direct guild share
    active_moisture_days_per_yr: float = 15.0  # Rain + dew resuscitation duty cycle (days/yr)

    # Atmospheric H2 & CO uptake reference envelope at reference GE (mg gas m^-2 yr^-1)
    # From lifeonair_4_energetics.py & Buzzard et al. (2022) / Cordero et al. (2019)
    h2_uptake_ref_mg_m2_yr_low: float = 5.0
    h2_uptake_ref_mg_m2_yr_point: float = 12.5
    h2_uptake_ref_mg_m2_yr_high: float = 25.0
    co_uptake_ref_mg_m2_yr_low: float = 20.0
    co_uptake_ref_mg_m2_yr_point: float = 50.0
    co_uptake_ref_mg_m2_yr_high: float = 100.0
    ref_h2_power_uw_m3: float = 28.0
    ref_co_power_uw_m3: float = 1.95

    # Stoichiometry & GWP100 (Warwick et al. 2023 ACP: 12 +/- 6; Paulot et al. 2024 ACP: 11.6)
    # 1 mol CO2 (44.01 g) via CBB requires ~8 mol H2 (16.13 g H2) or ~8 mol CO (224.08 g CO)
    mol_h2_per_mol_co2_cbb: float = 8.0
    anabolic_electron_fraction_low: float = 0.05
    anabolic_electron_fraction_point: float = 0.12
    anabolic_electron_fraction_high: float = 0.25
    h2_gwp100_low: float = 6.0
    h2_gwp100_point: float = 11.6
    h2_gwp100_high: float = 16.0

    # Vegetation NPP reference for full-canopy Haloxylon shrubland (t CO2 km^-2 yr^-1 at f=1.0, MAP=25 mm)
    veg_full_canopy_t_co2_km2_yr_low: float = 190.0
    veg_full_canopy_t_co2_km2_yr_point: float = 295.0
    veg_full_canopy_t_co2_km2_yr_high: float = 420.0

    # Pedogenic carbonate baseline (Schlesinger 1985 Mojave caliche: 0.44 - 1.54 t CO2 km^-2 yr^-1)
    sic_bare_ref_t_co2_km2_yr_low: float = 0.44
    sic_bare_ref_t_co2_km2_yr_point: float = 0.88
    sic_bare_ref_t_co2_km2_yr_high: float = 1.54

    # Persistent stock retention fractions (converting gross fixation/precipitation to net retained C stock)
    soc_retention_fraction_low: float = 0.10
    soc_retention_fraction_point: float = 0.15
    soc_retention_fraction_high: float = 0.25
    non_carbonate_ca_fraction_low: float = 0.25
    non_carbonate_ca_fraction_point: float = 0.50
    non_carbonate_ca_fraction_high: float = 0.75
    microbial_necromass_retention_low: float = 0.15
    microbial_necromass_retention_point: float = 0.25
    microbial_necromass_retention_high: float = 0.40


class CarbonBudgetEngine:
    """Computes dimensional CO2 and CO2-equivalent budgets across states, scenarios, and deserts."""

    # Published cross-desert empirical benchmarks for validation
    CROSS_DESERT_BENCHMARKS: List[Dict[str, Any]] = [
        {
            "desert": "Global 25-Desert Synthesis (Atacama, Namib, Gobi, Negev, Mojave, Kalahari, Antarctic)",
            "citation": "Leung et al. 2025 (bioRxiv 10.1101/2025.11.29.691269)",
            "metric": "cbb_rubisco_community_share",
            "unit": "fraction_of_cells_or_mapped_reads",
            "literature_low": 0.17,
            "literature_point": 0.21,
            "literature_high": 0.28,
        },
        {
            "desert": "Global 25-Desert Synthesis (Atacama, Namib, Gobi, Negev, Mojave, Antarctic)",
            "citation": "Leung et al. 2025 (bioRxiv 10.1101/2025.11.29.691269)",
            "metric": "dark_and_h2_stimulated_co2_fixation_nmol_gdw_d",
            "unit": "nmol_CO2_gdw^-1_d^-1",
            "literature_low": 0.01,
            "literature_point": 0.15,
            "literature_high": 1.00,
        },
        {
            "desert": "Namib Desert Hyperarid & Fog Transect (C14-2 to C14-18)",
            "citation": "Tribbia et al. 2026 (Appl. Environ. Microbiol. 10.1128/aem.00265-26)",
            "metric": "cell_specific_h2_oxidation_nmol_h_cell",
            "unit": "nmol_H2_h^-1_cell^-1",
            "literature_low": 1.36e-8,
            "literature_point": 2.82e-8,
            "literature_high": 5.31e-8,
        },
        {
            "desert": "Australian & Global Hyperarid Drylands",
            "citation": "Bay et al. 2021 (Nat. Microbiol. 10.1038/s41564-020-00811-w); Buzzard et al. 2022",
            "metric": "areal_h2_uptake_mg_m2_yr",
            "unit": "mg_H2_m^-2_yr^-1",
            "literature_low": 5.0,
            "literature_point": 12.5,
            "literature_high": 25.0,
        },
        {
            "desert": "Negev Desert Aridity Gradient (Hyperarid to Sub-humid)",
            "citation": "Bay et al. 2021 (ISME J 10.1038/s41396-021-01001-0)",
            "metric": "actinomycetota_rubisco_cocarriage_fraction",
            "unit": "fraction_of_mags",
            "literature_low": 0.17,
            "literature_point": 0.28,
            "literature_high": 0.35,
        },
        {
            "desert": "Mojave Desert Aridisols (Unirrigated Caliche)",
            "citation": "Schlesinger 1985 (Geochim. Cosmochim. Acta 49:57-66)",
            "metric": "bare_pedogenic_carbonate_t_co2_km2_yr",
            "unit": "t_CO2_km^-2_yr^-1",
            "literature_low": 0.44,
            "literature_point": 0.88,
            "literature_high": 1.54,
        },
        {
            "desert": "Taklamakan / Tarim & Gobi Desert Haloxylon Plantations & Alkaline Soils",
            "citation": "Ma et al. 2013 (Sci. Rep. 3:2025); Ma et al. 2021 (Glob. Ecol. Conserv. 28:e01607)",
            "metric": "revegetated_corridor_total_t_co2_km2_yr",
            "unit": "t_CO2_km^-2_yr^-1",
            "literature_low": 45.0,
            "literature_point": 95.0,
            "literature_high": 185.0,
        },
    ]

    def __init__(self, params: Optional[CarbonBudgetParams] = None) -> None:
        self.params = params or CarbonBudgetParams()
        self.fba = MetabolicCommunityEngine()

    @staticmethod
    def infer_baseline_shrub_cover(site_row: Dict[str, Any]) -> float:
        """Infer baseline fractional shrub cover from site NDVI or transect zone."""
        if "shrub_cover_fraction" in site_row and site_row["shrub_cover_fraction"] is not None:
            return float(np.clip(site_row["shrub_cover_fraction"], 0.0, 0.85))
        ndvi = float(site_row.get("ndvi_W60d", 0.091))
        map_mm = float(site_row.get("map_mm", 24.0))
        # Hyperarid bare sand NDVI baseline is ~0.075; interdune corridors have 1-5% shrub cover
        cover_from_ndvi = (ndvi - 0.075) / 0.45
        cover_from_map = 0.012 + 0.0009 * max(0.0, map_mm - 14.0)
        return float(np.clip(0.5 * cover_from_ndvi + 0.5 * cover_from_map, 0.005, 0.08))

    def compute_state_carbon_budget(
        self,
        site_row: Dict[str, Any],
        shrub_cover_fraction: Optional[float] = None,
        rak_pr_pgpr_dose: float = 0.0,
        rbcl_hyab_plasmid_dose: float = 0.0,
        ca_so4_to_nacl_ratio: float = 1.0,
        discharge_salinity_ms_cm: float = 1.8,
        effective_porosity: float = 0.40,
        vwc_theta: float = 0.018,
        hard_maintenance_feasible: bool = True,
        cbb_flux_multiplier: float = 1.0,
        biomass_density_multiplier: float = 1.0,
    ) -> Dict[str, Any]:
        """Compute the 4-pathway dimensional CO2 and CO2-eq budget for a single state or scenario."""
        p = self.params
        cover = (
            float(np.clip(shrub_cover_fraction, 0.0, 0.85))
            if shrub_cover_fraction is not None
            else self.infer_baseline_shrub_cover(site_row)
        )
        pgpr = float(np.clip(rak_pr_pgpr_dose, 0.0, 1.0))
        plasmid = float(np.clip(rbcl_hyab_plasmid_dose, 0.0, 1.0))

        h1l_ge = max(0.0, float(site_row.get("hydB_1l_ge", 0.75)))
        h1h_ge = max(0.0, float(site_row.get("hydB_1h_ge", 0.10)))
        h2a_ge = max(0.0, float(site_row.get("hydB_2a_ge", 0.03)))
        cox_ge = max(0.0, float(site_row.get("coxL_formI_ge", 0.10)))

        # Empirical CBB co-carriage share from the state cube (unclipped conditional ratio)
        share_1l = max(0.0, float(site_row.get("share_1l_mapped", 0.2535)))
        share_1l_cbb = max(0.0, float(site_row.get("share_1l_cbb_mapped", 0.0885)))
        share_1l_cbb_3g = max(0.0, float(site_row.get("share_1l_cbb_3gene_mapped", 0.0497)))
        share_cbb_total = max(0.0, float(site_row.get("share_cbb_mapped", 0.2235)))
        share_cbb_3g_total = max(0.0, float(site_row.get("share_cbb_3gene_mapped", 0.1095)))

        cbb_cocarriage_cond = (
            float(np.clip(share_1l_cbb / share_1l, 0.0, 1.0)) if share_1l > 1e-9 else 0.0
        )
        cbb_3gene_cocarriage_cond = (
            float(np.clip(share_1l_cbb_3g / share_1l, 0.0, 1.0)) if share_1l > 1e-9 else 0.0
        )

        # Millington-Quirk gas diffusivity relative to reference (eps=0.40, theta=0.018)
        eps_tot = float(np.clip(effective_porosity, 0.12, 0.55))
        theta = float(np.clip(vwc_theta, 0.001, eps_tot - 1e-4))
        eps_air = max(0.0, eps_tot - theta)
        mq_factor = (eps_air ** (10.0 / 3.0)) / (eps_tot**2)
        ref_mq = ((0.40 - 0.018) ** (10.0 / 3.0)) / (0.40**2)
        rel_diff = mq_factor / max(1e-9, ref_mq)
        moisture_gate = 0.75 + 0.65 * float(np.sin(np.pi * np.clip(theta / 0.40, 0.0, 1.0)))

        # Trace-gas oxidation power and scaled uptake (mg gas m^-2 yr^-1)
        h2_power_uw = (29.8 * h1l_ge + 8.5 * (h1h_ge + h2a_ge)) * rel_diff * moisture_gate
        co_power_uw = 1.95 * cox_ge * rel_diff * moisture_gate
        h2_scale = (h2_power_uw / max(1e-6, p.ref_h2_power_uw_m3)) * max(
            0.2, min(2.5, biomass_density_multiplier ** 0.35)
        )
        co_scale = (co_power_uw / max(1e-6, p.ref_co_power_uw_m3)) * max(
            0.2, min(2.5, biomass_density_multiplier ** 0.35)
        )

        # Lower-bound single-density (10^6 cells/gdw) thermodynamic power conversion vs empirical deposition envelope
        bare_thermo_h2_mg_m2_yr = (
            (h2_power_uw * 1e-6 * p.active_soil_depth_m * 31557600.0 / 203000.0) * 2.016 * 1000.0
        )

        h2_mg_m2_yr_low = p.h2_uptake_ref_mg_m2_yr_low * h2_scale
        h2_mg_m2_yr_point = p.h2_uptake_ref_mg_m2_yr_point * h2_scale
        h2_mg_m2_yr_high = p.h2_uptake_ref_mg_m2_yr_high * h2_scale

        co_mg_m2_yr_low = p.co_uptake_ref_mg_m2_yr_low * co_scale
        co_mg_m2_yr_point = p.co_uptake_ref_mg_m2_yr_point * co_scale
        co_mg_m2_yr_high = p.co_uptake_ref_mg_m2_yr_high * co_scale

        # Pathway 1a: Direct Microbial Dark RuBisCO + 3-HP/4-HB CO2 Fixation
        # Zeroed strictly when hard maintenance ATP requirement is violated (e.g., >200 mS/cm Na-Cl brine)
        soil_gdw_m2 = p.dry_bulk_density_g_m3 * p.active_soil_depth_m
        guild_cells_m2 = soil_gdw_m2 * p.hyperarid_cells_per_gdw * p.guild_cell_fraction
        guild_cell_h2_nmol_h_cell = (
            (h2_mg_m2_yr_point * 1e-3 / 2.016) * 1e9 / (8760.0 * max(1.0, guild_cells_m2))
        )

        if not hard_maintenance_feasible or (h2_power_uw + co_power_uw) <= 1e-9:
            rubisco_t_km2_low = 0.0
            rubisco_t_km2_point = 0.0
            rubisco_t_km2_high = 0.0
            rubisco_3g_t_km2_point = 0.0
            dark_co2_nmol_gdw_d_annual = 0.0
            dark_co2_nmol_gdw_d_wetted_active = 0.0
        else:
            # Convert mg H2 m^-2 yr^-1 (MW=2.016) and mg CO m^-2 yr^-1 (MW=28.01) to mmol electron-pair donors
            mmol_donor_low = (h2_mg_m2_yr_low / 2.016) + 0.70 * (co_mg_m2_yr_low / 28.01)
            mmol_donor_point = (h2_mg_m2_yr_point / 2.016) + 0.70 * (co_mg_m2_yr_point / 28.01)
            mmol_donor_high = (h2_mg_m2_yr_high / 2.016) + 0.70 * (co_mg_m2_yr_high / 28.01)

            # Plasmid inoculation and rhizosphere mixotrophy boost CBB electron allocation
            plasmid_cbb_boost = (1.0 + 0.65 * plasmid * max(0.15, cover)) * max(0.0, cbb_flux_multiplier)
            eff_cbb_share = min(0.85, cbb_cocarriage_cond * plasmid_cbb_boost)
            eff_cbb_3g_share = min(0.85, cbb_3gene_cocarriage_cond * plasmid_cbb_boost)

            # mmol CO2 fixed m^-2 yr^-1 -> t CO2 km^-2 yr^-1
            # 1 mmol CO2 m^-2 yr^-1 = 44.01e-3 g CO2 m^-2 yr^-1 = 0.04401 t CO2 km^-2 yr^-1
            mmol_co2_m2_yr_low = (
                (mmol_donor_low / p.mol_h2_per_mol_co2_cbb)
                * eff_cbb_share
                * p.anabolic_electron_fraction_low
            )
            mmol_co2_m2_yr_point = (
                (mmol_donor_point / p.mol_h2_per_mol_co2_cbb)
                * eff_cbb_share
                * p.anabolic_electron_fraction_point
            )
            mmol_co2_m2_yr_high = (
                (mmol_donor_high / p.mol_h2_per_mol_co2_cbb)
                * eff_cbb_share
                * p.anabolic_electron_fraction_high
            )
            mmol_co2_3g_point = (
                (mmol_donor_point / p.mol_h2_per_mol_co2_cbb)
                * eff_cbb_3g_share
                * p.anabolic_electron_fraction_point
            )

            rubisco_t_km2_low = 0.04401 * mmol_co2_m2_yr_low
            rubisco_t_km2_point = 0.04401 * mmol_co2_m2_yr_point
            rubisco_t_km2_high = 0.04401 * mmol_co2_m2_yr_high
            rubisco_3g_t_km2_point = 0.04401 * mmol_co2_3g_point

            # (i) 365-day annualized dry-column rate (0-10 cm, 1.6e5 gdw m^-2)
            dark_co2_nmol_gdw_d_annual = (mmol_co2_m2_yr_point * 1e6) / (soil_gdw_m2 * 365.0)
            # (ii) Active wetted-window rate in top 0-2 cm horizon (comparable to 24-72h wetted 14CO2 microcosms)
            top_gdw_m2 = p.dry_bulk_density_g_m3 * p.top_biocrust_assay_depth_m
            dark_co2_nmol_gdw_d_wetted_active = (mmol_co2_m2_yr_point * 1e6) / (
                top_gdw_m2 * p.active_moisture_days_per_yr
            )

        # Pathway 1b: Indirect Atmospheric H2 Scavenging Sink (t CO2-eq km^-2 yr^-1)
        # 1 mg H2 m^-2 yr^-1 = 1e-3 t H2 km^-2 yr^-1
        h2_co2eq_t_km2_low = (h2_mg_m2_yr_low * 1e-3) * p.h2_gwp100_low
        h2_co2eq_t_km2_point = (h2_mg_m2_yr_point * 1e-3) * p.h2_gwp100_point
        h2_co2eq_t_km2_high = (h2_mg_m2_yr_high * 1e-3) * p.h2_gwp100_high

        # Pathway 2: Photosynthetic Shrub + Rhizosphere Gross NPP (and net retained SOC)
        map_mm = float(site_row.get("map_mm", 25.0))
        precip_factor = float(np.clip((map_mm / 25.0) ** 0.40, 0.65, 1.55))
        inoc_factor = 1.0 + 0.25 * pgpr + 0.15 * plasmid
        salinity_penalty = (
            float(np.exp(-max(0.0, discharge_salinity_ms_cm - 35.0) / 60.0))
            if hard_maintenance_feasible
            else 0.05
        )
        veg_t_km2_low = p.veg_full_canopy_t_co2_km2_yr_low * cover * precip_factor * inoc_factor * salinity_penalty
        veg_t_km2_point = (
            p.veg_full_canopy_t_co2_km2_yr_point * cover * precip_factor * inoc_factor * salinity_penalty
        )
        veg_t_km2_high = p.veg_full_canopy_t_co2_km2_yr_high * cover * precip_factor * inoc_factor * salinity_penalty

        # Pathway 3: Pedogenic & Sabkha Carbonate Binding (CaCO3 / CaMg(CO3)2)
        ca_wt = max(0.1, float(site_row.get("Ca", 4.5)))
        mg_wt = max(0.05, float(site_row.get("Mg", 1.1)))
        ph_val = float(site_row.get("ph_mean", 8.15))
        alkaline_gate = float(np.clip((ph_val - 7.1) / 1.2, 0.25, 1.45))
        divalent_mineral_factor = float(np.clip(((ca_wt + 1.2 * mg_wt) / 6.0) ** 0.35, 0.60, 1.75))
        # Root/microbial respiration + Ca-SO4 irrigation provides HCO3^- and Ca2+ for secondary CaCO3
        biogenic_ca_so4_boost = 1.0 + 4.5 * cover * (1.0 + 0.45 * max(0.0, min(2.5, ca_so4_to_nacl_ratio - 0.5)))
        sic_t_km2_low = p.sic_bare_ref_t_co2_km2_yr_low * alkaline_gate * divalent_mineral_factor * biogenic_ca_so4_boost
        sic_t_km2_point = (
            p.sic_bare_ref_t_co2_km2_yr_point * alkaline_gate * divalent_mineral_factor * biogenic_ca_so4_boost
        )
        sic_t_km2_high = p.sic_bare_ref_t_co2_km2_yr_high * alkaline_gate * divalent_mineral_factor * biogenic_ca_so4_boost

        # Gross direct CO2 bound across the three carbon-fixation/mineralization pathways
        total_direct_co2_low = rubisco_t_km2_low + veg_t_km2_low + sic_t_km2_low
        total_direct_co2_point = rubisco_t_km2_point + veg_t_km2_point + sic_t_km2_point
        total_direct_co2_high = rubisco_t_km2_high + veg_t_km2_high + sic_t_km2_high

        # Estimated net persistent carbon storage (SOC/woody retention + non-carbonate Ca2+ SIC + necromass)
        net_persistent_co2_low = (
            p.microbial_necromass_retention_low * rubisco_t_km2_low
            + p.soc_retention_fraction_low * veg_t_km2_low
            + p.non_carbonate_ca_fraction_low * sic_t_km2_low
        )
        net_persistent_co2_point = (
            p.microbial_necromass_retention_point * rubisco_t_km2_point
            + p.soc_retention_fraction_point * veg_t_km2_point
            + p.non_carbonate_ca_fraction_point * sic_t_km2_point
        )
        net_persistent_co2_high = (
            p.microbial_necromass_retention_high * rubisco_t_km2_high
            + p.soc_retention_fraction_high * veg_t_km2_high
            + p.non_carbonate_ca_fraction_high * sic_t_km2_high
        )

        total_co2eq_low = total_direct_co2_low + h2_co2eq_t_km2_low
        total_co2eq_point = total_direct_co2_point + h2_co2eq_t_km2_point
        total_co2eq_high = total_direct_co2_high + h2_co2eq_t_km2_high

        return {
            "shrub_cover_fraction": round(cover, 4),
            "cbb_community_share_mapped": round(share_cbb_total, 4),
            "cbb_3gene_community_share_mapped": round(share_cbb_3g_total, 4),
            "cbb_1l_cocarriage_conditional": round(cbb_cocarriage_cond, 4),
            "cbb_1l_3gene_cocarriage_conditional": round(cbb_3gene_cocarriage_cond, 4),
            "bare_thermodynamic_1e6_h2_mg_m2_yr": round(bare_thermo_h2_mg_m2_yr, 4),
            "h2_uptake_mg_m2_yr": {
                "low": round(h2_mg_m2_yr_low, 4),
                "point": round(h2_mg_m2_yr_point, 4),
                "high": round(h2_mg_m2_yr_high, 4),
            },
            "co_uptake_mg_m2_yr": {
                "low": round(co_mg_m2_yr_low, 4),
                "point": round(co_mg_m2_yr_point, 4),
                "high": round(co_mg_m2_yr_high, 4),
            },
            "dark_co2_fixation_annual_avg_nmol_gdw_d": round(dark_co2_nmol_gdw_d_annual, 6),
            "dark_co2_fixation_wetted_active_nmol_gdw_d": round(dark_co2_nmol_gdw_d_wetted_active, 5),
            "guild_cell_specific_h2_oxidation_nmol_h_cell": float(f"{guild_cell_h2_nmol_h_cell:.3e}"),
            "direct_microbial_rubisco_t_co2_km2_yr": {
                "low": round(rubisco_t_km2_low, 5),
                "point": round(rubisco_t_km2_point, 5),
                "high": round(rubisco_t_km2_high, 5),
            },
            "direct_microbial_rubisco_3gene_strict_t_co2_km2_yr_point": round(rubisco_3g_t_km2_point, 5),
            "indirect_h2_sink_t_co2eq_km2_yr": {
                "low": round(h2_co2eq_t_km2_low, 4),
                "point": round(h2_co2eq_t_km2_point, 4),
                "high": round(h2_co2eq_t_km2_high, 4),
            },
            "vegetation_npp_t_co2_km2_yr": {
                "low": round(veg_t_km2_low, 4),
                "point": round(veg_t_km2_point, 4),
                "high": round(veg_t_km2_high, 4),
            },
            "pedogenic_carbonate_t_co2_km2_yr": {
                "low": round(sic_t_km2_low, 4),
                "point": round(sic_t_km2_point, 4),
                "high": round(sic_t_km2_high, 4),
            },
            "total_direct_co2_bound_t_km2_yr": {
                "low": round(total_direct_co2_low, 4),
                "point": round(total_direct_co2_point, 4),
                "high": round(total_direct_co2_high, 4),
            },
            "net_persistent_co2_sequestered_t_km2_yr": {
                "low": round(net_persistent_co2_low, 4),
                "point": round(net_persistent_co2_point, 4),
                "high": round(net_persistent_co2_high, 4),
            },
            "total_co2eq_bound_t_km2_yr": {
                "low": round(total_co2eq_low, 4),
                "point": round(total_co2eq_point, 4),
                "high": round(total_co2eq_high, 4),
            },
        }

    def _evaluate_policy_on_cell(
        self,
        base_site_row: Dict[str, Any],
        cell_map_mm: float,
        policy: Optional[Dict[str, Any]],
    ) -> Dict[str, Any]:
        """Evaluate baseline or counterfactual policy on a single supported corridor grid cell."""
        r_dict = dict(base_site_row)
        r_dict["map_mm"] = float(cell_map_mm)
        if policy is None:
            return self.compute_state_carbon_budget(r_dict)

        cov = float(policy.get("shrub_cover_fraction", 0.04))
        pgpr = float(policy.get("rak_pr_pgpr_dose", 0.0))
        plasmid = float(policy.get("site60_rbcl_hyab_plasmid_dose", 0.0))
        ec = float(policy.get("discharge_salinity_ms_cm", 15.0))
        q_ratio = float(policy.get("ca_so4_to_nacl_ratio", 1.0))
        ect_dose = float(policy.get("sphingomonas_ect_kdp_plasmid_dose", 0.0))
        betaine_mm = float(policy.get("exogenous_choline_betaine_mm", 0.0))
        camels = float(policy.get("camel_head_per_100km2", 2.0))
        grading = float(policy.get("cut_line_grading_fraction", 0.10))
        nodal = bool(policy.get("low_impact_nodal_mode", False))
        rot = bool(policy.get("rotational_water_point_management", False))

        camel_units = camels * 2.8 * (0.55 if rot else 1.0)
        seis_mult = 0.18 if nodal else 1.0
        eff_porosity = max(0.22, 0.40 - 0.0035 * camel_units - 0.11 * grading * seis_mult)

        effective_nacl_ec = ec / (1.0 + 0.55 * q_ratio)
        fba_res = self.fba.solve_community_fba(
            CommunityMetabolicInputs(
                h2_uptake_bound=1.6 * np.exp(-max(0.0, effective_nacl_ec - 12.0) / 65.0),
                co_uptake_bound=0.30 * np.exp(-max(0.0, effective_nacl_ec - 12.0) / 75.0),
                root_exudate_doc_bound=0.04 if effective_nacl_ec <= 25.0 else 0.01,
                exogenous_nh4_no3_bound=0.10,
                sulfate_so4_bound=0.45 * min(2.0, q_ratio),
                salinity_ms_cm=effective_nacl_ec,
                choline_betaine_bound=betaine_mm,
                rbcl_hyab_plasmid_inoc=plasmid if effective_nacl_ec < 35.0 else 0.05 * plasmid,
                ect_kdp_plasmid_inoc=ect_dose,
            )
        )
        hard_feasible = bool(fba_res["hard_maintenance_feasible"])

        root_bio_mult = 1.0 + 9.34 * (0.88 + 0.24 * pgpr + 0.12 * plasmid)
        mosaic_bio_mult = cov * root_bio_mult + (1.0 - cov) * 1.0
        cbb_mult = 1.0 + 0.65 * plasmid * cov

        base_1l = float(r_dict.get("hydB_1l_ge", 0.75))
        root_1l = base_1l * (1.0 - 0.793 * (0.95 + 0.05 * pgpr))
        r_dict["hydB_1l_ge"] = cov * root_1l + (1.0 - cov) * base_1l

        return self.compute_state_carbon_budget(
            r_dict,
            shrub_cover_fraction=cov,
            rak_pr_pgpr_dose=pgpr,
            rbcl_hyab_plasmid_dose=plasmid,
            ca_so4_to_nacl_ratio=q_ratio,
            discharge_salinity_ms_cm=ec,
            effective_porosity=eff_porosity,
            hard_maintenance_feasible=hard_feasible,
            cbb_flux_multiplier=cbb_mult,
            biomass_density_multiplier=mosaic_bio_mult,
        )

    def upscale_corridor_and_basin_budget(
        self,
        comp_df: pd.DataFrame,
        grid_df: pd.DataFrame,
        policy: Optional[Dict[str, Any]] = None,
        portfolio_by_code: Optional[Dict[str, Dict[str, Any]]] = None,
    ) -> Dict[str, Any]:
        """Upscale CO2 and CO2-eq fluxes across the 128 supported corridor cells and the Rub' al-Khali.

        Applies active management policies ONLY to the 128 supported corridor cells (92,160 km^2),
        evaluating each cell using its interpolated climate (`map_mm`) and nearest primary site
        geochemical/microbial state, while the remaining unmanaged desert interior remains at the
        bare/native hyperarid baseline (`unmanaged_reserve_delta_kt_yr == 0.0` identically).
        """
        p = self.params
        surf_df = comp_df[comp_df["compartment_norm"] == "Surface"].copy()
        site_lookup: Dict[int, Dict[str, Any]] = {
            int(row["site_num"]): row.to_dict() for _, row in surf_df.iterrows()
        }

        supp_cells = grid_df[grid_df["in_empirical_support"]].copy()
        n_supported_cells = len(supp_cells)
        supported_area_km2 = n_supported_cells * p.corridor_cell_area_km2
        unmanaged_desert_area_km2 = max(0.0, p.rub_al_khali_total_area_km2 - supported_area_km2)
        unmanaged_sand_sea_area_km2 = max(0.0, p.rub_al_khali_sand_sea_polygon_km2 - supported_area_km2)

        # 1. Compute bare hyperarid background rates across the unmanaged interior (f_shrub = 0.002)
        bg_budgets = [
            self.compute_state_carbon_budget(row.to_dict(), shrub_cover_fraction=0.002)
            for _, row in surf_df.iterrows()
        ]

        # 2. Compute baseline and policy budgets across the 128 supported corridor cells
        baseline_corridor_budgets: List[Dict[str, Any]] = []
        corridor_budgets: List[Dict[str, Any]] = []
        for _, cell in supp_cells.iterrows():
            nearest_site_num = int(cell.get("nearest_primary_site", 35))
            base_row = site_lookup.get(nearest_site_num, surf_df.iloc[0].to_dict())
            cell_map = float(cell.get("map_mm", base_row.get("map_mm", 24.0)))
            base_cell_b = self._evaluate_policy_on_cell(base_row, cell_map, policy=None)
            baseline_corridor_budgets.append(base_cell_b)
            if portfolio_by_code is not None:
                p_code = str(cell.get("prescribed_portfolio", "P2"))
                cell_pol = portfolio_by_code.get(p_code, policy)
                corridor_budgets.append(self._evaluate_policy_on_cell(base_row, cell_map, policy=cell_pol))
            elif policy is None:
                corridor_budgets.append(base_cell_b)
            else:
                corridor_budgets.append(self._evaluate_policy_on_cell(base_row, cell_map, policy=policy))

        def _agg_pathway(key: str) -> Dict[str, float]:
            corr_vals = {
                b: float(np.mean([x[key][b] for x in corridor_budgets]))
                for b in ("low", "point", "high")
            }
            base_corr_vals = {
                b: float(np.mean([x[key][b] for x in baseline_corridor_budgets]))
                for b in ("low", "point", "high")
            }
            bg_vals = {
                b: float(np.mean([x[key][b] for x in bg_budgets]))
                for b in ("low", "point", "high")
            }
            delta_t_km2 = corr_vals["point"] - base_corr_vals["point"]
            # Convert t km^-2 yr^-1 * km^2 -> kt yr^-1 (divide by 1000)
            return {
                "corridor_t_km2_yr_point": round(corr_vals["point"], 4),
                "delta_vs_baseline_t_km2_yr_point": round(delta_t_km2, 4),
                "corridor_kt_yr_low": round(corr_vals["low"] * supported_area_km2 / 1e3, 3),
                "corridor_kt_yr_point": round(corr_vals["point"] * supported_area_km2 / 1e3, 3),
                "corridor_kt_yr_high": round(corr_vals["high"] * supported_area_km2 / 1e3, 3),
                "delta_vs_baseline_corridor_kt_yr_low": round(
                    (corr_vals["low"] - base_corr_vals["low"]) * supported_area_km2 / 1e3, 3
                ),
                "delta_vs_baseline_corridor_kt_yr_point": round(
                    delta_t_km2 * supported_area_km2 / 1e3, 3
                ),
                "delta_vs_baseline_corridor_kt_yr_high": round(
                    (corr_vals["high"] - base_corr_vals["high"]) * supported_area_km2 / 1e3, 3
                ),
                "sand_sea_500k_kt_yr_point": round(
                    (corr_vals["point"] * supported_area_km2 + bg_vals["point"] * unmanaged_sand_sea_area_km2)
                    / 1e3,
                    3,
                ),
                "basin_kt_yr_low": round(
                    (corr_vals["low"] * supported_area_km2 + bg_vals["low"] * unmanaged_desert_area_km2)
                    / 1e3,
                    3,
                ),
                "basin_kt_yr_point": round(
                    (
                        corr_vals["point"] * supported_area_km2
                        + bg_vals["point"] * unmanaged_desert_area_km2
                    )
                    / 1e3,
                    3,
                ),
                "basin_kt_yr_high": round(
                    (
                        corr_vals["high"] * supported_area_km2
                        + bg_vals["high"] * unmanaged_desert_area_km2
                    )
                    / 1e3,
                    3,
                ),
            }

        return {
            "policy_name": policy["name"] if policy and "name" in policy else "Current_Baseline_State",
            "supported_corridor_area_km2": supported_area_km2,
            "unmanaged_reserve_area_km2": unmanaged_desert_area_km2,
            "sand_sea_polygon_area_km2": p.rub_al_khali_sand_sea_polygon_km2,
            "total_basin_area_km2": p.rub_al_khali_total_area_km2,
            "unmanaged_reserve_management_delta_kt_yr": 0.0,
            "direct_microbial_rubisco_co2": _agg_pathway("direct_microbial_rubisco_t_co2_km2_yr"),
            "indirect_h2_sink_co2eq": _agg_pathway("indirect_h2_sink_t_co2eq_km2_yr"),
            "vegetation_npp_co2": _agg_pathway("vegetation_npp_t_co2_km2_yr"),
            "pedogenic_carbonate_co2": _agg_pathway("pedogenic_carbonate_t_co2_km2_yr"),
            "total_direct_co2_bound": _agg_pathway("total_direct_co2_bound_t_km2_yr"),
            "net_persistent_co2_sequestered": _agg_pathway("net_persistent_co2_sequestered_t_km2_yr"),
            "total_co2eq_bound": _agg_pathway("total_co2eq_bound_t_km2_yr"),
        }
