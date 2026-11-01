"""Slot resolver, formatter and supervisor-written science-page texts for EQ-BioDT v0.8."""
from __future__ import annotations

import json
import math
from pathlib import Path
import re
from typing import Any

ROOT = Path(__file__).resolve().parents[2]

READABLE_INPUT_NAMES: dict[str, str] = {
    "k_h2_s": "the H₂ rate constant k_H2",
    "k_co_s": "the CO rate constant k_CO",
    "retention_b": "the retention exponent b",
    "psi_entry_mpa": "the air-entry potential ψ_e",
    "q10": "Q10",
    "wet_activation_mpa": "the activation potential",
    "dry_threshold_mpa": "the dry threshold",
    "activation_shape": "the activation shape exponent",
    "d_h2_m2_s": "the H₂ diffusivity",
    "d_co_m2_s": "the CO diffusivity",
    "porosity": "porosity",
    "drainage_fraction_day": "the drainage fraction",
    "field_capacity_fraction": "field capacity",
    "cloud_seeding_rain_increase_frac": "the seeding efficacy",
    "mulch_evap_reduction_frac": "the mulch evaporation factor",
    "shade_net_soil_temp_reduction_c": "the shade-net cooling",
    "solar_panel_evap_reduction_frac": "the panel evaporation reduction",
    "solar_panel_soil_temp_reduction_c": "the panel cooling",
    "panel_wash_l_m2_per_wash": "the washing volume",
    "evap_scaling_per_c": "the evaporation scaling",
    "biochar_porosity_increase_per_pct": "the biochar porosity response",
    "offroad_bd_max_increase_g_cm3": "the maximum compaction",
    "offroad_bd_passes_efold": "the compaction rate",
    "sa187_mock_ratio_fresh": "the digitised mock yield ratio",
    "sa187_mock_ratio_dry": "the digitised dry-mass yield ratio",
    "sa187_benefit_b": "the SA187 yield multiplier",
    "fao29_ece_per_ecw": "the FAO-29 leaching factor",
    "mh_alfalfa_slope_pct_per_dS_m": "the alfalfa salt-tolerance slope",
    "mh_alfalfa_threshold_dS_m": "the alfalfa salt-tolerance threshold",
    "knallgas_co2_per_h2": "the knallgas yield",
    "nacl_mol_l_per_dS_m": "the conductivity-to-NaCl conversion",
    # Extended consistently for all remaining registry inputs
    "air_temp_offset_C": "the air-temperature offset",
    "carbon_molar_mass_g_mol": "the molar mass of carbon",
    "crop_water_demand_alfalfa_mm_day": "the alfalfa water demand",
    "crop_water_demand_desert_mm_day": "the desert crop water demand",
    "h2_background_ppm": "the background H₂ mixing ratio",
    "mh_barley_slope_pct_per_dS_m": "the barley salt-tolerance slope",
    "mh_barley_threshold_dS_m": "the barley salt-tolerance threshold",
    "mh_date_palm_slope_pct_per_dS_m": "the date-palm salt-tolerance slope",
    "mh_date_palm_threshold_dS_m": "the date-palm salt-tolerance threshold",
    "mh_tomato_slope_pct_per_dS_m": "the tomato salt-tolerance slope",
    "mh_tomato_threshold_dS_m": "the tomato salt-tolerance threshold",
    "offroad_bd_initial_g_cm3": "the initial bulk density",
    "panel_wash_interval_days": "the panel-washing interval",
    "quartz_particle_density_g_cm3": "the quartz particle density",
    "rain_scaling_per_c": "the rainfall scaling with warming",
    "sa187_ecw_high_dS_m": "the high irrigation salinity",
    "sa187_ecw_low_dS_m": "the low irrigation salinity",
    "sa187_soil_ec_dS_m": "the site soil salinity",
    "shade_net_evap_reduction_frac": "the shade-net evaporation reduction",
    "fao29_ecsw_per_ece": "the soil-water-to-extract EC ratio",
    "crop_etc_mm_per_day": "the desert crop water demand",
    "alfalfa_etc_mm_per_day": "the alfalfa water demand",
    "molar_mass_c_g_mol": "the molar mass of carbon",
    "days_per_year": "the number of days per year",
    "biocrust_net_c_fixation_g_c_m2_yr": "the biocrust net carbon fixation rate",
    "h2_hydrogen_economy_global_ppm": "the global hydrogen-economy H₂ concentration",
    "particle_density_g_cm3": "the particle density",
    "soil_temp_per_air_temp_warming": "the soil-to-air warming ratio",
    "sa187_effect_fresh_low_s1": "the season-1 fresh-mass SA187 effect",
    "sa187_effect_fresh_low_s2": "the season-2 fresh-mass SA187 effect",
    "sa187_effect_dry_low_s1": "the season-1 dry-mass SA187 effect",
    "sa187_effect_dry_low_s2": "the season-2 dry-mass SA187 effect",
}

SCIENCE_TEMPLATES: dict[str, dict[str, str]] = {
    "rain-pattern": {
        "title_science": "Rainfall partitioning and 14-day H₂ uptake",
        "answer_science": (
            "Partitioning the same 12 mm of rain into smaller events reduced H₂ uptake in every parameter draw. "
            "Relative to a single 12 mm event, four 3 mm events reduced 14-day uptake by a median of "
            "{rain-schedule/2:h2_uptake_change_pct:median|absneg} (90 % interval "
            "{rain-schedule/2:h2_uptake_change_pct:ci90|absneg}; reference parameterisation "
            "{rain-schedule/2:h2_uptake_change_pct:central}), and daily drizzle reduced it by "
            "{rain-schedule/3:h2_uptake_change_pct:median|absneg} "
            "({rain-schedule/3:h2_uptake_change_pct:ci90|absneg}). Each small addition is largely removed by the "
            "daily evaporative demand, so the soil water potential spends fewer days above the activation threshold. "
            "The magnitude of the reduction depends mainly on the water-retention curve "
            "({rain-schedule/2:h2_uptake_change_pct:st1name}, S_T = {rain-schedule/2:h2_uptake_change_pct:st1|idx}; "
            "{rain-schedule/2:h2_uptake_change_pct:st2name}, S_T = {rain-schedule/2:h2_uptake_change_pct:st2|idx}); "
            "its direction holds across all sampled parameter values."
        ),
        "limits_science": (
            "The intervals combine uncertainty in the twin's retention, activation and kinetic parameters. "
            "The twin represents the top 10 cm of sand as one well-mixed layer at 25 °C with 1 mm d⁻¹ evaporative "
            "demand, so it resolves the timing of wetting but no vertical redistribution, runoff or vapour transport."
        ),
    },
    "drip-or-flood": {
        "title_science": "Irrigation frequency at constant irrigation volume",
        "answer_science": (
            "Delivering the same irrigation volume, on top of the 12 mm rain, daily, every three days or every ten "
            "days changed 14-day H₂ uptake by a median of {irrigation-interval/0:h2_uptake_change_pct:median}, "
            "{irrigation-interval/1:h2_uptake_change_pct:median} and "
            "{irrigation-interval/2:h2_uptake_change_pct:median}, respectively. The 90 % intervals of all three "
            "schedules include both decreases and increases ({irrigation-interval/0:h2_uptake_change_pct:ci90}, "
            "{irrigation-interval/1:h2_uptake_change_pct:ci90} and "
            "{irrigation-interval/2:h2_uptake_change_pct:ci90}). Irrigation acts in two directions in the twin: it "
            "lengthens the period above the activation threshold, and it fills pores, which slows H₂ diffusion to "
            "the reacting community. The balance is set by the retention curve; "
            "{irrigation-interval/2:h2_uptake_change_pct:st1name} and "
            "{irrigation-interval/2:h2_uptake_change_pct:st2name} explain most of the variance for the ten-day "
            "schedule (S_T = {irrigation-interval/2:h2_uptake_change_pct:st1|idx} and "
            "{irrigation-interval/2:h2_uptake_change_pct:st2|idx}). Ranking the schedules therefore requires a "
            "measured retention curve for Rub' al Khali sand."
        ),
        "limits_science": (
            "The irrigation volume follows the FAO-56 order-of-magnitude crop demand and is a scenario input. "
            "The twin uses one well-mixed 10 cm layer, so drip emitters, wetting fronts and deeper drainage are "
            "represented only through the daily water balance of that layer."
        ),
    },
    "cloud-seeding": {
        "title_science": "Additional rain from cloud seeding",
        "answer_science": (
            "With the seeding efficacy sampled between no effect and 20 % additional rain on rain days (central "
            "10 %), 14-day H₂ uptake rose by a median of {seeding-intensity/2:h2_uptake_change_pct:median|abspos} "
            "(90 % interval {seeding-intensity/2:h2_uptake_change_pct:ci90|abspos}; reference parameterisation "
            "{seeding-intensity/2:h2_uptake_change_pct:central}); doubling the seeding intensity gave "
            "{seeding-intensity/3:h2_uptake_change_pct:median|abspos} "
            "({seeding-intensity/3:h2_uptake_change_pct:ci90|abspos}). Every draw increased uptake because the extra "
            "water lengthens the period above the activation threshold. The efficacy itself is the largest single "
            "source of variance at the central intensity (S_T = {seeding-intensity/2:h2_uptake_change_pct:st1|idx} "
            "for {seeding-intensity/2:h2_uptake_change_pct:st1name}), followed by "
            "{seeding-intensity/2:h2_uptake_change_pct:st2name} (S_T = "
            "{seeding-intensity/2:h2_uptake_change_pct:st2|idx})."
        ),
        "limits_science": (
            "The efficacy range follows the WMO statement on weather modification, which reports no established "
            "quantitative rain enhancement; the result is conditional on the additional rain occurring. The baseline "
            "contains one rain day, so seeding acts on a single event within the 14-day window."
        ),
    },
    "salty-water": {
        "title_science": "Saline rain and daily brackish irrigation",
        "answer_science": (
            "A single 12 mm rain carrying NaCl at 2 dS/m or more concentrates as the soil dries until the pore-water "
            "NaCl exceeds 0.1 mol/L, the limit of the ideal-dilute osmotic model; the twin returns no estimate for "
            "any draw of this series. Daily brackish irrigation keeps the solution dilute (refusals "
            "{daily-irrigation/3:h2_uptake_change_pct:refusal} at 8 dS/m). Relative to the single fresh rain, daily "
            "irrigation at 8 dS/m changed uptake by a median of {daily-irrigation/3:h2_uptake_change_pct:median} "
            "(90 % interval {daily-irrigation/3:h2_uptake_change_pct:ci90}), and daily irrigation with near-fresh "
            "water (0.5 dS/m) by {daily-irrigation/0:h2_uptake_change_pct:median} "
            "({daily-irrigation/0:h2_uptake_change_pct:ci90}). Salinity therefore adds a few percentage points to an "
            "effect set by the irrigation schedule, whose sign depends on the retention curve. With the FAO-29 "
            "leaching relation and the Maas–Hoffman response, alfalfa relative yield at 8 dS/m is "
            "{daily-irrigation/3:relative_yield_alfalfa:median|frac2pct} "
            "({daily-irrigation/3:relative_yield_alfalfa:ci90|frac2pct}), a crop-model estimate driven by the twin's "
            "salt balance."
        ),
        "limits_science": (
            "The osmotic model is valid for dilute solutions; the refusal of the saline-rain series marks the "
            "boundary of the model, and estimates for concentrated pore water require an activity model. Crop yields "
            "come from the FAO-29 and FAO-56 salt-tolerance relations; the twin contributes the soil salt balance of "
            "bare sand."
        ),
    },
    "alfalfa-bacterium": {
        "title_science": "Alfalfa yield under saline irrigation with Enterobacter sp. SA187",
        "answer_science": (
            "We estimated alfalfa relative yield at the high irrigation salinity of the Hada Al-Sham field trial "
            "({const:sa187_ecw_high_dS_m|plain} dS/m) by combining the Maas–Hoffman response with three root-zone "
            "salinity models. With the leaching factor fitted to the trial's season-1 yield ratio, SA187-coated seed "
            "reached a median relative yield of {C_fitted_fresh/4:relative_yield_alfalfa_sa187:median|frac2pct} "
            "(90 % interval {C_fitted_fresh/4:relative_yield_alfalfa_sa187:ci90|frac2pct}) and uncoated seed "
            "{C_fitted_fresh/4:relative_yield_alfalfa:median|frac2pct}; the field-trial inputs carry this uncertainty "
            "({C_fitted_fresh/4:relative_yield_alfalfa_sa187:st1name}, S_T = "
            "{C_fitted_fresh/4:relative_yield_alfalfa_sa187:st1|idx}). The choice of root-zone model dominates the "
            "estimate: the FAO-29 fixed leaching factor gives "
            "{A_fao_fixed_leaching/4:relative_yield_alfalfa_sa187:median|frac2pct} "
            "({A_fao_fixed_leaching/4:relative_yield_alfalfa_sa187:ci90|frac2pct}) and the twin's bare-sand salt "
            "balance {B_twin_physics/4:relative_yield_alfalfa_sa187:median|frac2pct} "
            "({B_twin_physics/4:relative_yield_alfalfa_sa187:ci90|frac2pct}) for coated seed. The fitted leaching "
            "factor implies a root-zone salinity close to the soil salinity measured at the site "
            "({const:sa187_soil_ec_dS_m|plain} dS/m)."
        ),
        "limits_science": (
            "The fitted model reproduces the trial's mock yield ratio by construction, from two salinity levels in "
            "one season of one trial (de Zélicourt et al. 2018). The SA187 benefit enters as a constant multiplier; "
            "its dependence on salinity is unidentified by these data. Yields are crop-model estimates."
        ),
    },
    "hotter": {
        "title_science": "Soil warming and climate warming",
        "answer_science": (
            "Warming the soil at fixed moisture raised 14-day H₂ uptake in every draw, by a median of "
            "{warmer-soil-only/3:h2_uptake_change_pct:median|abspos} at +4 °C "
            "({warmer-soil-only/3:h2_uptake_change_pct:ci90|abspos}) and "
            "{warmer-soil-only/4:h2_uptake_change_pct:median|abspos} at +10 °C "
            "({warmer-soil-only/4:h2_uptake_change_pct:ci90|abspos}), with Q10 explaining most of the variance "
            "(S_T = {warmer-soil-only/4:h2_uptake_change_pct:st1|idx}). Climate warming of the same magnitude also "
            "raises the evaporative demand, which shortens the wet period. At +4 °C with unchanged rainfall the two "
            "effects nearly cancel: the reference parameterisation gives "
            "{warmer-climate/3:h2_uptake_change_pct:central}, the median is "
            "{warmer-climate/3:h2_uptake_change_pct:median}, and the 90 % interval includes decreases and increases "
            "({warmer-climate/3:h2_uptake_change_pct:ci90}). The sign of the climate response is set by the balance "
            "between {warmer-climate/3:h2_uptake_change_pct:st1name} (S_T = "
            "{warmer-climate/3:h2_uptake_change_pct:st1|idx}) and {warmer-climate/3:h2_uptake_change_pct:st2name} "
            "(S_T = {warmer-climate/3:h2_uptake_change_pct:st2|idx})."
        ),
        "limits_science": (
            "Soil temperature follows air temperature one to one, and evaporative demand scales with warming after "
            "Scheff and Frierson (2014). The kinetic temperature response is a single Q10 relation, valid within the "
            "0 to 50 °C range of the twin."
        ),
    },
    "mulch-shade": {
        "title_science": "Organic mulch and shade nets",
        "answer_science": (
            "Full organic mulch raised 14-day H₂ uptake by a median of "
            "{mulch-only/2:h2_uptake_change_pct:median|abspos} (90 % interval "
            "{mulch-only/2:h2_uptake_change_pct:ci90|abspos}) in every draw, because the reduced evaporation keeps "
            "the soil above the activation threshold for longer. A full shade net lowered uptake by "
            "{shade-only/2:h2_uptake_change_pct:median|absneg} ({shade-only/2:h2_uptake_change_pct:ci90|absneg}), "
            "also in every draw: the net cools the soil and slows the reaction through Q10 (S_T = "
            "{shade-only/2:h2_uptake_change_pct:st1|idx}). Combining full mulch and shade changed uptake by "
            "{mulch-and-shade/2:h2_uptake_change_pct:median}, with an interval that includes both directions "
            "({mulch-and-shade/2:h2_uptake_change_pct:ci90}). The mulch effect depends mostly on "
            "{mulch-only/2:h2_uptake_change_pct:st1name} (S_T = {mulch-only/2:h2_uptake_change_pct:st1|idx})."
        ),
        "limits_science": (
            "The mulch evaporation factor extrapolates the FAO-56 rule for partial organic mulch cover linearly to "
            "full cover, and the shade-net cooling is a constant soil-temperature offset; both are modeller "
            "assumptions."
        ),
    },
    "solar-park": {
        "title_science": "Microclimate under photovoltaic panels",
        "answer_science": (
            "Full panel cover raised 14-day H₂ uptake by a median of "
            "{panel-cover/2:h2_uptake_change_pct:median|abspos} (90 % interval "
            "{panel-cover/2:h2_uptake_change_pct:ci90|abspos}; reference parameterisation "
            "{panel-cover/2:h2_uptake_change_pct:central}), and every draw gave an increase. The panels act on the "
            "soil in three ways in the chain: they reduce evaporation and add washing water, both of which lengthen "
            "the wet period, and they cool the soil, which slows the reaction through Q10. At half cover "
            "{panel-cover/1:h2_uptake_change_pct:st1name} explains the largest share of the variance (S_T = "
            "{panel-cover/1:h2_uptake_change_pct:st1|idx}); at full cover "
            "{panel-cover/2:h2_uptake_change_pct:st1name} does (S_T = {panel-cover/2:h2_uptake_change_pct:st1|idx})."
        ),
        "limits_science": (
            "The soil cooling and evaporation reduction come from studies at temperate sites (Armstrong et al. 2016; "
            "Hassanpour Adeh et al. 2018), and the washing volume and interval are scenario values; all three are "
            "modeller assumptions for the Rub' al Khali."
        ),
    },
    "compaction": {
        "title_science": "Soil compaction by vehicle traffic",
        "answer_science": (
            "Fifty vehicle passes lowered porosity from 0.40 to {offroad-passes/3:porosity:central|plain} and "
            "changed 14-day H₂ uptake by a median of {offroad-passes/3:h2_uptake_change_pct:median}; the 90 % "
            "interval includes decreases and increases ({offroad-passes/3:h2_uptake_change_pct:ci90}). Lower "
            "porosity acts in two directions: it reduces the gas-filled pore space, which slows H₂ diffusion, and "
            "it raises the water saturation reached after the rain, which lengthens activation. The net response to "
            "compaction is therefore non-monotonic and depends on {offroad-passes/3:h2_uptake_change_pct:st1name} "
            "and {offroad-passes/3:h2_uptake_change_pct:st2name} (S_T = "
            "{offroad-passes/3:h2_uptake_change_pct:st1|idx} and {offroad-passes/3:h2_uptake_change_pct:st2|idx})."
        ),
        "limits_science": (
            "The conversion from vehicle passes to bulk density is an exponential approach to a maximum increase, "
            "fitted by the modeller to Mojave off-road studies (Webb 1983). Compaction changes porosity only; the "
            "retention curve keeps its parameters."
        ),
    },
    "biochar": {
        "title_science": "Biochar amendment",
        "answer_science": (
            "Adding 5 % biochar by weight raised porosity to {biochar-dose/3:porosity:central|plain} and changed "
            "14-day H₂ uptake by a median of {biochar-dose/3:h2_uptake_change_pct:median}, with a 90 % interval "
            "that includes both directions ({biochar-dose/3:h2_uptake_change_pct:ci90}). The additional pore space "
            "lowers the water saturation reached after the 12 mm rain and enlarges the gas-filled volume; the two "
            "effects nearly cancel. The variance comes mainly from {biochar-dose/3:h2_uptake_change_pct:st1name} "
            "(S_T = {biochar-dose/3:h2_uptake_change_pct:st1|idx}) and "
            "{biochar-dose/3:h2_uptake_change_pct:st2name} (S_T = {biochar-dose/3:h2_uptake_change_pct:st2|idx})."
        ),
        "limits_science": (
            "The porosity response per percent biochar is the modeller's conversion of a meta-analysis mean "
            "(Omondi et al. 2016). Biochar effects on water retention, sorption and microbial activity are outside "
            "the chain."
        ),
    },
    "missing-co-genes": {
        "title_science": "Deductive null effects of absent CO-oxidation capacity",
        "answer_science": (
            "In samples whose reasoner-entailed CO-carrier share is zero ({special:deduced_samples}), CO uptake is "
            "zero by deduction: the community factor γ_CO is zero, so uptake is zero for every parameter draw and "
            "the interval is [0, 0]. In samples with carriers ({special:open_samples}), 14-day CO uptake has a "
            "median of {co-gate-samples/6:co_uptake_umol_m2:median|um} (90 % interval "
            "{co-gate-samples/6:co_uptake_umol_m2:ci90|um}), with {co-gate-samples/6:co_uptake_umol_m2:st1name} "
            "explaining most of the variance (S_T = {co-gate-samples/6:co_uptake_umol_m2:st1|idx}). For "
            "{special:refused_samples} the read-mapping weights needed to compute carrier shares are missing, and "
            "the twin refuses the request."
        ),
        "limits_science": (
            "Absence of entailment reflects the genomic evidence in the sample; a carrier below the detection limit "
            "of the metagenome would change the deduction. The open-sample uptake uses the reference community factor."
        ),
    },
    "adding-microbes": {
        "title_science": "Inoculation as a constitutive intervention on the community",
        "answer_science": (
            "Inoculation changes the community: the reasoner memberships of the added genome are reused and the "
            "community factor γ_CO is recomputed for the hypothetical community. Adding a CO-carrier genome to "
            "51PRr3, whose deduced carrier share is zero, opens CO uptake. The amount depends on the coupling "
            "between carrier share and rate: under presence coupling, uptake is "
            "{51PRr3-gate/1:co_uptake_umol_m2:median|um} (90 % interval "
            "{51PRr3-gate/1:co_uptake_umol_m2:ci90|um}) at any dose, whereas under share-proportional coupling it "
            "rises from {51PRr3-proportional/0:co_uptake_umol_m2:median|um} at 0.1 % to "
            "{51PRr3-proportional/2:co_uptake_umol_m2:median|um} at 5 % inoculum. A non-carrier genome leaves uptake "
            "at zero under both couplings. Where carriers are already present (V53PRr2), proportional coupling gives "
            "gains of {V53PRr2-proportional/1:co_uptake_change_pct:median|abspos} at 1 % and "
            "{V53PRr2-proportional/2:co_uptake_change_pct:median|abspos} at 5 %."
        ),
        "limits_science": (
            "The inoculum share is a scenario value; survival and establishment of the added organism in desert sand "
            "are outside the model. Which coupling holds is untested, and the two couplings bracket the response."
        ),
    },
    "air-hydrogen": {
        "title_science": "Soil H₂ sink under elevated atmospheric H₂",
        "answer_science": (
            "Raising atmospheric H₂ to 1.5, 2 and 10 times the background mixing ratio "
            "({const:h2_background_ppm|plain} ppm; Novelli et al. 1999) increased 14-day soil uptake by "
            "{h2-multiplier/1:h2_uptake_change_pct:median|abspos}, "
            "{h2-multiplier/2:h2_uptake_change_pct:median|abspos} and "
            "{h2-multiplier/3:h2_uptake_change_pct:median|abspos}, identically in every parameter draw. "
            "Proportionality follows from the first-order kinetics of the twin; the percentage change has no "
            "variance, so sensitivity indices are not defined for it. The absolute additional sink at twice "
            "background is {h2-multiplier/2:h2_uptake_contrast_umol_m2:median|um} over 14 days (90 % interval "
            "{h2-multiplier/2:h2_uptake_contrast_umol_m2:ci90|um}), dominated by "
            "{h2-multiplier/2:h2_uptake_contrast_umol_m2:st1name} (S_T = "
            "{h2-multiplier/2:h2_uptake_contrast_umol_m2:st1|idx})."
        ),
        "limits_science": (
            "First-order kinetics holds at trace H₂ concentrations; at much higher concentrations, saturation of "
            "high-affinity hydrogenases would make the response sublinear. The concentration scenarios are inputs."
        ),
    },
    "carbon-fixation": {
        "title_science": "Upper bound on CO₂ fixation by H₂-oxidising autotrophs",
        "answer_science": (
            "We bounded CO₂ fixation by H₂-oxidising autotrophs by converting the twin's H₂ uptake with the "
            "knallgas yield ({const:knallgas_co2_per_h2|plain} mol CO₂ per mol H₂; Bongers 1970), assuming that all "
            "oxidised H₂ supports growth. In sample F14Dr1, the bound after a single 12 mm rain is "
            "{baseline-12mm/1:c_fixed_upper_bound_ug_m2:median|ug} over 14 days; repeating such a rain every two "
            "weeks gives {baseline-12mm/1:c_fixed_upper_bound_g_m2_yr:median|g} (90 % interval "
            "{baseline-12mm/1:c_fixed_upper_bound_g_m2_yr:ci90|g}), with "
            "{baseline-12mm/1:c_fixed_upper_bound_g_m2_yr:st1name} explaining most of the variance (S_T = "
            "{baseline-12mm/1:c_fixed_upper_bound_g_m2_yr:st1|idx}). In sample 14PRr2 the reasoner entails no "
            "H₂-autotroph candidate, and the bound is zero by deduction."
        ),
        "limits_science": (
            "The bound assumes that all energy from H₂ oxidation goes into growth; atmospheric H₂ oxidisers use "
            "much of it for persistence (Jordaan et al. 2020), so realised fixation lies below the bound. Net carbon "
            "removal requires measured stock changes."
        ),
    },
}

_SLOT_RE = re.compile(r"\{([^{}]+)\}")


def _load_registry_constants() -> dict[str, float]:
    reg: dict[str, float] = {}
    for rel in ("data/external/composition_constants.json", "data/external/twin_parameter_ranges.json"):
        doc = json.loads((ROOT / rel).read_text(encoding="utf-8"))
        for list_key in ("constants", "parameters"):
            for entry in doc.get(list_key, []):
                if isinstance(entry, dict) and "id" in entry and "value" in entry:
                    reg[str(entry["id"])] = float(entry["value"])
        for k, v in doc.items():
            if not k.startswith("_") and isinstance(v, dict) and "value" in v:
                reg[k] = float(v["value"])
    return reg


_REGISTRY_CONSTANTS: dict[str, float] | None = None


def get_registry_constant(const_id: str) -> float:
    global _REGISTRY_CONSTANTS
    if _REGISTRY_CONSTANTS is None:
        _REGISTRY_CONSTANTS = _load_registry_constants()
    if const_id not in _REGISTRY_CONSTANTS:
        raise KeyError(f"Unknown constant id {const_id!r}")
    return _REGISTRY_CONSTANTS[const_id]


def readable_input_name(name: str) -> str:
    if name not in READABLE_INPUT_NAMES:
        raise KeyError(f"Missing readable input name for {name!r}")
    return READABLE_INPUT_NAMES[name]


def _fmt_pct_mag(val: float) -> str:
    """Format an unsigned percentage magnitude: 2 sig figs below 10, integer at or above 10."""
    av = abs(float(val))
    if av < 1e-12:
        return "0"
    exp = math.floor(math.log10(av))
    decimals = max(0, 1 - exp)
    rounded_2sf = round(av, decimals)
    if rounded_2sf >= 10.0:
        return f"{round(av):d}"
    return f"{rounded_2sf:.{decimals}f}"


def _fmt_signed_pct(val: float) -> str:
    """Format a signed percentage with true minus sign (−), '+' on positive, no sign on 0."""
    v = float(val)
    mag = _fmt_pct_mag(v)
    if mag == "0":
        return "0 %"
    sign = "−" if v < 0 else "+"
    return f"{sign}{mag} %"


def _fmt_signed_pct_endpoint(val: float, span_zero: bool) -> str:
    v = float(val)
    mag = _fmt_pct_mag(v)
    if mag == "0":
        return "0"
    if v < 0:
        return f"−{mag}"
    return f"+{mag}" if span_zero else mag


def _fmt_sig(val: float, sig: int) -> str:
    """Format a scalar to `sig` significant digits without scientific notation."""
    v = float(val)
    av = abs(v)
    if av < 1e-15:
        return "0"
    sign = "−" if v < 0 else ""
    exp = math.floor(math.log10(av))
    decimals = sig - 1 - exp
    rounded = round(av, decimals)
    if rounded > 0:
        exp = math.floor(math.log10(rounded))
        decimals = sig - 1 - exp
        rounded = round(av, decimals)
    if decimals <= 0:
        return f"{sign}{int(rounded):d}"
    return f"{sign}{rounded:.{decimals}f}"


def _fmt_plain(val: float) -> str:
    """Format a number with up to 3 significant digits and no unit (strips trailing zeros)."""
    v = float(val)
    if abs(v) < 1e-15:
        return "0"
    sign = "−" if v < 0 else ""
    s = f"{abs(v):.3g}"
    if "e" in s or "E" in s:
        return _fmt_sig(v, 3)
    return f"{sign}{s}"


def _central_value(level_doc: dict[str, Any], output_key: str) -> float:
    central = level_doc.get("central") or {}
    if output_key in central and central[output_key] is not None:
        return float(central[output_key])
    if output_key in {"h2_uptake_change_pct", "co_uptake_change_pct"}:
        gas = "h2" if output_key.startswith("h2_") else "co"
        u_key = f"{gas}_uptake_umol_m2"
        c_key = f"{gas}_uptake_contrast_umol_m2"
        if u_key in central and c_key in central and central[u_key] is not None and central[c_key] is not None:
            uptake = float(central[u_key])
            contrast = float(central[c_key])
            base = uptake - contrast
            if abs(base) > 1e-12:
                return 100.0 * contrast / base
            if abs(contrast) <= 1e-12:
                return 0.0
    raise KeyError(f"Central value for {output_key!r} not available on level {level_doc.get('level')}")


def _ranked_st(out_doc: dict[str, Any], series_id: str, level_idx: int, output_key: str) -> list[tuple[str, float]]:
    if out_doc.get("mode") == "constant":
        raise ValueError(
            f"Sensitivity indices are undefined for constant output {series_id}/{level_idx}:{output_key}"
        )
    st = out_doc.get("ST") or {}
    if not st:
        raise ValueError(f"No ST indices on {series_id}/{level_idx}:{output_key}")
    return sorted(((k, float(v[0])) for k, v in st.items()), key=lambda kv: (-kv[1], kv[0]))


def _format_scalar(val: float, output_key: str, mod: str | None) -> str:
    if mod is None:
        if output_key.endswith("_pct"):
            return _fmt_signed_pct(val)
        return _fmt_plain(val)
    if mod == "absneg":
        return f"{_fmt_pct_mag(abs(val))} %"
    if mod == "abspos":
        return f"{_fmt_pct_mag(val)} %"
    if mod == "frac2pct":
        return f"{_fmt_pct_mag(100.0 * val)} %"
    if mod == "um":
        return f"{_fmt_sig(val, 3)} µmol m⁻²"
    if mod == "ug":
        return f"{_fmt_sig(val, 3)} µg C m⁻²"
    if mod == "g":
        return f"{_fmt_sig(val, 2)} g C m⁻² yr⁻¹"
    if mod == "idx":
        return f"{float(val):.2f}"
    if mod == "plain":
        return _fmt_plain(val)
    raise ValueError(f"Unknown slot modifier {mod!r}")


def _format_ci90(q05: float, q95: float, output_key: str, mod: str | None) -> str:
    lo = float(q05)
    hi = float(q95)
    if mod is None:
        if output_key.endswith("_pct"):
            span_zero = lo < 0.0 < hi
            lo_s = _fmt_signed_pct_endpoint(lo, span_zero)
            hi_s = _fmt_signed_pct_endpoint(hi, span_zero)
            return f"{lo_s} to {hi_s} %"
        return f"{_fmt_plain(lo)} to {_fmt_plain(hi)}"
    if mod == "absneg":
        mags = sorted([abs(lo), abs(hi)])
        return f"{_fmt_pct_mag(mags[0])} to {_fmt_pct_mag(mags[1])} %"
    if mod == "abspos":
        vals = sorted([lo, hi])
        return f"{_fmt_pct_mag(vals[0])} to {_fmt_pct_mag(vals[1])} %"
    if mod == "frac2pct":
        vals = sorted([100.0 * lo, 100.0 * hi])
        return f"{_fmt_pct_mag(vals[0])} to {_fmt_pct_mag(vals[1])} %"
    if mod == "um":
        return f"{_fmt_sig(lo, 3)} to {_fmt_sig(hi, 3)} µmol m⁻²"
    if mod == "ug":
        return f"{_fmt_sig(lo, 3)} to {_fmt_sig(hi, 3)} µg C m⁻²"
    if mod == "g":
        return f"{_fmt_sig(lo, 2)} to {_fmt_sig(hi, 2)} g C m⁻² yr⁻¹"
    raise ValueError(f"Unsupported modifier {mod!r} for ci90")


def _check_sign_guard(
    mod: str | None,
    level_doc: dict[str, Any],
    out_doc: dict[str, Any],
    output_key: str,
    expr: str,
) -> None:
    if mod not in {"absneg", "abspos"}:
        return
    c_val = _central_value(level_doc, output_key)
    med = out_doc.get("median")
    q05 = out_doc.get("q05")
    q95 = out_doc.get("q95")
    if med is None or q05 is None or q95 is None:
        raise ValueError(f"Missing quantile/median values for sign-checked slot {expr!r}")
    four = {"central": float(c_val), "median": float(med), "q05": float(q05), "q95": float(q95)}
    if mod == "absneg":
        bad = {k: v for k, v in four.items() if v > 0.0}
        if bad:
            raise ValueError(f"|absneg sign check failed for {expr!r}: positive values {bad}")
    elif mod == "abspos":
        bad = {k: v for k, v in four.items() if v < 0.0}
        if bad:
            raise ValueError(f"|abspos sign check failed for {expr!r}: negative values {bad}")


def _resolve_special(name: str, uq_doc: dict[str, Any]) -> str:
    levels = uq_doc["series"][0]["levels"]
    if name == "deduced_samples":
        samples = [
            str(lvl["x"])
            for lvl in levels
            if lvl.get("outputs", {}).get("co_uptake_umol_m2", {}).get("mode") == "deduced"
        ]
        if not samples:
            raise ValueError("No deduced samples found for {special:deduced_samples}")
        return ", ".join(samples)
    if name == "open_samples":
        samples = [
            str(lvl["x"])
            for lvl in levels
            if lvl.get("outputs", {}).get("co_uptake_umol_m2", {}).get("mode") == "variance"
        ]
        if not samples:
            raise ValueError("No open samples found for {special:open_samples}")
        return ", ".join(samples)
    if name == "refused_samples":
        samples = [
            str(lvl["x"])
            for lvl in levels
            if float(lvl.get("outputs", {}).get("co_uptake_umol_m2", {}).get("refusal_fraction", 0.0)) == 1.0
        ]
        if not samples:
            raise ValueError("No refused samples found for {special:refused_samples}")
        if len(samples) == 1:
            return samples[0]
        if len(samples) == 2:
            return f"{samples[0]} and {samples[1]}"
        return ", ".join(samples[:-1]) + f" and {samples[-1]}"
    raise ValueError(f"Unknown special slot {name!r}")


def resolve_slot(expr: str, uq_doc: dict[str, Any]) -> str:
    """Resolve a single `{...}` slot expression against `uq_doc` and registries."""
    if "|" in expr:
        core, mod = expr.split("|", 1)
    else:
        core, mod = expr, None

    if core.startswith("const:"):
        const_id = core.split(":", 1)[1]
        val = get_registry_constant(const_id)
        return _format_scalar(val, const_id, mod)

    if core.startswith("special:"):
        if mod is not None:
            raise ValueError(f"Modifier not supported on special slot {expr!r}")
        spec_name = core.split(":", 1)[1]
        return _resolve_special(spec_name, uq_doc)

    parts = core.split(":")
    if len(parts) != 3 or "/" not in parts[0]:
        raise ValueError(f"Malformed slot expression {expr!r}")
    series_level, output_key, stat = parts
    series_id, level_s = series_level.split("/", 1)
    level_idx = int(level_s)

    series_doc = next((s for s in uq_doc.get("series", []) if s.get("id") == series_id), None)
    if series_doc is None:
        raise KeyError(f"Series {series_id!r} not found in uq_doc for {uq_doc.get('slug')!r}")
    level_doc = next(
        (l for l in series_doc.get("levels", []) if int(l.get("level", -1)) == level_idx),
        None,
    )
    if level_doc is None:
        raise KeyError(f"Level {level_idx} not found in series {series_id!r}")
    out_doc = (level_doc.get("outputs") or {}).get(output_key)
    if out_doc is None and stat != "central":
        raise KeyError(f"Output {output_key!r} not found in {series_id}/{level_idx}")

    if mod in {"absneg", "abspos"}:
        if out_doc is None:
            raise KeyError(f"Output {output_key!r} required for sign check in {expr!r}")
        _check_sign_guard(mod, level_doc, out_doc, output_key, expr)

    if stat == "central":
        val = _central_value(level_doc, output_key)
        return _format_scalar(val, output_key, mod)
    if stat in {"median", "q05", "q95"}:
        raw = out_doc.get(stat)
        if raw is None:
            raise ValueError(f"Statistic {stat!r} is None for {expr!r}")
        return _format_scalar(float(raw), output_key, mod)
    if stat == "ci90":
        q05 = out_doc.get("q05")
        q95 = out_doc.get("q95")
        if q05 is None or q95 is None:
            raise ValueError(f"ci90 endpoints missing for {expr!r}")
        return _format_ci90(float(q05), float(q95), output_key, mod)
    if stat in {"st1", "st2"}:
        ranked = _ranked_st(out_doc, series_id, level_idx, output_key)
        idx = 0 if stat == "st1" else 1
        return _format_scalar(ranked[idx][1], output_key, mod or "idx")
    if stat in {"st1name", "st2name"}:
        if mod is not None:
            raise ValueError(f"Modifier not allowed on {stat} in {expr!r}")
        ranked = _ranked_st(out_doc, series_id, level_idx, output_key)
        idx = 0 if stat == "st1name" else 1
        return readable_input_name(ranked[idx][0])
    if stat == "assum":
        if out_doc.get("mode") == "constant":
            raise ValueError(f"group_ST is undefined for constant output in {expr!r}")
        gst = out_doc.get("group_ST") or {}
        if "modeller_assumption" not in gst:
            raise ValueError(f"modeller_assumption missing from group_ST in {expr!r}")
        return _format_scalar(float(gst["modeller_assumption"][0]), output_key, mod or "idx")
    if stat == "refusal":
        rf = float(out_doc.get("refusal_fraction", 0.0)) * 100.0
        return f"{_fmt_pct_mag(rf)} %"

    raise ValueError(f"Unknown slot statistic {stat!r} in {expr!r}")


def render_template(template: str, uq_doc: dict[str, Any]) -> str:
    """Replace all `{...}` slots in `template` using `uq_doc`."""
    def _repl(match: re.Match[str]) -> str:
        return resolve_slot(match.group(1), uq_doc)

    rendered = _SLOT_RE.sub(_repl, template)
    if "{" in rendered or "}" in rendered:
        raise ValueError(f"Unresolved braces remain in rendered text: {rendered!r}")
    return rendered


def render_science_texts(slug: str, uq_doc: dict[str, Any]) -> dict[str, str]:
    """Render `title_science`, `answer_science`, `limits_science` and `headline_science` for `slug`."""
    if slug not in SCIENCE_TEMPLATES:
        raise KeyError(f"Unknown question slug {slug!r}")
    tmpl = SCIENCE_TEMPLATES[slug]
    title_science = render_template(tmpl["title_science"], uq_doc)
    answer_science = render_template(tmpl["answer_science"], uq_doc)
    limits_science = render_template(tmpl["limits_science"], uq_doc)
    first_sentence = answer_science.split(". ", 1)[0]
    if not first_sentence.endswith("."):
        first_sentence += "."
    return {
        "title_science": title_science,
        "answer_science": answer_science,
        "limits_science": limits_science,
        "headline_science": first_sentence,
    }
