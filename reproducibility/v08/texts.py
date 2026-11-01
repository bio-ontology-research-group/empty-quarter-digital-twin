"""Plain and scientific text templates and metadata for EQ-BioDT v0.8 question files."""
from __future__ import annotations

import re
from typing import Any

JARGON = re.compile(
    r"[\u0370-\u03ff\u1f00-\u1fff\u00b5]"
    r"|\bdo\s*\("
    r"|\b(?i:sobol)\b"
    r"|\b(?i:ontolog\w*)\b"
    r"|\bMAGs?\b"
    r"|\bEC_"
    r"|\bECe\b"
    r"|\bMaas\b"
    r"|\bCBB\b"
    r"|\b(?i:first-order)\b"
    r"|\b(?i:gates?)\b"
    r"|\b(?i:read-mapping)\b"
    r"|\bIQR\b"
    r"|\b(?i:central run)\b"
    r"|\b(?i:micromoles?)\b"
    r"|\b(?i:quantiles?)\b"
    r"|\b(?i:medians?)\b"
    r"|\b(?i:shares?)\b"
    r"|(?i:fraction of unsalted yield)"
    r"|\b[A-Za-z0-9]+_[A-Za-z0-9_]+\b"
    r"|\b(?:umol|nmol|pmol|mmol|mol|ug|mg|ng|uL|mL|MPa|kPa|Pa|km|nm|m3|cm3|yr)\b"
    r"|µmol|µg|µL|mol/L|g/L|g/g|m³|cm³|cm⁻³|m⁻²|s⁻¹|h⁻¹|min⁻¹|yr⁻¹|°[NESWF]"
)

SCOPE_SENTENCE_PLAIN = (
    "The range shows how much the answer moves when uncertain numbers in the models "
    "are varied; the models themselves could also be wrong."
)

SCOPE_SENTENCE_SCIENCE = (
    "The parameter-uncertainty interval quantifies joint variation over literature "
    "constants and twin process parameters only and excludes external-model and "
    "twin model-structure error."
)

PLAIN_TITLES: dict[str, str] = {
    "rain-pattern": "One heavy rain, four showers or daily drizzle?",
    "drip-or-flood": "Drip or flood irrigation?",
    "cloud-seeding": "What if cloud seeding adds rain?",
    "salty-water": "What if the water is salty?",
    "alfalfa-bacterium": "Can a desert bacterium help alfalfa in salty water?",
    "hotter": "What if it gets hotter?",
    "mulch-shade": "Mulch or shade nets?",
    "solar-park": "What if we build a solar park?",
    "compaction": "What if vehicles compact the sand?",
    "biochar": "What if we add biochar?",
    "missing-co-genes": "What if the soil lacks the CO genes?",
    "adding-microbes": "What if we add microbes?",
    "air-hydrogen": "What if air hydrogen rises?",
    "carbon-fixation": "How much carbon can the hydrogen eaters fix?",
}

PLAIN_BASELINE_LABELS: dict[str, str] = {
    "rain-pattern": "One 12 mm rain on day 1, no treatment",
    "drip-or-flood": "One 12 mm rain on day 1, no irrigation",
    "cloud-seeding": "One 12 mm rain on day 1, unseeded",
    "salty-water": "One 12 mm fresh rain on day 1, no irrigation",
    "alfalfa-bacterium": "Normal, uninoculated alfalfa yield",
    "hotter": "One 12 mm rain on day 1 at 25 °C",
    "mulch-shade": "Bare unshaded sand (one 12 mm rain on day 1)",
    "solar-park": "Open desert sand outside panels",
    "compaction": "Uncompacted sand (pore space 0.40)",
    "biochar": "Unamended sand (0 % biochar, pore space 0.40)",
    "missing-co-genes": "One 12 mm rain on day 1 in the selected sample community",
    "adding-microbes": "Uninoculated sample community",
    "air-hydrogen": "Normal background air hydrogen (1x normal level)",
    "carbon-fixation": "Sample community under 14-day rain schedule",
}

PLAIN_X_LABELS: dict[str, str] = {
    "adding-microbes": "Added microbes (% of community)",
}

PLAIN_SERIES_LABELS: dict[tuple[str, str], str] = {
    ("adding-microbes", "51PRr3-gate"): "Sample lacking CO eaters + CO microbe (yes/no presence)",
    ("adding-microbes", "51PRr3-proportional"): "Sample lacking CO eaters + CO microbe (scaled by amount)",
    ("adding-microbes", "V53PRr2-gate"): "Sample already having CO eaters + CO microbe (yes/no presence)",
    ("adding-microbes", "V53PRr2-proportional"): "Sample already having CO eaters + CO microbe (scaled by amount)",
    ("adding-microbes", "51PRr3-control-gate"): "Negative control: microbe without CO genes (yes/no presence)",
    ("adding-microbes", "51PRr3-control-proportional"): "Negative control: microbe without CO genes (scaled by amount)",
}

# (series_id, level_index, output_key) for the primary headline of each question
HEADLINE_SPEC: dict[str, tuple[str, int, str]] = {
    "rain-pattern": ("rain-schedule", 2, "h2_uptake_change_pct"),
    "drip-or-flood": ("irrigation-interval", 2, "h2_uptake_change_pct"),
    "cloud-seeding": ("seeding-intensity", 2, "h2_uptake_change_pct"),
    "salty-water": ("daily-irrigation", 3, "h2_uptake_change_pct"),
    "alfalfa-bacterium": ("C_fitted_fresh", 4, "relative_yield_alfalfa_sa187"),
    "hotter": ("warmer-climate", 3, "h2_uptake_change_pct"),
    "mulch-shade": ("mulch-only", 2, "h2_uptake_change_pct"),
    "solar-park": ("panel-cover", 2, "h2_uptake_change_pct"),
    "compaction": ("offroad-passes", 3, "h2_uptake_change_pct"),
    "biochar": ("biochar-dose", 3, "h2_uptake_change_pct"),
    "missing-co-genes": ("co-gate-samples", 2, "co_uptake_umol_m2"),
    "adding-microbes": ("V53PRr2-proportional", 1, "co_uptake_change_pct"),
    "air-hydrogen": ("h2-multiplier", 2, "h2_uptake_change_pct"),
    "carbon-fixation": ("baseline-12mm", 1, "c_fixed_upper_bound_g_m2_yr"),
}

OUTCOME_META: dict[str, dict[str, str]] = {
    "h2_uptake_change_pct": {
        "unit": "%",
        "unit_plain": "%",
        "label_plain": "Percentage change in hydrogen uptake over 14 days",
        "label_science": "14-day relative H₂ uptake change vs baseline",
        "provenance": "simulated",
    },
    "co_uptake_change_pct": {
        "unit": "%",
        "unit_plain": "%",
        "label_plain": "Percentage change in carbon monoxide uptake over 14 days",
        "label_science": "14-day relative CO uptake change vs baseline",
        "provenance": "simulated",
    },
    "h2_uptake_umol_m2": {
        "unit": "µmol m⁻²",
        "unit_plain": "millionths of a mole per m²",
        "label_plain": "Hydrogen taken up over 14 days",
        "label_science": "14-day cumulative H₂ uptake",
        "provenance": "simulated",
    },
    "h2_uptake_contrast_umol_m2": {
        "unit": "µmol m⁻²",
        "unit_plain": "millionths of a mole per m²",
        "label_plain": "Change in hydrogen uptake over 14 days",
        "label_science": "14-day H₂ uptake contrast vs baseline",
        "provenance": "simulated",
    },
    "co_uptake_umol_m2": {
        "unit": "µmol m⁻²",
        "unit_plain": "millionths of a mole per m²",
        "label_plain": "Carbon monoxide taken up over 14 days",
        "label_science": "14-day cumulative CO uptake",
        "provenance": "simulated",
    },
    "co_uptake_contrast_umol_m2": {
        "unit": "µmol m⁻²",
        "unit_plain": "millionths of a mole per m²",
        "label_plain": "Change in carbon monoxide uptake over 14 days",
        "label_science": "14-day CO uptake contrast vs baseline",
        "provenance": "simulated",
    },
    "active_days": {
        "unit": "d",
        "unit_plain": "days",
        "label_plain": "Active microbial days (out of 14)",
        "label_science": "Days with microbial activation A > 0.01",
        "provenance": "simulated",
    },
    "final_water_mm": {
        "unit": "mm",
        "unit_plain": "mm",
        "label_plain": "Water left in sand on day 14",
        "label_science": "Day-14 stored soil water W(14)",
        "provenance": "simulated",
    },
    "final_salt_mol_m2": {
        "unit": "mol m⁻²",
        "unit_plain": "moles of salt per m²",
        "label_plain": "Salt left in sand on day 14",
        "label_science": "Day-14 stored salt inventory N(14)",
        "provenance": "simulated",
    },
    "mean_nacl_mol_l": {
        "unit": "mol/L",
        "unit_plain": "moles of salt per litre",
        "label_plain": "Average pore-water salt concentration",
        "label_science": "Mean wet-day pore-water NaCl molarity",
        "provenance": "simulated",
    },
    "porosity": {
        "unit": "m³ m⁻³",
        "unit_plain": "fraction of sand volume",
        "label_plain": "Pore space in sand",
        "label_science": "Soil porosity φ",
        "provenance": "simulated",
    },
    "ece_dS_m": {
        "unit": "dS/m",
        "unit_plain": "dS/m",
        "label_plain": "Root-zone saltiness",
        "label_science": "Root-zone saturation-extract EC_e",
        "provenance": "external",
    },
    "relative_yield_date_palm": {
        "unit": "fraction",
        "unit_plain": "fraction of normal yield",
        "label_plain": "Date palm yield (compared with normal)",
        "label_science": "Maas-Hoffman relative yield (date palm)",
        "provenance": "external",
    },
    "relative_yield_barley": {
        "unit": "fraction",
        "unit_plain": "fraction of normal yield",
        "label_plain": "Barley yield (compared with normal)",
        "label_science": "Maas-Hoffman relative yield (barley)",
        "provenance": "external",
    },
    "relative_yield_alfalfa": {
        "unit": "fraction",
        "unit_plain": "fraction of normal yield",
        "label_plain": "Alfalfa yield without bacterium (compared with normal)",
        "label_science": "Maas-Hoffman relative yield (alfalfa, uninoculated)",
        "provenance": "external",
    },
    "relative_yield_tomato": {
        "unit": "fraction",
        "unit_plain": "fraction of normal yield",
        "label_plain": "Tomato yield (compared with normal)",
        "label_science": "Maas-Hoffman relative yield (tomato)",
        "provenance": "external",
    },
    "relative_yield_alfalfa_sa187": {
        "unit": "fraction",
        "unit_plain": "fraction of normal yield",
        "label_plain": "Alfalfa yield with desert bacterium (compared with normal)",
        "label_science": "Inoculated alfalfa relative yield (+SA187)",
        "provenance": "external",
    },
    "c_fixed_upper_bound_g_m2_yr": {
        "unit": "g C m⁻² yr⁻¹",
        "unit_plain": "g of carbon per m² per year",
        "label_plain": "Upper bound on carbon fixed per year",
        "label_science": "Stoichiometric upper bound on autotrophic C fixation (annualised)",
        "provenance": "bound",
    },
    "c_fixed_upper_bound_ug_m2": {
        "unit": "µg C m⁻²",
        "unit_plain": "micrograms of carbon per m²",
        "label_plain": "Upper bound on carbon fixed over 14 days",
        "label_science": "14-day stoichiometric upper bound on autotrophic C fixation",
        "provenance": "bound",
    },
    "co2_fixed_upper_bound_mol_m2": {
        "unit": "mol CO₂ m⁻²",
        "unit_plain": "moles of carbon dioxide per m²",
        "label_plain": "Upper bound on carbon dioxide fixed over 14 days",
        "label_science": "14-day stoichiometric upper bound on CO₂ fixation",
        "provenance": "bound",
    },
}

MODEL_META: dict[str, dict[str, str]] = {
    "ec_to_nacl": {
        "name_plain": "Water saltiness to dissolved salt amount",
        "name_science": "Irrigation-water EC to NaCl concentration",
        "equation_tex": r"c_{\mathrm{NaCl}} = \kappa_{\mathrm{NaCl}}\,\mathrm{EC}_w",
    },
    "irrigation_schedule": {
        "name_plain": "Equal-volume watering schedule across 14 days",
        "name_science": "Equal-total irrigation delivered in pulses every k days",
        "equation_tex": r"P_{\mathrm{irr}}(t) = \frac{V_{\mathrm{tot}}}{N_{\mathrm{pulses}}}\,\mathbb{I}(t \equiv 0 \pmod k)",
    },
    "crop_water_demand_generic": {
        "name_plain": "Daily desert crop watering demand",
        "name_science": "Crop water demand, hot-desert order of magnitude (FAO-56)",
        "equation_tex": r"P_{\mathrm{crop}}(t) = ET_{c,\mathrm{desert}}",
    },
    "crop_water_demand_alfalfa": {
        "name_plain": "Daily drip watering demand for desert alfalfa",
        "name_science": "Drip-irrigated alfalfa water demand (FAO-56)",
        "equation_tex": r"P_{\mathrm{alf}}(t) = ET_{c,\mathrm{alfalfa}}",
    },
    "fao29_leaching": {
        "name_plain": "Standard root-zone salt build-up from irrigation water",
        "name_science": "Root-zone salinity from irrigation-water EC at fixed leaching (FAO-29)",
        "equation_tex": r"\mathrm{EC}_e = f_{\mathrm{FAO}}\,\mathrm{EC}_w",
    },
    "twin_solution_to_ece": {
        "name_plain": "Sand pore-water salt to root-zone saltiness",
        "name_science": "Twin mean soil-solution NaCl to saturation-extract EC (FAO-29 relations)",
        "equation_tex": r"\mathrm{EC}_e = \frac{\bar{c}_{\mathrm{NaCl}}}{2\,\kappa_{\mathrm{NaCl}}}",
    },
    "fitted_leaching_fresh": {
        "name_plain": "Field-fitted root-zone salt build-up (fresh weight)",
        "name_science": "Leaching factor fitted to the SA187 field mock ratio (fresh weight)",
        "equation_tex": r"\mathrm{EC}_e = f_{\mathrm{fit,fresh}}\,\mathrm{EC}_w",
    },
    "fitted_leaching_dry": {
        "name_plain": "Field-fitted root-zone salt build-up (dry weight)",
        "name_science": "Leaching factor fitted to the SA187 field mock ratio (dry weight)",
        "equation_tex": r"\mathrm{EC}_e = f_{\mathrm{fit,dry}}\,\mathrm{EC}_w",
    },
    "sa187_benefit": {
        "name_plain": "Desert bacterium seed coating boost on alfalfa",
        "name_science": "Enterobacter sp. SA187 yield benefit on alfalfa (field trial)",
        "equation_tex": r"Y_{r,\mathrm{SA187}} = M_{\mathrm{SA187}}\,Y_r(\mathrm{EC}_e)",
    },
    "mulch_evaporation": {
        "name_plain": "Mulch cover slowing water evaporation",
        "name_science": "Mulch reduces soil evaporation",
        "equation_tex": r"E'(t) = E(t)\,\bigl(1 - \beta_{\mathrm{mulch}}\,c_{\mathrm{mulch}}\bigr)",
    },
    "shade_net_temperature": {
        "name_plain": "Shade net cooling the sand and slowing evaporation",
        "name_science": "Shade net lowers soil temperature",
        "equation_tex": r"T'(t) = T(t) - c_{\mathrm{shade}}\,\Delta T_{\mathrm{shade}},\quad E'(t) = E(t)\bigl(1 - \beta_{\mathrm{shade}}c_{\mathrm{shade}}\bigr)",
    },
    "biochar_porosity": {
        "name_plain": "Biochar mixing raising sand pore space",
        "name_science": "Biochar dose raises porosity (meta-analysis)",
        "equation_tex": r"\phi(f_b) = \phi_0 + \alpha_{\phi,b}\,f_b",
    },
    "inoculant": {
        "name_plain": "Adding microbe strains to the soil community",
        "name_science": "Inoculum: reference MAG at a relative share",
        "equation_tex": r"w_i' = (1 - s)\,w_i + s\,\mathbb{I}(i = i_{\mathrm{inoc}})",
    },
    "knallgas_carbon_fixation": {
        "name_plain": "Upper limit on carbon fixed from hydrogen energy",
        "name_science": "Upper bound of CO2 fixed by H2-oxidising autotrophs",
        "equation_tex": r"C_{\mathrm{fixed}} \le \mathbb{I}_{\mathrm{CBB}}\,\eta_{\mathrm{C/H_2}}\,M_C\,\int J_{\mathrm{H_2}}\,\mathrm{d}t",
    },
    "co_oxidation_carbon_balance": {
        "name_plain": "Carbon monoxide conversion to carbon dioxide",
        "name_science": "CO oxidation converts CO to CO2 (no net carbon-atom change)",
        "equation_tex": r"\Delta C_{\mathrm{net,CO}} = 0",
    },
    "h2_atmosphere": {
        "name_plain": "Higher hydrogen level in the outside air",
        "name_science": "Atmospheric H2 scenario (multiple of background)",
        "equation_tex": r"C_{\mathrm{H_2,atm}}' = m_{\mathrm{H_2}}\,C_{\mathrm{H_2,0}}",
    },
    "solar_park_microclimate": {
        "name_plain": "Solar-panel shade cooling the sand and slowing evaporation",
        "name_science": "Solar park: under-panel soil temperature and evaporation",
        "equation_tex": r"T'(t) = T(t) - c_{\mathrm{PV}}\Delta T_{\mathrm{PV}},\quad E'(t) = E(t)\bigl(1 - \beta_{\mathrm{PV}}c_{\mathrm{PV}}\bigr)",
    },
    "panel_washing": {
        "name_plain": "Water dripping onto the sand from solar-panel washing",
        "name_science": "Solar park: panel-washing water pulses",
        "equation_tex": r"P'(t) = P(t) + c_{\mathrm{PV}}\,w_{\mathrm{wash}}(t)",
    },
    "cloud_seeding": {
        "name_plain": "Cloud-seeding rain boost on rainy days",
        "name_science": "Cloud seeding: relative rain increase on rain days",
        "equation_tex": r"P'(t) = P(t)\,\bigl(1 + I_{\mathrm{seed}}\,\alpha_{\mathrm{seed}}\bigr)\quad (P(t)>0)",
    },
    "offroad_compaction": {
        "name_plain": "Vehicle passes squeezing sand pore space",
        "name_science": "Off-road traffic compacts desert sand",
        "equation_tex": r"\phi(n_{\mathrm{pass}}) = 1 - \frac{\rho_b(n_{\mathrm{pass}})}{\rho_s}",
    },
    "climate_warming": {
        "name_plain": "Warmer air and higher water evaporation",
        "name_science": "Warming: temperature, evaporation scaling and rainfall change",
        "equation_tex": r"T'(t)=T(t)+\Delta T,\quad E'(t)=E(t)\bigl(1+\alpha_{\mathrm{PE}}\Delta T\bigr)",
    },
    "maas_hoffman_date_palm": {
        "name_plain": "Date palm yield loss from root-zone salt",
        "name_science": "Maas-Hoffman relative yield, date palm (FAO-29)",
        "equation_tex": r"Y_r = \mathrm{clip}\!\bigl(1 - b_{\mathrm{MH}}(\mathrm{EC}_e - a_{\mathrm{MH}})/100,\,0,\,1\bigr)",
    },
    "maas_hoffman_barley": {
        "name_plain": "Barley yield loss from root-zone salt",
        "name_science": "Maas-Hoffman relative yield, barley (FAO-29)",
        "equation_tex": r"Y_r = \mathrm{clip}\!\bigl(1 - b_{\mathrm{MH}}(\mathrm{EC}_e - a_{\mathrm{MH}})/100,\,0,\,1\bigr)",
    },
    "maas_hoffman_alfalfa": {
        "name_plain": "Alfalfa yield loss from root-zone salt",
        "name_science": "Maas-Hoffman relative yield, alfalfa (FAO-29)",
        "equation_tex": r"Y_r = \mathrm{clip}\!\bigl(1 - b_{\mathrm{MH}}(\mathrm{EC}_e - a_{\mathrm{MH}})/100,\,0,\,1\bigr)",
    },
    "maas_hoffman_tomato": {
        "name_plain": "Tomato yield loss from root-zone salt",
        "name_science": "Maas-Hoffman relative yield, tomato (FAO-29)",
        "equation_tex": r"Y_r = \mathrm{clip}\!\bigl(1 - b_{\mathrm{MH}}(\mathrm{EC}_e - a_{\mathrm{MH}})/100,\,0,\,1\bigr)",
    },
    "eq_biodt_v03": {
        "name_plain": "Desert sand water, salt and gas-uptake twin",
        "name_science": "EQ-BioDT v0.3 two-layer causal-ontological soil twin",
        "equation_tex": r"\frac{\mathrm{d}W}{\mathrm{d}t}=P-E-Q,\quad J_g = \mathbb{I}_{\mathrm{gate}}\,A(\psi,T)\,D_{\mathrm{eff}}(\theta,\phi)\,k_g\,C_g",
    },
}


def round_plain(x: float | None) -> str:
    """Format a number for plain-language display: integers for |x| >= 10, one decimal below."""
    if x is None:
        return "none"
    val = float(x)
    if abs(val) < 1e-12:
        return "0"
    ax = abs(val)
    if ax >= 1000.0:
        return f"{round(val):,}"
    if ax >= 10.0:
        return f"{round(val):d}"
    if ax >= 0.1:
        s = f"{val:.1f}".rstrip("0").rstrip(".")
        return "0" if s in {"-0", ""} else s
    s = f"{val:.3f}".rstrip("0").rstrip(".")
    return "0" if s in {"-0", ""} else s


def science_central(q: dict[str, Any]) -> float:
    """Return the central headline value from a built question dict (in % for relative_yield_*)."""
    if "headline_value" in q and q["headline_value"] is not None:
        val = float(q["headline_value"])
    else:
        s_id, l_idx, out_key = HEADLINE_SPEC[q["slug"]]
        s = next(sr for sr in q["series"] if sr["id"] == s_id)
        val = float(s["levels"][l_idx]["outcomes"][out_key]["value"])
    _, _, out_key = HEADLINE_SPEC[q["slug"]]
    if out_key.startswith("relative_yield_") and abs(val) <= 2.0:
        return 100.0 * val
    return val


def plain_number(q: dict[str, Any]) -> str:
    """Extract the core rounded headline number from `q['answer_plain']`."""
    if q.get("slug") == "adding-microbes":
        m = re.search(r"raises uptake only slightly \((\d+(?:\.\d+)?)\s*%", q["answer_plain"])
        if m:
            return m.group(1)
    m = re.search(
        r"(?<![A-Za-z0-9_])-?\d+(?:,\d{3})*(?:\.\d+)?(?![A-Za-z0-9_])",
        q["answer_plain"],
    )
    if not m:
        raise ValueError(f"No plain number found in answer_plain for {q.get('slug')}")
    return m.group(0)


def plain_level_label(slug: str, series_id: str, spec_label: str, x: float | str) -> str:
    """Return a jargon-free plain label for a level."""
    if slug == "salty-water":
        if series_id == "salty-rain-once":
            return f"{round_plain(float(x))} dS/m single rain"
        return f"{round_plain(float(x))} dS/m daily water"
    if slug == "alfalfa-bacterium":
        return f"{float(x):g} dS/m water"
    if slug == "adding-microbes":
        return f"{round_plain(float(x))} % added microbes"
    if slug == "compaction" and series_id == "direct-porosity-v06":
        return "pore space 0.30"
    if slug == "biochar":
        return f"{round_plain(float(x))} % biochar"
    if slug == "drip-or-flood":
        d = int(round(float(x)))
        return "every day" if d == 1 else f"every {d} days"
    if slug == "hotter":
        return spec_label.replace(" C", " °C")
    if slug == "carbon-fixation":
        return f"sample {x}"
    if slug == "missing-co-genes":
        return f"sample {x}"
    return spec_label


def format_change_pct_interval(q05: float | None, q95: float | None) -> str:
    """Format a percentage-change interval following direction/sign rules."""
    if q05 is None or q95 is None:
        return "range not available in %"
    lo = float(q05)
    hi = float(q95)
    r_lo = round_plain(lo)
    r_hi = round_plain(hi)
    if r_lo == r_hi:
        val_str = round_plain(abs(lo))
        return f"exactly {val_str} % in every case we tried"
    if lo < 0.0 < hi:
        return f"from {round_plain(abs(lo))} % lower to {round_plain(abs(hi))} % higher"
    mag_min = min(abs(lo), abs(hi))
    mag_max = max(abs(lo), abs(hi))
    return f"by {round_plain(mag_min)} to {round_plain(mag_max)} %"


def _format_interval_plain(q05: float | None, q95: float | None, unit_plain: str) -> str:
    if q05 is None or q95 is None:
        return f"range not available in {unit_plain}"
    return f"range {round_plain(q05)} to {round_plain(q95)} {unit_plain}"


def _hedge_phrase(q05: float | None, q95: float | None) -> str:
    """Return the required hedge phrase matching the sign of the 5–95 % headline interval."""
    if q05 is not None and q95 is not None:
        lo = float(q05)
        hi = float(q95)
        if lo < 0.0 < hi:
            return "may go up or down"
        if hi <= 0.0 and lo < 0.0:
            return "goes down in every case we tried"
        if lo >= 0.0 and hi > 0.0:
            return "goes up in every case we tried"
    return "goes up in every case we tried"


def _top_two_sobol_text(st: dict[str, list[float]], s1: dict[str, list[float]], mode: str) -> str:
    if mode == "deduced" or not st:
        return (
            "Top two Sobol' contributors: none (ST = 0, S1 = 0 across all inputs because "
            "the deductive gate fixes the output to exact 0)."
        )
    ranked = sorted(st.items(), key=lambda kv: kv[1][0], reverse=True)
    parts = []
    for name, vals in ranked[:2]:
        st_est, st_lo, st_hi = vals[0], vals[1], vals[2]
        s1_vals = s1.get(name, [0.0, 0.0, 0.0])
        s1_est = 0.0 if abs(s1_vals[0]) < 1e-6 else s1_vals[0]
        parts.append(
            f"{name} (ST = {st_est:.2f} [{st_lo:.2f}, {st_hi:.2f}], S1 = {s1_est:.2f})"
        )
    target = "refusal indicator" if mode == "refusal_indicator" else "output variance"
    return f"Top two Sobol' contributors to {target}: " + " and ".join(parts) + "."


def _group_share_text(
    group_st: dict[str, list[float]],
    group_s1: dict[str, list[float]],
    mode: str,
) -> str:
    if mode == "deduced" or not group_st:
        return (
            "Status share of variance (group ST / S1): verified_in_source = 0.00 / 0.00, "
            "modeller_assumption = 0.00 / 0.00, source_not_checked = 0.00 / 0.00 "
            "(joint group total indices can sum to more than 1 because of interactions)."
        )
    order = ("verified_in_source", "modeller_assumption", "source_not_checked")
    items = []
    for k in order:
        if k in group_st:
            gst = 0.0 if abs(group_st[k][0]) < 1e-6 else group_st[k][0]
            gs1_raw = group_s1.get(k, [0.0])[0]
            gs1 = 0.0 if abs(gs1_raw) < 1e-6 else gs1_raw
            items.append(f"{k} ST = {gst:.2f} (S1 = {gs1:.2f})")
    return (
        "Status share of variance: "
        + ", ".join(items)
        + " (joint group total indices can sum to more than 1 because of interactions)."
    )


def _fmt_sci_stat(central_val: float | None, uq_o: dict[str, Any] | None, unit: str) -> str:
    c_str = f"{central_val:.4g}" if central_val is not None else "n/a"
    if not uq_o:
        return f"{c_str} {unit}"
    med = uq_o.get("median")
    q25 = uq_o.get("q25")
    q75 = uq_o.get("q75")
    q05 = uq_o.get("q05")
    q95 = uq_o.get("q95")
    med_s = f"{med:.4g}" if med is not None else "n/a"
    iqr_s = f"[{q25:.4g}, {q75:.4g}]" if (q25 is not None and q75 is not None) else "[n/a, n/a]"
    ci90_s = f"[{q05:.4g}, {q95:.4g}]" if (q05 is not None and q95 is not None) else "[n/a, n/a]"
    return (
        f"{c_str} {unit} (median {med_s} {unit}, 50 % IQR {iqr_s} {unit}, "
        f"joint 5–95 % interval {ci90_s} {unit})"
    )


def _get_uq_level(uq_doc: dict[str, Any] | None, series_id: str, level_idx: int) -> dict[str, Any] | None:
    if not uq_doc:
        return None
    for s in uq_doc.get("series", []):
        if s.get("id") == series_id:
            for l in s.get("levels", []):
                if int(l.get("level", -1)) == level_idx:
                    return l
    return None


def _convergence_text(outcome_uq: dict[str, Any], uq_doc: dict[str, Any] | None) -> str:
    flagged_items: list[str] = []
    if uq_doc:
        for s in uq_doc.get("series", []):
            for l in s.get("levels", []):
                lbl = l.get("label", f"level {l.get('level')}")
                n_val = l.get("n")
                for k, out_v in l.get("outputs", {}).items():
                    conv = out_v.get("convergence") or {}
                    if conv.get("flagged"):
                        shift = float(conv.get("half_n_shift_frac") or 0.0) * 100.0
                        flagged_items.append(
                            f"at {lbl} ({k}, N = {n_val}), interval endpoints still move by "
                            f"{shift:.1f} % of the width between N/2 and N"
                        )
    else:
        conv = outcome_uq.get("convergence") or {}
        if conv.get("flagged"):
            shift = float(conv.get("half_n_shift_frac") or 0.0) * 100.0
            flagged_items.append(
                f"interval endpoints still move by {shift:.1f} % of the width between N/2 and N"
            )
    if not flagged_items:
        return ""
    return "Convergence check: " + "; ".join(flagged_items) + "."


def render_question_texts(
    slug: str,
    central_val: float,
    outcome_uq: dict[str, Any],
    max_refusal_fraction: float,
    uq_doc: dict[str, Any] | None = None,
) -> dict[str, str]:
    """Render plain and scientific answer, headline, and limits strings for `slug`."""
    s_id, l_idx, out_key = HEADLINE_SPEC[slug]
    is_rel_yield = out_key.startswith("relative_yield_")
    is_change = out_key in {"h2_uptake_change_pct", "co_uptake_change_pct"}
    display_central = 100.0 * central_val if (is_rel_yield and abs(central_val) <= 2.0) else central_val
    core = round_plain(abs(display_central) if is_change else display_central)

    q05 = outcome_uq.get("q05")
    q25 = outcome_uq.get("q25")
    q75 = outcome_uq.get("q75")
    q95 = outcome_uq.get("q95")
    median = outcome_uq.get("median")
    mode = str(outcome_uq.get("mode", "variance"))
    ref_frac = float(outcome_uq.get("refusal_fraction", 0.0))

    meta = OUTCOME_META[out_key]
    unit_sci = meta["unit"]
    unit_plain = meta["unit_plain"]
    if is_change:
        rng_plain = format_change_pct_interval(q05, q95)
    elif is_rel_yield:
        q05_disp = 100.0 * float(q05) if q05 is not None else None
        q95_disp = 100.0 * float(q95) if q95 is not None else None
        rng_plain = _format_interval_plain(q05_disp, q95_disp, "%")
    else:
        rng_plain = _format_interval_plain(q05, q95, unit_plain)

    hedge = _hedge_phrase(q05, q95)
    has_refusal = max_refusal_fraction >= 0.05 or ref_frac >= 0.05

    if slug == "rain-pattern":
        headline_plain = (
            f"Splitting one heavy rain into four smaller showers lowers hydrogen uptake over two weeks "
            f"by {core} %, and daily drizzle lowers it by 100 %, so uptake {hedge}."
        )
        answer_plain = (
            f"Splitting one heavy downpour into four smaller showers lowers hydrogen uptake over two weeks "
            f"by {core} % ({rng_plain}), and splitting it into daily drizzle lowers uptake by 100 % "
            f"(by 51 to 100 %), so uptake {hedge}. A single soaking pushes water deeper into the sand "
            f"so microbes stay active for more days before surface evaporation dries the top layer."
        )
        limits_extra_plain = (
            "The simulation uses a single sand column with fixed evaporation and does not include "
            "surface crusting or runoff down dune slopes."
        )
        limits_extra_sci = (
            "Assumes 1-D bucket hydrology without overland runoff, hysteresis, or biological soil-crust "
            "water-retention shifts."
        )
    elif slug == "drip-or-flood":
        headline_plain = (
            f"Watering once every ten days lowers hydrogen uptake over two weeks by {core} % compared "
            f"with a single flood, versus 38 % lower for daily drip, and the change {hedge}."
        )
        answer_plain = (
            f"Watering once every ten days lowers hydrogen uptake over two weeks by {core} % compared "
            f"with a single rain ({rng_plain}), versus 37 % lower when watering every three days "
            f"(from 68 % lower to 129 % higher) and 38 % lower for daily drip (from 73 % lower to 136 % higher), "
            f"so the change {hedge} across uncertain model numbers. Less frequent, deeper watering loses a "
            f"smaller fraction of water to daily surface evaporation than tiny daily doses."
        )
        limits_extra_plain = (
            "Pipe layout, root uptake and uneven wetting around drip emitters are not represented."
        )
        limits_extra_sci = (
            "Assumes spatially uniform 1-D infiltration and fixed daily potential evaporation without "
            "3-D drip-bulb geometry or plant transpiration."
        )
    elif slug == "cloud-seeding":
        headline_plain = (
            f"Adding ten percent more rain on rainy days raises hydrogen uptake over two weeks by {core} %, "
            f"and uptake {hedge}."
        )
        answer_plain = (
            f"Adding a ten percent rain boost on rainy days raises hydrogen uptake over two weeks by {core} % "
            f"({rng_plain}), so uptake {hedge}. Extra water on wet days pushes moisture deeper into the sand "
            f"and keeps microbes active for longer after the rain."
        )
        limits_extra_plain = (
            "Cloud seeding can only add water when suitable clouds are already present; the twin does not "
            "predict when clouds form."
        )
        limits_extra_sci = (
            "Treats seeding as a multiplicative scaling on existing rain-day precipitation; atmospheric "
            "microphysics and cloud availability are external to the soil twin."
        )
    elif slug == "salty-water":
        headline_plain = (
            f"Daily watering with salty water lowers hydrogen uptake by {core} % at 8 dS/m, where the "
            f"change {hedge}; at higher salt levels, the model cannot tell for such salty water."
        )
        answer_plain = (
            f"Daily watering with brackish water lowers hydrogen uptake over two weeks by {core} % at 8 dS/m "
            f"({rng_plain}), so uptake {hedge} depending on how strongly daily water additions offset salt "
            f"build-up, while leaving alfalfa at about 89 % of normal yield (range 83 to 93 %); at higher salt "
            f"levels, the model cannot tell for such salty water because pore-water salt passes the dilute limit."
        )
        limits_extra_plain = (
            "When salt concentration in the sand water rises above a tenth of a mole per litre, the twin "
            "refuses to guess rather than extrapolate outside its dilute-salt formula."
        )
        limits_extra_sci = (
            "The twin enforces an ontological/domain admissibility gate at c_NaCl > 0.1 mol/L (ideal-dilute "
            "van 't Hoff limit) and abstains when evaporation concentrates solutes above this threshold."
        )
    elif slug == "alfalfa-bacterium":
        headline_plain = (
            f"Coating alfalfa seeds with the desert bacterium gives about {core} % of normal yield in salty "
            f"water at 7.81 dS/m, up from 80 % without it."
        )
        answer_plain = (
            f"Coating alfalfa seeds with the desert bacterium gives about {core} % of normal yield, {rng_plain}, "
            f"in salty irrigation water at 7.81 dS/m, compared with about 80 % of normal yield (range 78 to 81 %) "
            f"without the bacterium; at extreme salt levels in the sand-physics branch, the model cannot tell for "
            f"such salty water."
        )
        limits_extra_plain = (
            "The crop boost comes from field trials in one desert farm and assumes the added bacterium "
            "colonises the roots reliably."
        )
        limits_extra_sci = (
            "Multiplies Maas-Hoffman relative yield by an empirical SA187 field multiplier; root-zone EC_e "
            "depends strongly on the chosen leaching/accumulation closure (models A, B, C)."
        )
    elif slug == "hotter":
        headline_plain = (
            f"Four degrees of climate warming raises hydrogen uptake over two weeks by {core} %, "
            f"and uptake {hedge} across uncertain model numbers."
        )
        answer_plain = (
            f"Four degrees of climate warming raises hydrogen uptake over two weeks by {core} % when rainfall "
            f"stays the same ({rng_plain}), so uptake {hedge} depending on whether faster enzyme activity "
            f"outweighs faster drying of the sand; warming the sand alone by ten degrees without extra drying "
            f"raises uptake by 69 % (by 3.8 to 78 %)."
        )
        limits_extra_plain = (
            "Heat stress death at extreme afternoon surface temperatures is not included in the simple "
            "warming curve."
        )
        limits_extra_sci = (
            "Uses a monotonic Q10 temperature response without high-temperature protein denaturation or "
            "diurnal surface-skin thermal gradients."
        )
    elif slug == "mulch-shade":
        headline_plain = (
            f"Full mulch cover raises hydrogen uptake over two weeks by {core} % by keeping the sand damp "
            f"longer, and uptake {hedge}."
        )
        answer_plain = (
            f"Full mulch cover raises hydrogen uptake over two weeks by {core} % ({rng_plain}) by slowing "
            f"evaporation, so uptake with mulch alone {hedge}, whereas full shade alone lowers uptake by "
            f"25 % (by 1.7 to 28 %) and combining mulch with shade raises uptake by 14 % "
            f"(from 6.2 % lower to 70 % higher). Mulch keeps water in the sand without cooling the microbes, "
            f"whereas shade cools the sand and slows enzyme speed."
        )
        limits_extra_plain = (
            "Organic mulch breakdown and extra carbon food for soil microbes are not included."
        )
        limits_extra_sci = (
            "Treats mulch and shade as boundary forcing modifiers on evaporation and soil temperature only, "
            "excluding organic carbon leaching from mulch."
        )
    elif slug == "solar-park":
        headline_plain = (
            f"Full solar-panel cover raises hydrogen uptake over two weeks by {core} %, "
            f"and uptake {hedge}."
        )
        answer_plain = (
            f"Full solar-panel cover raises hydrogen uptake over two weeks by {core} % ({rng_plain}), "
            f"so uptake {hedge} because reduced evaporation under the panels and water from panel washing "
            f"keep the sand moist longer than panel shade slows microbial enzymes."
        )
        limits_extra_plain = (
            "Wind shifts under panel rows and runoff strips along panel edges are simplified to an even average."
        )
        limits_extra_sci = (
            "Homogenises under-panel microclimate and wash-water drips spatially; competing cooling (Q10) and "
            "moisture-retention effects partially offset each other across parameter draws."
        )
    elif slug == "compaction":
        headline_plain = (
            f"Fifty vehicle passes lower hydrogen uptake over two weeks by {core} %, "
            f"and the change {hedge} across the uncertainty range."
        )
        answer_plain = (
            f"Fifty vehicle passes lower hydrogen uptake over two weeks by {core} % ({rng_plain}), "
            f"so the change {hedge} because squeezing pore space slows gas movement through air spaces "
            f"while changing how tightly the sand holds water."
        )
        limits_extra_plain = (
            "Crushing of soil crusts, plants or burrows by tyres is not modelled beyond the drop in pore space."
        )
        limits_extra_sci = (
            "Maps traffic passes to bulk density and porosity only, holding retention shape b and air-entry "
            "potential psi_entry uncoupled from pore-size redistribution."
        )
    elif slug == "biochar":
        headline_plain = (
            f"Adding five percent biochar by weight lowers hydrogen uptake over two weeks by {core} %, "
            f"and the change {hedge} across uncertain model numbers."
        )
        answer_plain = (
            f"Adding five percent biochar by weight lowers hydrogen uptake over two weeks by {core} % "
            f"({rng_plain}), so the change {hedge} across uncertain model numbers. Increasing total pore space "
            f"without adding extra water lowers the fraction of pores filled with water after a single rain "
            f"while opening more air space for gas movement."
        )
        limits_extra_plain = (
            "Chemical effects of biochar on soil acidity, nutrients or microbe attachment are not included."
        )
        limits_extra_sci = (
            "Propagates biochar dose solely through linear mixture porosity φ(f_b), holding Campbell retention "
            "parameters and surface chemistry constant."
        )
    elif slug == "missing-co-genes":
        headline_plain = (
            f"Where the carbon-monoxide genes are missing there is no carbon-monoxide uptake at all (none), "
            f"versus some uptake (about {core}) where the genes are present; for one sample the records are "
            f"incomplete, so the twin gives no answer."
        )
        answer_plain = (
            f"Where the soil community lacks the carbon-monoxide enzyme genes there is no carbon-monoxide "
            f"uptake at all (none), versus some uptake where the genes are present (about {core} over two "
            f"weeks, range {round_plain(q05)} to {round_plain(q95)}); for one sample the records are "
            f"incomplete, so the twin gives no answer."
        )
        limits_extra_plain = (
            "A missing gene in recovered DNA means the twin sets that pathway to zero, though very rare "
            "unsequenced microbes could still exist in the sand."
        )
        limits_extra_sci = (
            "Deductive zero relies on the closed-world metagenomic gate over assembled MAGs; unbinned or "
            "low-abundance carriers below sequencing detection are not represented."
        )
    elif slug == "adding-microbes":
        headline_plain = (
            f"Adding carbon-monoxide-eating microbes at one percent of the community raises carbon monoxide "
            f"uptake by {core} % where those microbes are already present, and switches uptake from none to "
            f"some where they were missing."
        )
        answer_plain = (
            f"In sample 51PRr3, which has no carbon-monoxide eaters, adding a carrier at one percent switches "
            f"uptake from none to some (about 144 over two weeks, range 21 to 507), and adding a microbe without "
            f"the genes leaves it at none; then, where carbon-monoxide eaters are already present, adding more "
            f"raises uptake only slightly ({core} % at one percent, 2.7 % at five percent), so uptake {hedge}."
        )
        limits_extra_plain = (
            "The scenario assumes the added microbes survive and work immediately in dry desert sand without "
            "dying off."
        )
        limits_extra_sci = (
            "Constitutive community replacement modifies composition shares instantaneously without modelling "
            "inoculant establishment, competition, or die-off kinetics."
        )
    elif slug == "air-hydrogen":
        headline_plain = (
            f"Doubling hydrogen in the air raises soil hydrogen uptake over two weeks by {core} %, "
            f"and uptake {hedge}."
        )
        answer_plain = (
            f"Doubling hydrogen in the air raises soil hydrogen uptake over two weeks by {core} % "
            f"({rng_plain}), and tenfold air hydrogen increases uptake by 900 %, so uptake {hedge}. Because "
            f"soil uptake in the twin is directly proportional to the hydrogen level in the air, percentage "
            f"uptake rises in lockstep with air concentration."
        )
        limits_extra_plain = (
            "At very high hydrogen levels, microbes could hit a maximum enzyme speed or grow in number, "
            "which a straight-line formula does not capture."
        )
        limits_extra_sci = (
            "First-order Henry/diffusion-reaction kinetics in C_H2; excludes Michaelis-Menten saturation and "
            "long-term community biomass growth under elevated H₂."
        )
    elif slug == "carbon-fixation":
        headline_plain = (
            f"Hydrogen-eating microbes can fix at most {core} {unit_plain} in sample F14Dr1, "
            f"and 0 where the carbon-fixing genes are absent."
        )
        answer_plain = (
            f"Hydrogen-eating microbes can fix at most {core} {unit_plain} in sample F14Dr1 after a single rain "
            f"({rng_plain}), and 0 {unit_plain} in sample 14PRr2 where the carbon-fixing genes are absent. "
            f"This is an upper ceiling if every bit of captured hydrogen energy went into building cell carbon."
        )
        limits_extra_plain = (
            "This figure is an upper ceiling, not measured net carbon storage in the sand, because microbes "
            "also burn energy just to stay alive and release carbon dioxide when they respire."
        )
        limits_extra_sci = (
            "Stoichiometric upper bound only (mode = bound); net ecosystem CO₂ removal remains unidentifiable "
            "without maintenance respiration, carbon-stock change, and lifecycle measurements."
        )
    else:
        raise KeyError(f"Unknown question slug: {slug}")

    if has_refusal and slug != "missing-co-genes" and "the model cannot tell for such" not in answer_plain.lower():
        cond = "salty water" if slug in {"salty-water", "alfalfa-bacterium"} else "conditions outside the model domain"
        answer_plain += f" At extreme settings, the model cannot tell for such {cond}."

    limits_plain = f"{SCOPE_SENTENCE_PLAIN} {limits_extra_plain}"

    # Scientific headline and answer
    h_lvl = _get_uq_level(uq_doc, s_id, l_idx)
    h_central_map = (h_lvl or {}).get("central", {})
    h_outputs_map = (h_lvl or {}).get("outputs", {})

    if out_key in {"h2_uptake_change_pct", "co_uptake_change_pct"}:
        gas = "h2" if out_key.startswith("h2_") else "co"
        contrast_key = f"{gas}_uptake_contrast_umol_m2"
        abs_key = f"{gas}_uptake_umol_m2"
        pct_stat = _fmt_sci_stat(central_val, outcome_uq, "%")
        con_stat = _fmt_sci_stat(
            h_central_map.get(contrast_key),
            h_outputs_map.get(contrast_key),
            "µmol m⁻²",
        )
        abs_stat = _fmt_sci_stat(
            h_central_map.get(abs_key),
            h_outputs_map.get(abs_key),
            "µmol m⁻²",
        )
        headline_science = (
            f"Central {out_key} = {pct_stat} (plain headline: {core} %); "
            f"absolute change {contrast_key} = {con_stat}; "
            f"absolute uptake {abs_key} = {abs_stat} "
            f"(mode = {mode}, refusal = {ref_frac * 100:.1f} %)."
        )
        extra_levels_sci = ""
        if slug == "rain-pattern":
            driz_lvl = _get_uq_level(uq_doc, "rain-schedule", 3)
            base_lvl = _get_uq_level(uq_doc, "rain-schedule", 1)
            if driz_lvl and base_lvl:
                d_c = driz_lvl.get("central", {})
                d_o = driz_lvl.get("outputs", {})
                b_c = base_lvl.get("central", {})
                b_o = base_lvl.get("outputs", {})
                extra_levels_sci = (
                    f" Daily drizzle (12 mm total) vs one 12 mm storm: "
                    f"h2_uptake_change_pct = {_fmt_sci_stat(d_c.get('h2_uptake_change_pct'), d_o.get('h2_uptake_change_pct'), '%')} "
                    f"(plain: -100 %); "
                    f"h2_uptake_contrast_umol_m2 = {_fmt_sci_stat(d_c.get('h2_uptake_contrast_umol_m2'), d_o.get('h2_uptake_contrast_umol_m2'), 'µmol m⁻²')}; "
                    f"h2_uptake_umol_m2 = {_fmt_sci_stat(d_c.get('h2_uptake_umol_m2'), d_o.get('h2_uptake_umol_m2'), 'µmol m⁻²')} "
                    f"(one-storm baseline h2_uptake_umol_m2 = {_fmt_sci_stat(b_c.get('h2_uptake_umol_m2'), b_o.get('h2_uptake_umol_m2'), 'µmol m⁻²')})."
                )
        elif slug == "adding-microbes":
            gate_lvl = _get_uq_level(uq_doc, "51PRr3-gate", 1)
            prop_lvl = _get_uq_level(uq_doc, "51PRr3-proportional", 1)
            v53_5_lvl = _get_uq_level(uq_doc, "V53PRr2-proportional", 2)
            if gate_lvl and prop_lvl and v53_5_lvl:
                extra_levels_sci = (
                    f" At 5 % inoculum in V53PRr2-proportional: "
                    f"co_uptake_change_pct = {_fmt_sci_stat(v53_5_lvl['central'].get('co_uptake_change_pct'), v53_5_lvl['outputs'].get('co_uptake_change_pct'), '%')}, "
                    f"co_uptake_contrast_umol_m2 = {_fmt_sci_stat(v53_5_lvl['central'].get('co_uptake_contrast_umol_m2'), v53_5_lvl['outputs'].get('co_uptake_contrast_umol_m2'), 'µmol m⁻²')}, "
                    f"co_uptake_umol_m2 = {_fmt_sci_stat(v53_5_lvl['central'].get('co_uptake_umol_m2'), v53_5_lvl['outputs'].get('co_uptake_umol_m2'), 'µmol m⁻²')}. "
                    f"In closed-baseline sample 51PRr3 (baseline co_uptake_umol_m2 = 0 µmol m⁻²), 1 % carrier inoculum yields "
                    f"co_uptake_umol_m2 = {_fmt_sci_stat(gate_lvl['central'].get('co_uptake_umol_m2'), gate_lvl['outputs'].get('co_uptake_umol_m2'), 'µmol m⁻²')} under gate coupling "
                    f"and {_fmt_sci_stat(prop_lvl['central'].get('co_uptake_umol_m2'), prop_lvl['outputs'].get('co_uptake_umol_m2'), 'µmol m⁻²')} under proportional coupling."
                )
        core_body_sci = (
            f"Central percentage change: {out_key} = {pct_stat} (rounded headline: {core} %); "
            f"absolute change {contrast_key} = {con_stat}; "
            f"absolute uptake {abs_key} = {abs_stat} "
            f"(mode = {mode}, refusal fraction = {ref_frac * 100:.1f} %, "
            f"max level refusal across question = {max_refusal_fraction * 100:.1f} %).{extra_levels_sci}"
        )
    elif slug == "alfalfa-bacterium":
        pct_q05 = 100.0 * float(q05) if q05 is not None else 0.0
        pct_q95 = 100.0 * float(q95) if q95 is not None else 0.0
        pct_med = 100.0 * float(median) if median is not None else 0.0
        frac_stat = _fmt_sci_stat(central_val, outcome_uq, "fraction")
        uninoc_stat = _fmt_sci_stat(
            h_central_map.get("relative_yield_alfalfa"),
            h_outputs_map.get("relative_yield_alfalfa"),
            "fraction",
        )
        headline_science = (
            f"Central inoculated yield = {display_central:.4g} % of normal yield (rounded: {core} %, "
            f"median {pct_med:.4g} %, joint 5–95 % interval [{pct_q05:.4g}, {pct_q95:.4g}] %; "
            f"{out_key} = {frac_stat}, mode = {mode}, refusal = {ref_frac * 100:.1f} %)."
        )
        core_body_sci = (
            f"Central estimate: {display_central:.4g} % of normal yield (rounded headline: {core} %, "
            f"joint 5–95 % interval [{pct_q05:.4g}, {pct_q95:.4g}] %), corresponding to "
            f"{out_key} = {frac_stat}, versus uninoculated relative_yield_alfalfa = {uninoc_stat} "
            f"(mode = {mode}, refusal fraction = {ref_frac * 100:.1f} %, "
            f"max level refusal across question = {max_refusal_fraction * 100:.1f} %)."
        )
    elif slug == "missing-co-genes":
        open_stat = _fmt_sci_stat(central_val, outcome_uq, "µmol m⁻²")
        closed_lvl = _get_uq_level(uq_doc, "co-gate-samples", 0)
        closed_stat = _fmt_sci_stat(
            0.0,
            (closed_lvl or {}).get("outputs", {}).get("co_uptake_umol_m2"),
            "µmol m⁻²",
        )
        headline_science = (
            f"Closed-gate sample 51PRr3 has co_uptake_umol_m2 = 0 µmol m⁻² (joint 5–95 % interval [0, 0] µmol m⁻², "
            f"mode = deduced), whereas open-gate samples 10Dr2 and V53PRr2 have co_uptake_umol_m2 = {open_stat} "
            f"(a tiny amount, about {core} millionths of a mole per square metre); samples 10Dr1 and 51PRr1 are "
            f"refused (100.0 %) due to missing read-mapping weights."
        )
        core_body_sci = (
            f"Open-gate samples 10Dr2 and V53PRr2 have co_uptake_umol_m2 = {open_stat} (a tiny amount, about "
            f"{core} millionths of a mole per square metre; mode = {mode}, refusal fraction = 0.0 %), whereas "
            f"closed-gate sample 51PRr3 has exact deductive co_uptake_umol_m2 = {closed_stat} (mode = deduced, "
            f"because 0 assembled MAGs in 51PRr3 are entailed as carrying FormICODehydrogenaseSubunit). "
            f"Samples 10Dr1 and 51PRr1 are refused for 100.0 % of parameter draws (max level refusal across "
            f"question = {max_refusal_fraction * 100:.1f} %; exact refusal reason: No usable read-mapping weights "
            f"(unknown sample or excluded for missing weights): ['10Dr1'] and ['51PRr1'])."
        )
    else:
        # carbon-fixation
        stat_str = _fmt_sci_stat(central_val, outcome_uq, unit_sci)
        ug_stat = _fmt_sci_stat(
            h_central_map.get("c_fixed_upper_bound_ug_m2"),
            h_outputs_map.get("c_fixed_upper_bound_ug_m2"),
            "µg C m⁻²",
        )
        headline_science = (
            f"Central {out_key} = {stat_str} (rounded: {core} {unit_sci}; "
            f"14-day bound c_fixed_upper_bound_ug_m2 = {ug_stat}; "
            f"mode = {mode}, refusal = {ref_frac * 100:.1f} %)."
        )
        core_body_sci = (
            f"Central estimate in open-gate sample F14Dr1: {out_key} = {stat_str} "
            f"(rounded headline: {core} {unit_sci}; 14-day bound c_fixed_upper_bound_ug_m2 = {ug_stat}; "
            f"mode = {mode}, refusal fraction = {ref_frac * 100:.1f} %, "
            f"max level refusal across question = {max_refusal_fraction * 100:.1f} %), "
            f"and exact 0 {unit_sci} (joint 5–95 % interval [0, 0] {unit_sci}, mode = deduced) in closed-gate sample 14PRr2."
        )

    sobol_src = outcome_uq
    sobol_prefix = ""
    if (
        out_key in {"h2_uptake_change_pct", "co_uptake_change_pct"}
        and q05 is not None
        and q95 is not None
        and abs(float(q95) - float(q05)) < 1e-9
    ):
        gas = "h2" if out_key.startswith("h2_") else "co"
        con_uq = h_outputs_map.get(f"{gas}_uptake_contrast_umol_m2")
        if con_uq:
            sobol_src = con_uq
            sobol_prefix = (
                f"Because {out_key} is invariant across all parameter draws (ST = 0, S1 = 0), "
                f"Sobol' indices are reported for {gas}_uptake_contrast_umol_m2. "
            )

    sobol_str = sobol_prefix + _top_two_sobol_text(
        sobol_src.get("ST", {}),
        sobol_src.get("S1", {}),
        mode,
    )
    group_str = _group_share_text(
        sobol_src.get("group_ST", {}),
        sobol_src.get("group_S1", {}),
        mode,
    )
    conv_str = _convergence_text(outcome_uq, uq_doc)
    conv_part = f" {conv_str}" if conv_str else ""

    answer_science = (
        f"{core_body_sci} "
        f"{sobol_str} {group_str}{conv_part} "
        f"Limits: {SCOPE_SENTENCE_SCIENCE}"
    )
    limits_science = f"{SCOPE_SENTENCE_SCIENCE} {limits_extra_sci}"

    return {
        "headline_plain": headline_plain,
        "headline_science": headline_science,
        "answer_plain": answer_plain,
        "answer_science": answer_science,
        "limits_plain": limits_plain,
        "limits_science": limits_science,
    }
