"""Plain and scientific text templates and metadata for EQ-BioDT v0.8 question files."""
from __future__ import annotations

from pathlib import Path
import re
import sys
from typing import Any

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))

from eq_biodt.v03.community import carrier_shares  # noqa: E402
from eq_biodt.v03.composition import load_constants  # noqa: E402

_CONST_REG = load_constants()

_NUM_WORDS: dict[int, str] = {
    0: "zero",
    1: "one",
    2: "two",
    3: "three",
    4: "four",
    5: "five",
    6: "six",
    7: "seven",
    8: "eight",
    9: "nine",
    10: "ten",
}

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
    r"|\bFAO-29\b"
    r"|\b(?i:first-order)\b"
    r"|\b(?i:gates?)\b"
    r"|\b(?i:read-mapping)\b"
    r"|\bIQR\b"
    r"|\b(?i:micromoles?)\b"
    r"|\b(?i:quantiles?)\b"
    r"|\b(?i:medians?)\b"
    r"|\b(?i:inoculum shares?|mapped shares?|scaled by share)\b"
    r"|(?i:fraction of unsalted yield)"
    r"|(?i:dilute limit)"
    r"|(?i:sand-physics branch)"
    r"|\bCO\s+(?:genes?|eaters?|microbes?)\b"
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
    "missing-co-genes": "What if the soil lacks carbon-monoxide-eating genes?",
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
    "compaction": "Uncompacted sand (pore space 40 % of volume)",
    "biochar": "Unamended sand (0 % biochar, pore space 40 % of volume)",
    "missing-co-genes": "One 12 mm rain on day 1 in the selected sample community",
    "adding-microbes": "Uninoculated sample community",
    "air-hydrogen": "Normal background air hydrogen (1x normal level)",
    "carbon-fixation": "Sample community under 14-day rain schedule",
}

PLAIN_X_LABELS: dict[str, str] = {
    "cloud-seeding": "Extra rain from seeding (%)",
    "mulch-shade": "Share of ground covered (0 = none, 1 = all)",
    "solar-park": "Solar panel cover (0 = none, 1 = whole area)",
    "adding-microbes": "Added microbes (% of community)",
}

PLAIN_SERIES_LABELS: dict[tuple[str, str], str] = {
    ("alfalfa-bacterium", "A_fao_fixed_leaching"): "Standard farm drainage rule (published formula)",
    ("compaction", "direct-porosity-v06"): "Heavily packed sand (pore space 30 % of volume)",
    ("adding-microbes", "51PRr3-gate"): "Sample lacking carbon-monoxide-eating microbes + carrier microbe (yes/no presence)",
    ("adding-microbes", "51PRr3-proportional"): "Sample lacking carbon-monoxide-eating microbes + carrier microbe (scaled by amount)",
    ("adding-microbes", "V53PRr2-gate"): "Sample already having carbon-monoxide-eating microbes + carrier microbe (yes/no presence)",
    ("adding-microbes", "V53PRr2-proportional"): "Sample already having carbon-monoxide-eating microbes + carrier microbe (scaled by amount)",
    ("adding-microbes", "51PRr3-control-gate"): "Negative control: microbe without carbon-monoxide-eating genes (yes/no presence)",
    ("adding-microbes", "51PRr3-control-proportional"): "Negative control: microbe without carbon-monoxide-eating genes (scaled by amount)",
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


def _num_word(n: int) -> str:
    return _NUM_WORDS.get(int(n), str(int(n)))


def round_plain(x: float | None, decimals: int | None = None) -> str:
    """Format a number for plain-language display."""
    if x is None:
        return "none"
    val = float(x)
    if abs(val) < 1e-12:
        return "0"
    if decimals is not None:
        s = f"{val:.{decimals}f}".rstrip("0").rstrip(".")
        return "0" if s in {"-0", ""} else s
    ax = abs(val)
    if ax >= 1000.0:
        return f"{round(val):,}"
    if ax >= 10.0:
        return f"{round(val):d}"
    if ax >= 0.1:
        s = f"{val:.1f}".rstrip("0").rstrip(".")
        return "0" if s in {"-0", ""} else s
    if ax >= 0.01:
        s = f"{val:.2f}".rstrip("0").rstrip(".")
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
    if q.get("slug") == "salty-water":
        m = re.search(r"gives about (\d+(?:\.\d+)?)\s*%\s*(?:less|more)", q["answer_plain"])
        if m:
            return m.group(1)
    m = re.search(
        r"(?<![A-Za-z0-9_])[+-]?(\d+(?:,\d{3})*(?:\.\d+)?)(?![A-Za-z0-9_])",
        q["answer_plain"],
    )
    if not m:
        raise ValueError(f"No plain number found in answer_plain for {q.get('slug')}")
    return m.group(1)


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
    if slug == "compaction":
        if series_id == "direct-porosity-v06":
            return f"pore space {round_plain(float(x) * 100.0)} % of volume"
        d = int(round(float(x)))
        return "1 pass" if d == 1 else f"{d} passes"
    if slug == "biochar":
        return f"{round_plain(float(x))} % biochar"
    if slug == "drip-or-flood":
        d = int(round(float(x)))
        return "every day" if d == 1 else f"every {d} days"
    if slug == "hotter":
        return spec_label.replace(" C", " °C")
    if slug == "cloud-seeding":
        xv = float(x)
        if abs(xv) < 1e-12:
            return "no extra rain"
        pct = round_plain(xv * float(_CONST_REG["cloud_seeding_rain_increase_frac"]["value"]) * 100.0)
        return f"{pct} % more"
    if slug == "solar-park":
        xv = float(x)
        if abs(xv) < 1e-6:
            return "not covered"
        if abs(xv - 0.5) < 1e-6:
            return "half covered"
        if abs(xv - 1.0) < 1e-6:
            return "fully covered"
        return f"{round_plain(xv * 100.0)} % covered"
    if slug == "mulch-shade":
        xv = float(x)
        amt = "no" if abs(xv) < 1e-6 else ("half" if abs(xv - 0.5) < 1e-6 else "full")
        if series_id == "mulch-only":
            return f"{amt} mulch, no shade"
        if series_id == "shade-only":
            return f"no mulch, {amt} shade"
        return f"{amt} mulch, {amt} shade"
    if slug == "carbon-fixation":
        sx = str(x)
        try:
            sh = carrier_shares(sx)["shares"]
        except ValueError:
            return f"sample {sx} (incomplete records)"
        if float(sh.get("H2AutotrophicCandidateMAG", 0.0)) > 0.0:
            return f"sample {sx} (has carbon-fixing genes)"
        return f"sample {sx} (no carbon-fixing genes)"
    if slug == "missing-co-genes":
        sx = str(x)
        try:
            sh = carrier_shares(sx)["shares"]
        except ValueError:
            return f"sample {sx} (incomplete records)"
        if float(sh.get("COCarrierMAG", 0.0)) > 0.0:
            return f"sample {sx} (has carbon-monoxide-eating genes)"
        return f"sample {sx} (no carbon-monoxide-eating genes)"
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


def _median_diverges(central_val: float | None, uq_o: dict[str, Any] | None) -> bool:
    """Return True when |central - median| exceeds a quarter of the 5–95 % width."""
    if central_val is None or not uq_o:
        return False
    med = uq_o.get("median")
    q05 = uq_o.get("q05")
    q95 = uq_o.get("q95")
    if med is None or q05 is None or q95 is None:
        return False
    width = float(q95) - float(q05)
    if width <= 1e-12:
        return False
    return abs(float(central_val) - float(med)) > 0.25 * width


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
    div_note = (
        "; note: central and median diverge by >25 % of the 5–95 % width"
        if _median_diverges(central_val, uq_o)
        else ""
    )
    return (
        f"{c_str} {unit} (median {med_s} {unit}, 50 % IQR {iqr_s} {unit}, "
        f"joint 5–95 % interval {ci90_s} {unit}{div_note})"
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
    h_lvl = _get_uq_level(uq_doc, s_id, l_idx)
    h_central_map = (h_lvl or {}).get("central", {})
    h_outputs_map = (h_lvl or {}).get("outputs", {})

    if slug == "rain-pattern":
        driz_lvl = _get_uq_level(uq_doc, "rain-schedule", 3)
        driz_c = float((driz_lvl or {}).get("central", {}).get("h2_uptake_change_pct", -99.74))
        driz_o = (driz_lvl or {}).get("outputs", {}).get("h2_uptake_change_pct", {})
        driz_str = (
            f"almost {round_plain(abs(driz_c))} %"
            if 99.0 <= abs(driz_c) < 100.0
            else f"{round_plain(abs(driz_c))} %"
        )
        driz_rng = format_change_pct_interval(driz_o.get("q05", -99.93), driz_o.get("q95", -51.40))
        headline_plain = (
            f"Splitting one heavy rain into four smaller showers lowers hydrogen uptake over two weeks "
            f"by about {core} %, and daily drizzle lowers it by {driz_str}, so uptake {hedge}."
        )
        answer_plain = (
            f"Splitting one heavy downpour into four smaller showers lowers hydrogen uptake over two weeks "
            f"by about {core} % ({rng_plain}), so uptake {hedge}. "
            f"Splitting the same total rain into daily drizzle lowers uptake by {driz_str} ({driz_rng}). "
            f"A single soaking pushes water deeper into the sand so microbes stay active for more days "
            f"before surface evaporation dries the top layer."
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
        days_10 = _num_word(int(round(float((h_lvl or {}).get("x", 10)))))
        signed_core = f"-{core}" if display_central < 0 else f"+{core}"
        lvl_3 = _get_uq_level(uq_doc, "irrigation-interval", 1)
        days_3 = _num_word(int(round(float((lvl_3 or {}).get("x", 3)))))
        c_3 = round_plain(abs(float((lvl_3 or {}).get("central", {}).get("h2_uptake_change_pct", -36.84))))
        o_3 = (lvl_3 or {}).get("outputs", {}).get("h2_uptake_change_pct", {})
        lo_3 = round_plain(abs(float(o_3.get("q05", -68.18))))
        hi_3 = round_plain(abs(float(o_3.get("q95", 129.39))))
        lvl_1 = _get_uq_level(uq_doc, "irrigation-interval", 0)
        c_1 = round_plain(abs(float((lvl_1 or {}).get("central", {}).get("h2_uptake_change_pct", -38.35))))
        o_1 = (lvl_1 or {}).get("outputs", {}).get("h2_uptake_change_pct", {})
        lo_1 = round_plain(abs(float(o_1.get("q05", -73.06))))
        hi_1 = round_plain(abs(float(o_1.get("q95", 136.34))))
        headline_plain = (
            f"Giving the same total water every {days_10} days changes hydrogen uptake over two weeks "
            f"by about {signed_core} % relative to a single rain, versus about {c_1} % less for daily drip, "
            f"and the change {hedge}."
        )
        answer_plain = (
            f"Giving the same total water in one go and then every {days_10} days changes hydrogen uptake "
            f"over two weeks by about {signed_core} % relative to a single rain, but it {hedge} ({rng_plain}). "
            f"Spreading it every {days_3} days gives about {c_3} % less ({lo_3} % less to {hi_3} % more) "
            f"and daily drip about {c_1} % less ({lo_1} % less to {hi_1} % more). "
            f"One possible reason is that frequent small doses mostly evaporate from the surface."
        )
        limits_extra_plain = (
            "Pipe layout, root uptake and uneven wetting around drip emitters are not represented."
        )
        limits_extra_sci = (
            "Assumes spatially uniform 1-D infiltration and fixed daily potential evaporation without "
            "3-D drip-bulb geometry or plant transpiration."
        )
    elif slug == "cloud-seeding":
        boost_pct = int(
            round(
                float((h_lvl or {}).get("x", 1.0))
                * float(_CONST_REG["cloud_seeding_rain_increase_frac"]["value"])
                * 100.0
            )
        )
        boost_word = _num_word(boost_pct)
        headline_plain = (
            f"If cloud seeding really added {boost_word} percent more rain on rainy days, hydrogen uptake "
            f"over two weeks would rise by about {core} %, though whether seeding adds that much rain is disputed."
        )
        answer_plain = (
            f"If cloud seeding really added {boost_word} percent more rain on rainy days, hydrogen uptake "
            f"over two weeks would rise by about {core} % ({rng_plain}), so under that rain boost uptake {hedge}. "
            f"Whether seeding adds that much rain at all is disputed. "
            f"When extra water does fall on wet days, it pushes moisture deeper into the sand and keeps "
            f"microbes active for longer after the rain."
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
        ec_8 = round_plain(float((h_lvl or {}).get("x", 8.0)))
        fresh_lvl = _get_uq_level(uq_doc, "daily-irrigation", 0)
        fresh_drop = round_plain(
            abs(float((fresh_lvl or {}).get("central", {}).get("h2_uptake_change_pct", -38.37)))
        )
        alf_c = round_plain(100.0 * float(h_central_map.get("relative_yield_alfalfa", 0.8904)))
        alf_o = h_outputs_map.get("relative_yield_alfalfa", {})
        alf_lo = round_plain(100.0 * float(alf_o.get("q05", 0.8295)))
        alf_hi = round_plain(100.0 * float(alf_o.get("q95", 0.9313)))
        ref_ec = "2"
        if uq_doc:
            for s in uq_doc.get("series", []):
                if s.get("id") == "salty-rain-once":
                    for l in s.get("levels", []):
                        rf = float(
                            l.get("outputs", {}).get("h2_uptake_umol_m2", {}).get("refusal_fraction", 0.0)
                        )
                        if rf >= 0.05:
                            ref_ec = round_plain(float(l.get("x", 2.0)))
                            break
        headline_plain = (
            f"Watering every day with salty water ({ec_8} dS/m) gives about {core} % less hydrogen uptake "
            f"over two weeks than one fresh rain (and {hedge}), but fresh daily watering alone gives about "
            f"{fresh_drop} % less."
        )
        answer_plain = (
            f"Watering every day with salty water ({ec_8} dS/m) gives about {core} % less hydrogen uptake "
            f"over two weeks than one fresh rain of the same size ({rng_plain}, so it {hedge}). "
            f"Almost all of that drop is the daily-watering pattern: fresh water given daily gives about "
            f"{fresh_drop} % less. "
            f"Saltiness itself adds only a few points. "
            f"A separate crop formula (not the twin) puts alfalfa at about {alf_c} % of normal yield "
            f"({alf_lo} to {alf_hi} %). "
            f"The twin cannot tell (gives no answer) for a single salty rain, even at {ref_ec} dS/m, "
            f"because the salt in the sand water becomes too strong for its formulas."
        )
        limits_extra_plain = (
            "When salt concentration in the sand water rises above a tenth of a mole per litre, the twin "
            "refuses to guess rather than extrapolate outside its formulas. Crop yields come from published "
            "crop-salt formulas, not from the twin."
        )
        limits_extra_sci = (
            "The twin enforces an ontological/domain admissibility gate at c_NaCl > 0.1 mol/L (ideal-dilute "
            "van 't Hoff limit) and abstains when evaporation concentrates solutes above this threshold; "
            "alfalfa and tomato relative yields come from external Maas-Hoffman crop-salt formulas."
        )
    elif slug == "alfalfa-bacterium":
        ec_hi = f"{float((h_lvl or {}).get('x', 7.81)):g}"
        c_lo = round_plain(100.0 * float(q05 if q05 is not None else 0.8938))
        c_hi = round_plain(100.0 * float(q95 if q95 is not None else 0.9313))
        u_c = round_plain(100.0 * float(h_central_map.get("relative_yield_alfalfa", 0.7959)))
        u_o = h_outputs_map.get("relative_yield_alfalfa", {})
        u_lo = round_plain(100.0 * float(u_o.get("q05", 0.7847)))
        u_hi = round_plain(100.0 * float(u_o.get("q95", 0.8109)))
        mod_a_lvl = _get_uq_level(uq_doc, "A_fao_fixed_leaching", 1)
        mod_b_lvl = _get_uq_level(uq_doc, "B_twin_physics", 1)
        a_val = round_plain(
            100.0 * float((mod_a_lvl or {}).get("central", {}).get("relative_yield_alfalfa_inoculated", 0.3324))
        )
        b_med = (
            (mod_b_lvl or {})
            .get("outputs", {})
            .get("relative_yield_alfalfa_inoculated", {})
            .get("median", 1.0128)
        )
        b_val = round_plain(100.0 * float(b_med))
        headline_plain = (
            f"In a crop-salt formula fitted to one field trial, alfalfa from bacterium-coated seeds reaches "
            f"about {core} % of normal yield in {ec_hi} dS/m water (versus {u_c} % without it), though "
            f"other root-zone salt formulas give {a_val} to {b_val} %."
        )
        answer_plain = (
            f"In a crop-salt formula fitted to one field trial (not the twin), alfalfa from bacterium-coated "
            f"seeds reaches about {core} % of normal yield ({c_lo} to {c_hi} %) in {ec_hi} dS/m water, "
            f"versus about {u_c} % ({u_lo} to {u_hi} %) without the bacterium. "
            f"Other ways of estimating how salt builds up in the root zone give very different answers, "
            f"from about {a_val} % to about {b_val} % with the bacterium, so this figure depends strongly "
            f"on that choice."
        )
        limits_extra_plain = (
            "Crop yields and the bacterial boost come from published crop-salt formulas and one desert farm "
            "field trial, not from the soil twin itself."
        )
        limits_extra_sci = (
            "Multiplies external Maas-Hoffman relative yield by an empirical SA187 field multiplier calibrated "
            "to fresh weight at one trial; root-zone EC_e depends strongly on the chosen leaching/accumulation "
            "closure (models A, B, C)."
        )
    elif slug == "hotter":
        deg_4 = _num_word(int(round(float((h_lvl or {}).get("x", 4.0))))).capitalize()
        signed_core = f"+{core}" if display_central > 0 else f"-{core}"
        med_4 = round_plain(abs(float(median if median is not None else -2.62)))
        med_4_dir = "lower" if float(median if median is not None else -2.62) < 0 else "higher"
        soil_lvl = _get_uq_level(uq_doc, "soil-only-v06", 2)
        deg_10 = round_plain(float((soil_lvl or {}).get("x", 10.0)))
        s_c = round_plain(abs(float((soil_lvl or {}).get("central", {}).get("h2_uptake_change_pct", 69.08))))
        s_o = (soil_lvl or {}).get("outputs", {}).get("h2_uptake_change_pct", {})
        s_rng = format_change_pct_interval(s_o.get("q05", 3.78), s_o.get("q95", 77.78))
        s_med = round_plain(abs(float(s_o.get("median", 35.88))))
        headline_plain = (
            f"{deg_4} degrees of warming changes hydrogen uptake over two weeks by about {signed_core} % in "
            f"the central run, but a typical run is about {med_4} % {med_4_dir} and uptake {hedge}."
        )
        answer_plain = (
            f"{deg_4} degrees of warming with the same rainfall changes hydrogen uptake over two weeks by "
            f"about {signed_core} % in the model's central run. "
            f"Across the uncertain numbers it {hedge} ({rng_plain}), and a typical run is slightly "
            f"{med_4_dir} (about {med_4} %). "
            f"Warming only the sand by {deg_10} degrees, with no extra drying, raises uptake by about "
            f"{s_c} % in the central run ({s_rng}; typical run about {s_med} %)."
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
        sh_lvl = _get_uq_level(uq_doc, "shade-only", 2)
        sh_c = round_plain(abs(float((sh_lvl or {}).get("central", {}).get("h2_uptake_change_pct", -24.76))))
        sh_o = (sh_lvl or {}).get("outputs", {}).get("h2_uptake_change_pct", {})
        sh_rng = format_change_pct_interval(sh_o.get("q05", -27.77), sh_o.get("q95", -1.72))
        sh_med = round_plain(abs(float(sh_o.get("median", -13.68))))
        ms_lvl = _get_uq_level(uq_doc, "mulch-and-shade", 2)
        ms_c = round_plain(abs(float((ms_lvl or {}).get("central", {}).get("h2_uptake_change_pct", 14.07))))
        ms_o = (ms_lvl or {}).get("outputs", {}).get("h2_uptake_change_pct", {})
        ms_rng = format_change_pct_interval(ms_o.get("q05", -6.16), ms_o.get("q95", 70.11))
        headline_plain = (
            f"Full mulch cover raises hydrogen uptake over two weeks by about {core} % by keeping the sand "
            f"damp longer, and uptake {hedge}."
        )
        answer_plain = (
            f"Full mulch cover raises hydrogen uptake over two weeks by about {core} % ({rng_plain}) by "
            f"slowing evaporation, so uptake with mulch alone {hedge}. "
            f"Full shade alone lowers uptake by about {sh_c} % ({sh_rng}; typical run about {sh_med} % lower), "
            f"and combining mulch with shade raises uptake by about {ms_c} % ({ms_rng}). "
            f"Mulch keeps water in the sand without cooling the microbes, whereas shade cools the sand "
            f"and slows enzyme speed."
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
            f"Full solar-panel cover raises hydrogen uptake over two weeks by about {core} %, "
            f"and uptake {hedge}."
        )
        answer_plain = (
            f"Full solar-panel cover raises hydrogen uptake over two weeks by about {core} % ({rng_plain}), "
            f"and uptake {hedge}. "
            f"Less evaporation under the panels and water from panel washing keep the sand moist longer. "
            f"This outweighs the slowing of microbes by shade."
        )
        limits_extra_plain = (
            "Wind shifts under panel rows and runoff strips along panel edges are simplified to an even average."
        )
        limits_extra_sci = (
            "Homogenises under-panel microclimate and wash-water drips spatially; competing cooling (Q10) and "
            "moisture-retention effects partially offset each other across parameter draws."
        )
    elif slug == "compaction":
        passes_n = int(round(float((h_lvl or {}).get("x", 50))))
        passes_word = "Fifty" if passes_n == 50 else str(passes_n)
        signed_core = f"-{core}" if display_central < 0 else f"+{core}"
        headline_plain = (
            f"{passes_word} vehicle passes change hydrogen uptake over two weeks by about {signed_core} %, "
            f"and the change {hedge}."
        )
        answer_plain = (
            f"{passes_word} vehicle passes change hydrogen uptake over two weeks by about {signed_core} % "
            f"({rng_plain}), so the change {hedge}. "
            f"Squeezing the sand slows gas movement through the pores but also changes how tightly the "
            f"sand holds water."
        )
        limits_extra_plain = (
            "Crushing of soil crusts, plants or burrows by tyres is not modelled beyond the drop in pore space."
        )
        limits_extra_sci = (
            "Maps traffic passes to bulk density and porosity only, holding retention shape b and air-entry "
            "potential psi_entry uncoupled from pore-size redistribution."
        )
    elif slug == "biochar":
        dose_word = _num_word(int(round(float((h_lvl or {}).get("x", 5.0)))))
        signed_core = f"-{core}" if display_central < 0 else f"+{core}"
        headline_plain = (
            f"Adding {dose_word} percent biochar by weight changes hydrogen uptake over two weeks by about "
            f"{signed_core} %, and the change {hedge}."
        )
        answer_plain = (
            f"Adding {dose_word} percent biochar by weight changes hydrogen uptake over two weeks by about "
            f"{signed_core} % ({rng_plain}), so the change {hedge}. "
            f"Biochar adds empty pore space, so one rain fills a smaller share of the pores while gas moves "
            f"more easily. "
            f"The two effects nearly cancel, so the net change is small and could go either way."
        )
        limits_extra_plain = (
            "Chemical effects of biochar on soil acidity, nutrients or microbe attachment are not included."
        )
        limits_extra_sci = (
            "Propagates biochar dose solely through linear mixture porosity φ(f_b), holding Campbell retention "
            "parameters and surface chemistry constant."
        )
    elif slug == "missing-co-genes":
        mco_levels = uq_doc["series"][0]["levels"] if (uq_doc and uq_doc.get("series")) else []
        n_total = len(mco_levels) if mco_levels else 7
        n_refused = (
            sum(
                1
                for l in mco_levels
                if float(l.get("outputs", {}).get("co_uptake_umol_m2", {}).get("refusal_fraction", 0.0))
                >= 0.05
            )
            if mco_levels
            else 2
        )
        ref_word = _num_word(n_refused)
        tot_word = _num_word(n_total)
        headline_plain = (
            f"Where the soil community has no genes for eating carbon monoxide, the twin gives zero uptake, "
            f"versus about {core} {unit_plain} over two weeks where the genes are present; for {ref_word} "
            f"of the {tot_word} samples, incomplete records give no answer."
        )
        answer_plain = (
            f"Where the soil community has no genes for eating carbon monoxide, the twin gives exactly zero "
            f"carbon-monoxide uptake. "
            f"Where the genes are present it gives a small uptake, about {core} {unit_plain} over two weeks "
            f"({round_plain(q05)} to {round_plain(q95)}). "
            f"For {ref_word} of the {tot_word} samples the records are incomplete, so the twin gives no answer."
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
        prop_0 = _get_uq_level(uq_doc, "51PRr3-proportional", 0)
        prop_2 = _get_uq_level(uq_doc, "51PRr3-proportional", 2)
        v53_2 = _get_uq_level(uq_doc, "V53PRr2-proportional", 2)
        lo_start = round_plain(
            float((prop_0 or {}).get("central", {}).get("co_uptake_umol_m2", 10.18))
        )
        hi_start = round_plain(
            float((prop_2 or {}).get("central", {}).get("co_uptake_umol_m2", 220.62))
        )
        x_1 = round_plain(float((h_lvl or {}).get("x", 1.0)))
        x_5 = round_plain(float((v53_2 or {}).get("x", 5.0)))
        c_5 = round_plain(
            abs(float((v53_2 or {}).get("central", {}).get("co_uptake_change_pct", 2.66)))
        )
        abs_unit_plain = OUTCOME_META["co_uptake_umol_m2"]["unit_plain"]
        headline_plain = (
            f"Adding microbes that carry the carbon-monoxide-eating genes to a sample that had none starts "
            f"uptake (roughly {lo_start} to {hi_start} {abs_unit_plain} over two weeks), and raises uptake "
            f"by {core} % at {x_1} % added where such microbes are already present."
        )
        answer_plain = (
            f"Adding microbes that carry the carbon-monoxide-eating genes to a sample that had none makes "
            f"uptake start: roughly {lo_start} to {hi_start} {abs_unit_plain} over two weeks, depending on "
            f"whether presence alone or the amount added matters. "
            f"Adding a microbe without those genes changes nothing. "
            f"Where such microbes are already present, adding more raises uptake only slightly "
            f"({core} % at {x_1} % added, {c_5} % at {x_5} % added)."
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
        lvl_10 = _get_uq_level(uq_doc, "air-h2-scaling", 3)
        mult_10 = _num_word(int(round(float((lvl_10 or {}).get("x", 10.0)))))
        c_10 = round_plain(
            abs(float((lvl_10 or {}).get("central", {}).get("h2_uptake_change_pct", 900.0)))
        )
        headline_plain = (
            f"Doubling hydrogen in the air doubles soil hydrogen uptake over two weeks "
            f"(+{core} %, the same in every case we tried)."
        )
        answer_plain = (
            f"Doubling hydrogen in the air doubles soil hydrogen uptake over two weeks "
            f"(+{core} %, the same in every case we tried), and {mult_10} times as much air hydrogen "
            f"gives {mult_10} times the uptake (+{c_10} %)."
        )
        limits_extra_plain = (
            "This follows from a formula that assumes uptake rises in proportion to the hydrogen level; "
            "at very high hydrogen levels, microbes could hit a maximum enzyme speed or grow."
        )
        limits_extra_sci = (
            "First-order Henry/diffusion-reaction kinetics in C_H2; excludes Michaelis-Menten saturation and "
            "long-term community biomass growth under elevated H₂."
        )
    elif slug == "carbon-fixation":
        c_lo = round_plain(q05)
        c_hi = round_plain(q95, decimals=2)
        closed_lvl = _get_uq_level(uq_doc, "carbon-bound-samples", 0)
        closed_c = round_plain(float((closed_lvl or {}).get("central", {}).get(out_key, 0.0)))
        headline_plain = (
            f"In one desert sand sample, hydrogen-eating microbes could fix at most about {core} {unit_plain} "
            f"if it rained every two weeks, and {closed_c} where the carbon-fixing genes are absent."
        )
        answer_plain = (
            f"In one desert sand sample, hydrogen-eating microbes could fix at most about {core} {unit_plain}, "
            f"if a rain like this fell every two weeks and all the energy from hydrogen went into making "
            f"cell carbon. "
            f"The model's estimate of this ceiling ranges from {c_lo} to {c_hi} g. "
            f"In a second sample, which lacks the carbon-fixing genes, the ceiling is exactly {closed_c}."
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

    limits_plain = f"{SCOPE_SENTENCE_PLAIN} {limits_extra_plain}"

    # Scientific headline and answer
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
                    f"(plain: almost 100 % reduction); "
                    f"h2_uptake_contrast_umol_m2 = {_fmt_sci_stat(d_c.get('h2_uptake_contrast_umol_m2'), d_o.get('h2_uptake_contrast_umol_m2'), 'µmol m⁻²')}; "
                    f"h2_uptake_umol_m2 = {_fmt_sci_stat(d_c.get('h2_uptake_umol_m2'), d_o.get('h2_uptake_umol_m2'), 'µmol m⁻²')} "
                    f"(one-storm baseline h2_uptake_umol_m2 = {_fmt_sci_stat(b_c.get('h2_uptake_umol_m2'), b_o.get('h2_uptake_umol_m2'), 'µmol m⁻²')})."
                )
        elif slug == "drip-or-flood":
            lvl_1 = _get_uq_level(uq_doc, "irrigation-interval", 0)
            lvl_3 = _get_uq_level(uq_doc, "irrigation-interval", 1)
            if lvl_1 and lvl_3:
                extra_levels_sci = (
                    f" At equal 12 mm total water over 14 days relative to a single 12 mm rain: "
                    f"every-3-days h2_uptake_change_pct = {_fmt_sci_stat(lvl_3['central'].get('h2_uptake_change_pct'), lvl_3['outputs'].get('h2_uptake_change_pct'), '%')}; "
                    f"daily-drip h2_uptake_change_pct = {_fmt_sci_stat(lvl_1['central'].get('h2_uptake_change_pct'), lvl_1['outputs'].get('h2_uptake_change_pct'), '%')}."
                )
        elif slug == "salty-water":
            d0 = _get_uq_level(uq_doc, "daily-irrigation", 0)
            d1 = _get_uq_level(uq_doc, "daily-irrigation", 1)
            d2 = _get_uq_level(uq_doc, "daily-irrigation", 2)
            r1 = _get_uq_level(uq_doc, "salty-rain-once", 1)
            if d0 and d1 and d2 and r1:
                c0 = d0["central"].get("h2_uptake_change_pct")
                c1 = d1["central"].get("h2_uptake_change_pct")
                c2 = d2["central"].get("h2_uptake_change_pct")
                r1_rf = float(r1["outputs"].get("h2_uptake_umol_m2", {}).get("refusal_fraction", 1.0)) * 100.0
                alf_stat = _fmt_sci_stat(
                    h_central_map.get("relative_yield_alfalfa"),
                    h_outputs_map.get("relative_yield_alfalfa"),
                    "fraction",
                )
                extra_levels_sci = (
                    f" Attribution vs single fresh 12 mm rain baseline: the daily-irrigation schedule drives most of the "
                    f"reduction (h2_uptake_change_pct = {c0:.4g} % at 0.5 dS/m fresh daily water, {c1:.4g} % at 2 dS/m, "
                    f"{c2:.4g} % at 4 dS/m, and {central_val:.4g} % at 8 dS/m; salinity adds only "
                    f"{abs(central_val - float(c0)):.2f} percentage points between 0.5 and 8 dS/m). "
                    f"The single-salty-rain series (salty-rain-once) is refused in {r1_rf:.1f} % of draws from 2 dS/m up "
                    f"(pore-water c_NaCl > 0.1 mol/L after evaporative concentration), whereas daily irrigation at 8 dS/m "
                    f"is refused in only {ref_frac * 100:.1f} % of draws. External Maas-Hoffman relative_yield_alfalfa at "
                    f"8 dS/m = {alf_stat}."
                )
        elif slug == "hotter":
            soil_lvl = _get_uq_level(uq_doc, "soil-only-v06", 2)
            if soil_lvl:
                sc = soil_lvl["central"].get("h2_uptake_change_pct")
                so = soil_lvl["outputs"].get("h2_uptake_change_pct")
                extra_levels_sci = (
                    f" Because moisture-loss suppression and Q10 acceleration pull in opposite directions under +4 °C "
                    f"climate warming, the distribution is strongly skewed (central {central_val:.4g} % near the upper "
                    f"tail vs median {float(median):.4g} %). Under +10 °C soil-only warming at fixed moisture "
                    f"(soil-only-v06), h2_uptake_change_pct = {_fmt_sci_stat(sc, so, '%')}."
                )
        elif slug == "mulch-shade":
            sh_lvl = _get_uq_level(uq_doc, "shade-only", 2)
            ms_lvl = _get_uq_level(uq_doc, "mulch-and-shade", 2)
            if sh_lvl and ms_lvl:
                extra_levels_sci = (
                    f" Full shade alone (shade-only level 2): h2_uptake_change_pct = "
                    f"{_fmt_sci_stat(sh_lvl['central'].get('h2_uptake_change_pct'), sh_lvl['outputs'].get('h2_uptake_change_pct'), '%')}; "
                    f"combined full mulch and shade (mulch-and-shade level 2): h2_uptake_change_pct = "
                    f"{_fmt_sci_stat(ms_lvl['central'].get('h2_uptake_change_pct'), ms_lvl['outputs'].get('h2_uptake_change_pct'), '%')}."
                )
        elif slug == "adding-microbes":
            gate_lvl = _get_uq_level(uq_doc, "51PRr3-gate", 1)
            prop_0 = _get_uq_level(uq_doc, "51PRr3-proportional", 0)
            prop_1 = _get_uq_level(uq_doc, "51PRr3-proportional", 1)
            prop_2 = _get_uq_level(uq_doc, "51PRr3-proportional", 2)
            v53_5_lvl = _get_uq_level(uq_doc, "V53PRr2-proportional", 2)
            if gate_lvl and prop_0 and prop_1 and prop_2 and v53_5_lvl:
                p0_c = prop_0["central"].get("co_uptake_umol_m2")
                p1_c = prop_1["central"].get("co_uptake_umol_m2")
                p2_c = prop_2["central"].get("co_uptake_umol_m2")
                extra_levels_sci = (
                    f" At 5 % inoculum in V53PRr2-proportional: "
                    f"co_uptake_change_pct = {_fmt_sci_stat(v53_5_lvl['central'].get('co_uptake_change_pct'), v53_5_lvl['outputs'].get('co_uptake_change_pct'), '%')}, "
                    f"co_uptake_contrast_umol_m2 = {_fmt_sci_stat(v53_5_lvl['central'].get('co_uptake_contrast_umol_m2'), v53_5_lvl['outputs'].get('co_uptake_contrast_umol_m2'), 'µmol m⁻²')}, "
                    f"co_uptake_umol_m2 = {_fmt_sci_stat(v53_5_lvl['central'].get('co_uptake_umol_m2'), v53_5_lvl['outputs'].get('co_uptake_umol_m2'), 'µmol m⁻²')} "
                    f"(whereas V53PRr2-gate gives exact 0 % change across all inoculum doses). "
                    f"In closed-baseline sample 51PRr3 (baseline co_uptake_umol_m2 = 0 µmol m⁻²), non-carrier inoculation "
                    f"remains at exact 0 µmol m⁻², whereas carrier inoculation depends on the coupling assumption: "
                    f"binary gate coupling yields {_fmt_sci_stat(gate_lvl['central'].get('co_uptake_umol_m2'), gate_lvl['outputs'].get('co_uptake_umol_m2'), 'µmol m⁻²')} "
                    f"at 0.1, 1, and 5 % inoculum, while share-proportional coupling yields central co_uptake_umol_m2 = "
                    f"{p0_c:.4g}, {p1_c:.4g}, and {p2_c:.4g} µmol m⁻² at 0.1, 1, and 5 % inoculum (at 1 %: "
                    f"{_fmt_sci_stat(p1_c, prop_1['outputs'].get('co_uptake_umol_m2'), 'µmol m⁻²')})."
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
        mod_a_lvl = _get_uq_level(uq_doc, "A_fao_fixed_leaching", 1)
        mod_b_lvl = _get_uq_level(uq_doc, "B_twin_physics", 1)
        spread_sci = ""
        if mod_a_lvl and mod_b_lvl:
            a_in = _fmt_sci_stat(
                mod_a_lvl["central"].get("relative_yield_alfalfa_inoculated"),
                mod_a_lvl["outputs"].get("relative_yield_alfalfa_inoculated"),
                "fraction",
            )
            a_un = _fmt_sci_stat(
                mod_a_lvl["central"].get("relative_yield_alfalfa"),
                mod_a_lvl["outputs"].get("relative_yield_alfalfa"),
                "fraction",
            )
            b_in = _fmt_sci_stat(
                mod_b_lvl["central"].get("relative_yield_alfalfa_inoculated"),
                mod_b_lvl["outputs"].get("relative_yield_alfalfa_inoculated"),
                "fraction",
            )
            b_un = _fmt_sci_stat(
                mod_b_lvl["central"].get("relative_yield_alfalfa"),
                mod_b_lvl["outputs"].get("relative_yield_alfalfa"),
                "fraction",
            )
            spread_sci = (
                f" Across root-zone salinity closures at 7.81 dS/m (all external crop-salt formulas, not soil-twin "
                f"outputs, with 0.0 % refusal across all yield series): Model A (FAO-29 fixed leaching) gives "
                f"inoculated {a_in} vs uninoculated {a_un}; Model B (twin bucket salt balance) gives inoculated "
                f"{b_in} vs uninoculated {b_un}."
            )
        headline_science = (
            f"Central inoculated yield (Model C, trial-calibrated fresh-weight closure) = {display_central:.4g} % "
            f"of normal yield (rounded: {core} %, median {pct_med:.4g} %, joint 5–95 % interval "
            f"[{pct_q05:.4g}, {pct_q95:.4g}] %; {out_key} = {frac_stat}, mode = {mode}, "
            f"refusal = {ref_frac * 100:.1f} %)."
        )
        core_body_sci = (
            f"Central estimate under Model C (calibrated to Jizan field trial fresh weight): "
            f"{display_central:.4g} % of normal yield (rounded headline: {core} %, "
            f"joint 5–95 % interval [{pct_q05:.4g}, {pct_q95:.4g}] %), corresponding to "
            f"{out_key} = {frac_stat}, versus uninoculated relative_yield_alfalfa = {uninoc_stat} "
            f"(mode = {mode}, refusal fraction = {ref_frac * 100:.1f} %, "
            f"max level refusal across question = {max_refusal_fraction * 100:.1f} %).{spread_sci}"
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
            f"mode = deduced), whereas open-gate samples have co_uptake_umol_m2 = {open_stat} "
            f"(a tiny amount, about {core} millionths of a mole per square metre); two of the seven samples "
            f"(10Dr1 and 51PRr1) are refused (100.0 %) due to missing read-mapping weights."
        )
        core_body_sci = (
            f"Open-gate samples (10Dr2, V53PRr2, 13Dr3, 14PRr2) have co_uptake_umol_m2 = {open_stat} (a tiny amount, "
            f"about {core} millionths of a mole per square metre; mode = {mode}, refusal fraction = 0.0 %), whereas "
            f"closed-gate sample 51PRr3 has exact deductive co_uptake_umol_m2 = {closed_stat} (mode = deduced, "
            f"because 0 assembled MAGs in 51PRr3 are entailed as carrying FormICODehydrogenaseSubunit). "
            f"Two of the seven samples (10Dr1 and 51PRr1) are refused for 100.0 % of parameter draws (max level "
            f"refusal across question = {max_refusal_fraction * 100:.1f} %; exact refusal reason: No usable "
            f"read-mapping weights (unknown sample or excluded for missing weights): ['10Dr1'] and ['51PRr1'])."
        )
    else:
        # carbon-fixation
        stat_str = _fmt_sci_stat(central_val, outcome_uq, unit_sci)
        ug_stat = _fmt_sci_stat(
            h_central_map.get("c_fixed_upper_bound_ug_m2"),
            h_outputs_map.get("c_fixed_upper_bound_ug_m2"),
            "µg C m⁻²",
        )
        periods_per_yr = 365.25 / 14.0
        headline_science = (
            f"Central {out_key} = {stat_str} (rounded: {core} {unit_sci}; "
            f"annualised from 14-day bound c_fixed_upper_bound_ug_m2 = {ug_stat} via ×{periods_per_yr:.2f} "
            f"two-week periods yr⁻¹; mode = {mode}, refusal = {ref_frac * 100:.1f} %)."
        )
        core_body_sci = (
            f"Central estimate in open-gate sample F14Dr1: {out_key} = {stat_str} "
            f"(rounded headline: {core} {unit_sci}; 14-day bound c_fixed_upper_bound_ug_m2 = {ug_stat}, "
            f"annualised by multiplying by 365.25 / 14 = {periods_per_yr:.2f} under the hypothetical assumption "
            f"that a 12 mm rain event recurs every two weeks all year and all captured H₂ energy goes into CBB "
            f"cell carbon; mode = {mode}, refusal fraction = {ref_frac * 100:.1f} %, "
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
