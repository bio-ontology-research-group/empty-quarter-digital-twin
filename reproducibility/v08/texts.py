"""Plain and scientific text templates and metadata for EQ-BioDT v0.8 question files."""
from __future__ import annotations

from pathlib import Path
import re
import sys
from typing import Any

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))

from eq_biodt.v03.community import carrier_shares  # noqa: E402
from eq_biodt.v03.composition import load_constants  # noqa: E402
from reproducibility.v08.science_texts import render_science_texts  # noqa: E402

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
    "missing-co-genes": ("co-gate-samples", 6, "co_uptake_umol_m2"),
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
    if ax >= 1.0:
        s = f"{val:.1f}".rstrip("0").rstrip(".")
        return "0" if s in {"-0", ""} else s
    if ax >= 0.01:
        s = f"{val:.2f}".rstrip("0").rstrip(".")
        return "0" if s in {"-0", ""} else s
    s = f"{val:.3f}".rstrip("0").rstrip(".")
    return "0" if s in {"-0", ""} else s


def science_median(q: dict[str, Any]) -> float:
    """Return the core (median) headline value from a built question dict (in % for relative_yield_*)."""
    if "headline_value" in q and q["headline_value"] is not None:
        val = float(q["headline_value"])
    elif "uncertainty" in q and q["uncertainty"].get("median") is not None:
        val = float(q["uncertainty"]["median"])
    else:
        s_id, l_idx, out_key = HEADLINE_SPEC[q["slug"]]
        s = next(sr for sr in q["series"] if sr["id"] == s_id)
        val = float(s["levels"][l_idx]["outcomes"][out_key]["value"])
    _, _, out_key = HEADLINE_SPEC[q["slug"]]
    if out_key.startswith("relative_yield_") and abs(val) <= 2.0:
        return 100.0 * val
    return val


science_central = science_median


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
            xv = float(x)
            pct = xv * 100.0 if xv <= 1.0 else xv
            return f"pore space {round_plain(pct)} % of volume"
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


def _hedge_phrase(
    q05: float | None,
    q95: float | None,
    frac_neg: float | None = None,
    frac_pos: float | None = None,
) -> str:
    """Return the required hedge phrase matching direction probabilities and the 5–95 % headline interval."""
    if q05 is not None and q95 is not None:
        lo = float(q05)
        hi = float(q95)
        if lo < 0.0 < hi:
            return "may go up or down"
        if hi <= 0.0 and lo < 0.0:
            if frac_neg is not None and float(frac_neg) < 1.0:
                n_runs = min(99, max(1, round(100.0 * float(frac_neg))))
                return f"goes down in about {n_runs} of every 100 model runs"
            return "goes down in every case we tried"
        if lo >= 0.0 and hi > 0.0:
            if frac_pos is not None and float(frac_pos) < 1.0:
                n_runs = min(99, max(1, round(100.0 * float(frac_pos))))
                return f"goes up in about {n_runs} of every 100 model runs"
            return "goes up in every case we tried"
    return "goes up in every case we tried"


def _get_uq_level(uq_doc: dict[str, Any] | None, series_id: str, level_idx: int) -> dict[str, Any] | None:
    if not uq_doc:
        return None
    for s in uq_doc.get("series", []):
        if s.get("id") == series_id:
            for l in s.get("levels", []):
                if int(l.get("level", -1)) == level_idx:
                    return l
    return None


def _level_median(lvl: dict[str, Any] | None, key: str, default: float) -> float:
    if lvl:
        out = (lvl.get("outputs") or {}).get(key) or {}
        if out.get("median") is not None:
            return float(out["median"])
        central = lvl.get("central") or {}
        if central.get(key) is not None:
            return float(central[key])
    return float(default)


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
    q05 = outcome_uq.get("q05")
    q25 = outcome_uq.get("q25")
    q75 = outcome_uq.get("q75")
    q95 = outcome_uq.get("q95")
    median = outcome_uq.get("median")
    core_val = float(median) if median is not None else float(central_val)
    display_core = 100.0 * core_val if (is_rel_yield and abs(core_val) <= 2.0) else core_val
    core = round_plain(abs(display_core) if is_change else display_core)

    meta = OUTCOME_META[out_key]
    unit_plain = meta["unit_plain"]
    if is_change:
        rng_plain = format_change_pct_interval(q05, q95)
    elif is_rel_yield:
        q05_disp = 100.0 * float(q05) if q05 is not None else None
        q95_disp = 100.0 * float(q95) if q95 is not None else None
        rng_plain = _format_interval_plain(q05_disp, q95_disp, "%")
    else:
        rng_plain = _format_interval_plain(q05, q95, unit_plain)

    hedge = _hedge_phrase(
        q05,
        q95,
        outcome_uq.get("frac_neg"),
        outcome_uq.get("frac_pos"),
    )
    h_lvl = _get_uq_level(uq_doc, s_id, l_idx)
    h_outputs_map = (h_lvl or {}).get("outputs", {})

    if slug == "rain-pattern":
        driz_lvl = _get_uq_level(uq_doc, "rain-schedule", 3)
        driz_m = _level_median(driz_lvl, "h2_uptake_change_pct", -99.93)
        driz_o = (driz_lvl or {}).get("outputs", {}).get("h2_uptake_change_pct", {})
        driz_str = (
            f"almost {round_plain(abs(driz_m))} %"
            if 99.0 <= abs(driz_m) < 100.0
            else f"{round_plain(abs(driz_m))} %"
        )
        driz_rng = format_change_pct_interval(driz_o.get("q05", -100.0), driz_o.get("q95", -51.32))
        headline_plain = (
            f"Splitting one heavy rain into four smaller showers lowers hydrogen uptake over two weeks "
            f"by about {core} %, and daily drizzle lowers it by {driz_str}, so uptake {hedge}."
        )
        answer_plain = (
            f"Splitting one heavy downpour into four smaller showers lowers hydrogen uptake over two weeks "
            f"by about {core} % ({rng_plain}), so uptake {hedge}. "
            f"Splitting the same total rain into daily drizzle lowers uptake by {driz_str} ({driz_rng}). "
            f"One soaking keeps the sand wet enough for the microbes for several days, while small daily amounts "
            f"mostly evaporate the same day."
        )
        limits_extra_plain = (
            "The simulation uses a single sand column with fixed evaporation and does not include "
            "surface crusting or runoff down dune slopes."
        )
    elif slug == "drip-or-flood":
        days_10 = _num_word(int(round(float((h_lvl or {}).get("x", 10)))))
        signed_core = f"-{core}" if display_core < 0 else f"+{core}"
        lvl_3 = _get_uq_level(uq_doc, "irrigation-interval", 1)
        days_3 = _num_word(int(round(float((lvl_3 or {}).get("x", 3)))))
        c_3 = round_plain(abs(_level_median(lvl_3, "h2_uptake_change_pct", -26.64)))
        o_3 = (lvl_3 or {}).get("outputs", {}).get("h2_uptake_change_pct", {})
        lo_3 = round_plain(abs(float(o_3.get("q05", -68.34))))
        hi_3 = round_plain(abs(float(o_3.get("q95", 128.57))))
        lvl_1 = _get_uq_level(uq_doc, "irrigation-interval", 0)
        c_1 = round_plain(abs(_level_median(lvl_1, "h2_uptake_change_pct", -28.93)))
        o_1 = (lvl_1 or {}).get("outputs", {}).get("h2_uptake_change_pct", {})
        lo_1 = round_plain(abs(float(o_1.get("q05", -72.56))))
        hi_1 = round_plain(abs(float(o_1.get("q95", 136.47))))
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
            f"When extra rain does fall, the sand stays wet for longer and the microbes stay active for longer."
        )
        limits_extra_plain = (
            "Cloud seeding can only add water when suitable clouds are already present; the twin does not "
            "predict when clouds form."
        )
    elif slug == "salty-water":
        ec_8 = round_plain(float((h_lvl or {}).get("x", 8.0)))
        fresh_lvl = _get_uq_level(uq_doc, "daily-irrigation", 0)
        fresh_drop = round_plain(abs(_level_median(fresh_lvl, "h2_uptake_change_pct", -29.28)))
        alf_c = round_plain(100.0 * _level_median(h_lvl, "relative_yield_alfalfa", 0.8866))
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
            f"Watering every day with salty water gives about {core} % less hydrogen uptake "
            f"over two weeks than one fresh rain (at {ec_8} dS/m, and it {hedge}), while fresh daily "
            f"watering alone gives about {fresh_drop} % less."
        )
        answer_plain = (
            f"Watering every day with salty water gives about {core} % less hydrogen uptake "
            f"over two weeks than one fresh rain of the same size (at {ec_8} dS/m; {rng_plain}, so it {hedge}). "
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
    elif slug == "alfalfa-bacterium":
        ec_hi = f"{float((h_lvl or {}).get('x', 7.81)):g}"
        c_lo = round_plain(100.0 * float(q05 if q05 is not None else 0.8874))
        c_hi = round_plain(100.0 * float(q95 if q95 is not None else 0.9310))
        u_c = round_plain(100.0 * _level_median(h_lvl, "relative_yield_alfalfa", 0.7946))
        u_o = h_outputs_map.get("relative_yield_alfalfa", {})
        u_lo = round_plain(100.0 * float(u_o.get("q05", 0.7847)))
        u_hi = round_plain(100.0 * float(u_o.get("q95", 0.8109)))
        mod_a_lvl = _get_uq_level(uq_doc, "A_fao_fixed_leaching", l_idx)
        mod_b_lvl = _get_uq_level(uq_doc, "B_twin_physics", l_idx)
        a_val = round_plain(100.0 * _level_median(mod_a_lvl, "relative_yield_alfalfa_sa187", 0.3204))
        b_val = round_plain(100.0 * _level_median(mod_b_lvl, "relative_yield_alfalfa_sa187", 1.0128))
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
    elif slug == "hotter":
        deg_4 = _num_word(int(round(float((h_lvl or {}).get("x", 4.0))))).capitalize()
        core_dir = "downwards" if display_core < 0 else "upwards"
        soil_lvl = _get_uq_level(uq_doc, "warmer-soil-only", 4)
        deg_10 = round_plain(float((soil_lvl or {}).get("x", 10.0)))
        s_m = round_plain(abs(_level_median(soil_lvl, "h2_uptake_change_pct", 35.53)))
        s_o = (soil_lvl or {}).get("outputs", {}).get("h2_uptake_change_pct", {})
        s_rng = format_change_pct_interval(s_o.get("q05", 3.79), s_o.get("q95", 77.83))
        headline_plain = (
            f"{deg_4} degrees of warming with the same rainfall changes hydrogen uptake over two weeks "
            f"by about {core} % {core_dir} in a typical case, and uptake {hedge}."
        )
        answer_plain = (
            f"{deg_4} degrees of warming with the same rainfall changes hydrogen uptake over two weeks "
            f"by about {core} % {core_dir} in a typical case, but it {hedge} ({rng_plain}). "
            f"Warming only the sand by {deg_10} degrees, with no extra drying, raises uptake by about "
            f"{s_m} % ({s_rng})."
        )
        limits_extra_plain = (
            "Heat stress death at extreme afternoon surface temperatures is not included in the simple "
            "warming curve."
        )
    elif slug == "mulch-shade":
        sh_lvl = _get_uq_level(uq_doc, "shade-only", 2)
        sh_m = round_plain(abs(_level_median(sh_lvl, "h2_uptake_change_pct", -13.71)))
        sh_o = (sh_lvl or {}).get("outputs", {}).get("h2_uptake_change_pct", {})
        sh_rng = format_change_pct_interval(sh_o.get("q05", -27.87), sh_o.get("q95", -1.74))
        ms_lvl = _get_uq_level(uq_doc, "mulch-and-shade", 2)
        ms_m = round_plain(abs(_level_median(ms_lvl, "h2_uptake_change_pct", 28.17)))
        ms_o = (ms_lvl or {}).get("outputs", {}).get("h2_uptake_change_pct", {})
        ms_rng = format_change_pct_interval(ms_o.get("q05", -6.22), ms_o.get("q95", 70.17))
        headline_plain = (
            f"Full mulch cover raises hydrogen uptake over two weeks by about {core} % by keeping the sand "
            f"damp longer, and uptake {hedge}."
        )
        answer_plain = (
            f"Full mulch cover raises hydrogen uptake over two weeks by about {core} % ({rng_plain}) by "
            f"slowing evaporation, so uptake with mulch alone {hedge}. "
            f"Full shade alone lowers uptake by about {sh_m} % ({sh_rng}), "
            f"and combining mulch with shade raises uptake by about {ms_m} % ({ms_rng}). "
            f"Mulch keeps water in the sand without cooling the microbes, whereas shade cools the sand "
            f"and slows enzyme speed."
        )
        limits_extra_plain = (
            "Organic mulch breakdown and extra carbon food for soil microbes are not included."
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
    elif slug == "compaction":
        passes_n = int(round(float((h_lvl or {}).get("x", 50))))
        passes_word = "Fifty" if passes_n == 50 else str(passes_n)
        signed_core = f"-{core}" if display_core < 0 else f"+{core}"
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
    elif slug == "biochar":
        dose_word = _num_word(int(round(float((h_lvl or {}).get("x", 5.0)))))
        core_dir = "downwards" if display_core < 0 else "upwards"
        headline_plain = (
            f"Adding {dose_word} percent biochar by weight changes hydrogen uptake over two weeks by about "
            f"{core} % {core_dir}, and the change {hedge}."
        )
        answer_plain = (
            f"Adding {dose_word} percent biochar by weight changes hydrogen uptake over two weeks by about "
            f"{core} % {core_dir} ({rng_plain}), so the change {hedge}. "
            f"Biochar adds empty pore space, so one rain fills a smaller share of the pores while gas moves "
            f"more easily. "
            f"The two effects nearly cancel, so the net change is small and could go either way."
        )
        limits_extra_plain = (
            "Chemical effects of biochar on soil acidity, nutrients or microbe attachment are not included."
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
            else 0
        )
        ref_word = _num_word(n_refused)
        tot_word = _num_word(n_total)
        if n_refused > 0:
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
        else:
            headline_plain = (
                f"Where the soil community has no genes for eating carbon monoxide, the twin gives zero uptake, "
                f"versus about {core} {unit_plain} over two weeks where the genes are present."
            )
            answer_plain = (
                f"Where the soil community has no genes for eating carbon monoxide, the twin gives exactly zero "
                f"carbon-monoxide uptake. "
                f"Where the genes are present it gives a small uptake, about {core} {unit_plain} over two weeks "
                f"({round_plain(q05)} to {round_plain(q95)})."
            )
        limits_extra_plain = (
            "A missing gene in recovered DNA means the twin sets that pathway to zero, though very rare "
            "unsequenced microbes could still exist in the sand."
        )
    elif slug == "adding-microbes":
        prop_0 = _get_uq_level(uq_doc, "51PRr3-proportional", 0)
        prop_2 = _get_uq_level(uq_doc, "51PRr3-proportional", 2)
        v53_2 = _get_uq_level(uq_doc, "V53PRr2-proportional", 2)
        lo_start = round_plain(_level_median(prop_0, "co_uptake_umol_m2", 8.75), decimals=2)
        hi_start = round_plain(_level_median(prop_2, "co_uptake_umol_m2", 199.38))
        x_1 = round_plain(float((h_lvl or {}).get("x", 1.0)))
        x_5 = round_plain(float((v53_2 or {}).get("x", 5.0)))
        c_5 = round_plain(abs(_level_median(v53_2, "co_uptake_change_pct", 2.72)))
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
    elif slug == "air-hydrogen":
        lvl_10 = _get_uq_level(uq_doc, "h2-multiplier", 3)
        mult_10 = _num_word(int(round(float((lvl_10 or {}).get("x", 10.0)))))
        c_10 = round_plain(abs(_level_median(lvl_10, "h2_uptake_change_pct", 900.0)))
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
    elif slug == "carbon-fixation":
        c_lo = round_plain(q05)
        c_hi = round_plain(q95, decimals=2)
        closed_lvl = _get_uq_level(uq_doc, "carbon-bound-samples", 0)
        closed_c = round_plain(_level_median(closed_lvl, out_key, 0.0))
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
    else:
        raise KeyError(f"Unknown question slug: {slug}")

    limits_plain = f"{SCOPE_SENTENCE_PLAIN} {limits_extra_plain}"

    if uq_doc is None:
        raise ValueError(f"uq_doc is required to render science texts for {slug!r}")
    sci = render_science_texts(slug, uq_doc)

    return {
        "title_science": sci["title_science"],
        "headline_plain": headline_plain,
        "headline_science": sci["headline_science"],
        "answer_plain": answer_plain,
        "answer_science": sci["answer_science"],
        "limits_plain": limits_plain,
        "limits_science": sci["limits_science"],
    }
