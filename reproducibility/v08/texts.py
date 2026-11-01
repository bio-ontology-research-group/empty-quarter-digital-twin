"""Plain and scientific text templates and metadata for EQ-BioDT v0.8 question files."""
from __future__ import annotations

import re
from typing import Any

JARGON = re.compile(
    r"[\u0370-\u03ff\u1f00-\u1fff\u00b5]"
    r"|\bdo\s*\("
    r"|\bSobol\b"
    r"|\bontolog\w*\b"
    r"|\bMAGs?\b"
    r"|\bEC_"
    r"|\bECe\b"
    r"|\bMaas\b"
    r"|\bCBB\b"
    r"|\b[A-Za-z0-9]+_[A-Za-z0-9_]+\b"
    r"|\b(?:umol|nmol|pmol|mmol|mol|ug|mg|ng|uL|mL|MPa|kPa|Pa|cm|km|nm|m3|cm3|yr)\b"
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
    "alfalfa-bacterium": "Unsalted, uninoculated alfalfa yield",
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

# (series_id, level_index, output_key) for the primary headline of each question
HEADLINE_SPEC: dict[str, tuple[str, int, str]] = {
    "rain-pattern": ("rain-schedule", 1, "h2_uptake_umol_m2"),
    "drip-or-flood": ("irrigation-interval", 2, "h2_uptake_umol_m2"),
    "cloud-seeding": ("seeding-intensity", 2, "h2_uptake_contrast_umol_m2"),
    "salty-water": ("daily-irrigation", 3, "h2_uptake_contrast_umol_m2"),
    "alfalfa-bacterium": ("C_fitted_fresh", 4, "relative_yield_alfalfa_sa187"),
    "hotter": ("warmer-climate", 3, "h2_uptake_contrast_umol_m2"),
    "mulch-shade": ("mulch-only", 2, "h2_uptake_contrast_umol_m2"),
    "solar-park": ("panel-cover", 2, "h2_uptake_contrast_umol_m2"),
    "compaction": ("offroad-passes", 3, "h2_uptake_contrast_umol_m2"),
    "biochar": ("biochar-dose", 3, "h2_uptake_contrast_umol_m2"),
    "missing-co-genes": ("co-gate-samples", 0, "co_uptake_umol_m2"),
    "adding-microbes": ("51PRr3-gate", 1, "co_uptake_contrast_umol_m2"),
    "air-hydrogen": ("h2-multiplier", 2, "h2_uptake_contrast_umol_m2"),
    "carbon-fixation": ("baseline-12mm", 1, "c_fixed_upper_bound_g_m2_yr"),
}

OUTCOME_META: dict[str, dict[str, str]] = {
    "h2_uptake_umol_m2": {
        "unit": "µmol m⁻²",
        "unit_plain": "micromoles per m²",
        "label_plain": "Hydrogen taken up over 14 days",
        "label_science": "14-day cumulative H₂ uptake",
        "provenance": "simulated",
    },
    "h2_uptake_contrast_umol_m2": {
        "unit": "µmol m⁻²",
        "unit_plain": "micromoles per m²",
        "label_plain": "Change in hydrogen uptake over 14 days",
        "label_science": "14-day H₂ uptake contrast vs baseline",
        "provenance": "simulated",
    },
    "co_uptake_umol_m2": {
        "unit": "µmol m⁻²",
        "unit_plain": "micromoles per m²",
        "label_plain": "Carbon monoxide taken up over 14 days",
        "label_science": "14-day cumulative CO uptake",
        "provenance": "simulated",
    },
    "co_uptake_contrast_umol_m2": {
        "unit": "µmol m⁻²",
        "unit_plain": "micromoles per m²",
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
        "unit": "1",
        "unit_plain": "fraction of unsalted yield",
        "label_plain": "Date palm yield (share of unsalted)",
        "label_science": "Maas-Hoffman relative yield (date palm)",
        "provenance": "external",
    },
    "relative_yield_barley": {
        "unit": "1",
        "unit_plain": "fraction of unsalted yield",
        "label_plain": "Barley yield (share of unsalted)",
        "label_science": "Maas-Hoffman relative yield (barley)",
        "provenance": "external",
    },
    "relative_yield_alfalfa": {
        "unit": "1",
        "unit_plain": "fraction of unsalted yield",
        "label_plain": "Alfalfa yield without bacterium (share of unsalted)",
        "label_science": "Maas-Hoffman relative yield (alfalfa, uninoculated)",
        "provenance": "external",
    },
    "relative_yield_tomato": {
        "unit": "1",
        "unit_plain": "fraction of unsalted yield",
        "label_plain": "Tomato yield (share of unsalted)",
        "label_science": "Maas-Hoffman relative yield (tomato)",
        "provenance": "external",
    },
    "relative_yield_alfalfa_sa187": {
        "unit": "1",
        "unit_plain": "fraction of unsalted yield",
        "label_plain": "Alfalfa yield with desert bacterium (share of unsalted)",
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
        "name_plain": "Adding microbe strains at a chosen community share",
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
    """Format a number for plain-language display and cross-mode headline matching."""
    if x is None:
        return "none"
    val = float(x)
    if abs(val) < 1e-12:
        return "0"
    ax = abs(val)
    if ax >= 1000.0:
        return f"{round(val):,}"
    if ax >= 100.0:
        return f"{round(val):d}"
    if ax >= 10.0:
        s = f"{val:.1f}".rstrip("0").rstrip(".")
        return "0" if s in {"-0", ""} else s
    if ax >= 1.0:
        s = f"{val:.2f}".rstrip("0").rstrip(".")
        return "0" if s in {"-0", ""} else s
    s = f"{val:.3f}".rstrip("0").rstrip(".")
    return "0" if s in {"-0", ""} else s


def science_central(q: dict[str, Any]) -> float:
    """Return the exact central headline value from a built question dict."""
    if "headline_value" in q and q["headline_value"] is not None:
        return float(q["headline_value"])
    s_id, l_idx, out_key = HEADLINE_SPEC[q["slug"]]
    s = next(sr for sr in q["series"] if sr["id"] == s_id)
    return float(s["levels"][l_idx]["outcomes"][out_key]["value"])


def plain_number(q: dict[str, Any]) -> str:
    """Extract the core rounded headline number from `q['answer_plain']`."""
    m = re.search(
        r"(?<![A-Za-z0-9_])-?\d[\d,]*(?:\.\d+)?(?![A-Za-z0-9_])",
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
        return f"{round_plain(float(x))} dS/m water"
    if slug == "adding-microbes":
        return f"{round_plain(float(x))} % added share"
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


def _format_interval_plain(q05: float | None, q95: float | None, unit_plain: str) -> str:
    if q05 is None or q95 is None:
        return f"range not available in {unit_plain}"
    return f"range {round_plain(q05)} to {round_plain(q95)} {unit_plain}"


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
        s1_est = s1_vals[0]
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
            "modeller_assumption = 0.00 / 0.00, source_not_checked = 0.00 / 0.00."
        )
    order = ("verified_in_source", "modeller_assumption", "source_not_checked")
    items = []
    for k in order:
        if k in group_st:
            gst = group_st[k][0]
            gs1 = group_s1.get(k, [0.0])[0]
            items.append(f"{k} ST = {gst:.2f} (S1 = {gs1:.2f})")
    return "Status share of variance: " + ", ".join(items) + "."


def render_question_texts(
    slug: str,
    central_val: float,
    outcome_uq: dict[str, Any],
    max_refusal_fraction: float,
) -> dict[str, str]:
    """Render plain and scientific answer, headline, and limits strings for `slug`."""
    core = round_plain(central_val)
    q05 = outcome_uq.get("q05")
    q25 = outcome_uq.get("q25")
    q75 = outcome_uq.get("q75")
    q95 = outcome_uq.get("q95")
    median = outcome_uq.get("median")
    mode = str(outcome_uq.get("mode", "variance"))
    ref_frac = float(outcome_uq.get("refusal_fraction", 0.0))

    _, _, out_key = HEADLINE_SPEC[slug]
    meta = OUTCOME_META[out_key]
    unit_sci = meta["unit"]
    unit_plain = meta["unit_plain"]
    rng_plain = _format_interval_plain(q05, q95, unit_plain)

    spans_zero = q05 is not None and q95 is not None and (float(q05) < 0.0 < float(q95))
    has_refusal = max_refusal_fraction >= 0.05 or ref_frac >= 0.05

    if slug == "rain-pattern":
        headline_plain = (
            f"One heavy rain gives {core} {unit_plain} of hydrogen uptake over two weeks, "
            f"versus 405 for four smaller showers and 2.87 for daily drizzle."
        )
        answer_plain = (
            f"One heavy downpour gives {core} {unit_plain} of hydrogen uptake over two weeks "
            f"({rng_plain}), compared with 405 for four smaller showers and 2.87 for daily drizzle. "
            f"A single soaking pushes water deeper into the sand so microbes stay active for more days "
            f"before evaporation dries the surface."
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
            f"Watering once every ten days gives {core} {unit_plain} of hydrogen uptake over two weeks, "
            f"compared with 680 for daily drip."
        )
        answer_plain = (
            f"Watering once every ten days gives {core} {unit_plain} of hydrogen uptake over two weeks "
            f"({rng_plain}), compared with 692 when watering every three days and 680 when watering every day. "
            f"Infrequent deeper watering loses a smaller share to daily surface evaporation than tiny daily doses."
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
            f"Adding ten percent more rain on rainy days changes hydrogen uptake by {core} {unit_plain} "
            f"over two weeks."
        )
        answer_plain = (
            f"Adding a ten percent rain boost on rainy days changes hydrogen uptake by {core} {unit_plain} "
            f"over two weeks ({rng_plain}), raising total uptake from 1,102 to 1,242 {unit_plain}. "
            f"Extra water on wet days keeps the sand moist past the end of the first week."
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
            f"Daily watering with salty water changes hydrogen uptake by {core} {unit_plain} at 8 dS/m; "
            f"the model cannot tell for such salty water when pore salt passes the dilute limit."
        )
        answer_plain = (
            f"Daily watering with brackish water changes hydrogen uptake by {core} {unit_plain} "
            f"over two weeks at 8 dS/m ({rng_plain}) and leaves alfalfa at 88.7 % of its unsalted yield; at high salt "
            f"levels, the model cannot tell for such salty water because the pore water passes the dilute-salt limit."
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
            f"Coating alfalfa seeds with the desert bacterium gives {core} of unsalted crop yield "
            f"at 7.81 dS/m salty water, up from 0.796 without it."
        )
        answer_plain = (
            f"Coating alfalfa seeds with the desert bacterium gives {core} of unsalted crop yield "
            f"({rng_plain}) in salty irrigation water at 7.81 dS/m, compared with 0.796 without the bacterium. "
            f"Across the three salt-accumulation formulas, the bacterium lifts yield when salt stress is moderate "
            f"to high, though at extreme salt levels the model cannot tell for such salty water in the twin-physics branch."
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
            f"Four degrees of climate warming changes hydrogen uptake by {core} {unit_plain} "
            f"over two weeks when rain stays the same."
        )
        answer_plain = (
            f"Four degrees of climate warming changes hydrogen uptake by {core} {unit_plain} over two weeks "
            f"({rng_plain}), while warming the sand alone by ten degrees raises uptake by 762 {unit_plain}. "
            f"Warmer sand speeds up microbial enzymes while water is present, but higher air temperature also "
            f"dries the sand faster."
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
            f"Full mulch cover increases hydrogen uptake by {core} {unit_plain} over two weeks "
            f"by keeping the sand damp longer."
        )
        answer_plain = (
            f"Full mulch cover increases hydrogen uptake by {core} {unit_plain} over two weeks "
            f"({rng_plain}) by slowing evaporation, whereas full shade alone changes uptake by -273 {unit_plain} "
            f"and combining mulch with shade gives 155 {unit_plain}. Mulch keeps water in the sand without cooling "
            f"the microbes, whereas shade cools the sand and slows enzyme speed."
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
            f"Full solar-panel cover changes hydrogen uptake by {core} {unit_plain} in the central run, "
            f"and uptake may go up or down across uncertain model numbers."
        )
        answer_plain = (
            f"Full solar-panel cover changes hydrogen uptake by {core} {unit_plain} in the central run "
            f"({rng_plain}), and the change may go up or down depending on how strongly panel shade cuts "
            f"evaporation and cools the sand versus how much panel-washing water drips onto the ground."
        )
        limits_extra_plain = (
            "Wind shifts under panel rows and runoff strips along panel edges are simplified to an even average."
        )
        limits_extra_sci = (
            "Homogenises under-panel microclimate and wash-water drips spatially; competing cooling (Q10) and "
            "moisture-retention effects yield sign ambiguity across parameter draws."
        )
    elif slug == "compaction":
        headline_plain = (
            f"Fifty vehicle passes change hydrogen uptake by {core} {unit_plain} in the central run, "
            f"and uptake may go up or down across the uncertainty range."
        )
        answer_plain = (
            f"Fifty vehicle passes change hydrogen uptake by {core} {unit_plain} in the central run "
            f"({rng_plain}), and the change may go up or down because squeezing pore space slows gas movement "
            f"through air spaces while changing how tightly the sand holds water."
        )
        limits_extra_plain = (
            "Crushing of soil crusts, plants or burrows by tyres is not modelled beyond the drop in pore space."
        )
        limits_extra_sci = (
            "Maps traffic passes to bulk density and porosity only, holding retention shape b and air-entry "
            "potentialpsi_entry uncoupled from pore-size redistribution."
        )
    elif slug == "biochar":
        headline_plain = (
            f"Adding five percent biochar by weight changes hydrogen uptake by {core} {unit_plain} "
            f"over two weeks."
        )
        answer_plain = (
            f"Adding five percent biochar by weight changes hydrogen uptake by {core} {unit_plain} over two weeks "
            f"({rng_plain}). Increasing total pore space without adding extra water lowers the fraction of pores "
            f"filled with water after a fixed 12 mm rain."
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
            f"Carbon monoxide uptake is {core} {unit_plain} when the soil community lacks the "
            f"carbon-monoxide genes; the model cannot tell for such incomplete sample records."
        )
        answer_plain = (
            f"Carbon monoxide uptake is {core} {unit_plain} ({rng_plain}) in sample 51PRr3 where the "
            f"carbon-monoxide enzyme genes are absent, compared with 144 {unit_plain} in samples that carry them; "
            f"when read-mapping weights are missing, the model cannot tell for such incomplete sample records."
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
            f"Adding carrier microbes increases carbon monoxide uptake by {core} {unit_plain} at 1 % share "
            f"under the on-off gene gate, or 77.4 when speed scales with share."
        )
        answer_plain = (
            f"Adding carrier microbes increases carbon monoxide uptake by {core} {unit_plain} "
            f"over two weeks at 1 % share ({rng_plain}) under the on-off gene gate, or by 77.4 {unit_plain} when reaction "
            f"speed scales with carrier abundance. Adding a microbe that lacks the carbon-monoxide genes leaves "
            f"uptake at 0 {unit_plain}."
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
            f"Doubling hydrogen in the air increases soil hydrogen uptake by {core} {unit_plain} "
            f"over two weeks."
        )
        answer_plain = (
            f"Doubling hydrogen in the air increases soil hydrogen uptake by {core} {unit_plain} over two weeks "
            f"({rng_plain}), doubling total uptake from 1,102 to 2,205 {unit_plain}. Because the reaction in the "
            f"twin is first-order in air hydrogen, uptake scales in lockstep with the air concentration."
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
            f"Hydrogen-eating microbes can fix at most {core} {unit_plain} in sample F14Dr1 after a 12 mm rain "
            f"({rng_plain}), and 0 {unit_plain} in sample 14PRr2 where the carbon-fixing genes are absent. "
            f"This is an upper limit if every bit of captured hydrogen energy went into building cell carbon."
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

    if spans_zero and "may go up or down" not in answer_plain.lower():
        answer_plain += " Across the uncertainty range, the answer may go up or down."
    if has_refusal and "the model cannot tell for such" not in answer_plain.lower():
        cond = "salty water" if slug in {"salty-water", "alfalfa-bacterium"} else "conditions outside the model domain"
        answer_plain += f" At extreme settings, the model cannot tell for such {cond}."

    limits_plain = f"{SCOPE_SENTENCE_PLAIN} {limits_extra_plain}"

    # Scientific headline and answer
    if q05 is not None and q95 is not None:
        ci90_str = f"[{q05:.4g}, {q95:.4g}]"
        iqr_str = f"[{q25:.4g}, {q75:.4g}]" if (q25 is not None and q75 is not None) else "n/a"
        med_str = f"{median:.4g}" if median is not None else "n/a"
    else:
        ci90_str = "[n/a, n/a]"
        iqr_str = "[n/a, n/a]"
        med_str = "n/a"

    headline_science = (
        f"Central {out_key} = {central_val:.4g} {unit_sci} (rounded: {core}; "
        f"median {med_str}, 90 % interval {ci90_str} {unit_sci}, mode = {mode}, "
        f"refusal = {ref_frac * 100:.1f} %)."
    )

    sobol_str = _top_two_sobol_text(
        outcome_uq.get("ST", {}),
        outcome_uq.get("S1", {}),
        mode,
    )
    group_str = _group_share_text(
        outcome_uq.get("group_ST", {}),
        outcome_uq.get("group_S1", {}),
        mode,
    )

    answer_science = (
        f"Central estimate: {out_key} = {central_val:.6g} {unit_sci} (rounded headline: {core} {unit_sci}; "
        f"median = {med_str} {unit_sci}, 50 % IQR = {iqr_str} {unit_sci}, joint 90 % interval = {ci90_str} {unit_sci}; "
        f"mode = {mode}, refusal fraction = {ref_frac * 100:.1f} %, max level refusal across question = {max_refusal_fraction * 100:.1f} %). "
        f"{sobol_str} {group_str} "
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
