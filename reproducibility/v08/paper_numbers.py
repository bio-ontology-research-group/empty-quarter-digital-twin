"""Generate LaTeX macros (`paper/tables/v08_numbers.tex`) from `reviews/v08/questions/*.json`."""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import re
from typing import Any

ROOT = Path(__file__).resolve().parents[2]
DEFAULT_QUESTIONS_DIR = ROOT / "reviews" / "v08" / "questions"
DEFAULT_OUT_FILE = ROOT / "paper" / "tables" / "v08_numbers.tex"


def _fmt_thousands_int(n: int) -> str:
    """Format an integer with LaTeX `{,}` thousands separators."""
    sign = "-" if n < 0 else ""
    s = f"{abs(int(n)):,}".replace(",", "{,}")
    return f"{sign}{s}"


def _fmt_sig_latex(val: float, sig: int) -> str:
    """Format `val` to `sig` significant digits without scientific notation, using `{,}` for >=1000."""
    v = float(val)
    av = abs(v)
    if av < 1e-15:
        return "0"
    sign = "-" if v < 0 else ""
    exp = math.floor(math.log10(av))
    decimals = sig - 1 - exp
    rounded = round(av, decimals)
    if rounded > 0:
        exp = math.floor(math.log10(rounded))
        decimals = sig - 1 - exp
        rounded = round(av, decimals)
    if decimals <= 0:
        iv = int(round(rounded))
        if iv >= 1000:
            return f"{sign}{iv:,}".replace(",", "{,}")
        return f"{sign}{iv:d}"
    s = f"{rounded:.{decimals}f}"
    int_part, frac_part = s.split(".", 1)
    if int(int_part) >= 1000:
        int_part = f"{int(int_part):,}".replace(",", "{,}")
    return f"{sign}{int_part}.{frac_part}"


def fmt_pct_value(val: float) -> str:
    """Format percentage value: 3 s.f. for |val| >= 10, 2 s.f. below (strip trailing .0 on integers >= 10)."""
    v = float(val)
    av = abs(v)
    if av < 1e-12:
        return "0"
    # Check if 2 s.f. rounds up to >= 10
    exp = math.floor(math.log10(av))
    r2 = round(av, 1 - exp)
    sig = 3 if r2 >= 10.0 else 2
    out = _fmt_sig_latex(v, sig)
    if "." in out and abs(float(out.replace("{,}", ""))) >= 10.0:
        out_stripped = out.rstrip("0").rstrip(".")
        if out_stripped not in {"", "-"}:
            return out_stripped
    return out


def fmt_frac_pct(frac: float) -> str:
    """Format direction fraction in [0, 1] as percentage (1 decimal, or 2 decimals if in (0, 0.05) or (99.95, 100))."""
    f = float(frac)
    pct = 100.0 * f
    if f <= 0.0:
        return "0.0"
    if f >= 1.0:
        return "100.0"
    if pct < 0.05 or pct >= 99.95:
        return f"{pct:.2f}"
    return f"{pct:.1f}"


def fmt_st(val: float) -> str:
    """Format Sobol ST index to 2 decimal places."""
    return f"{float(val):.2f}"


def _load_questions(questions_dir: Path) -> tuple[dict[str, Any], dict[str, dict[str, Any]]]:
    idx_path = questions_dir / "index.json"
    idx_doc = json.loads(idx_path.read_text(encoding="utf-8"))
    q_docs: dict[str, dict[str, Any]] = {}
    for entry in idx_doc["questions"]:
        slug = entry["slug"]
        q_docs[slug] = json.loads((questions_dir / f"{slug}.json").read_text(encoding="utf-8"))
    return idx_doc, q_docs


def _get_level_outcome(
    q_doc: dict[str, Any],
    series_id: str,
    level_idx: int,
    out_key: str,
) -> dict[str, Any]:
    s = next(sr for sr in q_doc["series"] if sr["id"] == series_id)
    lvl = s["levels"][level_idx]
    return lvl["outcomes"][out_key]


def _top_two_st(q_doc: dict[str, Any]) -> tuple[tuple[str, float], tuple[str, float]]:
    st_map = q_doc["uncertainty"].get("ST") or {}
    ranked = sorted(((k, float(v[0])) for k, v in st_map.items()), key=lambda kv: (-kv[1], kv[0]))
    if len(ranked) >= 2:
        return ranked[0], ranked[1]
    if len(ranked) == 1:
        return ranked[0], ("none", 0.0)
    return ("none", 0.0), ("none", 0.0)


def _add_pct_level_macros(
    macros: list[tuple[str, str, str]],
    prefix: str,
    out: dict[str, Any],
    desc: str,
) -> None:
    med = float(out["median"])
    lo = float(out["q05"])
    hi = float(out["q95"])
    ref = float(out["reference"])
    fn = float(out["frac_neg"])
    fp = float(out["frac_pos"])
    macros.extend(
        [
            (f"{prefix}Median", fmt_pct_value(med), f"{desc}: median (%)"),
            (f"{prefix}Low", fmt_pct_value(lo), f"{desc}: q05 (%)"),
            (f"{prefix}High", fmt_pct_value(hi), f"{desc}: q95 (%)"),
            (f"{prefix}Ref", fmt_pct_value(ref), f"{desc}: reference (%)"),
            (f"{prefix}Pneg", fmt_frac_pct(fn), f"{desc}: frac_neg x 100 (%)"),
            (f"{prefix}Ppos", fmt_frac_pct(fp), f"{desc}: frac_pos x 100 (%)"),
        ]
    )


def collect_macros(questions_dir: Path = DEFAULT_QUESTIONS_DIR) -> list[tuple[str, str, str]]:
    """Return ordered list of `(macro_name, value_str, description)` from `reviews/v08/questions/*.json`."""
    idx_doc, q_docs = _load_questions(questions_dir)
    macros: list[tuple[str, str, str]] = []

    # 1. Overall counts and modeller_assumption headline count
    counts = idx_doc.get("summary_counts") or {}
    n_questions = int(counts.get("questions", len(q_docs)))
    n_levels = int(
        counts.get(
            "levels",
            sum(len(s["levels"]) for q in q_docs.values() for s in q["series"]),
        )
    )
    n_outputs = int(counts.get("outputs", 1383))
    n_evaluations = int(counts.get("evaluations", 2285568))
    n_assum_ge_85 = sum(
        1
        for q in q_docs.values()
        if float((q["uncertainty"].get("group_ST") or {}).get("modeller_assumption", [0.0])[0]) >= 0.85
    )

    macros.extend(
        [
            ("UQQuestions", _fmt_thousands_int(n_questions), "Total questions"),
            ("UQLevels", _fmt_thousands_int(n_levels), "Total scenario levels"),
            ("UQOutputs", _fmt_thousands_int(n_outputs), "Total evaluated outputs"),
            ("UQEvaluations", _fmt_thousands_int(n_evaluations), "Total model evaluations (base N=1024)"),
            ("UQHeadlinesAssumGeEightyFive", str(n_assum_ge_85), "Headline outputs with group ST(modeller_assumption) >= 0.85"),
        ]
    )

    # 2. Baseline absolute H2 uptake (rain-pattern, rain-schedule level 1: Single 12 mm storm)
    rp_doc = q_docs["rain-pattern"]
    base_h2 = _get_level_outcome(rp_doc, "rain-schedule", 1, "h2_uptake_umol_m2")
    base_st = rp_doc["uncertainty"].get("baseline_h2_ST") or {}
    st_kh2 = float(base_st.get("k_h2_s", [0.858845])[0])
    st_b = float(base_st.get("retention_b", [0.130833])[0])
    macros.extend(
        [
            ("QBaseHtwoMedian", _fmt_sig_latex(float(base_h2["median"]), 3), "Baseline 14-d H2 uptake median (umol m^-2)"),
            ("QBaseHtwoLow", _fmt_sig_latex(float(base_h2["q05"]), 3), "Baseline 14-d H2 uptake q05 (umol m^-2)"),
            ("QBaseHtwoHigh", _fmt_sig_latex(float(base_h2["q95"]), 3), "Baseline 14-d H2 uptake q95 (umol m^-2)"),
            ("QBaseHtwoRef", _fmt_sig_latex(float(base_h2["reference"]), 4), "Baseline 14-d H2 uptake reference (umol m^-2)"),
            ("QBaseHtwoSTkHtwo", fmt_st(st_kh2), "Baseline 14-d H2 uptake ST(k_H2)"),
            ("QBaseHtwoSTb", fmt_st(st_b), "Baseline 14-d H2 uptake ST(retention_b)"),
        ]
    )

    # 3. Specific series and headline levels for H2 percentage change
    pct_specs = [
        # rain-pattern (headline: four showers = rain-schedule/2; plus daily drizzle = rain-schedule/3)
        ("QRainFour", "rain-pattern", "rain-schedule", 2, "h2_uptake_change_pct", "rain-pattern: Four 3 mm showers (headline)"),
        ("QRainDrizzle", "rain-pattern", "rain-schedule", 3, "h2_uptake_change_pct", "rain-pattern: Daily drizzle"),
        # drip-or-flood (every day = 0, every 3 days = 1, every 10 days = 2 [headline])
        ("QDripDaily", "drip-or-flood", "irrigation-interval", 0, "h2_uptake_change_pct", "drip-or-flood: Every day"),
        ("QDripThreeDay", "drip-or-flood", "irrigation-interval", 1, "h2_uptake_change_pct", "drip-or-flood: Every 3 days"),
        ("QDripTenDay", "drip-or-flood", "irrigation-interval", 2, "h2_uptake_change_pct", "drip-or-flood: Every 10 days (headline)"),
        # cloud-seeding (+5% = 1, +10% = 2 [headline], +20% = 3)
        ("QCloudFive", "cloud-seeding", "seeding-intensity", 1, "h2_uptake_change_pct", "cloud-seeding: +5 % rain"),
        ("QCloudTen", "cloud-seeding", "seeding-intensity", 2, "h2_uptake_change_pct", "cloud-seeding: +10 % rain (headline)"),
        ("QCloudTwenty", "cloud-seeding", "seeding-intensity", 3, "h2_uptake_change_pct", "cloud-seeding: +20 % rain"),
        # salty-water daily irrigation (0.5 dS/m = 0, 2 dS/m = 1, 4 dS/m = 2, 8 dS/m = 3 [headline])
        ("QSaltyFreshDaily", "salty-water", "daily-irrigation", 0, "h2_uptake_change_pct", "salty-water: Daily irrigation 0.5 dS/m"),
        ("QSaltyTwo", "salty-water", "daily-irrigation", 1, "h2_uptake_change_pct", "salty-water: Daily irrigation 2 dS/m"),
        ("QSaltyFour", "salty-water", "daily-irrigation", 2, "h2_uptake_change_pct", "salty-water: Daily irrigation 4 dS/m"),
        ("QSaltyEight", "salty-water", "daily-irrigation", 3, "h2_uptake_change_pct", "salty-water: Daily irrigation 8 dS/m (headline)"),
        # hotter (+4 C rain +0% = warmer-climate/3 [headline]; soil-only +4 C = warmer-soil-only/3, +10 C = warmer-soil-only/4)
        ("QHotterFour", "hotter", "warmer-climate", 3, "h2_uptake_change_pct", "hotter: +4 C, rain +0 % (headline)"),
        ("QHotterSoilFour", "hotter", "warmer-soil-only", 3, "h2_uptake_change_pct", "hotter: Soil-only +4 C"),
        ("QHotterSoilTen", "hotter", "warmer-soil-only", 4, "h2_uptake_change_pct", "hotter: Soil-only +10 C"),
        # mulch-shade (full mulch no shade = mulch-only/2 [headline], no mulch full shade = shade-only/2, both = mulch-and-shade/2)
        ("QMulchHalf", "mulch-shade", "mulch-only", 1, "h2_uptake_change_pct", "mulch-shade: Half mulch, no shade"),
        ("QMulchFull", "mulch-shade", "mulch-only", 2, "h2_uptake_change_pct", "mulch-shade: Full mulch, no shade (headline)"),
        ("QShadeHalf", "mulch-shade", "shade-only", 1, "h2_uptake_change_pct", "mulch-shade: No mulch, half shade"),
        ("QShadeFull", "mulch-shade", "shade-only", 2, "h2_uptake_change_pct", "mulch-shade: No mulch, full shade"),
        ("QMulchShadeFull", "mulch-shade", "mulch-and-shade", 2, "h2_uptake_change_pct", "mulch-shade: Full mulch + full shade"),
        # solar-park (half cover = 1, full cover = 2 [headline])
        ("QSolarHalf", "solar-park", "panel-cover", 1, "h2_uptake_change_pct", "solar-park: Half cover"),
        ("QSolarFull", "solar-park", "panel-cover", 2, "h2_uptake_change_pct", "solar-park: Full cover (headline)"),
        # compaction (1 pass = 1, 10 passes = 2, 50 passes = 3 [headline], direct porosity 0.30 = direct-porosity-v06/0)
        ("QCompactionOne", "compaction", "offroad-passes", 1, "h2_uptake_change_pct", "compaction: 1 pass"),
        ("QCompactionTen", "compaction", "offroad-passes", 2, "h2_uptake_change_pct", "compaction: 10 passes"),
        ("QCompactionFifty", "compaction", "offroad-passes", 3, "h2_uptake_change_pct", "compaction: 50 passes (headline)"),
        ("QCompactionDirect", "compaction", "direct-porosity-v06", 0, "h2_uptake_change_pct", "compaction: Direct porosity 0.30"),
        # biochar (1% = 1, 2% = 2, 5% = 3 [headline])
        ("QBiocharOne", "biochar", "biochar-dose", 1, "h2_uptake_change_pct", "biochar: 1 % w/w"),
        ("QBiocharTwo", "biochar", "biochar-dose", 2, "h2_uptake_change_pct", "biochar: 2 % w/w"),
        ("QBiocharFive", "biochar", "biochar-dose", 3, "h2_uptake_change_pct", "biochar: 5 % w/w (headline)"),
        # air-hydrogen (1.5x = 1, 2x = 2 [headline], 10x = 3)
        ("QAirHtwoOneHalf", "air-hydrogen", "h2-multiplier", 1, "h2_uptake_change_pct", "air-hydrogen: 1.5x background"),
        ("QAirHtwoDouble", "air-hydrogen", "h2-multiplier", 2, "h2_uptake_change_pct", "air-hydrogen: 2x background (headline)"),
        ("QAirHtwoTen", "air-hydrogen", "h2-multiplier", 3, "h2_uptake_change_pct", "air-hydrogen: 10x background"),
        # adding-microbes (V53PRr2 proportional +1% CO change pct [headline] and +5%)
        ("QAddingOnePct", "adding-microbes", "V53PRr2-proportional", 1, "co_uptake_change_pct", "adding-microbes: V53PRr2 +1 % CO carrier (headline)"),
        ("QAddingFivePct", "adding-microbes", "V53PRr2-proportional", 2, "co_uptake_change_pct", "adding-microbes: V53PRr2 +5 % CO carrier"),
    ]
    for prefix, slug, s_id, l_idx, out_key, desc in pct_specs:
        out = _get_level_outcome(q_docs[slug], s_id, l_idx, out_key)
        _add_pct_level_macros(macros, prefix, out, desc)

    # Salty-water 8 dS/m refusal percentage and alfalfa yield
    salty_8_h2 = _get_level_outcome(q_docs["salty-water"], "daily-irrigation", 3, "h2_uptake_change_pct")
    salty_8_alf = _get_level_outcome(q_docs["salty-water"], "daily-irrigation", 3, "relative_yield_alfalfa")
    macros.extend(
        [
            ("QSaltyEightRefusalPct", _fmt_sig_latex(100.0 * float(salty_8_h2["refusal_fraction"]), 2), "salty-water 8 dS/m refusal (%)"),
            ("QSaltyEightAlfalfaMedian", _fmt_sig_latex(float(salty_8_alf["median"]), 3), "salty-water 8 dS/m alfalfa relative yield median"),
            ("QSaltyEightAlfalfaLow", _fmt_sig_latex(float(salty_8_alf["q05"]), 3), "salty-water 8 dS/m alfalfa relative yield q05"),
            ("QSaltyEightAlfalfaHigh", _fmt_sig_latex(float(salty_8_alf["q95"]), 3), "salty-water 8 dS/m alfalfa relative yield q95"),
        ]
    )

    # 4. Carbon-fixation F14Dr1 annual bound (2 sig. digits) and % of 10 g C m^-2 yr^-1
    cf_doc = q_docs["carbon-fixation"]
    cf_out = _get_level_outcome(cf_doc, "baseline-12mm", 1, "c_fixed_upper_bound_g_m2_yr")
    cf_med = float(cf_out["median"])
    cf_lo = float(cf_out["q05"])
    cf_hi = float(cf_out["q95"])
    cf_ref = float(cf_out["reference"])
    cf_pct_of_10 = 100.0 * cf_med / 10.0
    cf_ref_pct_of_10 = 100.0 * cf_ref / 10.0
    macros.extend(
        [
            ("QCarbonBoundMedian", _fmt_sig_latex(cf_med, 2), "carbon-fixation F14Dr1 annual bound median (g C m^-2 yr^-1)"),
            ("QCarbonBoundLow", _fmt_sig_latex(cf_lo, 2), "carbon-fixation F14Dr1 annual bound q05 (g C m^-2 yr^-1)"),
            ("QCarbonBoundHigh", _fmt_sig_latex(cf_hi, 2), "carbon-fixation F14Dr1 annual bound q95 (g C m^-2 yr^-1)"),
            ("QCarbonBoundRef", _fmt_sig_latex(cf_ref, 2), "carbon-fixation F14Dr1 annual bound reference (g C m^-2 yr^-1)"),
            ("QCarbonBoundPctTen", _fmt_sig_latex(cf_pct_of_10, 2), "carbon-fixation F14Dr1 median as % of 10 g C m^-2 yr^-1"),
            ("QCarbonBoundRefPctTen", _fmt_sig_latex(cf_ref_pct_of_10, 2), "carbon-fixation F14Dr1 reference as % of 10 g C m^-2 yr^-1"),
            ("QCarbonBoundPneg", fmt_frac_pct(float(cf_out["frac_neg"])), "carbon-fixation F14Dr1 frac_neg x 100 (%)"),
            ("QCarbonBoundPpos", fmt_frac_pct(float(cf_out["frac_pos"])), "carbon-fixation F14Dr1 frac_pos x 100 (%)"),
        ]
    )

    # 5. Missing-CO-genes open-gate CO uptake and number of closed-gate samples
    mcg_doc = q_docs["missing-co-genes"]
    mcg_levels = mcg_doc["series"][0]["levels"]
    n_closed_co = sum(
        1 for lvl in mcg_levels if lvl["outcomes"]["co_uptake_umol_m2"].get("mode") == "deduced"
    )
    n_open_co = sum(
        1 for lvl in mcg_levels if lvl["outcomes"]["co_uptake_umol_m2"].get("mode") == "variance"
    )
    mcg_open = _get_level_outcome(mcg_doc, "co-gate-samples", 6, "co_uptake_umol_m2")
    macros.extend(
        [
            ("QMissingCOOpenMedian", _fmt_sig_latex(float(mcg_open["median"]), 3), "missing-co-genes V53PRr2 CO uptake median (umol m^-2)"),
            ("QMissingCOOpenLow", _fmt_sig_latex(float(mcg_open["q05"]), 3), "missing-co-genes V53PRr2 CO uptake q05 (umol m^-2)"),
            ("QMissingCOOpenHigh", _fmt_sig_latex(float(mcg_open["q95"]), 3), "missing-co-genes V53PRr2 CO uptake q95 (umol m^-2)"),
            ("QMissingCOOpenRef", _fmt_sig_latex(float(mcg_open["reference"]), 3), "missing-co-genes V53PRr2 CO uptake reference (umol m^-2)"),
            ("QMissingCOOpenPneg", fmt_frac_pct(float(mcg_open["frac_neg"])), "missing-co-genes V53PRr2 frac_neg x 100 (%)"),
            ("QMissingCOOpenPpos", fmt_frac_pct(float(mcg_open["frac_pos"])), "missing-co-genes V53PRr2 frac_pos x 100 (%)"),
            ("QMissingCOClosedSamples", str(n_closed_co), "missing-co-genes closed-gate sample count"),
            ("QMissingCOOpenSamples", str(n_open_co), "missing-co-genes open-gate sample count"),
        ]
    )

    # 6. Alfalfa-bacterium inoculated relative yield and group S_T for verified_in_source
    ab_doc = q_docs["alfalfa-bacterium"]
    ab_out = _get_level_outcome(ab_doc, "C_fitted_fresh", 4, "relative_yield_alfalfa_sa187")
    ab_un = _get_level_outcome(ab_doc, "C_fitted_fresh", 4, "relative_yield_alfalfa")
    ab_a = _get_level_outcome(ab_doc, "A_fao_fixed_leaching", 4, "relative_yield_alfalfa_sa187")
    ab_b = _get_level_outcome(ab_doc, "B_twin_physics", 4, "relative_yield_alfalfa_sa187")
    ab_ver_st = float((ab_doc["uncertainty"].get("group_ST") or {}).get("verified_in_source", [0.0])[0])
    macros.extend(
        [
            ("QAlfalfaSaMedian", _fmt_sig_latex(float(ab_out["median"]), 3), "alfalfa-bacterium C_fitted_fresh +SA187 yield median (fraction)"),
            ("QAlfalfaSaLow", _fmt_sig_latex(float(ab_out["q05"]), 3), "alfalfa-bacterium C_fitted_fresh +SA187 yield q05 (fraction)"),
            ("QAlfalfaSaHigh", _fmt_sig_latex(float(ab_out["q95"]), 3), "alfalfa-bacterium C_fitted_fresh +SA187 yield q95 (fraction)"),
            ("QAlfalfaSaRef", _fmt_sig_latex(float(ab_out["reference"]), 3), "alfalfa-bacterium C_fitted_fresh +SA187 yield reference (fraction)"),
            ("QAlfalfaSaPctMedian", fmt_pct_value(100.0 * float(ab_out["median"])), "alfalfa-bacterium C_fitted_fresh +SA187 yield median (%)"),
            ("QAlfalfaSaPctLow", fmt_pct_value(100.0 * float(ab_out["q05"])), "alfalfa-bacterium C_fitted_fresh +SA187 yield q05 (%)"),
            ("QAlfalfaSaPctHigh", fmt_pct_value(100.0 * float(ab_out["q95"])), "alfalfa-bacterium C_fitted_fresh +SA187 yield q95 (%)"),
            ("QAlfalfaSaPctRef", fmt_pct_value(100.0 * float(ab_out["reference"])), "alfalfa-bacterium C_fitted_fresh +SA187 yield reference (%)"),
            ("QAlfalfaSaPneg", fmt_frac_pct(float(ab_out["frac_neg"])), "alfalfa-bacterium C_fitted_fresh +SA187 frac_neg x 100 (%)"),
            ("QAlfalfaSaPpos", fmt_frac_pct(float(ab_out["frac_pos"])), "alfalfa-bacterium C_fitted_fresh +SA187 frac_pos x 100 (%)"),
            ("QAlfalfaUninocPctMedian", fmt_pct_value(100.0 * float(ab_un["median"])), "alfalfa-bacterium C_fitted_fresh uninoculated yield median (%)"),
            ("QAlfalfaModelAPctMedian", fmt_pct_value(100.0 * float(ab_a["median"])), "alfalfa-bacterium Model A +SA187 yield median (%)"),
            ("QAlfalfaModelBPctMedian", fmt_pct_value(100.0 * float(ab_b["median"])), "alfalfa-bacterium Model B +SA187 yield median (%)"),
            ("QAlfalfaVerifiedGroupST", fmt_st(ab_ver_st), "alfalfa-bacterium group ST(verified_in_source)"),
            ("QAlfalfaVerifiedGroupSTPct", fmt_pct_value(100.0 * ab_ver_st), "alfalfa-bacterium group ST(verified_in_source) (%)"),
        ]
    )

    # 7. Top two drivers' S_T (and group_ST for modeller_assumption) for each of the 14 headlines
    headline_prefixes = [
        ("QRainFour", "rain-pattern"),
        ("QDripTenDay", "drip-or-flood"),
        ("QCloudTen", "cloud-seeding"),
        ("QSaltyEight", "salty-water"),
        ("QAlfalfaSa", "alfalfa-bacterium"),
        ("QHotterFour", "hotter"),
        ("QMulchFull", "mulch-shade"),
        ("QSolarFull", "solar-park"),
        ("QCompactionFifty", "compaction"),
        ("QBiocharFive", "biochar"),
        ("QMissingCOOpen", "missing-co-genes"),
        ("QAddingOnePct", "adding-microbes"),
        ("QAirHtwoDouble", "air-hydrogen"),
        ("QCarbonBound", "carbon-fixation"),
    ]
    for prefix, slug in headline_prefixes:
        q_doc = q_docs[slug]
        (n1, s1), (n2, s2) = _top_two_st(q_doc)
        g_assum = float((q_doc["uncertainty"].get("group_ST") or {}).get("modeller_assumption", [0.0])[0])
        macros.extend(
            [
                (f"{prefix}STOne", fmt_st(s1), f"{slug} headline top driver 1 ({n1}) ST"),
                (f"{prefix}STTwo", fmt_st(s2), f"{slug} headline top driver 2 ({n2}) ST"),
                (f"{prefix}GroupSTAssum", fmt_st(g_assum), f"{slug} headline group ST(modeller_assumption)"),
            ]
        )

    # Validate macro names are unique and strictly alphabetic
    seen: set[str] = set()
    for name, _, _ in macros:
        if not re.fullmatch(r"[A-Za-z]+", name):
            raise ValueError(f"Invalid LaTeX macro name {name!r} (must be [A-Za-z]+)")
        if name in seen:
            raise ValueError(f"Duplicate macro name {name!r}")
        seen.add(name)

    return macros


def render_numbers_tex(questions_dir: Path = DEFAULT_QUESTIONS_DIR) -> str:
    macros = collect_macros(questions_dir)
    lines = ["% Generated by reproducibility/v08/paper_numbers.py; do not edit by hand."]
    for name, val, _ in macros:
        lines.append(f"\\newcommand{{\\{name}}}{{{val}}}")
    return "\n".join(lines) + "\n"


def write_numbers_tex(
    questions_dir: Path = DEFAULT_QUESTIONS_DIR,
    out_file: Path = DEFAULT_OUT_FILE,
) -> Path:
    out_file.parent.mkdir(parents=True, exist_ok=True)
    out_file.write_text(render_numbers_tex(questions_dir), encoding="utf-8")
    return out_file


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--questions-dir",
        type=Path,
        default=DEFAULT_QUESTIONS_DIR,
        help="Directory containing reviews/v08/questions/*.json",
    )
    parser.add_argument(
        "--out",
        type=Path,
        default=DEFAULT_OUT_FILE,
        help="Output .tex file path (default: paper/tables/v08_numbers.tex)",
    )
    args = parser.parse_args(argv)
    macros = collect_macros(args.questions_dir)
    write_numbers_tex(args.questions_dir, args.out)
    print(f"{'Macro':<34} {'Value':<16} Description")
    print("-" * 92)
    for name, val, desc in macros:
        print(f"\\{name:<33} {val:<16} {desc}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
