"""Fourteen question specifications for EQ-BioDT v0.8 as replayable composition chains.

Unifies the six v0.6 process-model scenarios and the twelve v0.7 external-model composition
experiments into 14 public/scientific questions across five themes (water, salt, heat, soil,
microbes). Every level is a LevelSpec evaluated by `evaluate()`, which calls `composition.run_chain`
with `envelope=False`.
"""
from __future__ import annotations

from copy import deepcopy
from dataclasses import dataclass, field
from typing import Literal

from .community import DATA
from .composition import load_constants, run_chain, twin_request
from .external_models import MODELS
from .process import Parameters, State

Theme = Literal["water", "salt", "heat", "soil", "microbes"]
LevelKind = Literal["numeric", "categorical", "sample"]

_REG = load_constants()
_NACL_PER_DSM = float(_REG["nacl_mol_l_per_dS_m"]["value"])

DEFAULT_FORCING: list[dict[str, float]] = [
    {"water_mm": 12.0 if d == 0 else 0.0, "nacl_mol_l": 0.0, "evaporation_mm": 1.0, "temp_c": 25.0}
    for d in range(14)
]

SCHEDULES_V06: dict[str, list[float]] = {
    "dry": [0.0] * 14,
    "single": [12.0 if d == 0 else 0.0 for d in range(14)],
    "four": [3.0 if d in (0, 3, 6, 9) else 0.0 for d in range(14)],
    "daily": [12.0 / 14.0] * 14,
}


@dataclass(frozen=True)
class LevelSpec:
    label: str
    x: float | str
    action: dict
    upstream: list[str]
    twin: dict
    downstream: list[str]
    intervened_params: frozenset[str] = field(default_factory=frozenset)
    sample: str | None = None


@dataclass(frozen=True)
class Series:
    id: str
    label_plain: str
    label_science: str
    levels: list[LevelSpec]


@dataclass(frozen=True)
class Question:
    slug: str
    theme: Theme
    title_plain: str
    title_science: str
    level_kind: LevelKind
    x_label_plain: str
    x_label_science: str
    series: list[Series]
    headline: list[str]
    baseline_label: str


def _forcing_water_iv(waters: list[float]) -> dict:
    return {"forcing": {str(d): {"water_mm": float(w)} for d, w in enumerate(waters)}}


def _forcing_temp_iv(temp_c: float) -> dict:
    return {"forcing": {str(d): {"temp_c": float(temp_c)} for d in range(14)}}


def _build_questions() -> dict[str, Question]:
    q_list: list[Question] = [
        # 1. rain-pattern (water) — v0.6 schedules
        Question(
            slug="rain-pattern",
            theme="water",
            title_plain="One heavy rain, four showers or daily drizzle?",
            title_science="Rainfall schedule contrast at fixed 12 mm total depth",
            level_kind="categorical",
            x_label_plain="Rain pattern",
            x_label_science="14-day rainfall schedule",
            baseline_label="One 12 mm rain on day 1, no treatment",
            headline=["h2_uptake_umol_m2", "h2_uptake_contrast_umol_m2"],
            series=[
                Series(
                    id="rain-schedule",
                    label_plain="Rain pattern",
                    label_science="Reference 14-day rainfall schedules (v0.6)",
                    levels=[
                        LevelSpec("No rain (0 mm)", "dry", {}, [], {"intervention": _forcing_water_iv(SCHEDULES_V06["dry"])}, []),
                        LevelSpec("One 12 mm storm", "single", {}, [], {"intervention": _forcing_water_iv(SCHEDULES_V06["single"])}, []),
                        LevelSpec("Four 3 mm showers", "four", {}, [], {"intervention": _forcing_water_iv(SCHEDULES_V06["four"])}, []),
                        LevelSpec("Daily drizzle (12 mm total)", "daily", {}, [], {"intervention": _forcing_water_iv(SCHEDULES_V06["daily"])}, []),
                    ],
                )
            ],
        ),
        # 2. drip-or-flood (water) — v0.7 exp02_drip_vs_flood
        Question(
            slug="drip-or-flood",
            theme="water",
            title_plain="Drip or flood irrigation?",
            title_science="Equal-total irrigation split across delivery intervals (exp02)",
            level_kind="numeric",
            x_label_plain="Days between waterings",
            x_label_science="Irrigation interval (days)",
            baseline_label="One 12 mm rain on day 1, no irrigation",
            headline=["h2_uptake_umol_m2", "h2_uptake_contrast_umol_m2"],
            series=[
                Series(
                    id="irrigation-interval",
                    label_plain="Days between waterings",
                    label_science="Equal seasonal water demand split by interval",
                    levels=[
                        LevelSpec(
                            label=f"every {k} d",
                            x=float(k),
                            action={"irrigation_interval_days": k, "irrigation_nacl_mol_l": 0.0},
                            upstream=["crop_water_demand_generic", "irrigation_schedule"],
                            twin={},
                            downstream=[],
                        )
                        for k in (1, 3, 10)
                    ],
                )
            ],
        ),
        # 3. cloud-seeding (water) — v0.7 exp09_cloud_seeding
        Question(
            slug="cloud-seeding",
            theme="water",
            title_plain="What if cloud seeding adds rain?",
            title_science="Relative rainfall enhancement on existing rain days (exp09)",
            level_kind="numeric",
            x_label_plain="Extra rain from seeding (%)",
            x_label_science="Seeding intensity (1 = central +10% rain)",
            baseline_label="One 12 mm rain on day 1, unseeded",
            headline=["h2_uptake_umol_m2", "h2_uptake_contrast_umol_m2"],
            series=[
                Series(
                    id="seeding-intensity",
                    label_plain="Cloud seeding",
                    label_science="Seeding intensity on existing rain events",
                    levels=[
                        LevelSpec(
                            label=f"intensity {s}",
                            x=float(s),
                            action={"seeding_intensity": s},
                            upstream=["cloud_seeding"],
                            twin={},
                            downstream=[],
                        )
                        for s in (0, 0.5, 1, 2)
                    ],
                )
            ],
        ),
        # 4. salty-water (salt) — v0.6 salty rain + v0.7 exp01_brackish_irrigation
        Question(
            slug="salty-water",
            theme="salt",
            title_plain="What if the water is salty?",
            title_science="Saline rain and daily brackish irrigation across water conductivities (v0.6 + exp01)",
            level_kind="numeric",
            x_label_plain="Water saltiness (dS/m)",
            x_label_science="Irrigation or rain electrical conductivity EC_w (dS m^-1)",
            baseline_label="One 12 mm fresh rain on day 1, no irrigation",
            headline=[
                "h2_uptake_contrast_umol_m2",
                "relative_yield_date_palm",
                "relative_yield_barley",
                "relative_yield_alfalfa",
                "relative_yield_tomato",
            ],
            series=[
                Series(
                    id="salty-rain-once",
                    label_plain="One salty rain (12 mm)",
                    label_science="Single 12 mm saline rain pulse (twin only, v0.6)",
                    levels=[
                        LevelSpec(
                            label=f"{nacl_m:g} M NaCl ({nacl_m / _NACL_PER_DSM:g} dS/m)",
                            x=float(nacl_m / _NACL_PER_DSM),
                            action={},
                            upstream=[],
                            twin={"intervention": {"forcing": {"0": {"nacl_mol_l": float(nacl_m)}}}},
                            downstream=[],
                        )
                        for nacl_m in (0.0, 0.02, 0.05, 0.08, 0.12)
                    ],
                ),
                Series(
                    id="daily-irrigation",
                    label_plain="Daily salty irrigation",
                    label_science="Daily brackish irrigation + FAO-29 / Maas-Hoffman crop yields (exp01)",
                    levels=[
                        LevelSpec(
                            label=f"EC_w={ec} dS/m",
                            x=float(ec),
                            action={"ec_w_dS_m": ec, "irrigation_interval_days": 1},
                            upstream=["ec_to_nacl", "crop_water_demand_generic", "irrigation_schedule"],
                            twin={},
                            downstream=[
                                "twin_solution_to_ece",
                                "maas_hoffman_date_palm",
                                "maas_hoffman_barley",
                                "maas_hoffman_alfalfa",
                                "maas_hoffman_tomato",
                            ],
                        )
                        for ec in (0.5, 2, 4, 8)
                    ],
                ),
            ],
        ),
        # 5. alfalfa-bacterium (salt) — v0.7 exp12_inoculated_alfalfa_salinity
        Question(
            slug="alfalfa-bacterium",
            theme="salt",
            title_plain="Can a desert bacterium help alfalfa in salty water?",
            title_science="Enterobacter sp. SA187 alfalfa salinity response across three root-zone salinity models (exp12)",
            level_kind="numeric",
            x_label_plain="Water saltiness (dS/m)",
            x_label_science="Irrigation water electrical conductivity EC_w (dS m^-1)",
            baseline_label="Normal, uninoculated alfalfa yield",
            headline=["relative_yield_alfalfa_sa187", "relative_yield_alfalfa"],
            series=[
                Series(
                    id=mname,
                    label_plain=plain_lbl,
                    label_science=sci_lbl,
                    levels=[
                        LevelSpec(
                            label=f"{mname}; EC_w={ec}",
                            x=float(ec),
                            action={"ec_w_dS_m": ec, **({"irrigation_interval_days": 1} if ups else {})},
                            upstream=list(ups),
                            twin={},
                            downstream=list(downs),
                        )
                        for ec in (0.5, 2, 3.12, 4, 7.81, 8, 12)
                    ],
                )
                for mname, plain_lbl, sci_lbl, ups, downs in (
                    (
                        "A_fao_fixed_leaching",
                        "Standard farm drainage rule (published formula)",
                        "Model A: FAO-29 fixed leaching factor (1.5x EC_w)",
                        [],
                        ["fao29_leaching", "maas_hoffman_alfalfa", "sa187_benefit"],
                    ),
                    (
                        "B_twin_physics",
                        "Desert sand water-and-salt balance (twin)",
                        "Model B: Twin bare-sand salt balance + FAO-29 conversion",
                        ["ec_to_nacl", "crop_water_demand_alfalfa", "irrigation_schedule"],
                        ["twin_solution_to_ece", "maas_hoffman_alfalfa", "sa187_benefit"],
                    ),
                    (
                        "C_fitted_fresh",
                        "Calibrated to field trial (fresh weight)",
                        "Model C (fresh): Leaching factor fitted to Hada Al-Sham fresh-weight mock ratio",
                        [],
                        ["fitted_leaching_fresh", "maas_hoffman_alfalfa", "sa187_benefit"],
                    ),
                    (
                        "C_fitted_dry",
                        "Calibrated to field trial (dry weight)",
                        "Model C (dry): Leaching factor fitted to Hada Al-Sham dry-weight mock ratio",
                        [],
                        ["fitted_leaching_dry", "maas_hoffman_alfalfa", "sa187_benefit"],
                    ),
                )
            ],
        ),
        # 6. hotter (heat) — v0.6 35 C + v0.7 exp11_warming
        Question(
            slug="hotter",
            theme="heat",
            title_plain="What if it gets hotter?",
            title_science="Soil-only Q10 warming vs coupled climate warming with evaporation scaling (v0.6 + exp11)",
            level_kind="numeric",
            x_label_plain="Warming (°C)",
            x_label_science="Temperature increase ΔT (°C)",
            baseline_label="One 12 mm rain on day 1 at 25 °C",
            headline=["h2_uptake_umol_m2", "h2_uptake_contrast_umol_m2"],
            series=[
                Series(
                    id="warmer-soil-only",
                    label_plain="Warmer soil only",
                    label_science="Direct soil temperature increase at fixed evaporation (v0.6)",
                    levels=[
                        LevelSpec(
                            label=f"+{dt:g} C soil ({25.0 + dt:g} C)",
                            x=float(dt),
                            action={},
                            upstream=[],
                            twin={"intervention": _forcing_temp_iv(25.0 + dt)},
                            downstream=[],
                        )
                        for dt in (0.0, 1.5, 2.0, 4.0, 10.0)
                    ],
                ),
                Series(
                    id="warmer-climate",
                    label_plain="Warmer air and faster drying",
                    label_science="Coupled soil warming and evaporation scaling at unchanged rain (exp11)",
                    levels=[
                        LevelSpec(
                            label=f"+{w} C, rain +0%",
                            x=float(w),
                            action={"warming_c": w, "rain_change_frac": 0},
                            upstream=["climate_warming"],
                            twin={},
                            downstream=[],
                        )
                        for w in (0, 1.5, 2, 4)
                    ],
                ),
            ],
        ),
        # 7. mulch-shade (heat) — v0.7 exp03_mulch_shade
        Question(
            slug="mulch-shade",
            theme="heat",
            title_plain="Mulch or shade nets?",
            title_science="Mulch evaporation reduction and shade-net soil cooling (exp03)",
            level_kind="numeric",
            x_label_plain="Share of ground covered (0 = none, 1 = all)",
            x_label_science="Mulch or shade-net cover fraction (0–1)",
            baseline_label="Bare unshaded sand (one 12 mm rain on day 1)",
            headline=["h2_uptake_umol_m2", "h2_uptake_contrast_umol_m2"],
            series=[
                Series(
                    id="mulch-only",
                    label_plain="Mulch (no shade)",
                    label_science="Mulch cover fraction at zero shade",
                    levels=[
                        LevelSpec(
                            label=f"mulch {m}, shade 0",
                            x=float(m),
                            action={"mulch_cover_fraction": m, "shade_net_cover_fraction": 0},
                            upstream=["mulch_evaporation", "shade_net_temperature"],
                            twin={},
                            downstream=[],
                        )
                        for m in (0, 0.5, 1)
                    ],
                ),
                Series(
                    id="shade-only",
                    label_plain="Shade net (no mulch)",
                    label_science="Shade-net cover fraction at zero mulch",
                    levels=[
                        LevelSpec(
                            label=f"mulch 0, shade {s}",
                            x=float(s),
                            action={"mulch_cover_fraction": 0, "shade_net_cover_fraction": s},
                            upstream=["mulch_evaporation", "shade_net_temperature"],
                            twin={},
                            downstream=[],
                        )
                        for s in (0, 0.5, 1)
                    ],
                ),
                Series(
                    id="mulch-and-shade",
                    label_plain="Both mulch and shade net",
                    label_science="Combined equal mulch and shade-net cover fraction",
                    levels=[
                        LevelSpec(
                            label=f"mulch {c}, shade {c}",
                            x=float(c),
                            action={"mulch_cover_fraction": c, "shade_net_cover_fraction": c},
                            upstream=["mulch_evaporation", "shade_net_temperature"],
                            twin={},
                            downstream=[],
                        )
                        for c in (0, 0.5, 1)
                    ],
                ),
            ],
        ),
        # 8. solar-park (heat) — v0.7 exp08_solar_park
        Question(
            slug="solar-park",
            theme="heat",
            title_plain="What if we build a solar park?",
            title_science="Under-panel soil cooling, evaporation reduction, and panel-washing water pulses (exp08)",
            level_kind="numeric",
            x_label_plain="Solar panel cover (0 = none, 1 = whole area)",
            x_label_science="Panel cover fraction (0–1)",
            baseline_label="Open desert sand outside panels",
            headline=["h2_uptake_umol_m2", "h2_uptake_contrast_umol_m2"],
            series=[
                Series(
                    id="panel-cover",
                    label_plain="Under solar panels",
                    label_science="Panel cover fraction",
                    levels=[
                        LevelSpec(
                            label=f"panel cover {f}",
                            x=float(f),
                            action={"panel_cover_fraction": f},
                            upstream=["solar_park_microclimate", "panel_washing"],
                            twin={},
                            downstream=[],
                        )
                        for f in (0, 0.5, 1)
                    ],
                )
            ],
        ),
        # 9. compaction (soil) — v0.6 porosity 0.3 + v0.7 exp10_compaction
        Question(
            slug="compaction",
            theme="soil",
            title_plain="What if vehicles compact the sand?",
            title_science="Off-road vehicle compaction curve and direct porosity clamp (v0.6 + exp10)",
            level_kind="numeric",
            x_label_plain="Vehicle passes",
            x_label_science="Off-road vehicle passes (count)",
            baseline_label="Uncompacted sand (pore space 40 % of volume)",
            headline=["h2_uptake_umol_m2", "h2_uptake_contrast_umol_m2", "porosity"],
            series=[
                Series(
                    id="offroad-passes",
                    label_plain="Off-road vehicle passes",
                    label_science="Bulk-density compaction curve (exp10)",
                    levels=[
                        LevelSpec(
                            label="1 pass" if n == 1 else f"{n} passes",
                            x=float(n),
                            action={"offroad_passes": n},
                            upstream=["offroad_compaction"],
                            twin={},
                            downstream=[],
                            intervened_params=frozenset({"porosity"}),
                        )
                        for n in (0, 1, 10, 50)
                    ],
                ),
                Series(
                    id="direct-porosity-v06",
                    label_plain="Heavily packed sand (pore space 30 % of volume)",
                    label_science="Direct porosity reduction to 0.30 (v0.6 scenario)",
                    levels=[
                        LevelSpec(
                            label="porosity 0.30 (v0.6)",
                            x=50.0,
                            action={},
                            upstream=[],
                            twin={"params": {"porosity": 0.3}, "intervention": {"parameters": {"porosity": 0.3}}},
                            downstream=[],
                            intervened_params=frozenset({"porosity"}),
                        )
                    ],
                ),
            ],
        ),
        # 10. biochar (soil) — v0.7 exp04_biochar
        Question(
            slug="biochar",
            theme="soil",
            title_plain="What if we add biochar?",
            title_science="Biochar dose porosity increase and soil gas-uptake response (exp04)",
            level_kind="numeric",
            x_label_plain="Biochar added (% by weight)",
            x_label_science="Biochar dose (% w/w)",
            baseline_label="Unamended sand (0 % biochar, pore space 40 % of volume)",
            headline=["h2_uptake_umol_m2", "h2_uptake_contrast_umol_m2", "porosity"],
            series=[
                Series(
                    id="biochar-dose",
                    label_plain="Biochar dose",
                    label_science="Biochar amendment (% w/w)",
                    levels=[
                        LevelSpec(
                            label=f"{d} % w/w",
                            x=float(d),
                            action={"biochar_dose_pct_w_w": d},
                            upstream=["biochar_porosity"],
                            twin={},
                            downstream=[],
                            intervened_params=frozenset({"porosity"}),
                        )
                        for d in (0, 1, 2, 5)
                    ],
                )
            ],
        ),
        # 11. missing-co-genes (microbes) — v0.6 closed-gate vs open-gate samples
        Question(
            slug="missing-co-genes",
            theme="microbes",
            title_plain="What if the soil lacks carbon-monoxide-eating genes?",
            title_science="Deductive zero CO uptake under closed COCarrierMAG gate across field samples (v0.6)",
            level_kind="sample",
            x_label_plain="Soil sample",
            x_label_science="Metagenomic sample ID",
            baseline_label="One 12 mm rain on day 1 in the selected sample community",
            headline=["co_uptake_umol_m2", "h2_uptake_umol_m2"],
            series=[
                Series(
                    id="co-gate-samples",
                    label_plain="Desert soil samples",
                    label_science="Closed-gate vs open-gate metagenomic communities (gate coupling)",
                    levels=[
                        LevelSpec(
                            label=f"{s} ({kind})",
                            x=s,
                            action={},
                            upstream=[],
                            twin={"sample": s, "coupling": "gate"},
                            downstream=[],
                            sample=s,
                        )
                        for s, kind in (
                            ("51PRr3", "CO gate closed"),
                            ("10Dr1", "CO gate closed"),
                            ("13Dr3", "CO gate closed"),
                            ("14PRr2", "CO gate closed"),
                            ("51PRr1", "CO gate closed"),
                            ("10Dr2", "CO gate open"),
                            ("V53PRr2", "CO gate open"),
                        )
                    ],
                )
            ],
        ),
        # 12. adding-microbes (microbes) — v0.7 exp05_inoculation
        Question(
            slug="adding-microbes",
            theme="microbes",
            title_plain="What if we add microbes?",
            title_science="Inoculation of closed-gate and open-gate communities under gate vs proportional coupling (exp05)",
            level_kind="numeric",
            x_label_plain="Added microbes (% of community)",
            x_label_science="Inoculum share (%)",
            baseline_label="Uninoculated sample community",
            headline=["co_uptake_umol_m2", "co_uptake_contrast_umol_m2"],
            series=[
                Series(
                    id=sid,
                    label_plain=plain_lbl,
                    label_science=sci_lbl,
                    levels=[
                        LevelSpec(
                            label=f"{case_name}; {mode}; share {s}",
                            x=float(100.0 * s),
                            action={"inoculum_mag_id": mag, "inoculum_share": s},
                            upstream=["inoculant"],
                            twin={"sample": sample, "coupling": mode},
                            downstream=[],
                            sample=sample,
                        )
                        for s in (0.001, 0.01, 0.05)
                    ],
                )
                for sid, plain_lbl, sci_lbl, case_name, sample, mag, mode in (
                    (
                        "51PRr3-gate",
                        "Sample lacking carbon-monoxide-eating microbes + carrier microbe (yes/no presence)",
                        "51PRr3 (CO gate closed) + COCarrierMAG, gate coupling",
                        "51PRr3 (CO gate closed) + CO carrier",
                        "51PRr3",
                        "11Dr3__single_SemiBin_57",
                        "gate",
                    ),
                    (
                        "51PRr3-proportional",
                        "Sample lacking carbon-monoxide-eating microbes + carrier microbe (scaled by amount)",
                        "51PRr3 (CO gate closed) + COCarrierMAG, proportional coupling",
                        "51PRr3 (CO gate closed) + CO carrier",
                        "51PRr3",
                        "11Dr3__single_SemiBin_57",
                        "proportional",
                    ),
                    (
                        "V53PRr2-gate",
                        "Sample already having carbon-monoxide-eating microbes + carrier microbe (yes/no presence)",
                        "V53PRr2 (CO gate open) + COCarrierMAG, gate coupling",
                        "V53PRr2 (CO gate already open) + CO carrier",
                        "V53PRr2",
                        "11Dr3__single_SemiBin_57",
                        "gate",
                    ),
                    (
                        "V53PRr2-proportional",
                        "Sample already having carbon-monoxide-eating microbes + carrier microbe (scaled by amount)",
                        "V53PRr2 (CO gate open) + COCarrierMAG, proportional coupling",
                        "V53PRr2 (CO gate already open) + CO carrier",
                        "V53PRr2",
                        "11Dr3__single_SemiBin_57",
                        "proportional",
                    ),
                    (
                        "51PRr3-control-gate",
                        "Negative control: microbe without carbon-monoxide-eating genes (yes/no presence)",
                        "51PRr3 + non-carrier MAG (negative control), gate coupling",
                        "51PRr3 + non-carrier MAG (negative control)",
                        "51PRr3",
                        "10Dr2__single_SemiBin_0",
                        "gate",
                    ),
                    (
                        "51PRr3-control-proportional",
                        "Negative control: microbe without carbon-monoxide-eating genes (scaled by amount)",
                        "51PRr3 + non-carrier MAG (negative control), proportional coupling",
                        "51PRr3 + non-carrier MAG (negative control)",
                        "51PRr3",
                        "10Dr2__single_SemiBin_0",
                        "proportional",
                    ),
                )
            ],
        ),
        # 13. air-hydrogen (microbes) — v0.7 exp07_h2_leakage
        Question(
            slug="air-hydrogen",
            theme="microbes",
            title_plain="What if air hydrogen rises?",
            title_science="Soil H2 sink scaling under elevated atmospheric H2 mixing ratios (exp07)",
            level_kind="numeric",
            x_label_plain="Air hydrogen (multiple of background)",
            x_label_science="Atmospheric H2 multiplier (× background)",
            baseline_label="Normal background air hydrogen (1x normal level)",
            headline=["h2_uptake_umol_m2", "h2_uptake_contrast_umol_m2"],
            series=[
                Series(
                    id="h2-multiplier",
                    label_plain="Air hydrogen multiple",
                    label_science="Multiple of background atmospheric H2 concentration",
                    levels=[
                        LevelSpec(
                            label=f"{m}x background",
                            x=float(m),
                            action={"h2_multiplier": m},
                            upstream=["h2_atmosphere"],
                            twin={},
                            downstream=[],
                        )
                        for m in (1, 1.5, 2, 10)
                    ],
                )
            ],
        ),
        # 14. carbon-fixation (microbes) — v0.7 exp06_carbon_fixation_bound
        Question(
            slug="carbon-fixation",
            theme="microbes",
            title_plain="How much carbon can the hydrogen eaters fix?",
            title_science="Stoichiometric upper bound on knallgas autotrophic carbon fixation (exp06)",
            level_kind="sample",
            x_label_plain="Soil sample",
            x_label_science="Metagenomic sample ID",
            baseline_label="Sample community under 14-day rain schedule",
            headline=["c_fixed_upper_bound_g_m2_yr", "c_fixed_upper_bound_ug_m2", "co2_fixed_upper_bound_mol_m2"],
            series=[
                Series(
                    id="baseline-12mm",
                    label_plain="Normal storm (12 mm)",
                    label_science="Baseline 12 mm single-pulse regime",
                    levels=[
                        LevelSpec(
                            label=f"baseline (12 mm pulse); {sample}",
                            x=sample,
                            action={},
                            upstream=[],
                            twin={"sample": sample, "coupling": "gate"},
                            downstream=["knallgas_carbon_fixation", "co_oxidation_carbon_balance"],
                            sample=sample,
                        )
                        for sample in ("14PRr2", "F14Dr1")
                    ],
                ),
                Series(
                    id="wet-36mm",
                    label_plain="Heavy storm (36 mm)",
                    label_science="Wet 36 mm single-pulse regime",
                    levels=[
                        LevelSpec(
                            label=f"wet (36 mm pulse); {sample}",
                            x=sample,
                            action={},
                            upstream=[],
                            twin={
                                "sample": sample,
                                "coupling": "gate",
                                "forcing": [{**f, "water_mm": 36.0 if d == 0 else 0.0} for d, f in enumerate(DEFAULT_FORCING)],
                            },
                            downstream=["knallgas_carbon_fixation", "co_oxidation_carbon_balance"],
                            sample=sample,
                        )
                        for sample in ("14PRr2", "F14Dr1")
                    ],
                ),
            ],
        ),
    ]
    return {q.slug: q for q in q_list}


QUESTIONS: dict[str, Question] = _build_questions()


def evaluate(
    spec: LevelSpec,
    constants: dict[str, float] | None = None,
    params: dict[str, float] | None = None,
    *,
    keep_twin_result: bool = False,
) -> dict:
    """Evaluate one LevelSpec at the given constants and twin parameters (with envelope=False).

    Returns a flat dictionary mapping output port names (plus derived contrasts and annualised bounds)
    to float | None, along with `"__refused__": str | None` and `"__deduced__": list[str]`.
    """
    base_kw: dict[str, float] = {"k_h2_s": 0.001, "k_co_s": 0.001}
    if params:
        for k, v in params.items():
            if k not in spec.intervened_params:
                base_kw[k] = float(v)

    downstream_ports = [o.name for d in spec.downstream for o in MODELS[d].outputs]
    extra_ports = (
        ["c_fixed_upper_bound_ug_m2_yr", "c_fixed_upper_bound_g_m2_yr"]
        if "knallgas_carbon_fixation" in spec.downstream
        else []
    )
    all_output_keys = [
        "h2_uptake_umol_m2",
        "co_uptake_umol_m2",
        "h2_uptake_contrast_umol_m2",
        "co_uptake_contrast_umol_m2",
        "active_days",
        "final_water_mm",
        "final_salt_mol_m2",
        "mean_nacl_mol_l",
        "porosity",
        *downstream_ports,
        *extra_ports,
    ]

    try:
        p_obj = Parameters(**base_kw)
    except ValueError as exc:
        out = {k: None for k in all_output_keys}
        out["__refused__"] = str(exc)
        out["__deduced__"] = []
        return out

    sample = spec.sample or spec.twin.get("sample")
    comm = spec.twin.get("community") or ({"sample": sample} if sample else None)
    coupling = spec.twin.get("coupling")
    iv = deepcopy(spec.twin.get("intervention", {}))
    if "params" in spec.twin:
        iv.setdefault("parameters", {}).update(spec.twin["params"])
    forcing = deepcopy(spec.twin.get("forcing", DEFAULT_FORCING))

    t_req = twin_request(
        forcing,
        initial=State(1.8, 0.0),
        params=p_obj,
        community=comm,
        coupling=coupling,
        intervention=iv if iv else None,
        data=DATA,
    )
    r = run_chain(
        spec.action,
        [MODELS[u] for u in spec.upstream],
        t_req,
        [MODELS[d] for d in spec.downstream],
        scenario="central",
        constants=constants,
        envelope=False,
        keep_twin_result=keep_twin_result,
        registry=_REG,
    )

    out: dict = {k: None for k in all_output_keys}
    deduced: list[str] = []

    if r["twin"] is not None:
        t_out = r["twin"]["outputs"]
        t_base = r["twin"]["baseline_outputs"] or {}
        h2_i = t_out["h2_uptake_mol_m2"]["value"]
        co_i = t_out["co_uptake_mol_m2"]["value"]
        h2_b = t_base.get("h2_uptake_mol_m2")
        co_b = t_base.get("co_uptake_mol_m2")

        out["h2_uptake_umol_m2"] = None if h2_i is None else 1e6 * float(h2_i)
        out["co_uptake_umol_m2"] = None if co_i is None else 1e6 * float(co_i)
        out["h2_uptake_contrast_umol_m2"] = (
            None if (h2_i is None or h2_b is None) else 1e6 * (float(h2_i) - float(h2_b))
        )
        out["co_uptake_contrast_umol_m2"] = (
            None if (co_i is None or co_b is None) else 1e6 * (float(co_i) - float(co_b))
        )
        out["active_days"] = (
            float(t_out["active_days"]["value"]) if t_out["active_days"]["value"] is not None else None
        )
        out["final_water_mm"] = t_out["final_water_mm"]["value"]
        out["final_salt_mol_m2"] = t_out["final_salt_mol_m2"]["value"]
        out["mean_nacl_mol_l"] = t_out["mean_nacl_mol_l"]["value"]
        iv_por = r["twin"]["intervention"].get("parameters", {}).get("porosity")
        out["porosity"] = float(iv_por) if iv_por is not None else float(p_obj.porosity)

        for td in r["twin"].get("twin_deduced", []):
            gas = td["gas"]
            for gk in (f"{gas}_uptake_umol_m2", f"{gas}_uptake_contrast_umol_m2"):
                if gk not in deduced:
                    deduced.append(gk)
                if out[gk] is not None:
                    out[gk] = 0.0

    for d_id in spec.downstream:
        d_map = r["downstream"].get(d_id, {})
        for o in MODELS[d_id].outputs:
            w = d_map.get(o.name)
            if w is not None:
                out[o.name] = w["value"]

    for dn in r.get("deduced_nulls", []):
        port = dn["port"]
        if port not in deduced:
            deduced.append(port)
        out[port] = 0.0

    if "knallgas_carbon_fixation" in spec.downstream:
        days_per_year = float((constants or {}).get("days_per_year", _REG["days_per_year"]["value"]))
        c_ug = out.get("c_fixed_upper_bound_ug_m2")
        if "c_fixed_upper_bound_ug_m2" in deduced:
            out["c_fixed_upper_bound_ug_m2_yr"] = 0.0
            out["c_fixed_upper_bound_g_m2_yr"] = 0.0
            for k in ("c_fixed_upper_bound_ug_m2_yr", "c_fixed_upper_bound_g_m2_yr"):
                if k not in deduced:
                    deduced.append(k)
        elif c_ug is not None:
            c_ug_yr = float(c_ug) * days_per_year / 14.0
            out["c_fixed_upper_bound_ug_m2_yr"] = c_ug_yr
            out["c_fixed_upper_bound_g_m2_yr"] = c_ug_yr / 1e6

    if r["refused"] is not None:
        out["__refused__"] = r["refused"]["reason"]
    elif r["twin"] is not None and r["twin"]["outputs"]["h2_uptake_mol_m2"]["value"] is None:
        out["__refused__"] = "twin abstains outside its ideal-dilute NaCl domain (> 0.1 mol/L)"
    else:
        out["__refused__"] = None

    out["__deduced__"] = deduced
    if keep_twin_result:
        out["__twin__"] = r["twin"]
        out["__chain__"] = r
    return out
