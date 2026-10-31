"""External models composed with the twin (v0.7). Pure functions over typed ports.

Every literature constant is read by id from data/external/composition_constants.json (passed in as the
``constants`` argument of ``fn``); none appears in this file. The validity intervals below are part of each model's
definition (the domain in which the model is applied, a design bound of this release, not a fitted constant); an
input outside the interval is refused by composition.check_validity.

Models are collected in ``MODELS`` (id -> ExternalModel). Naming of ports: a port name is the wire; producer and
consumer must agree on (term, unit). Series ports (baseline_*, water_mm, temp_c, evaporation_mm, nacl_mol_l) are
per-day lists of the baseline length; scalar ports are floats. Twin-input ports are defined in composition.py
(TWIN_INPUT_PORTS) and reused here.
"""
from __future__ import annotations

import json
from math import exp
from pathlib import Path

from .composition import (BASELINE_PORTS, TWIN_INPUT_PORTS, TWIN_OUTPUT_PORTS, CompositionError, ExternalModel, Port,
                          load_constants)

# ---- ports (composition-layer ports use the rak: namespace; ChEBI terms are those already used by the twin)
EC_W = Port("ec_w_dS_m", "rak:IrrigationWaterElectricalConductivity", "dS/m")
IRR_NACL = Port("irrigation_nacl_mol_l", "CHEBI:26710", "mol/L")
IRR_RATE = Port("irrigation_mm_per_day", "rak:IrrigationDepthRate", "mm/day")
IRR_INTERVAL = Port("irrigation_interval_days", "rak:IrrigationInterval", "day")
ECE = Port("ece_dS_m", "rak:SaturationExtractElectricalConductivity", "dS/m")
FITTED_F = Port("fitted_leaching_factor", "rak:LeachingFactor", "-")
MULCH = Port("mulch_cover_fraction", "rak:MulchCover", "-")
SHADE = Port("shade_net_cover_fraction", "rak:ShadeNetCover", "-")
BIOCHAR = Port("biochar_dose_pct_w_w", "rak:BiocharDose", "% w/w")
MAG = Port("inoculum_mag_id", "rak:InoculumMAG", "-")
SHARE = Port("inoculum_share", "rak:InoculumShare", "-")
H2_MULT = Port("h2_multiplier", "rak:AtmosphericH2Multiplier", "-")
H2_PPM = Port("h2_mixing_ratio_ppm", "CHEBI:18276", "ppm")
PANEL = Port("panel_cover_fraction", "rak:SolarPanelCover", "-")
SEED = Port("seeding_intensity", "rak:CloudSeedingIntensity", "-")
PASSES = Port("offroad_passes", "rak:OffRoadVehiclePasses", "-")
WARM = Port("warming_c", "rak:AirTemperatureWarming", "degC")
RAIN_CHG = Port("rain_change_frac", "rak:RainfallRelativeChange", "-")
CO2_UB = Port("co2_fixed_upper_bound_mol_m2", "CHEBI:16526", "mol/m2")
C_UB = Port("c_fixed_upper_bound_ug_m2", "rak:FixedCarbonMass", "ug/m2")
CO2_FROM_CO = Port("co2_from_co_mol_m2", "CHEBI:16526", "mol/m2")
NET_C = Port("net_atmospheric_c_change_mol_m2", "rak:NetAtmosphericCarbonChange", "mol/m2")
RY = {c: Port(f"relative_yield_{c}", "rak:CropRelativeYield", "-") for c in ("date_palm", "barley", "alfalfa", "tomato")}
RY_SA187 = Port("relative_yield_alfalfa_sa187", "rak:CropRelativeYield", "-")
H2_UPTAKE, CO_UPTAKE = TWIN_OUTPUT_PORTS["h2_uptake_mol_m2"], TWIN_OUTPUT_PORTS["co_uptake_mol_m2"]
MEAN_NACL = TWIN_OUTPUT_PORTS["mean_nacl_mol_l"]
tp = TWIN_INPUT_PORTS
bp = BASELINE_PORTS


# ---- helpers
def maas_hoffman(ece: float, threshold: float, slope_pct: float) -> float:
    """Maas-Hoffman relative yield: 1 - slope * max(0, ECe - threshold), clipped to [0, 1] (slope in % per dS/m)."""
    return max(0.0, min(1.0, 1.0 - slope_pct / 100.0 * max(0.0, ece - threshold)))


def _pulse_days(n: int, interval: int, offset: int = 0) -> list:
    interval = max(1, int(round(interval)))
    return [d for d in range(n) if (d + offset) % interval == 0]


# ---- model functions (inputs by port name, constants by id)
def _ec_to_nacl(i, c):
    return {"irrigation_nacl_mol_l": i["ec_w_dS_m"] * c["nacl_mol_l_per_dS_m"]}


def _irrigation_schedule(i, c):
    """Add a fixed seasonal total (rate x days) to the baseline water series as equal pulses every `interval` days; the
    applied water has NaCl concentration irrigation_nacl_mol_l (mixed with baseline water by mass of salt)."""
    water, conc = i["baseline_water_mm"], i["baseline_nacl_mol_l"]
    n = len(water)
    days = _pulse_days(n, i["irrigation_interval_days"])
    amount = i["irrigation_mm_per_day"] * n / len(days)
    w, s = list(water), [water[d] * conc[d] for d in range(n)]
    for d in days:
        w[d] += amount
        s[d] += amount * i["irrigation_nacl_mol_l"]
    return {"water_mm": w, "nacl_mol_l": [s[d] / w[d] if w[d] else 0.0 for d in range(n)]}


def _water_demand(const_id):
    return lambda i, c: {"irrigation_mm_per_day": c[const_id]}


def _fao29_leaching(i, c):
    return {"ece_dS_m": c["fao29_ece_per_ecw"] * i["ec_w_dS_m"]}


def _mh(crop):
    return lambda i, c: {f"relative_yield_{crop}": maas_hoffman(i["ece_dS_m"], c[f"mh_{crop}_threshold_dS_m"],
                                                                 c[f"mh_{crop}_slope_pct_per_dS_m"])}


def _twin_solution_to_ece(i, c):
    """EC_sw = NaCl / (mol/L per dS/m); ECe = EC_sw / (EC_sw per ECe)."""
    return {"ece_dS_m": i["mean_nacl_mol_l"] / c["nacl_mol_l_per_dS_m"] / c["fao29_ecsw_per_ece"]}


def fit_leaching_factor(c: dict, trait: str, crop: str = "alfalfa") -> float:
    """Factor f in ECe = f * EC_w such that the Maas-Hoffman ratio RY(EC_high)/RY(EC_low) equals the mock ratio of the
    digitised season-1 yields of `trait` ('fresh' or 'dry'). Bisection on the strictly decreasing branch; deterministic."""
    thr, slope = c[f"mh_{crop}_threshold_dS_m"], c[f"mh_{crop}_slope_pct_per_dS_m"]
    lo_ec, hi_ec = c["sa187_ecw_low_dS_m"], c["sa187_ecw_high_dS_m"]
    target = c[f"sa187_mock_ratio_{trait}"]
    if not 0 < target < 1:
        raise CompositionError("fit needs a mock yield ratio strictly between 0 and 1")
    a, b = thr / hi_ec, (100.0 / slope + thr) / hi_ec        # ratio = 1 at a, 0 at b
    ratio = lambda f: maas_hoffman(f * hi_ec, thr, slope) / maas_hoffman(f * lo_ec, thr, slope)
    for _ in range(200):
        m = (a + b) / 2
        if ratio(m) > target:
            a = m
        else:
            b = m
    return (a + b) / 2


def _fitted(trait):
    def fn(i, c):
        f = fit_leaching_factor(c, trait)
        return {"ece_dS_m": f * i["ec_w_dS_m"], "fitted_leaching_factor": f}
    return fn


def _fitted_constants(trait):
    return ("mh_alfalfa_threshold_dS_m", "mh_alfalfa_slope_pct_per_dS_m", "sa187_ecw_low_dS_m", "sa187_ecw_high_dS_m",
            f"sa187_mock_ratio_{trait}")


DIGITISED_FIG1 = Path(__file__).resolve().parents[3] / "data" / "external" / "de_zelicourt_2018" / "fig1_digitised.json"


def digitised_trial(path: Path | None = None) -> dict:
    """Pixel digitisation of de Zelicourt et al. 2018 Fig 1 (reproducibility/v07/digitise_fig1.py), read from disk."""
    return json.loads(Path(path or DIGITISED_FIG1).read_text())


def digitisation_error(doc: dict | None = None) -> dict:
    """Error of a digitised bar height: absolute = largest y-axis calibration residual (t/ha); relative = largest
    disagreement between the digitised SA187 effects and the percentages printed on the figure."""
    d = doc or digitised_trial()
    a = max(v for k, v in d["values"].items() if k.endswith("calibration_max_residual_t_ha"))
    return {"absolute_t_ha": a, "relative": d["max_abs_diff"]}


def mock_yield(trait: str, level: str, doc: dict | None = None) -> float:
    """Season-1 mock yield (t/ha) at 'low' (3.12) or 'high' (7.81 dS/m) irrigation salinity, from the digitised figure."""
    d = doc or digitised_trial()
    return d["values"][f"{'low_3.12' if level == 'low' else 'high_7.81'}/{trait}_s1/mock"]


def observed_mock_ratio(trait: str, doc: dict | None = None) -> dict:
    """Observed season-1 mock yield ratio (high/low salinity) from the pixel digitisation; low/high propagate the
    digitisation error (absolute calibration residual plus relative disagreement with the printed percentages) as the
    extreme combinations. Mirrored in the registry as sa187_mock_ratio_<trait> (checked by a test)."""
    d = doc or digitised_trial()
    lo, hi, e = mock_yield(trait, "low", d), mock_yield(trait, "high", d), digitisation_error(d)
    a, r = e["absolute_t_ha"], e["relative"]
    return {"value": hi / lo, "low": (hi * (1 - r) - a) / (lo * (1 + r) + a), "high": (hi * (1 + r) + a) / (lo * (1 - r) - a),
            "unit": "-", "source": "external:sa187_field_trial", "error": e, "yields_t_ha": {"low": lo, "high": hi}}


def _sa187(i, c):
    return {"relative_yield_alfalfa_sa187": i["relative_yield_alfalfa"] * (1.0 + c["sa187_benefit_b"])}


def _mulch(i, c):
    k = 1.0 - c["mulch_evap_reduction_frac"] * i["mulch_cover_fraction"]
    return {"evaporation_mm": [e * k for e in i["baseline_evaporation_mm"]]}


def _shade(i, c):
    d = c["shade_net_soil_temp_reduction_c"] * i["shade_net_cover_fraction"]
    return {"temp_c": [max(0.0, t - d) for t in i["baseline_temp_c"]]}


def _biochar(i, c):
    return {"porosity": i["baseline_porosity"] + c["biochar_porosity_increase_per_pct"] * i["biochar_dose_pct_w_w"]}


def _inoculant(i, c):
    return {"inoculation": {"mag_id": i["inoculum_mag_id"], "share": i["inoculum_share"]}}


def _knallgas(i, c):
    co2 = i["h2_uptake_mol_m2"] * c["knallgas_co2_per_h2"]
    return {"co2_fixed_upper_bound_mol_m2": co2, "c_fixed_upper_bound_ug_m2": co2 * c["molar_mass_c_g_mol"] * 1e6}


def _co_balance(i, c):
    """CO + 1/2 O2 -> CO2: one C atom per CO, so atmospheric carbon atoms are unchanged (CO becomes CO2)."""
    return {"co2_from_co_mol_m2": i["co_uptake_mol_m2"], "net_atmospheric_c_change_mol_m2": 0.0}


def _h2_atmosphere(i, c):
    return {"h2_mol_m3": i["h2_multiplier"] * i["baseline_h2_mol_m3"], "h2_mixing_ratio_ppm": i["h2_multiplier"] * c["h2_background_ppm"]}


def _solar_microclimate(i, c):
    f = i["panel_cover_fraction"]
    return {"temp_c": [max(0.0, t - c["solar_panel_soil_temp_reduction_c"] * f) for t in i["baseline_temp_c"]],
            "evaporation_mm": [e * (1.0 - c["solar_panel_evap_reduction_frac"] * f) for e in i["baseline_evaporation_mm"]]}


def _panel_washing(i, c):
    """Wash water (L per m2 of panel = mm over the panel footprint) reaches the soil beneath, on every k-th day."""
    w = list(i["baseline_water_mm"])
    for d in _pulse_days(len(w), c["panel_wash_interval_days"], offset=1):
        w[d] += c["panel_wash_l_m2_per_wash"] * i["panel_cover_fraction"]
    return {"water_mm": w}


def _cloud_seeding(i, c):
    k = 1.0 + i["seeding_intensity"] * c["cloud_seeding_rain_increase_frac"]
    return {"water_mm": [w * k for w in i["baseline_water_mm"]]}


def _offroad(i, c):
    rho_s = c["particle_density_g_cm3"]
    rho0 = (1.0 - i["baseline_porosity"]) * rho_s
    d = c["offroad_bd_max_increase_g_cm3"] * (1.0 - exp(-i["offroad_passes"] / c["offroad_bd_passes_efold"]))
    return {"porosity": 1.0 - (rho0 + d) / rho_s}


def _warming(i, c):
    w = i["warming_c"]
    return {"temp_c": [t + c["soil_temp_per_air_temp_warming"] * w for t in i["baseline_temp_c"]],
            "evaporation_mm": [e * (1.0 + c["evap_scaling_per_c"] * w) for e in i["baseline_evaporation_mm"]],
            "water_mm": [x * (1.0 + i["rain_change_frac"]) for x in i["baseline_water_mm"]]}


def annualise(value: float, period_days: float, registry: dict | None = None) -> float:
    """Scale a total over `period_days` to one year (linear; a labelled extrapolation, not a forecast)."""
    r = registry or load_constants()
    return value * r["days_per_year"]["value"] / period_days


def _m(id, name, kind, inputs, outputs, validity, constants, fn, gate=None):
    return ExternalModel(id, name, kind, tuple(inputs), tuple(outputs), validity, tuple(constants), fn, gate)


_CROPS = {"date_palm": "date palm", "barley": "barley", "alfalfa": "alfalfa", "tomato": "tomato"}
MODELS: dict = {}
for _m_ in (
    _m("ec_to_nacl", "Irrigation-water EC to NaCl concentration", "upstream", [EC_W], [IRR_NACL], {"ec_w_dS_m": (0.0, 16.0)},
       ["nacl_mol_l_per_dS_m"], _ec_to_nacl),
    _m("irrigation_schedule", "Equal-total irrigation delivered in pulses every k days", "upstream",
       [bp["baseline_water_mm"], bp["baseline_nacl_mol_l"], IRR_RATE, IRR_INTERVAL, IRR_NACL], [tp["water_mm"], tp["nacl_mol_l"]],
       {"irrigation_mm_per_day": (0.0, 50.0), "irrigation_interval_days": (1.0, 60.0), "irrigation_nacl_mol_l": (0.0, 0.2)},
       [], _irrigation_schedule),
    _m("crop_water_demand_generic", "Crop water demand, hot-desert order of magnitude (FAO-56)", "upstream", [], [IRR_RATE], {},
       ["crop_etc_mm_per_day"], _water_demand("crop_etc_mm_per_day")),
    _m("crop_water_demand_alfalfa", "Drip-irrigated alfalfa water demand (FAO-56)", "upstream", [], [IRR_RATE], {},
       ["alfalfa_etc_mm_per_day"], _water_demand("alfalfa_etc_mm_per_day")),
    _m("fao29_leaching", "Root-zone salinity from irrigation-water EC at fixed leaching (FAO-29)", "downstream", [EC_W], [ECE],
       {"ec_w_dS_m": (0.0, 16.0)}, ["fao29_ece_per_ecw"], _fao29_leaching),
    _m("twin_solution_to_ece", "Twin mean soil-solution NaCl to saturation-extract EC (FAO-29 relations)", "downstream", [MEAN_NACL], [ECE],
       {"mean_nacl_mol_l": (0.0, 0.2)}, ["nacl_mol_l_per_dS_m", "fao29_ecsw_per_ece"], _twin_solution_to_ece),
    _m("fitted_leaching_fresh", "Leaching factor fitted to the SA187 field mock ratio (fresh weight)", "downstream", [EC_W], [ECE, FITTED_F],
       {"ec_w_dS_m": (0.0, 16.0)}, _fitted_constants("fresh"), _fitted("fresh")),
    _m("fitted_leaching_dry", "Leaching factor fitted to the SA187 field mock ratio (dry weight)", "downstream", [EC_W], [ECE, FITTED_F],
       {"ec_w_dS_m": (0.0, 16.0)}, _fitted_constants("dry"), _fitted("dry")),
    _m("sa187_benefit", "Enterobacter sp. SA187 yield benefit on alfalfa (field trial)", "downstream", [RY["alfalfa"]], [RY_SA187],
       {"relative_yield_alfalfa": (0.0, 1.0)}, ["sa187_benefit_b"], _sa187),
    _m("mulch_evaporation", "Mulch reduces soil evaporation", "upstream", [MULCH, bp["baseline_evaporation_mm"]], [tp["evaporation_mm"]],
       {"mulch_cover_fraction": (0.0, 1.0)}, ["mulch_evap_reduction_frac"], _mulch),
    _m("shade_net_temperature", "Shade net lowers soil temperature", "upstream", [SHADE, bp["baseline_temp_c"]], [tp["temp_c"]],
       {"shade_net_cover_fraction": (0.0, 1.0)}, ["shade_net_soil_temp_reduction_c"], _shade),
    _m("biochar_porosity", "Biochar dose raises porosity (meta-analysis)", "upstream", [BIOCHAR, bp["baseline_porosity"]], [tp["porosity"]],
       {"biochar_dose_pct_w_w": (0.0, 10.0)}, ["biochar_porosity_increase_per_pct"], _biochar),
    _m("inoculant", "Inoculum: reference MAG at a relative share", "upstream", [MAG, SHARE], [tp["inoculation"]],
       {"inoculum_share": (0.0, 0.2)}, [], _inoculant),
    _m("knallgas_carbon_fixation", "Upper bound of CO2 fixed by H2-oxidising autotrophs", "downstream", [H2_UPTAKE], [CO2_UB, C_UB],
       {"h2_uptake_mol_m2": (0.0, 1e3)}, ["knallgas_co2_per_h2", "molar_mass_c_g_mol"], _knallgas, gate="H2AutotrophicCandidateMAG"),
    _m("co_oxidation_carbon_balance", "CO oxidation converts CO to CO2 (no net carbon-atom change)", "downstream", [CO_UPTAKE],
       [CO2_FROM_CO, NET_C], {"co_uptake_mol_m2": (0.0, 1e3)}, [], _co_balance),
    _m("h2_atmosphere", "Atmospheric H2 scenario (multiple of background)", "upstream", [H2_MULT, bp["baseline_h2_mol_m3"]],
       [tp["h2_mol_m3"], H2_PPM], {"h2_multiplier": (0.5, 20.0)}, ["h2_background_ppm"], _h2_atmosphere),
    _m("solar_park_microclimate", "Solar park: under-panel soil temperature and evaporation", "upstream",
       [PANEL, bp["baseline_temp_c"], bp["baseline_evaporation_mm"]], [tp["temp_c"], tp["evaporation_mm"]],
       {"panel_cover_fraction": (0.0, 1.0)}, ["solar_panel_soil_temp_reduction_c", "solar_panel_evap_reduction_frac"], _solar_microclimate),
    _m("panel_washing", "Solar park: panel-washing water pulses", "upstream", [PANEL, bp["baseline_water_mm"]], [tp["water_mm"]],
       {"panel_cover_fraction": (0.0, 1.0)}, ["panel_wash_l_m2_per_wash", "panel_wash_interval_days"], _panel_washing),
    _m("cloud_seeding", "Cloud seeding: relative rain increase on rain days", "upstream", [SEED, bp["baseline_water_mm"]], [tp["water_mm"]],
       {"seeding_intensity": (0.0, 3.0)}, ["cloud_seeding_rain_increase_frac"], _cloud_seeding),
    _m("offroad_compaction", "Off-road traffic compacts desert sand", "upstream", [PASSES, bp["baseline_porosity"]], [tp["porosity"]],
       {"offroad_passes": (0.0, 200.0)}, ["particle_density_g_cm3", "offroad_bd_max_increase_g_cm3", "offroad_bd_passes_efold"], _offroad),
    _m("climate_warming", "Warming: temperature, evaporation scaling and rainfall change", "upstream",
       [WARM, RAIN_CHG, bp["baseline_temp_c"], bp["baseline_evaporation_mm"], bp["baseline_water_mm"]],
       [tp["temp_c"], tp["evaporation_mm"], tp["water_mm"]], {"warming_c": (0.0, 6.0), "rain_change_frac": (-0.5, 0.5)},
       ["soil_temp_per_air_temp_warming", "evap_scaling_per_c"], _warming),
):
    MODELS[_m_.id] = _m_
for _c, _name in _CROPS.items():
    MODELS[f"maas_hoffman_{_c}"] = _m(f"maas_hoffman_{_c}", f"Maas-Hoffman relative yield, {_name} (FAO-29)", "downstream", [ECE], [RY[_c]],
                                      {"ece_dS_m": (0.0, 40.0)}, [f"mh_{_c}_threshold_dS_m", f"mh_{_c}_slope_pct_per_dS_m"], _mh(_c))


def get(model_id: str) -> ExternalModel:
    if model_id not in MODELS:
        raise CompositionError(f"unknown external model {model_id!r}; available: {sorted(MODELS)}")
    return MODELS[model_id]
