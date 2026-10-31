"""v0.7 composition experiments: the twin chained with other published models ("chained what-ifs").

Twelve experiments, one function each (the docstring is the question in plain words). Results go to
reviews/v07/composition/<id>.json and summary.json. Deterministic: sorted keys, no timestamps, floats rounded to six
significant digits. Every number carries its provenance: 'twin' (the process model), 'external:<model>' (another published
model) or 'deduced' (a zero deduced from the ontology). Baseline: the v0.6 primary scenario (reproducibility/v06/analyse.py,
reference_schedules()[1] 'single_pulse', Parameters(k_h2_s=.001, k_co_s=.001), initial water 1.8 mm).
"""
from __future__ import annotations

import importlib.util
import inspect
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
from eq_biodt.v03 import community as comm                                   # noqa: E402
from eq_biodt.v03 import external_models as xm                               # noqa: E402
from eq_biodt.v03.composition import (constant_values, constants_sha256, load_constants, run_chain,   # noqa: E402
                                      twin_request)
from eq_biodt.v03.process import Parameters, State                           # noqa: E402

OUT = ROOT / "reviews/v07/composition"
M = xm.MODELS
REG = load_constants()
UM = 1e6          # mol -> umol

_spec = importlib.util.spec_from_file_location("v06_analyse", ROOT / "reproducibility/v06/analyse.py")
_an = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_an)
FORCING = _an.reference_schedules()[1]["forcing"]       # single_pulse: 14 d, 12 mm on day 0, 1 mm/d evaporation, 25 C
DAYS = len(FORCING)
BASELINE = {"source": "reproducibility/v06/analyse.py reference_schedules()[1] ('single_pulse')",
            "forcing": "14 days; 12 mm rain on day 0, none after; evaporation 1 mm/day; soil temperature 25 C; no NaCl in rain",
            "initial_water_mm": 1.8, "parameters": "Parameters(k_h2_s=0.001, k_co_s=0.001); other parameters at their defaults "
            "(Rub' al Khali sand: porosity 0.4)", "community": "none (coupling 'reference': gamma = 1) unless the experiment names one",
            "period_days": DAYS}
TWIN_KEYS = ("h2_uptake_mol_m2", "co_uptake_mol_m2", "final_water_mm", "final_salt_mol_m2", "mean_nacl_mol_l", "active_days")
TRACE_KEYS = ("ontology_sha256", "evidence_sha256", "process_code_sha256", "contract_sha256", "graph_sha256", "targets_replaced")
SRC = "Rub' al Khali soil and community; not a prediction for other soils"


# ---------------------------------------------------------------- generic helpers
def rnd(o):
    """Round floats to six significant digits; reject non-finite numbers; tuples become lists."""
    if isinstance(o, bool) or o is None or isinstance(o, (str, int)):
        return o
    if isinstance(o, float):
        if o != o or o in (float("inf"), float("-inf")):
            raise ValueError("non-finite number in result")
        return float(f"{o:.6g}")
    if isinstance(o, dict):
        return {str(k): rnd(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [rnd(v) for v in o]
    raise TypeError(type(o))


def dump(path: Path, obj) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(rnd(obj), indent=2, sort_keys=True, allow_nan=False) + "\n")


def req(forcing=None, **kw):
    return twin_request(forcing or FORCING, **kw)


def compact(w: dict) -> dict:
    """Series-valued wrapped numbers (per-day forcing) are summarised by their total over the period."""
    if not isinstance(w["value"], list):
        return w
    tot = lambda s: None if s is None else sum(s)
    return {"series_total": {"value": tot(w["value"]), "low": tot(w["low"]), "high": tot(w["high"])}, "days": len(w["value"]),
            "unit": w["unit"], "provenance": w["provenance"], "constants": w["constants"]}


def diff_vs_baseline(r: dict) -> dict:
    base = r["twin"]["baseline_outputs"]
    out = {}
    for k in TWIN_KEYS:
        w, b = r["twin"]["outputs"][k], base[k]
        if w["value"] is None or b is None:
            out[k] = {"value": None, "low": None, "high": None, "unit": w["unit"]}
            continue
        d = {"value": w["value"] - b, "unit": w["unit"], "baseline": b}
        d["low"] = None if w["low"] is None else w["low"] - b
        d["high"] = None if w["high"] is None else w["high"] - b
        if b:
            d["percent"] = {"value": 100 * d["value"] / b, "low": None if d["low"] is None else 100 * d["low"] / b,
                            "high": None if d["high"] is None else 100 * d["high"] / b}
        out[k] = d
    return out


def level(label, action, ups, twin_req, downs, **kw) -> tuple:
    r = run_chain(action, [M[u] for u in ups], twin_req, [M[d] for d in downs], **kw)
    rec = {"label": label, "action": action, "refused": r["refused"], "deduced_nulls": r["deduced_nulls"],
           "envelope_refusals": r["envelope_refusals"]}
    if r["twin"] is not None:
        rec["upstream"] = {m: {p: compact(w) for p, w in d.items()} for m, d in r["upstream"].items()}
        rec["twin"] = r["twin"]["outputs"]
        rec["twin_baseline"] = r["twin"]["baseline_outputs"]
        rec["difference_vs_baseline"] = diff_vs_baseline(r)
        rec["downstream"] = {m: {p: compact(w) for p, w in d.items()} for m, d in r["downstream"].items()}
        rec["twin_checks_passed"] = r["twin"]["checks_passed"]
        rec["trace"] = {k: r["twin"]["trace"][k] for k in TRACE_KEYS}
        if r["twin"]["community"]["community"]:
            rec["community"] = {k: r["twin"]["community"][k] for k in ("mode", "community", "gamma", "samples")}
        if r["twin"]["inoculation"]:
            rec["inoculation"] = {k: r["twin"]["inoculation"][k] for k in ("spec", "mag_classes", "baseline_gamma", "gamma", "gate_opened")}
    return rec, r


def pct(rec, key="h2_uptake_mol_m2"):
    p = rec["difference_vs_baseline"][key].get("percent")
    return p


def head(label, w, unit, provenance, detail=""):
    return {"label": label, "value": w["value"], "low": w["low"], "high": w["high"], "unit": unit, "provenance": provenance,
            "detail": detail}


def h2_head(label, rec, detail=""):
    p = pct(rec)
    return head(label, p, "% of baseline H2 uptake", "twin", detail) if p else None


def finish(fn, name, chain, levels, results, headline, plot, plain, caveats, **extra):
    doc = {"id": name, "question": inspect.getdoc(fn).split("\n\n")[0].replace("\n", " "), "chain": chain, "baseline": BASELINE,
           "levels": levels, "refusals": [{"level": l["label"], **l["refused"]} for l in levels if l["refused"]],
           "deduced_nulls": [{"level": l["label"], **n} for l in levels for n in l["deduced_nulls"]],
           "plain_language": plain, "caveats": caveats, "headline": headline, "figure": plot,
           "constants_used": {k: v for r in results for k, v in r["constants_used"].items()},
           "registry_sha256": constants_sha256(), **extra}
    dump(OUT / f"{name}.json", doc)
    return doc


def fmt(x, d=3):
    return f"{x:,.0f}" if abs(x) >= 1000 else f"{x:.{d}g}"


def um(x):
    return None if x is None else x * UM


def uptake_series(levels, xs, key="h2_uptake_mol_m2", label="twin H2 uptake"):
    ys = [l["twin"][key] if "twin" in l else None for l in levels]
    return {"label": label, "x": xs, "y": [um(y["value"]) if y and y["value"] is not None else None for y in ys],
            "low": [um(y["low"]) if y and y["low"] is not None else None for y in ys],
            "high": [um(y["high"]) if y and y["high"] is not None else None for y in ys]}


def plot_spec(title, xlabel, ylabel, series, baseline_y=None):
    return {"title": title, "xlabel": xlabel, "ylabel": ylabel, "series": series, "baseline_y": baseline_y}


def base_h2_umol():
    return um(run_chain({}, [], req(), [], envelope=False)["twin"]["baseline_outputs"]["h2_uptake_mol_m2"])


# ---------------------------------------------------------------- the experiments
CROPS = ("date_palm", "barley", "alfalfa", "tomato")


def exp01_brackish_irrigation():
    """What if we irrigate with brackish groundwater or treated wastewater?

    The water-demand model sets how much water is applied (daily); its salt content follows from its conductivity. The twin
    tracks the salt and water in the soil; the FAO-29 relations and Maas-Hoffman curves (other models) turn the soil salinity into a crop's relative yield."""
    ups = ["ec_to_nacl", "crop_water_demand_generic", "irrigation_schedule"]
    downs = ["twin_solution_to_ece"] + [f"maas_hoffman_{c}" for c in CROPS]
    levels, runs = [], []
    ecs = (0.5, 2, 4, 8)
    for ec in ecs:
        rec, r = level(f"EC_w={ec} dS/m", {"ec_w_dS_m": ec, "irrigation_interval_days": 1}, ups, req(), downs)
        levels.append(rec); runs.append(r)
    ry = {c: [l["downstream"][f"maas_hoffman_{c}"][f"relative_yield_{c}"] if not l["refused"] else None for l in levels] for c in CROPS}
    series = [{"label": c.replace("_", " "), "x": list(ecs), "y": [w and w["value"] for w in ry[c]], "low": [w and w["low"] for w in ry[c]],
               "high": [w and w["high"] for w in ry[c]]} for c in CROPS]
    a = next(l for l in levels if l["action"]["ec_w_dS_m"] == 8)
    alf4 = a["downstream"]["maas_hoffman_alfalfa"]["relative_yield_alfalfa"]
    h = head("Alfalfa relative yield, EC_w = 8 dS/m (daily irrigation)", alf4, "fraction of unsalted yield", "external:maas_hoffman_alfalfa",
             "soil salinity from the twin's salt balance (14-day window) via FAO-29")
    e8 = levels[-1]
    plain = (f"We irrigated a sandy desert soil every day with water of increasing saltiness. The twin shows how much salt builds up in the soil "
             f"(mean {fmt(e8['twin']['mean_nacl_mol_l']['value'] * 1000, 3)} mmol/L at 8 dS/m) and that the soil microbes' hydrogen uptake falls "
             f"{fmt(-e8['difference_vs_baseline']['h2_uptake_mol_m2']['percent']['value'], 3)} % compared with the unirrigated dry soil, mostly because the soil stays wet. "
             f"How much a crop would suffer (e.g. alfalfa relative yield {fmt(alf4['value'], 3)} at 8 dS/m; none of the crops loses yield below 4 dS/m here, because the twin's bucket leaches most of the salt) comes from published crop-salinity curves, "
             "not from the twin; the twin does not model plants.")
    cav = ["The 14-day window with a 12 mm rain pulse is the v0.6 baseline; the 'season' is these 14 days, so salt accumulation is short-term.",
           "Irrigation amount = FAO-56 generic crop water demand (registry) every day; leaching is whatever the twin's bucket gives (drainage 50 % per day above field capacity), far higher than the 15-20 % assumed in FAO-29.",
           "Relative yields use Maas-Hoffman thresholds and slopes (steady-state tolerance), applied to a 14-day mean soil solution salinity converted with registry factors.",
           "The twin abstains when the soil solution leaves its ideal-dilute NaCl domain (0.1 mol/L); low/high are None where any envelope run abstains (see envelope_refusals). " + SRC]
    return finish(exp01_brackish_irrigation, "brackish_irrigation", ["ec_to_nacl", "crop_water_demand_generic", "irrigation_schedule", "twin"] + downs,
                  levels, runs, [h], plot_spec("1 Brackish irrigation", "irrigation water EC (dS/m)", "relative yield (Maas-Hoffman)", series),
                  plain, cav)


def exp02_drip_vs_flood():
    """Same water, delivered as daily drip or occasional flooding?

    A fixed total of irrigation water (crop-water-demand model times the period) is split into equal pulses every 1, 3 or 10 days; the twin computes the soil's hydrogen and CO uptake."""
    ups = ["crop_water_demand_generic", "irrigation_schedule"]
    levels, runs, ks = [], [], (1, 3, 10)
    for k in ks:
        rec, r = level(f"every {k} d", {"irrigation_interval_days": k, "irrigation_nacl_mol_l": 0.0}, ups, req(), [])
        levels.append(rec); runs.append(r)
    u = {k: l["twin"]["h2_uptake_mol_m2"]["value"] for k, l in zip(ks, levels)}
    contrast = {"h2_uptake_10d_over_1d": u[10] / u[1], "h2_uptake_3d_over_1d": u[3] / u[1],
                "v06_finding_reproduced_daily_lowest": u[1] < u[3] < u[10] or (u[1] < u[10] and u[1] < u[3])}
    h = h2_head("H2 uptake change vs unirrigated baseline, flood-type (every 10 days)", levels[2], "equal seasonal total, 14-day window")
    plain = (f"We gave the soil the same total amount of water either every day, every 3 days or only every 10 days. The twin says hydrogen uptake is "
             f"{fmt(u[1] * UM, 3)}, {fmt(u[3] * UM, 3)} and {fmt(u[10] * UM, 3)} micromole per square metre over two weeks, so "
             f"{'rare big irrigations keep more soil activity than daily small ones' if u[10] > u[1] else 'delivery pattern matters little here'}. "
             "The total water amount comes from a published crop-water-demand figure (another model); everything about the soil response is the twin's.")
    cav = ["Contrasts with no-irrigation baseline are large because 8-fold more water than the 12 mm rain keeps the soil near saturation; compare the three levels with each other.",
           "Total = FAO-56 generic demand x 14 days (registry range 4-9 mm/day) and is applied on top of the 12 mm rain pulse.",
           "Consistent with the v0.6 finding that daily pulses give the lowest uptake (wet soil limits gas diffusion).", SRC]
    return finish(exp02_drip_vs_flood, "drip_vs_flood", ups + ["twin"], levels, runs, [h],
                  plot_spec("2 Drip vs flood", "irrigation interval (days)", "H2 uptake (umol m$^{-2}$)", [uptake_series(levels, list(ks))], base_h2_umol()),
                  plain, cav, contrast=contrast)


def exp03_mulch_shade():
    """What if we mulch or put up shade nets?

    Mulch cuts evaporation from the soil and a shade net lowers soil temperature (both from published ranges); the twin recomputes the soil's gas uptake."""
    ups = ["mulch_evaporation", "shade_net_temperature"]
    combos = [(0, 0), (0.5, 0), (1, 0), (0, 0.5), (0, 1), (1, 1)]
    levels, runs = [], []
    for m, s in combos:
        rec, r = level(f"mulch {m}, shade {s}", {"mulch_cover_fraction": m, "shade_net_cover_fraction": s}, ups, req(), [])
        levels.append(rec); runs.append(r)
    ser = [uptake_series([l for l, c in zip(levels, combos) if c[1] == 0], [c[0] for c in combos if c[1] == 0], label="mulch (no shade)"),
           uptake_series([l for l, c in zip(levels, combos) if c[0] == 0], [c[1] for c in combos if c[0] == 0], label="shade net (no mulch)")]
    mul = levels[2]; sh = levels[4]
    h = h2_head("H2 uptake change, full mulch cover", mul, "evaporation reduced by registry fraction")
    plain = (f"Full mulch changes the soil's hydrogen uptake by {fmt(pct(mul)['value'], 3)} % and a full shade net by {fmt(pct(sh)['value'], 3)} % over two weeks "
             "(twin). How much mulch cuts evaporation and how much a net cools the soil are published ranges (other models), not twin results; "
             "neither changes how much water the soil receives.")
    cav = ["Mulch effect is bounded by the single 12 mm pulse of the baseline: little water to save in 14 days.",
           "Shade-net cooling is applied as a constant offset on soil temperature; the twin's Q10 domain is 0-50 C.", SRC]
    return finish(exp03_mulch_shade, "mulch_shade", ups + ["twin"], levels, runs, [h],
                  plot_spec("3 Mulch / shade net", "cover fraction", "H2 uptake (umol m$^{-2}$)", ser, base_h2_umol()), plain, cav)


def exp04_biochar():
    """What if we add biochar?

    Biochar dose (percent by weight) raises soil porosity by a meta-analysis rate; retention is left unchanged. The twin recomputes gas uptake with the new porosity."""
    ups = ["biochar_porosity"]
    ds = (0, 1, 2, 5)
    levels, runs = [], []
    for d in ds:
        rec, r = level(f"{d} % w/w", {"biochar_dose_pct_w_w": d}, ups, req(), [])
        levels.append(rec); runs.append(r)
    h = h2_head("H2 uptake change, 5 % w/w biochar", levels[-1])
    plain = (f"Adding 5 % biochar by weight changes the soil's hydrogen uptake by {fmt(pct(levels[-1])['value'], 3)} % over two weeks (twin). "
             "How much the biochar opens up the soil pores is a published average (another model); the twin only converts the new pore space into gas uptake.")
    cav = ["Only porosity is changed (registry: absolute porosity per % w/w); water retention and biochar's own chemistry/microbiology are not modelled.",
           "Higher porosity also raises the water capacity in the twin, which can change moisture and activation.", SRC]
    return finish(exp04_biochar, "biochar", ups + ["twin"], levels, runs, [h],
                  plot_spec("4 Biochar", "dose (% w/w)", "H2 uptake (umol m$^{-2}$)", [uptake_series(levels, list(ds))], base_h2_umol()), plain, cav)


def exp05_inoculation():
    """What if we add microbes that oxidise CO or H2?

    A reference genome (MAG) is added to a sample's community at 0.1, 1 or 5 percent; the ontology decides which gas classes it belongs to and the twin recomputes the community factors."""
    tab = comm._load()["table"]
    classes = comm._load()["classes"]
    co_mag = next(m for m in sorted(classes) if "COCarrierMAG" in classes[m])
    null_mag = next(m for m in sorted(classes) if not classes[m])
    closed_sample = "51PRr3"
    open_sample = tab["COCarrierMAG"].idxmax()
    cases = [("51PRr3 (CO gate closed) + CO carrier", closed_sample, co_mag),
             (f"{open_sample} (CO gate already open) + CO carrier", open_sample, co_mag),
             ("51PRr3 + non-carrier MAG (negative control)", closed_sample, null_mag)]
    shares = (0.001, 0.01, 0.05)
    levels, runs = [], []
    for name, sample, mag in cases:
        for mode in ("gate", "proportional"):
            for s in shares:
                rec, r = level(f"{name}; {mode}; share {s}", {"inoculum_mag_id": mag, "inoculum_share": s}, ["inoculant"],
                               req(community={"sample": sample}, coupling=mode), [])
                rec["case"] = name; rec["coupling"] = mode; rec["share"] = s; rec["sample"] = sample; rec["mag_id"] = mag
                levels.append(rec); runs.append(r)
    pick = lambda c, m, s: next(l for l in levels if l["case"].startswith(c) and l["coupling"] == m and l["share"] == s)
    g = pick("51PRr3 (CO", "gate", 0.01)
    p5 = pick("51PRr3 (CO", "proportional", 0.05)
    ctrl = pick("51PRr3 + non", "gate", 0.05)
    ctrlp = pick("51PRr3 + non", "proportional", 0.05)
    opn = pick(open_sample, "gate", 0.01)
    gate_summary = {"closed_sample": closed_sample, "co_gate_opened_by_carrier": g["inoculation"]["gate_opened"]["co"],
                    "co_uptake_before_umol_m2": um(g["twin_baseline"]["co_uptake_mol_m2"]), "co_uptake_after_umol_m2": um(g["twin"]["co_uptake_mol_m2"]["value"]),
                    "open_sample": open_sample, "gate_already_open": not opn["inoculation"]["gate_opened"]["co"],
                    "open_sample_co_change_gate_umol_m2": um(opn["difference_vs_baseline"]["co_uptake_mol_m2"]["value"]),
                    "negative_control_co_change_gate_umol_m2": um(ctrl["difference_vs_baseline"]["co_uptake_mol_m2"]["value"]),
                    "negative_control_h2_change_gate_umol_m2": um(ctrl["difference_vs_baseline"]["h2_uptake_mol_m2"]["value"]),
                    "negative_control_h2_change_proportional_percent": pct(ctrlp)["value"] if pct(ctrlp) else None,
                    "negative_control_note": "gate: exactly no change. proportional: the inoculum displaces the existing community (shares rescaled by 1-share), so uptake falls by about the share."}
    series = []
    for name, sample, mag in cases[:2]:
        for mode in ("gate", "proportional"):
            ls = [l for l in levels if l["case"] == name and l["coupling"] == mode]
            series.append(uptake_series(ls, [100 * x for x in shares], "co_uptake_mol_m2", f"{sample} {mode}"))
    h = head("CO uptake after 1 % CO-carrier inoculation of 51PRr3 (gate)", g["twin"]["co_uptake_mol_m2"], "mol m-2 over 14 d", "twin",
             "baseline CO uptake of this sample is exactly 0 (gate closed)")
    h = {**h, "value": um(h["value"]), "low": um(h["low"]), "high": um(h["high"]), "unit": "umol m-2 over 14 d"}
    plain = ("In the sample 51PRr3 the ontology finds no carbon-monoxide-oxidising organism, so the twin predicts no CO uptake at all. Adding even a tiny share "
             f"of a genome known to oxidise CO opens that process (CO uptake {fmt(gate_summary['co_uptake_after_umol_m2'], 3)} micromole per square metre over two weeks). "
             "The 'gate' result says only whether the process can run at all; the 'proportional' result says how much, assuming the rate scales with the share of such organisms. "
             "Adding an organism that is not a carrier changes nothing under the gate. Whether the added microbes would survive and establish in desert soil is not modelled: the share is an assumption of the scenario.")
    cav = ["Share is an assumption of the scenario; survival, establishment, competition and activity of the inoculum in desert soil are not modelled.",
           "Gate mode: any share > 0 sets gamma = 1, so the dose has no graded effect; the dose-response is the proportional mode, gamma = q/median(q) (assumes rate scales with carrier share).",
           "Proportional mode: a non-carrier inoculum displaces the existing community and lowers uptake by about its share (a model artefact of renormalisation, not a biological prediction).",
           "Memberships come from reasoner entailment (K and E unchanged); absence of entailment is not biological absence.", SRC]
    return finish(exp05_inoculation, "inoculation", ["inoculant", "twin"], levels, runs, [h],
                  plot_spec("5 Inoculation (CO uptake)", "inoculum share (%)", "CO uptake (umol m$^{-2}$)", series), plain, cav,
                  gate=gate_summary, inoculum_mags={"co_carrier": co_mag, "non_carrier": null_mag})


def exp06_carbon_fixation_bound():
    """How much CO2 could the H2 oxidisers fix?

    An upper bound: CO2 fixed cannot exceed the twin's H2 uptake times the CO2-per-H2 yield of knallgas autotrophy, and only if the community contains an organism entailed to be an H2-autotroph candidate."""
    tab = comm._load()["table"]
    zero = sorted(tab.index[tab["H2AutotrophicCandidateMAG"] == 0])[0]
    med = comm.medians()["H2AutotrophicCandidateMAG"]
    pos_pool = tab[tab["H2AutotrophicCandidateMAG"] > 0]["H2AutotrophicCandidateMAG"]
    pos = (pos_pool - med).abs().sort_values().index[0]
    scen = {"baseline (12 mm pulse)": FORCING,
            "wet (36 mm pulse)": [{**f, "water_mm": 36. if d == 0 else 0.} for d, f in enumerate(FORCING)]}
    downs = ["knallgas_carbon_fixation", "co_oxidation_carbon_balance"]
    levels, runs = [], []
    for sname, fc in scen.items():
        for sample in (zero, pos):
            rec, r = level(f"{sname}; {sample}", {}, [], req(fc, community={"sample": sample}, coupling="gate"), downs)
            k = rec["downstream"]["knallgas_carbon_fixation"]
            ann = lambda w: None if w is None else xm.annualise(w, DAYS)
            c_ug = k["c_fixed_upper_bound_ug_m2"]
            co2 = k["co2_fixed_upper_bound_mol_m2"]
            rec["sample"] = sample; rec["scenario"] = sname
            rec["h2_autotroph_share"] = comm.carrier_shares(sample)["shares"]["H2AutotrophicCandidateMAG"]
            rec["co2_fixed_upper_bound"] = {"per_period_umol_m2": {x: um(co2[x]) for x in ("value", "low", "high")},
                                            "per_year_umol_m2": {x: um(ann(co2[x])) for x in ("value", "low", "high")},
                                            "provenance": co2["provenance"]}
            rec["c_fixed_upper_bound"] = {"per_period_ug_c_m2": {x: c_ug[x] for x in ("value", "low", "high")},
                                          "per_year_ug_c_m2": {x: ann(c_ug[x]) for x in ("value", "low", "high")}, "provenance": c_ug["provenance"]}
            ctx = REG["biocrust_net_c_fixation_g_c_m2_yr"]
            rec["context"] = {"label": "context (not a twin result): desert biocrust net C fixation", "value_ug_c_m2_yr": ctx["value"] * 1e6,
                              "low": ctx["low"] * 1e6, "high": ctx["high"] * 1e6, "unit": "ug C m-2 yr-1", "constant": ctx["id"]}
            rec["bound_over_context"] = (rec["c_fixed_upper_bound"]["per_year_ug_c_m2"]["value"] / (ctx["value"] * 1e6)
                                         if rec["c_fixed_upper_bound"]["per_year_ug_c_m2"]["value"] is not None else None)
            rec["co_oxidation"] = {"co2_from_co_umol_m2": um(rec["downstream"]["co_oxidation_carbon_balance"]["co2_from_co_mol_m2"]["value"]),
                                   "net_atmospheric_c_change_mol_m2": rec["downstream"]["co_oxidation_carbon_balance"]["net_atmospheric_c_change_mol_m2"]["value"],
                                   "note": "CO + 1/2 O2 -> CO2: each CO carbon atom reappears as CO2; no net change in atmospheric carbon atoms."}
            levels.append(rec); runs.append(r)
    pb = next(l for l in levels if l["scenario"].startswith("baseline") and l["sample"] == pos)
    zb = next(l for l in levels if l["scenario"].startswith("baseline") and l["sample"] == zero)
    wet = next(l for l in levels if l["scenario"].startswith("wet") and l["sample"] == pos)
    cy = pb["c_fixed_upper_bound"]["per_year_ug_c_m2"]
    h = head(f"Upper bound of C fixed by H2 oxidisers, sample {pos}", cy, "ug C m-2 yr-1 (annualised)", "external:knallgas_carbon_fixation",
             "bound = H2 uptake (twin) x yield; annualised linearly from the 14-day period")
    series = [{"label": s, "x": [i], "y": [l["c_fixed_upper_bound"]["per_year_ug_c_m2"]["value"]],
               "low": [l["c_fixed_upper_bound"]["per_year_ug_c_m2"]["low"]], "high": [l["c_fixed_upper_bound"]["per_year_ug_c_m2"]["high"]]}
              for i, (s, l) in enumerate([("baseline", pb), ("wet", wet)])]
    plain = (f"If the soil's hydrogen eaters used all of that hydrogen to build biomass, the most carbon they could fix is about {fmt(cy['value'] / 1e6, 3)} "
             f"grams per square metre per year, in a typical sample (the twin supplies the hydrogen uptake; the carbon yield is another model's number). "
             f"In a sample where the ontology finds no such organism ({zero}) the answer is exactly zero by deduction. "
             f"For scale, desert biological soil crusts are reported to fix roughly {fmt(ctx['value'], 2)} g C per square metre per year (context only), and CO oxidation just turns CO into CO2 without changing the atmosphere's carbon.")
    cav = ["An upper bound, not a prediction or net carbon removal: autotroph growth, respiration, turnover and community composition are not modelled; yield 1/6 is a registry value with range 0.1-0.25.",
           "Annualisation is linear scaling of the 14-day total (days_per_year registry): a labelled extrapolation, not a forecast; real activity depends on the year's wetting events.",
           "Zero for the sample with no H2AutotrophicCandidateMAG share is deduced from the ontology (absence of entailment is not biological absence).",
           "Context value is a literature figure for biocrusts, unverified in the registry.", SRC]
    return finish(exp06_carbon_fixation_bound, "carbon_fixation_bound", ["twin"] + downs, levels, runs, [h],
                  plot_spec("6 Carbon-fixation bound", "scenario (0 baseline, 1 wet)", "C bound (ug C m$^{-2}$ yr$^{-1}$)", series), plain, cav,
                  gate={"zero_share_sample": zero, "zero_sample_bound_provenance": zb["downstream"]["knallgas_carbon_fixation"]["co2_fixed_upper_bound_mol_m2"]["provenance"],
                        "positive_share_sample": pos, "positive_sample_share": comm.carrier_shares(pos)["shares"]["H2AutotrophicCandidateMAG"]})


def exp07_h2_leakage():
    """What if a hydrogen economy raises atmospheric H2?

    The atmospheric H2 mixing ratio is set to 1, 1.5, 2 or 10 times the background; the twin's soil sink (H2 uptake per square metre) is recomputed. The 10x level is a hypothetical near-source plume."""
    ms = (1, 1.5, 2, 10)
    levels, runs, o2 = [], [], []
    p = Parameters(k_h2_s=.001, k_co_s=.001)
    for m in ms:
        rec, r = level(f"{m}x background", {"h2_multiplier": m}, ["h2_atmosphere"], req(), [], keep_twin_result=True)
        rows = r["twin"]["result"]["intervened"]["trajectory"]
        use = max((x["gas"]["oxygen_used_mol_m2_s"] for x in rows if x["gas"]), default=0.)
        o2.append(use / p.oxygen_supply_mol_m2_s)
        rec["oxygen_used_over_supply_max"] = use / p.oxygen_supply_mol_m2_s
        rec["upstream_ppm"] = rec["upstream"]["h2_atmosphere"]["h2_mixing_ratio_ppm"]
        levels.append(rec); runs.append(r)
    u = [l["twin"]["h2_uptake_mol_m2"]["value"] for l in levels]
    per_mult = [x / m for x, m in zip(u, ms)]
    lin = {"uptake_per_multiplier_mol_m2": per_mult, "max_relative_deviation_from_linear": max(abs(x / per_mult[0] - 1) for x in per_mult),
           "linear": max(abs(x / per_mult[0] - 1) for x in per_mult) < 1e-6, "max_oxygen_used_over_supply": max(o2),
           "oxygen_limited": max(o2) >= 1 - 1e-9,
           "note": "Deposition flux is linear in the atmospheric concentration; oxygen limits only if H2+CO demand exceeds the oxygen supply."}
    for l in levels:
        l.pop("twin_baseline", None)
    ctx = REG["h2_hydrogen_economy_global_ppm"]
    h = h2_head("H2 sink change at 2x background", levels[2], "sink per m2 of soil, 14 days")
    o2_text = "oxygen does not run short" if not lin["oxygen_limited"] else "oxygen becomes limiting"
    plain = (f"If the air's hydrogen doubles, the soil's hydrogen uptake doubles too ({fmt(u[2] * UM, 3)} instead of {fmt(u[0] * UM, 3)} micromole per square metre over two weeks): the twin's soil sink grows in proportion to the concentration. "
             f"At ten times background the uptake is {fmt(u[3] * UM, 3)} and {o2_text}. "
             f"The atmospheric concentrations are scenarios from the hydrogen-economy literature (another source, about {fmt(ctx['value'], 2)} ppm for a global hydrogen economy), not twin outputs.")
    cav = ["Linear response follows from the first-order uptake law at fixed activity; a larger sink does not change the atmosphere (no feedback, no mass balance with the atmosphere).",
           "Only soil uptake is modelled: no leakage rates, no atmospheric chemistry, no other sinks.", "10x is a hypothetical near-source plume.", SRC]
    return finish(exp07_h2_leakage, "h2_leakage", ["h2_atmosphere", "twin"], levels, runs, [h],
                  plot_spec("7 H$_2$ leakage", "H$_2$ multiple of background", "H2 uptake (umol m$^{-2}$)", [uptake_series(levels, list(ms))], None),
                  plain, cav, linearity=lin)


def exp08_solar_park():
    """What if we build a solar park?

    Under panels the soil is cooler and evaporates less (published ranges) and panel washing adds water; the twin recomputes the gas uptake under the panels."""
    ups = ["solar_park_microclimate", "panel_washing"]
    fs = (0, 0.5, 1)
    levels, runs = [], []
    for f in fs:
        rec, r = level(f"panel cover {f}", {"panel_cover_fraction": f}, ups, req(), [])
        levels.append(rec); runs.append(r)
    h = h2_head("H2 uptake change, full panel cover", levels[-1], "cooling, shading of evaporation and one washing pulse")
    plain = (f"Under a full solar-panel cover the twin's hydrogen uptake changes by {fmt(pct(levels[-1])['value'], 3)} % over two weeks. "
             "How much cooler and less evaporative it is under panels, and how much water panel washing adds, are published ranges (other models); the soil's response to them is the twin's.")
    cav = ["Washing is one wash per 14 days (registry) so at most one pulse falls in the window; its water is tiny against the 12 mm rain.",
           "Cover fraction scales effects linearly; panel shading of rain, dust and soil disturbance during construction are not modelled.", SRC]
    return finish(exp08_solar_park, "solar_park", ups + ["twin"], levels, runs, [h],
                  plot_spec("8 Solar park", "panel cover fraction", "H2 uptake (umol m$^{-2}$)", [uptake_series(levels, list(fs))], base_h2_umol()), plain, cav)


def exp09_cloud_seeding():
    """What if cloud seeding adds rain?

    Rain on existing rain days is increased by 0, 5, 10 or 20 percent (published efficacy is uncertain: 0-20 percent for full-intensity seeding); the twin recomputes the soil's gas uptake."""
    ints = (0, 0.5, 1, 2)
    levels, runs = [], []
    for s in ints:
        rec, r = level(f"intensity {s}", {"seeding_intensity": s}, ["cloud_seeding"], req(), [])
        levels.append(rec); runs.append(r)
    h = h2_head("H2 uptake change, full-intensity seeding (+10 % rain on rain days)", levels[2], "baseline has one rain day (12 mm) in 14")
    plain = (f"Adding 10 % to the rain that already falls changes the soil's hydrogen uptake by {fmt(pct(levels[2])['value'], 3)} % over two weeks (twin). "
             "In this baseline there is a single rain day, so the honest answer is a small effect; how much extra rain seeding really yields is uncertain and comes from the weather-modification literature, not the twin.")
    cav = ["Seeding cannot create rain on dry days; with one rain day in the baseline the possible effect is small by construction.",
           "Efficacy range 0-20 % (registry; WMO statement says evidence is weak); low end is no effect.", SRC]
    return finish(exp09_cloud_seeding, "cloud_seeding", ["cloud_seeding", "twin"], levels, runs, [h],
                  plot_spec("9 Cloud seeding", "seeding intensity (1 = +10 % rain)", "H2 uptake (umol m$^{-2}$)", [uptake_series(levels, list(ints))], base_h2_umol()), plain, cav)


def exp10_compaction():
    """What if off-road traffic compacts the soil?

    Vehicle passes (0, 1, 10, 50) raise bulk density following a published saturating curve; porosity falls accordingly and the twin recomputes the gas uptake."""
    ps = (0, 1, 10, 50)
    levels, runs = [], []
    for n in ps:
        rec, r = level(f"{n} passes", {"offroad_passes": n}, ["offroad_compaction"], req(), [])
        rec["porosity_after"] = rec["upstream"]["offroad_compaction"]["porosity"]
        levels.append(rec); runs.append(r)
    h = h2_head("H2 uptake change, 50 passes", levels[-1])
    plain = (f"After 50 vehicle passes the soil packs down and its hydrogen uptake changes by {fmt(pct(levels[-1])['value'], 3)} % over two weeks (twin). "
             "How much a pass compacts desert sand is a published relationship (another model); the twin converts the smaller pore space into gas uptake.")
    cav = ["Compaction is represented only by a porosity decrease with unchanged retention parameters.", "Crust destruction and surface disturbance (which may matter more) are not modelled.",
           "The response is not monotone: in the dry baseline a little compaction raises activation (the same water fills less pore space, so the soil is wetter) while diffusion falls; uptake peaks near porosity 0.36 and falls beyond (checked by direct porosity scan). Biochar (higher porosity) moves the other way.", SRC]
    return finish(exp10_compaction, "compaction", ["offroad_compaction", "twin"], levels, runs, [h],
                  plot_spec("10 Compaction", "vehicle passes", "H2 uptake (umol m$^{-2}$)", [uptake_series(levels, list(ps))], base_h2_umol()), plain, cav)


def exp11_warming():
    """What if the climate warms?

    Soil temperature rises with air warming (+1.5, +2, +4 C), evaporation scales with temperature and rain changes by -20, 0 or +20 percent; the twin recomputes the gas uptake."""
    ws, rs = (0, 1.5, 2, 4), (-0.2, 0, 0.2)
    levels, runs = [], []
    for w in ws:
        for rc in rs:
            rec, r = level(f"+{w} C, rain {rc:+.0%}", {"warming_c": w, "rain_change_frac": rc}, ["climate_warming"], req(), [])
            rec["warming_c"] = w; rec["rain_change_frac"] = rc
            levels.append(rec); runs.append(r)
    ser = [uptake_series([l for l in levels if l["rain_change_frac"] == rc], list(ws), label=f"rain {rc:+.0%}") for rc in rs]
    sel = next(l for l in levels if l["warming_c"] == 4 and l["rain_change_frac"] == 0)
    h = h2_head("H2 uptake change, +4 C, rain unchanged", sel)
    plain = (f"With 4 degrees of warming and the same rain, the twin's hydrogen uptake changes by {fmt(pct(sel)['value'], 3)} % over two weeks, because warm soil microbes work faster but the soil also dries sooner. "
             "How much the soil warms and how much more it evaporates are published rules of thumb (other models); the microbial response is the twin's.")
    cav = ["Warming is a uniform offset over the 14-day window; changes in rain timing, extremes and seasonality are not modelled.",
           "Temperature must stay in the twin's 0-50 C domain.", "Rain change scales the existing rain pulse only.", SRC]
    return finish(exp11_warming, "warming", ["climate_warming", "twin"], levels, runs, [h],
                  plot_spec("11 Warming", "air warming (C)", "H2 uptake (umol m$^{-2}$)", ser, base_h2_umol()), plain, cav)


def exp12_inoculated_alfalfa_salinity():
    """What if we irrigate alfalfa with brackish water and coat the seed with the desert bacterium Enterobacter sp. SA187?

    Three alternative models of the root-zone salinity (FAO fixed leaching, the twin's salt balance, and one factor fitted to the field trial's yield drop) feed the Maas-Hoffman alfalfa curve; the SA187 benefit is the field trial's (de Zelicourt et al. 2018, Hada Al-Sham, KSA)."""
    ecs = (0.5, 2, 3.12, 4, 7.81, 8, 12)
    doc = xm.digitised_trial()
    err = xm.digitisation_error(doc)
    fcons = constant_values(sorted(REG))
    models = {"A_fao_fixed_leaching": ([], ["fao29_leaching", "maas_hoffman_alfalfa", "sa187_benefit"]),
              "B_twin_physics": (["ec_to_nacl", "crop_water_demand_alfalfa", "irrigation_schedule"],
                                 ["twin_solution_to_ece", "maas_hoffman_alfalfa", "sa187_benefit"]),
              "C_fitted_fresh": ([], ["fitted_leaching_fresh", "maas_hoffman_alfalfa", "sa187_benefit"]),
              "C_fitted_dry": ([], ["fitted_leaching_dry", "maas_hoffman_alfalfa", "sa187_benefit"])}
    levels, runs, curves = [], [], {}
    for mname, (ups, downs) in models.items():
        for ec in ecs:
            act = {"ec_w_dS_m": ec, **({"irrigation_interval_days": 1} if ups else {})}
            rec, r = level(f"{mname}; EC_w={ec}", act, ups, req(), downs)
            rec["model"] = mname; rec["ec_w_dS_m"] = ec
            dn = rec.get("downstream", {})
            rec["relative_yield_mock"] = dn.get("maas_hoffman_alfalfa", {}).get("relative_yield_alfalfa")
            rec["relative_yield_sa187"] = dn.get("sa187_benefit", {}).get("relative_yield_alfalfa_sa187")
            levels.append(rec); runs.append(r)
            curves.setdefault(mname, {})[ec] = rec

    def ratio_env(mname):
        """Predicted RY(7.81)/RY(3.12): central, and the min/max over one-at-a-time constant ranges."""
        ups, downs = models[mname]
        ids = sorted({c for m in ups + downs for c in M[m].constants})
        vals = []
        for cid in [None] + ids:
            for end in (("value",) if cid is None else ("low", "high")):
                ov = {} if cid is None else {cid: REG[cid][end]}
                ry = []
                for ec in (3.12, 7.81):
                    r = run_chain({"ec_w_dS_m": ec, **({"irrigation_interval_days": 1} if ups else {})}, [M[u] for u in ups], req(),
                                  [M[d] for d in downs], constants=ov, envelope=False)
                    ry.append(None if r["refused"] else r["outputs"]["maas_hoffman_alfalfa:relative_yield_alfalfa"]["value"])
                if None not in ry and ry[0] > 0:
                    vals.append(ry[1] / ry[0])
        return {"value": vals[0], "low": min(vals), "high": max(vals)}

    obs = {t: xm.observed_mock_ratio(t, doc) for t in ("fresh", "dry")}
    pred = {m: ratio_env(m) for m in models}
    fit = {}
    for t in ("fresh", "dry"):
        o = REG[f"sa187_mock_ratio_{t}"]
        fv = {k: xm.fit_leaching_factor({**fcons, f"sa187_mock_ratio_{t}": o[k]}, t) for k in ("value", "low", "high")}
        fit[t] = {"f_central": fv["value"], "f_range_from_digitisation_error": sorted([fv["low"], fv["high"]]),
                  "ece_at_3_12_dS_m": fv["value"] * 3.12, "ece_at_3_12_range": sorted([fv["low"] * 3.12, fv["high"] * 3.12]),
                  "observed_ratio": obs[t]["value"], "observed_ratio_range": [obs[t]["low"], obs[t]["high"]]}
    soil = REG["sa187_soil_ec_dS_m"]["value"]
    fit["comparison_with_measured_soil_ec"] = {
        "measured_soil_ec_dS_m": soil, "note": "Independent measurement at the site (which EC definition is not stated in the registry); compare with f x 3.12.",
        "A_1.5x3.12": 1.5 * 3.12, "fitted_fresh_f_x_3.12": fit["fresh"]["ece_at_3_12_dS_m"], "fitted_dry_f_x_3.12": fit["dry"]["ece_at_3_12_dS_m"],
        "B_twin_ece_at_3.12": curves["B_twin_physics"][3.12]["downstream"]["twin_solution_to_ece"]["ece_dS_m"]["value"]
        if "downstream" in curves["B_twin_physics"][3.12] and not curves["B_twin_physics"][3.12]["refused"] else None}
    validation = {m: {"predicted_ratio": pred[m], **{f"observed_{t}": obs[t]["value"] for t in ("fresh", "dry")},
                      "gap_to_observed_fresh": pred[m]["value"] - obs["fresh"]["value"], "gap_to_observed_dry": pred[m]["value"] - obs["dry"]["value"]}
                  for m in ("A_fao_fixed_leaching", "B_twin_physics")}
    ratio_note = ("Predicted and observed are mock (uninoculated) yield ratios, high (7.81) over low (3.12 dS/m), season 1 only. Fitted models match "
                  "the observed ratio by construction (not a validation).")
    mock_yields = {t: {"low_3.12_t_ha": xm.mock_yield(t, "low", doc), "high_7.81_t_ha": xm.mock_yield(t, "high", doc)} for t in ("fresh", "dry")}
    err_note = (f"Digitisation error used: absolute {err['absolute_t_ha']} t/ha (largest y-axis calibration residual of the digitiser) plus relative "
                f"{err['relative']} (largest disagreement between digitised SA187 effects and the percentages printed on the figure); the eye-read +-0.15 t/ha is dropped. "
                "Ratio range = extreme combinations.")
    h0 = curves["C_fitted_fresh"][7.81]["relative_yield_sa187"]
    h = head("Inoculated alfalfa relative yield at EC_w 7.81 dS/m (salinity fitted to the field trial, fresh weight)", h0, "fraction of unsalted, uninoculated yield",
             "external:sa187_benefit", "fitted leaching model; relative yield is the crop model's and the field trial's, not the twin's")
    alt = {m: curves[m][7.81]["relative_yield_sa187"] for m in models}
    series = []
    for m in models:
        ys = [curves[m][e]["relative_yield_sa187"] for e in ecs]
        series.append({"label": m.split("_")[0] + " " + " ".join(m.split("_")[1:]), "x": list(ecs), "y": [w and w["value"] for w in ys],
                       "low": [w and w["low"] for w in ys], "high": [w and w["high"] for w in ys]})
    plain = ("At the field trial's two salinities the seed-coated bacterium gave roughly 12 to 17 percent more alfalfa than uncoated seed. "
             f"The trial's own yields fall by about {fmt(100 * (1 - obs['fresh']['value']), 2)} percent (fresh weight) from 3.12 to 7.81 dS/m irrigation water. "
             "The prediction of how salty the root zone gets, and how a crop reacts, comes from three alternative models; the twin contributes only one of them, a salt balance of a bare sandy desert soil, "
             "and plant growth is the crop model's and the field trial's, not the twin's. Whether the bacterium's benefit grows with salinity cannot be told from these data, which come from one trial at two salinities.")
    cav = ["Benefit b = 0.15 (0.12-0.17) is applied as a constant multiplier; the data (one trial, two salinities, two seasons) cannot show whether it depends on salinity.",
           "Season 2 is confounded by exceptional rain (authors) and is not used; mock yields are pixel-digitised bar heights, not a published table.",
           err_note, ratio_note,
           "Model B uses the twin's Rub' al Khali sand bucket (bare soil: evaporation 1 mm/day, 14 days, daily irrigation of registry alfalfa demand), not Hada Al-Sham soil; leaching is therefore much higher than FAO's 15-20 %, so B underestimates root-zone salinity relative to a real stand with high evapotranspiration. A 14-day mean is used for the season.",
           "Where the twin leaves its ideal-dilute NaCl domain, model B reports no value (refusal) and 12 dS/m may be refused.",
           "Maas-Hoffman alfalfa parameters (threshold 2.0, slope 7.3) and FAO-29 factors are unverified registry constants."]
    return finish(exp12_inoculated_alfalfa_salinity, "inoculated_alfalfa_salinity",
                  ["(A) fao29_leaching | (B) ec_to_nacl, crop_water_demand_alfalfa, irrigation_schedule, twin, twin_solution_to_ece | (C) fitted_leaching_{fresh,dry}",
                   "maas_hoffman_alfalfa", "sa187_benefit"], levels, runs, [h],
                  plot_spec("12 Inoculated alfalfa", "irrigation water EC (dS/m)", "inoculated relative yield", series), plain, cav,
                  fitted_leaching=fit, validation_mock_ratio=validation, observed=obs, mock_yields_digitised=mock_yields,
                  digitisation_error=err, headline_other_models=alt)


EXPERIMENTS = [exp01_brackish_irrigation, exp02_drip_vs_flood, exp03_mulch_shade, exp04_biochar, exp05_inoculation,
               exp06_carbon_fixation_bound, exp07_h2_leakage, exp08_solar_park, exp09_cloud_seeding, exp10_compaction,
               exp11_warming, exp12_inoculated_alfalfa_salinity]


def main():
    rows = []
    for i, fn in enumerate(EXPERIMENTS, 1):
        doc = fn()
        for hd in doc["headline"]:
            rows.append({"n": i, "id": doc["id"], "question": doc["question"], "headline": hd, "refused_levels": len(doc["refusals"]),
                         "deduced_nulls": len(doc["deduced_nulls"])})
        print(f"{i:2d} {doc['id']}")
    dump(OUT / "summary.json", {"note": "One row per experiment. Provenance 'twin' = process model; 'external:<id>' = another published model; 'deduced' = ontology.",
                                "registry_sha256": constants_sha256(), "experiments": rows})


if __name__ == "__main__":
    main()
