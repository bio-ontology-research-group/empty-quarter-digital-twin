"""Browser-side driver for the public Twin page (copied to py/eq_biodt/live.py in the web bundle).

Pure standard library plus eq_biodt.v03.process and eq_biodt.v03.theory, the released model code.
It evaluates a *paired regime*: a baseline arm and an intervened arm that share the initial state,
the non-targeted inputs, the parameters and the community factors gamma. Selecting a community and a
coupling mode conditions the twin on evidence; it is not an intervention.

Entry points (JSON string in, JSON string out): ``run`` and ``admissible``.
"""
from __future__ import annotations

import json
from dataclasses import replace
from itertools import product

from eq_biodt.v03 import theory
from eq_biodt.v03.process import Parameters, State, simulate

HORIZON_DAYS = 14
INITIAL_WATER_MM = 1.8
EVAPORATION_MM = 1.0
TEMPERATURE_C = 25.0
MODES = ("reference", "gate", "proportional")
GAS_CLASS = {"h2": "H2CarrierMAG", "co": "COCarrierMAG"}
VERSIONS = ("none", "single_pulse", "four_pulses", "daily_pulses", "custom_day")
# Request keys -> graph node ids for the admissibility verdict of each do() target.
TARGET_NODE = {"water": "I", "nacl": "c", "evaporation": "Evap", "temperature": "T", "porosity": "phi"}


def base_parameters() -> Parameters:
    return Parameters(k_h2_s=.001, k_co_s=.001)


def parameter_members(base: Parameters) -> list[Parameters]:
    """The 24-member sensitivity grid of eq_biodt.v03.planning.parameter_ensemble (kept identical by a test)."""
    return [replace(base, retention_b=b, dry_threshold_mpa=dry, activation_shape=shape,
                    k_h2_s=base.k_h2_s * rate, k_co_s=base.k_co_s * rate)
            for b, dry, shape, rate in product((3., 5.), (-70., -100.), (.5, 2.), (.1, 1., 10.))]


def pulse_days(version: str, custom_day: int = 0, horizon: int = HORIZON_DAYS) -> list[int]:
    if version == "none":
        return []
    if version == "single_pulse":
        return [0]
    if version == "four_pulses":
        return sorted({int(i * horizon / 4) for i in range(4)})
    if version == "daily_pulses":
        return list(range(horizon))
    if version == "custom_day":
        if not isinstance(custom_day, int) or not 0 <= custom_day < horizon:
            raise ValueError(f"custom day must be an integer in [0, {horizon - 1}]")
        return [custom_day]
    raise ValueError(f"unknown schedule version {version!r}; choose from {list(VERSIONS)}")


def make_forcing(total_mm: float, version: str, custom_day: int, nacl: float, temp: float, evap: float) -> list[dict]:
    days = pulse_days(version, custom_day)
    if total_mm < 0:
        raise ValueError("added water must be non-negative")
    if total_mm > 0 and not days:
        raise ValueError("a positive water total needs a schedule version")
    each = total_mm / len(days) if days else 0.
    day_set = set(days)
    return [{"water_mm": each if d in day_set else 0., "nacl_mol_l": float(nacl),
             "evaporation_mm": float(evap), "temp_c": float(temp)} for d in range(HORIZON_DAYS)]


def gamma(mode: str, shares: dict | None, medians: dict | None) -> dict:
    """gamma_h2, gamma_co (same rule as eq_biodt.v03.community.gamma)."""
    if mode not in MODES:
        raise ValueError(f"Unknown coupling mode {mode!r}; choose from {list(MODES)}")
    if mode == "reference":
        return {"h2": 1., "co": 1.}
    if not shares:
        raise ValueError(f"coupling mode {mode!r} needs a community (a sample or a site visit)")
    out = {}
    for gas, cls in GAS_CLASS.items():
        q = float(shares[cls])
        if q < 0 or q != q:
            raise ValueError(f"Invalid share for {cls}")
        if mode == "gate":
            out[gas] = 1. if q > 0 else 0.
        else:
            med = float(medians[cls])
            if med <= 0:
                raise ValueError(f"Median share of {cls} is zero; proportional coupling undefined")
            out[gas] = q / med
    return out


def _summarise(result: dict) -> dict:
    rows = result["trajectory"]
    out = {"day": [], "water_mm": [], "salt_mol_m2": [], "nacl_mol_l": [], "h2_umol_m2_day": [],
           "co_umol_m2_day": [], "chemistry_status": [], "water_input_mm": []}
    for r in rows:
        out["day"].append(r["day"] + 1)
        out["water_mm"].append(r["water_mm"])
        out["salt_mol_m2"].append(r["salt_mol_m2"])
        out["nacl_mol_l"].append(r["nacl_mol_l"])
        out["water_input_mm"].append(r["water_input_mm"])
        out["chemistry_status"].append(r["chemistry_status"])
        g = r["gas"]
        out["h2_umol_m2_day"].append(None if g is None else g["h2_mol_m2_s"] * 86400e6)
        out["co_umol_m2_day"].append(None if g is None else g["co_mol_m2_s"] * 86400e6)
    valid = all(g is not None for g in out["h2_umol_m2_day"])
    out["h2_total_umol_m2"] = sum(out["h2_umol_m2_day"]) if valid else None
    out["co_total_umol_m2"] = sum(out["co_umol_m2_day"]) if valid else None
    out["abstained_days"] = [d for d, g in zip(out["day"], out["h2_umol_m2_day"]) if g is None]
    out["final_water_mm"] = rows[-1]["water_mm"]
    out["final_salt_mol_m2"] = rows[-1]["salt_mol_m2"]
    return out


def _state_override(spec: dict, temp_c: float, params: Parameters):
    """State -> State for one do() on the causal state (W, N) or on a constituted target (theta, psi)."""
    target = spec["target"]
    value = float(spec["value"])
    record: dict = {}

    def apply(state: State) -> State:
        water, salt = state.water_mm, state.salt_mol_m2
        if target == "psi_mpa":
            t = theory.translate_constituted("psi", value, state, params, spec.get("version"), temp_c)
            water, salt = t["water_mm"], t["salt_mol_m2"]
            record.update(handle=t["handle"], version=t["version"])
        elif target == "theta":
            t = theory.translate_constituted("theta", value, state, params)
            water = t["water_mm"]
            record.update(handle="W", version="unique")
        elif target == "water_mm":
            water = value
            record.update(handle="W", version="direct")
        elif target == "nacl_mol_m2":
            salt = value
            record.update(handle="N", version="direct")
        else:
            raise ValueError(f"unknown state target {target!r}")
        record.update(state_before={"water_mm": state.water_mm, "salt_mol_m2": state.salt_mol_m2},
                      state_after={"water_mm": water, "salt_mol_m2": salt})
        return State(water, salt)

    apply.record = record
    return apply


def _simulate_pair(params: Parameters, base_forcing: list[dict], int_forcing: list[dict],
                   int_params: Parameters, state_spec: dict | None) -> tuple[dict, dict, dict | None]:
    base = simulate(State(INITIAL_WATER_MM), base_forcing, params)
    override = None
    overrides = None
    if state_spec:
        day = int(state_spec.get("day", 0))
        if not 0 <= day < HORIZON_DAYS:
            raise ValueError(f"state intervention day must be in [0, {HORIZON_DAYS - 1}]")
        override = _state_override(state_spec, int_forcing[day]["temp_c"], int_params)
        overrides = {day: override}
    inter = simulate(State(INITIAL_WATER_MM), int_forcing, int_params, overrides)
    return base, inter, (override.record if override else None)


def verdict(target: str) -> dict:
    """Admissibility verdict of one graph node (or alias) as a do() target."""
    try:
        node_id = theory.resolve(target)
    except ValueError as e:
        return {"target": target, "node": None, "verdict": "refused", "zone": None, "message": str(e)}
    n = theory.node(node_id)
    try:
        theory.check_admissible([node_id])
        return {"target": target, "node": node_id, "verdict": "admissible", "zone": n["zone"],
                "message": n["physical_action"] or "", "label": n["label"]}
    except ValueError as e:
        return {"target": target, "node": node_id, "verdict": "refused", "zone": n["zone"],
                "message": str(e), "label": n["label"]}


def admissible(payload: str) -> str:
    targets = json.loads(payload)
    return json.dumps([verdict(t) for t in targets])


def _relevance(target: str) -> dict:
    out = {}
    for outcome in ("A", "D_h2", "J_h2", "J_co"):
        r = theory.version_relevance(target, outcome)
        out[outcome] = {k: r[k] for k in ("handles", "relevant", "bypass_paths", "reading")}
    return out


def run(payload: str) -> str:
    """Paired-regime evaluation. See the web page for the request fields."""
    req = json.loads(payload)
    mode = req.get("mode", "reference")
    g = gamma(mode, req.get("shares"), req.get("medians"))
    base = replace(base_parameters(), gamma_h2=g["h2"], gamma_co=g["co"])
    b = req.get("baseline", {})
    iv = req.get("intervention", {})
    b_temp, b_evap, b_nacl = (float(b.get("temp_c", TEMPERATURE_C)), float(b.get("evaporation_mm", EVAPORATION_MM)),
                              float(b.get("nacl_mol_l", 0.)))
    if "porosity" in b:
        base = replace(base, porosity=float(b["porosity"]))
    base_forcing = make_forcing(float(b.get("water_mm", 0.)), b.get("version", "none"), int(b.get("custom_day", 0)),
                                b_nacl, b_temp, b_evap)

    int_forcing = [dict(f) for f in base_forcing]
    int_params = base
    verdicts, targets_applied = [], []
    if "water" in iv:
        w = iv["water"]
        fresh = make_forcing(float(w["total_mm"]), w.get("version", "single_pulse"), int(w.get("custom_day", 0)),
                             0., b_temp, b_evap)
        for f, n in zip(int_forcing, fresh):
            f["water_mm"] = n["water_mm"]
        verdicts.append(verdict("water_mm"))
        targets_applied.append({"target": "water_mm", "value": w["total_mm"], "version": w.get("version", "single_pulse"),
                                "custom_day": w.get("custom_day")})
    for key, lever in (("nacl", "nacl_mol_l"), ("temperature", "temp_c"), ("evaporation", "evaporation_mm")):
        if key in iv:
            for f in int_forcing:
                f[lever] = float(iv[key])
            verdicts.append(verdict(lever))
            targets_applied.append({"target": lever, "value": iv[key]})
    if "porosity" in iv:
        int_params = replace(base, porosity=float(iv["porosity"]))
        verdicts.append(verdict("porosity"))
        targets_applied.append({"target": "porosity", "value": iv["porosity"]})

    state_spec = iv.get("state")
    relevance = None
    if state_spec:
        t = state_spec["target"]
        if t == "psi_mpa":
            relevance = _relevance("psi")
            if state_spec.get("version") not in ("matric", "osmotic"):
                raise ValueError("psi is an ambiguous manipulation; specify version 'matric' (solve W) or 'osmotic' "
                                 "(solve N). The graph shows why: a path W -> theta -> D_g -> J_g bypasses psi.")
        v = verdict({"psi_mpa": "psi", "theta": "theta", "water_mm": "W", "nacl_mol_m2": "N"}[t])
        if v["verdict"] == "refused" and v["zone"] == "constitutive":
            handle = "W or N, by the stated version" if t == "psi_mpa" else "W"
            v = {**v, "verdict": "translated",
                 "message": f"{v['node']} is constitutive and cannot be clamped directly; do({v['node']}) is translated through its definition to a do() on {handle}."}
        verdicts.append(v)
        targets_applied.append(dict(state_spec))
    if not targets_applied:
        raise ValueError("choose at least one do() target")

    b_res, i_res, state_record = _simulate_pair(base, base_forcing, int_forcing, int_params, state_spec)
    baseline, intervened = _summarise(b_res), _summarise(i_res)

    def diff(key):
        a, c = baseline[key], intervened[key]
        return None if a is None or c is None else c - a

    effect = {"h2_umol_m2": diff("h2_total_umol_m2"), "co_umol_m2": diff("co_total_umol_m2"),
              "h2_daily": [None if x is None or y is None else y - x
                           for x, y in zip(baseline["h2_umol_m2_day"], intervened["h2_umol_m2_day"])],
              "co_daily": [None if x is None or y is None else y - x
                           for x, y in zip(baseline["co_umol_m2_day"], intervened["co_umol_m2_day"])]}
    nulls = {}
    for gas, cls in GAS_CLASS.items():
        if g[gas] == 0.:
            ok = effect[f"{gas}_umol_m2"] == 0. and all(v == 0. for v in intervened[f"{gas}_umol_m2_day"] if v is not None)
            nulls[gas] = {"zero_by_deduction": True, "numerically_confirmed": bool(ok),
                          "derivation": (f"No mapped {cls} genome in this community, so gamma_{gas} = 0. "
                                         f"k_{gas} = k_ref * gamma * A * Q10^(..) = 0 in every regime, hence J_{gas} = 0 and "
                                         f"every interventional contrast on {gas.upper()} uptake is zero. Deduced quantities "
                                         "have no causal ancestors, so no do() on water, salt, temperature or pores can change them.")}

    members = None
    if req.get("members", True):
        members = []
        for i, m in enumerate(parameter_members(base)):
            m_int = replace(m, porosity=int_params.porosity)
            try:
                mb, mi, _ = _simulate_pair(m, base_forcing, int_forcing, m_int, state_spec)
                sb, si = _summarise(mb), _summarise(mi)
                members.append({"member": i, "retention_b": m.retention_b, "dry_threshold_mpa": m.dry_threshold_mpa,
                                "activation_shape": m.activation_shape, "rate_multiplier": m.k_h2_s / base_parameters().k_h2_s,
                                "h2_baseline": sb["h2_total_umol_m2"], "h2_intervened": si["h2_total_umol_m2"],
                                "co_baseline": sb["co_total_umol_m2"], "co_intervened": si["co_total_umol_m2"]})
            except ValueError as e:
                members.append({"member": i, "error": str(e)})

    psi_versions = None
    if state_spec and state_spec["target"] == "psi_mpa":
        psi_versions = {}
        for ver in ("matric", "osmotic"):
            spec = dict(state_spec, version=ver)
            try:
                _, ir, rec = _simulate_pair(base, base_forcing, int_forcing, int_params, spec)
                s = _summarise(ir)
                psi_versions[ver] = {"handle": rec["handle"], "state_after": rec["state_after"],
                                     "h2_total_umol_m2": s["h2_total_umol_m2"], "co_total_umol_m2": s["co_total_umol_m2"],
                                     "abstained_days": s["abstained_days"], "nacl_mol_l_after": (
                                         rec["state_after"]["salt_mol_m2"] / rec["state_after"]["water_mm"]
                                         if rec["state_after"]["water_mm"] else None)}
            except ValueError as e:
                psi_versions[ver] = {"error": str(e)}

    return json.dumps({
        "mode": mode, "gamma": g, "parameters": {"k_h2_s": base.k_h2_s, "k_co_s": base.k_co_s, "porosity_baseline": base.porosity,
                                                 "porosity_intervened": int_params.porosity, "horizon_days": HORIZON_DAYS,
                                                 "initial_water_mm": INITIAL_WATER_MM},
        "targets": targets_applied, "verdicts": verdicts, "baseline": baseline, "intervened": intervened,
        "effect": effect, "nulls": nulls, "state_record": state_record, "relevance": relevance,
        "psi_versions": psi_versions, "members": members, "model_version": __import__("eq_biodt.v03", fromlist=["VERSION"]).VERSION,
    }, allow_nan=False)
