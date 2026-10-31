"""Stratified causal-ontological graph of the twin, as data (one time slice).

See docs/causal-ontological-theory.md. Zones: ontological (K), evidence (E), deduced
(D_K(E), q, gamma), exogenous_manipulable / exogenous_fixed (V_ex), constitutive
(V_const, definitions) and causal (V_caus, mechanisms). Temporal edges (W_t -> W_{t+1},
N_t -> N_{t+1}) carry lag=1 and are excluded from the within-slice acyclicity check.
"""
from __future__ import annotations

from functools import lru_cache
from math import isfinite

ZONES = ("ontological", "evidence", "deduced", "exogenous_manipulable",
         "exogenous_fixed", "constitutive", "causal")


def _n(id, label, zone, unit, equation, parents=(), manipulable=False, action=None, lag=()):
    return {"id": id, "label": label, "zone": zone, "unit": unit, "equation": equation,
            "parents": list(parents), "lag_parents": list(lag),
            "manipulable": manipulable, "physical_action": action}


def _nodes() -> list[dict]:
    ex = "exogenous_manipulable"
    fx = "exogenous_fixed"
    return [
        _n("K", "OWL 2 TBox (KEGG step classes, carrier classes)", "ontological", "-", "K"),
        _n("E", "Evidence: MAG KOs, marker guilds, mapped weights w_sg", "evidence", "-", "E"),
        _n("D_K", "Deduced class memberships g in C", "deduced", "-", "D_K(E) = {g : K, E |= g in C}", ["K", "E"]),
        _n("q_h2", "H2 carrier share q_s(H2CarrierMAG)", "deduced", "-",
           "q_s(C) = sum_{g in C} w_sg / sum_g w_sg", ["D_K", "E"]),
        _n("q_co", "CO carrier share q_s(COCarrierMAG)", "deduced", "-",
           "q_s(C) = sum_{g in C} w_sg / sum_g w_sg", ["D_K", "E"]),
        _n("gamma_h2", "Community factor for H2", "deduced", "-",
           "gamma = 1 | 1[q>0] | q/median(q)", ["q_h2"]),
        _n("gamma_co", "Community factor for CO", "deduced", "-",
           "gamma = 1 | 1[q>0] | q/median(q)", ["q_co"]),
        _n("I", "Added water", ex, "mm day-1", "I_t", manipulable=True,
           action="apply or withhold water (irrigation, rain event)"),
        _n("c", "NaCl concentration of added water", ex, "mol L-1", "c_t", manipulable=True,
           action="change the salinity of the applied water"),
        _n("Evap", "Prescribed evaporative demand", ex, "mm day-1", "E_t", manipulable=True,
           action="shade or mulch the surface"),
        _n("T", "Soil temperature", ex, "C", "T_t", manipulable=True,
           action="shade or warm the soil surface"),
        _n("phi", "Porosity", ex, "-", "phi", manipulable=True,
           action="compact or loosen the soil (measured post-disturbance value)"),
        _n("psi_e", "Air-entry potential", ex, "MPa", "psi_e", manipulable=True,
           action="change the soil texture class (measured retention)"),
        _n("b", "Campbell retention exponent", ex, "-", "b", manipulable=True,
           action="change the soil texture class (measured retention)"),
        _n("C_h2", "Atmospheric H2 concentration (boundary condition)", ex, "mol m-3", "C_h2,t (default C_a,h2)", ["C_a"],
           manipulable=True, action="raise or lower the ambient H2 concentration (e.g. a hydrogen-economy leakage scenario)"),
        _n("C_co", "Atmospheric CO concentration (boundary condition)", ex, "mol m-3", "C_co,t (default C_a,co)", ["C_a"],
           manipulable=True, action="raise or lower the ambient CO concentration"),
        _n("L", "Layer depth", fx, "m", "L"),
        _n("C_a", "Default atmospheric gas concentrations", fx, "mol m-3", "C_a,g"),
        _n("D0", "Free-air diffusivities", fx, "m2 s-1", "D_0,g"),
        _n("k_ref", "Reference reaction coefficients", fx, "s-1", "k_g,ref"),
        _n("Q10", "Temperature response", fx, "-", "Q10, T_ref"),
        _n("act_par", "Activation thresholds and shape", fx, "MPa", "psi_dry, psi_wet, shape"),
        _n("drain_par", "Drainage and field-capacity constants", fx, "-", "f_fc, f_d"),
        _n("osm_par", "Osmotic coefficient, gas constant", fx, "-", "varphi, R"),
        _n("O2", "Oxygen supply limit", fx, "mol m-2 s-1", "S_O2"),
        _n("capture", "Energy capture fraction", fx, "-", "eta"),
        _n("W", "Water storage", "causal", "mm",
           "W_{t+1} = clip(W_t + I_t - drain - evap, 0, W_max)", ["I", "Evap", "L", "phi", "W_max", "drain_par"],
           manipulable=True, action="wet or dry the soil to a set-point (state intervention)", lag=["W"]),
        _n("N", "Salt inventory", "causal", "mol m-2",
           "N_{t+1} = N_t + I_t c_t - leached salt", ["I", "c", "W", "drain_par"],
           manipulable=True, action="add or leach salt to a set-point (state intervention)", lag=["N"]),
        _n("W_max", "Pore storage", "constitutive", "mm", "W_max = 1000 L phi", ["L", "phi"]),
        _n("theta", "Volumetric water content", "constitutive", "m3 m-3", "theta = W/(1000 L)", ["W", "L"]),
        _n("psi_m", "Matric potential", "constitutive", "MPa", "psi_m = psi_e (theta/phi)^(-b)",
           ["theta", "phi", "psi_e", "b"]),
        _n("psi_o", "Osmotic potential", "constitutive", "MPa", "psi_o = -2 varphi (N/W) R T 1e-3",
           ["N", "W", "T", "osm_par"]),
        _n("psi", "Total water potential", "constitutive", "MPa", "psi = psi_m + psi_o", ["psi_m", "psi_o"]),
        _n("D_h2", "Effective H2 diffusivity", "constitutive", "m2 s-1",
           "D_g = D_0,g (phi - theta)^(10/3) / phi^2", ["D0", "phi", "theta"]),
        _n("D_co", "Effective CO diffusivity", "constitutive", "m2 s-1",
           "D_g = D_0,g (phi - theta)^(10/3) / phi^2", ["D0", "phi", "theta"]),
        _n("A", "Activation", "causal", "-", "A_t = A(psi_t)", ["psi", "act_par"]),
        _n("k_h2", "H2 reaction coefficient", "causal", "s-1",
           "k_g,t = k_g,ref gamma_g A_t Q10^((T_t - T_ref)/10)", ["k_ref", "gamma_h2", "A", "Q10", "T"]),
        _n("k_co", "CO reaction coefficient", "causal", "s-1",
           "k_g,t = k_g,ref gamma_g A_t Q10^((T_t - T_ref)/10)", ["k_ref", "gamma_co", "A", "Q10", "T"]),
        _n("J_h2", "H2 uptake", "causal", "mol m-2 s-1",
           "J_g = C_a sqrt(D_g k_g) tanh(L sqrt(k_g/D_g)), proportional oxygen limit (shared demand couples both gases)", ["C_h2", "D_h2", "k_h2", "D_co", "k_co", "L", "O2"]),
        _n("J_co", "CO uptake", "causal", "mol m-2 s-1",
           "J_g = C_a sqrt(D_g k_g) tanh(L sqrt(k_g/D_g)), proportional oxygen limit (shared demand couples both gases)", ["C_co", "D_co", "k_co", "D_h2", "k_h2", "L", "O2"]),
        _n("P", "Chemical power", "constitutive", "W m-2", "P = 202000 J_h2 + 243000 J_co", ["J_h2", "J_co"]),
        _n("P_cap", "Captured power", "constitutive", "W m-2", "P_cap = eta P", ["P", "capture"]),
    ]


@lru_cache(maxsize=1)
def _graph() -> dict:
    nodes = _nodes()
    kind = {"causal": "mechanism", "constitutive": "definition", "deduced": "deduction"}
    ids = {n["id"] for n in nodes}
    edges = []
    for n in nodes:
        for p in n["parents"]:
            edges.append({"source": p, "target": n["id"], "kind": kind.get(n["zone"], "definition"), "lag": 0})
        for p in n["lag_parents"]:
            edges.append({"source": p, "target": n["id"], "kind": "lag", "lag": 1})
    missing = [e for e in edges if e["source"] not in ids]
    if missing:
        raise ValueError(f"Graph edge from unknown node: {missing}")
    return {"nodes": nodes, "edges": edges, "zones": list(ZONES),
            "note": "One time-slice template; edges with kind 'lag' connect slice t to t+1."}


def graph() -> dict:
    import copy
    return copy.deepcopy(_graph())


def node(id: str) -> dict:
    for n in _graph()["nodes"]:
        if n["id"] == id:
            return n
    raise ValueError(f"Unknown node: {id}")


# Request keys -> node ids.
ALIASES = {"water_mm": "I", "nacl_mol_l": "c", "evaporation_mm": "Evap", "temp_c": "T",
           "porosity": "phi", "psi_entry_mpa": "psi_e", "retention_b": "b",
           "psi_mpa": "psi", "h2_mol_m3": "C_h2", "co_mol_m3": "C_co", "carrier_share": "q_h2", "membership": "D_K", "memberships": "D_K",
           "h2_uptake_mol_m2": "J_h2", "co_uptake_mol_m2": "J_co", "activation": "A"}


def resolve(target: str) -> str:
    target = ALIASES.get(target, target)
    node(target)
    return target


def admissible_targets() -> list[str]:
    return [n["id"] for n in _graph()["nodes"] if n["manipulable"]]


def _handles(id: str) -> list[str]:
    return [h for h in sorted(_ancestors(id) | {id}) if node(h)["zone"] == "causal" and node(h)["manipulable"] and h != id]


def check_admissible(targets) -> list[str]:
    """Return node ids if every target is a manipulable input or physical state; raise otherwise."""
    out = []
    for t in targets:
        try:
            i = resolve(t)
        except ValueError:
            raise ValueError(f"{t!r} is not a node of the causal graph; admissible targets: {admissible_targets()}") from None
        n = node(i)
        if n["manipulable"]:
            out.append(i)
            continue
        z = n["zone"]
        if z == "constitutive":
            hint = {"theta": "intervene on water storage W instead, or use a theta state target which is translated",
                    "psi": "use a psi state target with an explicit version ('matric' or 'osmotic'), which is translated to W or N"}.get(
                i, "intervene on its constituents " + ", ".join(_handles(i)))
            raise ValueError(f"{i} is constitutive ({n['equation']}); {hint}")
        if z == "causal":
            raise ValueError(f"{i} is a causal outcome mechanism without physical handle; intervene on its inputs")
        if z == "deduced":
            raise ValueError(f"{i} is in the deduced zone: carrier share is deduced from evidence and ontology; "
                             "change evidence by an observation update, not an intervention")
        if z in ("ontological", "evidence"):
            raise ValueError(f"{i} is in the {z} zone; change it by a model revision or an observation update, not an intervention")
        raise ValueError(f"{i} is an exogenous fixed constant of this release (zone {z}); it has no physical handle in the intervention interface")
    return out


def _ancestors(id: str, include_lag: bool = False) -> set[str]:
    g = _graph()
    seen: set[str] = set()
    stack = [id]
    while stack:
        cur = stack.pop()
        for e in g["edges"]:
            if e["target"] == cur and e["source"] not in seen and (include_lag or e["kind"] != "lag"):
                seen.add(e["source"])
                stack.append(e["source"])
    return seen


def slice_is_acyclic() -> bool:
    g = _graph()
    indeg = {n["id"]: 0 for n in g["nodes"]}
    out: dict[str, list[str]] = {i: [] for i in indeg}
    for e in g["edges"]:
        if e["kind"] != "lag":
            indeg[e["target"]] += 1
            out[e["source"]].append(e["target"])
    queue = [i for i, d in indeg.items() if d == 0]
    done = 0
    while queue:
        cur = queue.pop()
        done += 1
        for t in out[cur]:
            indeg[t] -= 1
            if indeg[t] == 0:
                queue.append(t)
    return done == len(indeg)


def version_relevance(target: str, outcome: str) -> dict:
    """Is the effect of do(target) on outcome version-dependent?

    Versions are the different physical handles realising a constituted target. They matter
    for an outcome iff some handle reaches the outcome by a directed path that avoids target.
    """
    t, o = resolve(target), resolve(outcome)
    handles = _handles(t)
    if not handles:
        handles = [t]
    adj: dict[str, list[str]] = {}
    for e in _graph()["edges"]:
        if e["kind"] != "lag":
            adj.setdefault(e["source"], []).append(e["target"])
    paths: list[list[str]] = []

    def walk(cur, trail):
        if cur == o:
            paths.append(list(trail))
            return
        for nxt in adj.get(cur, []):
            if nxt != t and nxt not in trail:
                walk(nxt, trail + [nxt])

    for h in handles:
        if h != o:
            walk(h, [h])
        else:
            paths.append([h])
    unique = len(handles) == 1
    return {"target": t, "outcome": o, "handles": handles,
            "relevant": bool(paths) and not unique, "bypass_paths": paths,
            "unique_translation": unique,
            "reading": ("version-independent: every route from the handles to the outcome passes through the target"
                        if not paths else
                        "unique translation, versions coincide" if unique else
                        "version-dependent: handles reach the outcome without passing through the target")}


# ---- translation of constituted targets ----
R_GAS = 8.314462618


def _osmotic_coefficient(p, temp_c: float) -> float:
    """psi_o = -coef * (N/W); MPa per (mol m-2 / mm)."""
    return 2 * p.osmotic_coefficient * 1000 * R_GAS * (temp_c + 273.15) / 1e6


def total_potential(water_mm: float, salt: float, p, temp_c: float) -> float:
    theta = water_mm / (1000 * p.depth_m)
    psi_m = p.psi_entry_mpa * (theta / p.porosity) ** (-p.retention_b)
    return psi_m - _osmotic_coefficient(p, temp_c) * (salt / water_mm)


def translate_constituted(target: str, value: float, state, params, version: str | None = None,
                          temp_c: float | None = None) -> dict:
    """Translate do(theta=v) or do(psi=v, version) into a do on W or N."""
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not isfinite(value):
        raise ValueError("Constituted target value must be a finite number")
    t = ALIASES.get(target, target)
    p = params
    temp = p.reference_temp_c if temp_c is None else temp_c
    capacity = 1000 * p.depth_m * p.porosity
    if t == "theta":
        if not 0 <= value <= p.porosity:
            raise ValueError(f"Infeasible theta={value}: require 0 <= theta <= porosity={p.porosity} (theta > phi violates W <= W_max)")
        return {"target": "theta", "value": value, "handle": "W", "water_mm": 1000 * p.depth_m * value,
                "salt_mol_m2": state.salt_mol_m2, "version": "unique"}
    if t != "psi":
        raise ValueError(f"{target!r} is not a constituted target (theta, psi_mpa)")
    if version not in ("matric", "osmotic"):
        raise ValueError("psi is an ambiguous manipulation; specify version 'matric' (solve W) or 'osmotic' (solve N)")
    if value >= 0:
        raise ValueError(f"Infeasible psi={value}: total potential must be negative")
    coef = _osmotic_coefficient(p, temp)
    if version == "matric":
        n = state.salt_mol_m2
        hi = total_potential(capacity, n, p, temp)
        if value >= hi:
            raise ValueError(f"Infeasible psi={value} MPa: matric potential would reach or exceed the air-entry potential "
                             f"({p.psi_entry_mpa} MPa) at saturation; maximum attainable potential at N={n} is {hi}")
        lo, up = capacity * 1e-12, capacity
        for _ in range(200):               # psi(W; N) is strictly increasing in W
            mid = (lo + up) / 2
            if total_potential(mid, n, p, temp) < value:
                lo = mid
            else:
                up = mid
        water = (lo + up) / 2
        return {"target": "psi", "value": value, "handle": "W", "water_mm": water, "salt_mol_m2": n, "version": "matric"}
    w = state.water_mm
    if w <= 0:
        raise ValueError("Osmotic version needs positive water storage")
    theta = w / (1000 * p.depth_m)
    psi_m = p.psi_entry_mpa * (theta / p.porosity) ** (-p.retention_b)
    psi_o = value - psi_m
    if psi_o > 0:
        raise ValueError(f"Infeasible psi={value} MPa: above the matric potential {psi_m} at current W; salt cannot raise potential")
    return {"target": "psi", "value": value, "handle": "N", "water_mm": w, "salt_mol_m2": -psi_o / coef * w, "version": "osmotic"}
