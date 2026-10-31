"""Composition of the twin with external models ("chained what-ifs").

A chain is  action -> upstream models -> twin -> downstream models.  Upstream models translate a management
action into twin inputs (forcing series, prescribed soil parameters, state, atmospheric gas concentration,
inoculation); downstream models translate twin outputs into an endpoint. Every model is a pure function over typed
ports whose literature constants are read by id from data/external/composition_constants.json.

Guarantees (see docs/causal-ontological-theory.md section 7 and the v0.7 design spec):

* typed ports: a value is passed from a producer to a consumer only if name, ontology term and unit are equal
  (otherwise CompositionError);
* validity domains are admissibility conditions of the external model: an input outside its closed interval gives a
  refusal object (a result, not an exception);
* the twin's own admissibility checks (theory.check_admissible via interventions.intervene) apply unchanged;
  a twin refusal is returned as a refusal object too;
* deductive gates: a downstream model with ``gate=C`` returns 0 with provenance 'deduced' (and does not call ``fn``)
  when the reasoner-entailed share of the class C in the chosen community is 0;
* the composed graph must be acyclic (no feedback from downstream to upstream) else CompositionError;
* every number is wrapped {value, low, high, unit, provenance, constants}; low/high are the envelope of a
  one-at-a-time rerun of the whole chain with each uncertain constant at its low and its high value (no sampling).

Nothing here is a prediction of the twin beyond its own outputs: provenance 'external:<id>' marks numbers that come
from another model.
"""
from __future__ import annotations

import hashlib
import json
from copy import deepcopy
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Literal

from .interventions import fingerprint, intervene
from .process import Parameters, State

REGISTRY_PATH = Path(__file__).resolve().parents[3] / "data" / "external" / "composition_constants.json"
DATA = Path(__file__).resolve().parents[3] / "data" / "v03"


class CompositionError(ValueError):
    """A structural error in a chain: port mismatch, missing input, duplicate producer, cycle, bad model output."""


@dataclass(frozen=True)
class Port:
    """A typed variable. Compatibility is equality of (term, unit); the name identifies the wire."""
    name: str
    term: str       # ontology CURIE (ChEBI/UO/ENVO/PATO/SIO) or rak: term
    unit: str       # UO CURIE or plain SI string, compared exactly


@dataclass(frozen=True)
class ExternalModel:
    id: str
    name: str
    kind: Literal["upstream", "downstream"]
    inputs: tuple
    outputs: tuple
    validity: dict           # input port name -> closed interval (lo, hi); outside -> refusal
    constants: tuple         # ids into composition_constants.json
    fn: Callable[[dict, dict], dict]   # (inputs by port name, constants by id) -> outputs by port name
    gate: str | None = None  # ontology class whose reasoner-entailed share must be > 0 (downstream only)

    def __post_init__(self):
        if self.kind not in ("upstream", "downstream"):
            raise CompositionError(f"{self.id}: kind must be 'upstream' or 'downstream'")
        names = {p.name for p in self.inputs}
        if set(self.validity) - names:
            raise CompositionError(f"{self.id}: validity given for unknown input ports {sorted(set(self.validity) - names)}")
        if self.gate is not None and self.kind != "downstream":
            raise CompositionError(f"{self.id}: only downstream models can carry a gate")


# ---------------------------------------------------------------- constants registry

_REG_CACHE: dict = {}


def load_constants(path: Path | None = None) -> dict:
    """Registry entries by id. Validates keys and low <= value <= high."""
    path = Path(path or REGISTRY_PATH)
    key = str(path)
    if key not in _REG_CACHE:
        doc = json.loads(path.read_text())
        out = {}
        for e in doc["constants"]:
            missing = {"id", "value", "low", "high", "unit", "meaning", "source", "verified"} - set(e)
            if missing:
                raise CompositionError(f"registry entry {e.get('id')!r} lacks {sorted(missing)}")
            if e["id"] in out:
                raise CompositionError(f"duplicate registry id {e['id']!r}")
            if not e["low"] <= e["value"] <= e["high"]:
                raise CompositionError(f"registry entry {e['id']!r}: need low <= value <= high")
            out[e["id"]] = e
        _REG_CACHE[key] = out
    return _REG_CACHE[key]


def constants_sha256(path: Path | None = None) -> str:
    return hashlib.sha256(Path(path or REGISTRY_PATH).read_bytes()).hexdigest()


def constant_values(ids, scenario: str = "central", overrides: dict | None = None, registry: dict | None = None) -> dict:
    """{id: float} for scenario in central|low|high, with optional per-id overrides."""
    reg = registry or load_constants()
    key = {"central": "value", "low": "low", "high": "high"}[scenario]
    out = {}
    for i in ids:
        if i not in reg:
            raise CompositionError(f"unknown constant id {i!r}")
        out[i] = (overrides or {}).get(i, reg[i][key])
    return out


# ---------------------------------------------------------------- twin ports

_T = lambda n, term, unit: Port(n, term, unit)
# twin inputs: port name == intervention key; binding says where the value goes in the intervention object
TWIN_INPUTS = {
    "water_mm": (_T("water_mm", "rak:AddedWater", "mm/day"), ("forcing", "water_mm")),
    "nacl_mol_l": (_T("nacl_mol_l", "CHEBI:26710", "mol/L"), ("forcing", "nacl_mol_l")),
    "evaporation_mm": (_T("evaporation_mm", "rak:PrescribedEvaporation", "mm/day"), ("forcing", "evaporation_mm")),
    "temp_c": (_T("temp_c", "rak:SoilTemperature", "degC"), ("forcing", "temp_c")),
    "h2_mol_m3": (_T("h2_mol_m3", "CHEBI:18276", "mol/m3"), ("forcing", "h2_mol_m3")),
    "co_mol_m3": (_T("co_mol_m3", "CHEBI:17245", "mol/m3"), ("forcing", "co_mol_m3")),
    "porosity": (_T("porosity", "rak:SoilPorosity", "-"), ("parameters", "porosity")),
    "psi_entry_mpa": (_T("psi_entry_mpa", "rak:AirEntryPotential", "MPa"), ("parameters", "psi_entry_mpa")),
    "retention_b": (_T("retention_b", "rak:CampbellRetentionExponent", "-"), ("parameters", "retention_b")),
    "state_water_mm": (_T("state_water_mm", "rak:SoilWaterStorage", "mm"), ("state", "water_mm")),
    "state_nacl_mol_m2": (_T("state_nacl_mol_m2", "CHEBI:26710", "mol/m2"), ("state", "nacl_mol_m2")),
    "inoculation": (_T("inoculation", "rak:Inoculum", "-"), ("inoculation", None)),
}
TWIN_INPUT_PORTS = {k: v[0] for k, v in TWIN_INPUTS.items()}
# twin outputs (intervened arm); also available to downstream models
TWIN_OUTPUT_PORTS = {
    "h2_uptake_mol_m2": _T("h2_uptake_mol_m2", "CHEBI:18276", "mol/m2"),
    "co_uptake_mol_m2": _T("co_uptake_mol_m2", "CHEBI:17245", "mol/m2"),
    "final_water_mm": _T("final_water_mm", "rak:SoilWaterStorage", "mm"),
    "final_salt_mol_m2": _T("final_salt_mol_m2", "CHEBI:26710", "mol/m2"),
    "mean_nacl_mol_l": _T("mean_nacl_mol_l", "CHEBI:26710", "mol/L"),
    "active_days": _T("active_days", "rak:ActiveDays", "day"),
}
# baseline values of the twin request that upstream models may read (name 'baseline_<twin input>')
BASELINE_PORTS = {f"baseline_{k}": Port(f"baseline_{k}", v.term, v.unit) for k, v in TWIN_INPUT_PORTS.items() if k != "inoculation"}


def twin_request(forcing: list, *, initial: State | None = None, params: Parameters | None = None, community: dict | None = None,
                 coupling: str | None = None, data: Path | None = None, intervention: dict | None = None) -> dict:
    """The twin side of a chain. Defaults reproduce the v0.6 analysis baseline (reproducibility/v06/analyse.py):
    Parameters(k_h2_s=.001, k_co_s=.001), State(1.8 mm, no salt). `intervention` is a base intervention (explicit
    do()s) merged under the values bound from upstream models."""
    iv = deepcopy(intervention or {})
    if community is not None:
        iv["community"] = community
    if coupling is not None:
        iv["coupling"] = coupling
    return {"initial": initial or State(1.8), "forcing": forcing, "params": params or Parameters(k_h2_s=.001, k_co_s=.001),
            "data": Path(data or DATA), "intervention": iv}


def _baseline_pool(req: dict) -> dict:
    forcing, p = req["forcing"], req["params"]
    vals = {"water_mm": [f.get("water_mm", 0.) for f in forcing], "nacl_mol_l": [f.get("nacl_mol_l", 0.) for f in forcing],
            "evaporation_mm": [f["evaporation_mm"] for f in forcing], "temp_c": [f["temp_c"] for f in forcing],
            "h2_mol_m3": p.h2_mol_m3, "co_mol_m3": p.co_mol_m3, "porosity": p.porosity,
            "psi_entry_mpa": p.psi_entry_mpa, "retention_b": p.retention_b,
            "state_water_mm": req["initial"].water_mm, "state_nacl_mol_m2": req["initial"].salt_mol_m2}
    return {f"baseline_{k}": v for k, v in vals.items()}


# ---------------------------------------------------------------- validity, gates, graph

def check_validity(model: ExternalModel, inputs: dict) -> dict | None:
    """Refusal object if an input lies outside its closed validity interval (or is missing because the twin abstained)."""
    for port, (lo, hi) in model.validity.items():
        v = inputs[port]
        for x in (v if isinstance(v, (list, tuple)) else [v]):
            if x is None:
                return {"refused": True, "model": model.id, "port": port, "value": None, "interval": [lo, hi],
                        "reason": "upstream value unavailable (the twin abstains outside its ideal-dilute NaCl domain)"}
            if not lo <= x <= hi:
                return {"refused": True, "model": model.id, "port": port, "value": x, "interval": [lo, hi],
                        "reason": f"{port}={x} lies outside the validity domain [{lo}, {hi}] of {model.id}"}
    return None


def gate_derivation(model: ExternalModel, shares: dict) -> str:
    return (f"q_s({model.gate}) = 0 under the reasoner memberships of the chosen community, so the gate of {model.id} is closed and its "
            "output is 0 in every admissible regime; absence of entailment is not read as biological absence.")


def check_dag(upstream: list, downstream: list, bound: set | None = None) -> None:
    """Composed graph upstream -> twin -> downstream must be a DAG. Edge A -> B when an output port name of A is an input
    port name of B; twin consumes upstream outputs named like a twin input and produces TWIN_OUTPUT_PORTS.
    Raises CompositionError on a cycle (feedback from downstream to upstream is out of scope for v0.7)."""
    nodes = {m.id: m for m in upstream + downstream}
    if len(nodes) != len(upstream) + len(downstream):
        raise CompositionError("duplicate model id in chain")
    edges = {n: set() for n in nodes} | {"twin": set()}
    producers: dict = {}
    for m in upstream + downstream:
        for o in m.outputs:
            if o.name in producers or o.name in TWIN_OUTPUT_PORTS:
                raise CompositionError(f"port {o.name!r} has more than one producer ({producers.get(o.name, 'twin')} and {m.id})")
            producers[o.name] = m.id
    for o in TWIN_OUTPUT_PORTS:
        producers[o] = "twin"
    for m in upstream + downstream:
        for i in m.inputs:
            src = producers.get(i.name)
            if src and src != m.id:
                edges[src].add(m.id)
            elif src == m.id:
                raise CompositionError(f"{m.id} consumes its own output {i.name!r}")
    for m in upstream:
        if any(o.name in TWIN_INPUT_PORTS for o in m.outputs):
            edges[m.id].add("twin")
    indeg = {n: 0 for n in edges}
    for a, bs in edges.items():
        for b in bs:
            indeg[b] += 1
    queue, done = [n for n, d in indeg.items() if d == 0], 0
    while queue:
        n = queue.pop()
        done += 1
        for b in edges[n]:
            indeg[b] -= 1
            if indeg[b] == 0:
                queue.append(b)
    if done != len(edges):
        raise CompositionError("the composed graph has a cycle (feedback from a downstream output to an upstream input is "
                               "out of scope): " + ", ".join(sorted(n for n, d in indeg.items() if d > 0)))


# ---------------------------------------------------------------- one evaluation of the chain

def _entry(value, port, provenance, consts, source):
    return {"value": value, "port": port, "provenance": provenance, "constants": frozenset(consts), "source": source}


def _gather(model: ExternalModel, pool: dict) -> tuple:
    inputs, lineage = {}, set()
    for port in model.inputs:
        e = pool.get(port.name)
        if e is None:
            raise CompositionError(f"{model.id}: input port {port.name!r} is not produced by the action, the twin baseline, "
                                   "an upstream model or the twin")
        if e["port"] is None:                       # action value: the first consumer fixes its type
            e["port"] = port
        elif (e["port"].term, e["port"].unit) != (port.term, port.unit):
            raise CompositionError(f"{model.id}: port {port.name!r} expects ({port.term}, {port.unit}) but "
                                   f"{e['source']} provides ({e['port'].term}, {e['port'].unit})")
        inputs[port.name] = e["value"]
        lineage |= e["constants"]
    return inputs, lineage


def _twin_outputs(res: dict, arm: str) -> dict:
    rows = res[arm]["trajectory"]
    gas = res["gas_uptake"]
    conc = [r["nacl_mol_l"] for r in rows if r["nacl_mol_l"] is not None]
    return {"h2_uptake_mol_m2": gas["h2_mol_m2"]["baseline" if arm == "baseline" else "intervened"],
            "co_uptake_mol_m2": gas["co_mol_m2"]["baseline" if arm == "baseline" else "intervened"],
            "final_water_mm": res[arm]["summary"]["final_water_mm"],
            "final_salt_mol_m2": res[arm]["summary"]["final_salt_mol_m2"],
            "mean_nacl_mol_l": (sum(conc) / len(conc)) if conc else None,
            "active_days": sum(1 for r in rows if r["gas"] is not None and r["gas"]["activation"] > 0)}


def _build_intervention(req: dict, bound: dict) -> dict:
    iv = deepcopy(req["intervention"])
    n = len(req["forcing"])
    for name, (value, _port) in bound.items():
        kind, key = TWIN_INPUTS[name][1]
        if kind == "forcing":
            series = list(value) if isinstance(value, (list, tuple)) else [value] * n
            if len(series) != n:
                raise CompositionError(f"{name}: series of length {len(series)} for a {n}-day baseline")
            for d, v in enumerate(series):
                iv.setdefault("forcing", {}).setdefault(str(d), {})[key] = v
        elif kind == "parameters":
            iv.setdefault("parameters", {})[key] = value
        elif kind == "state":
            iv.setdefault("state", {}).setdefault("0", {})[key] = value
        else:
            iv["inoculation"] = dict(value)
    return iv


def _evaluate(action: dict, upstream: list, req: dict, downstream: list, consts: dict, cache: dict) -> dict:
    """One pass of the chain with fixed constants. Returns {'refused', 'records', 'twin'}."""
    pool = {k: _entry(v, None, "action", (), "action") for k, v in action.items()}
    for k, v in _baseline_pool(req).items():
        if k in pool:
            raise CompositionError(f"{k!r} is reserved for the twin baseline")
        pool[k] = _entry(v, BASELINE_PORTS[k], "twin-baseline", (), "twin-baseline")
    records: dict = {}

    def run(model):
        inputs, lineage = _gather(model, pool)
        used = {c: consts[c] for c in model.constants}
        return inputs, lineage | set(model.constants), used

    for m in upstream:
        inputs, lineage, used = run(m)
        ref = check_validity(m, inputs)
        if ref:
            return {"refused": {**ref, "stage": "upstream"}, "records": records, "twin": None}
        out = m.fn(deepcopy(inputs), used)
        if set(out) != {o.name for o in m.outputs}:
            raise CompositionError(f"{m.id}: fn returned {sorted(out)} but declares {[o.name for o in m.outputs]}")
        for o in m.outputs:
            pool[o.name] = _entry(out[o.name], o, f"external:{m.id}", lineage, m.id)
            records[f"{m.id}:{o.name}"] = {"stage": m.id, "port": o.name, "value": out[o.name], "unit": o.unit,
                                          "provenance": f"external:{m.id}", "constants": sorted(lineage)}
    bound = {}
    for name, e in pool.items():
        if name in TWIN_INPUTS and e["provenance"] != "twin-baseline":
            tp = TWIN_INPUT_PORTS[name]
            if e["port"] is None:
                e["port"] = tp
            if (e["port"].term, e["port"].unit) != (tp.term, tp.unit):
                raise CompositionError(f"twin input {name!r} expects ({tp.term}, {tp.unit}) but {e['source']} provides "
                                       f"({e['port'].term}, {e['port'].unit})")
            bound[name] = (e["value"], e)
    twin_lineage = set().union(*[e["constants"] for _, e in bound.values()]) if bound else set()
    iv = _build_intervention(req, bound)
    key = fingerprint(iv)
    if key not in cache:
        try:
            cache[key] = intervene(req["initial"], req["forcing"], req["params"], iv, req["data"])
        except ValueError as e:
            cache[key] = {"refused": str(e)}
    res = cache[key]
    if "refused" in res and len(res) == 1:
        return {"refused": {"refused": True, "model": "twin", "port": None, "value": None, "interval": None,
                            "reason": res["refused"], "stage": "twin"}, "records": records, "twin": None}
    out_i, out_b = _twin_outputs(res, "intervened"), _twin_outputs(res, "baseline")
    for name, port in TWIN_OUTPUT_PORTS.items():
        pool[name] = _entry(out_i[name], port, "twin", twin_lineage, "twin")
        records[f"twin:{name}"] = {"stage": "twin", "port": name, "value": out_i[name], "unit": port.unit,
                                  "provenance": "twin", "constants": sorted(twin_lineage)}
    if res["inoculation"]:
        shares = res["inoculation"]["shares"]
    else:
        shares = res["community"]["shares"]
    for m in downstream:
        inputs, lineage, used = run(m)
        if m.gate is not None:
            if shares is None:
                raise CompositionError(f"{m.id} has gate {m.gate} but the twin request has no community to evaluate it in")
            if shares[m.gate] == 0:
                der = gate_derivation(m, shares)
                for o in m.outputs:
                    pool[o.name] = _entry(0.0, o, "deduced", (), m.id)
                    records[f"{m.id}:{o.name}"] = {"stage": m.id, "port": o.name, "value": 0.0, "unit": o.unit,
                                                  "provenance": "deduced", "constants": [], "derivation": der}
                continue
        ref = check_validity(m, inputs)
        if ref:
            return {"refused": {**ref, "stage": "downstream"}, "records": records, "twin": res}
        out = m.fn(deepcopy(inputs), used)
        if set(out) != {o.name for o in m.outputs}:
            raise CompositionError(f"{m.id}: fn returned {sorted(out)} but declares {[o.name for o in m.outputs]}")
        for o in m.outputs:
            pool[o.name] = _entry(out[o.name], o, f"external:{m.id}", lineage, m.id)
            records[f"{m.id}:{o.name}"] = {"stage": m.id, "port": o.name, "value": out[o.name], "unit": o.unit,
                                          "provenance": f"external:{m.id}", "constants": sorted(lineage)}
    return {"refused": None, "records": records, "twin": res, "baseline_outputs": out_b}


# ---------------------------------------------------------------- wrapping and envelope

def _span(values: list):
    """(low, high) over runs; elementwise for series; None if any run abstains or the value is not numeric."""
    if any(v is None for v in values) or any(isinstance(v, (dict, str, bool)) for v in values):
        return None, None
    if isinstance(values[0], (list, tuple)):
        return [min(c) for c in zip(*values)], [max(c) for c in zip(*values)]
    return min(values), max(values)


def wrap(record: dict, runs: list) -> dict:
    """{value, low, high, unit, provenance, constants} (+derivation) from the central record and all envelope runs."""
    lo, hi = _span([record["value"]] + runs)
    out = {"value": record["value"], "low": lo, "high": hi, "unit": record["unit"], "provenance": record["provenance"],
           "constants": record["constants"]}
    if "derivation" in record:
        out["derivation"] = record["derivation"]
    return out


def run_chain(action: dict, upstream: list, twin_req: dict, downstream: list, scenario: str = "central", *,
              constants: dict | None = None, envelope: bool = True, keep_twin_result: bool = False,
              registry: dict | None = None) -> dict:
    """Run action -> upstream -> twin -> downstream.

    action          {port name: value}; ports are typed by the first model that consumes them.
    upstream        list of ExternalModel(kind='upstream'), in execution order.
    twin_req        from twin_request(...).
    downstream      list of ExternalModel(kind='downstream').
    scenario        'central' | 'low' | 'high': every constant at its central / low / high value.
    constants       optional {id: value} overrides of registry values (reported under 'constants_used').
    envelope        with scenario 'central', rerun the chain one constant at a time at low and high and report
                    the min/max as low/high of every number (no sampling). Reruns that are refused are listed.

    Returns {'refused': None | refusal object, 'scenario', 'chain', 'upstream', 'twin', 'downstream', 'outputs',
    'deduced_nulls', 'envelope_refusals', 'constants_used', 'registry_sha256'}. Raises CompositionError for structural
    errors (port mismatch, missing input, duplicate producer, cycle).
    """
    if scenario not in ("central", "low", "high"):
        raise CompositionError("scenario must be central, low or high")
    if any(m.kind != "upstream" for m in upstream) or any(m.kind != "downstream" for m in downstream):
        raise CompositionError("upstream/downstream lists contain a model of the wrong kind")
    check_dag(upstream, downstream)
    reg = registry or load_constants()
    ids = sorted({c for m in upstream + downstream for c in m.constants})
    consts = constant_values(ids, scenario, constants, reg)
    cache: dict = {}
    central = _evaluate(action, upstream, twin_req, downstream, consts, cache)
    runs: dict = {}
    env_ref = []
    if envelope and scenario == "central" and central["refused"] is None:
        for cid in ids:
            if reg[cid]["low"] == reg[cid]["high"]:
                continue
            for end in ("low", "high"):
                c2 = {**consts, cid: reg[cid][end]}
                r = _evaluate(action, upstream, twin_req, downstream, c2, cache)
                if r["refused"]:
                    env_ref.append({"constant": cid, "end": end, "refusal": r["refused"]})
                    continue
                for k, rec in r["records"].items():
                    runs.setdefault(k, []).append(rec["value"])
    wrapped = {k: wrap(rec, runs.get(k, [])) for k, rec in central["records"].items()}
    by_stage: dict = {}
    for k, w in wrapped.items():
        stage, port = k.split(":", 1)
        by_stage.setdefault(stage, {})[port] = w
    twin = central["twin"]
    out = {"refused": central["refused"], "scenario": scenario,
           "chain": [{"id": m.id, "kind": m.kind, "name": m.name} for m in upstream] + [{"id": "twin", "kind": "twin", "name": "EQ-BioDT v0.7 process model"}]
                    + [{"id": m.id, "kind": m.kind, "name": m.name} for m in downstream],
           "upstream": {m.id: by_stage.get(m.id, {}) for m in upstream},
           "twin": None, "downstream": {m.id: by_stage.get(m.id, {}) for m in downstream},
           "outputs": {k: w for k, w in wrapped.items() if k.split(":", 1)[0] == "twin" or k.split(":", 1)[0] in {m.id for m in downstream}},
           "deduced_nulls": [{"model": w_k.split(":", 1)[0], "port": w_k.split(":", 1)[1], "derivation": w["derivation"]}
                             for w_k, w in wrapped.items() if w["provenance"] == "deduced"],
           "envelope_refusals": env_ref,
           "constants_used": {i: {"value": consts[i], "low": reg[i]["low"], "high": reg[i]["high"], "unit": reg[i]["unit"],
                                  "verified": reg[i]["verified"], "status": reg[i].get("status", "source_not_checked")} for i in ids},
           "registry_sha256": constants_sha256()}
    if twin is not None:
        out["twin"] = {"outputs": by_stage.get("twin", {}), "baseline_outputs": central.get("baseline_outputs"),
                       "intervention": twin["intervention"], "trace": twin["trace"],
                       "checks_passed": all(twin["checks"].values()), "twin_deduced": twin["deduced"],
                       "community": twin["community"], "inoculation": twin["inoculation"],
                       "gas_uptake": twin["gas_uptake"]}
        if keep_twin_result:
            out["twin"]["result"] = twin
    return out
