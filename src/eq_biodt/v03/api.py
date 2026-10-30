"""Current public API. Historical heuristic interventions are explicitly retired."""
from pathlib import Path
from . import VERSION
import json
from functools import lru_cache

from fastapi import FastAPI, HTTPException, Query
from fastapi.responses import HTMLResponse, JSONResponse
import pandas as pd
from pydantic import BaseModel, ConfigDict, Field

from . import community as comm
from . import theory
from .contract import contract
from .interventions import intervene
from .evidence import CBB_STEPS
from .ontology import class_annotations, load_classification
from .planning import compare, experiment, parameter_ensemble
from .process import Parameters, State, simulate
from .resources import resource_ceiling, carbon_account

PACKAGE = Path(__file__).parent
DATA = PACKAGE/"data" if (PACKAGE/"data").exists() else PACKAGE.parents[2]/"data"/"v03"
app = FastAPI(title="Empty Quarter Digital Twin", version=VERSION,
              description="Metagenomic evidence and conditional, resource-conserving scenarios. See /v3/contract.")


class ScenarioRequest(BaseModel):
    model_config = ConfigDict(extra="forbid")
    parameters: dict
    initial: dict[str, float]
    forcing: list[dict[str, float]] = Field(min_length=1, max_length=366)


class InterventionRequest(ScenarioRequest):
    intervention: dict
    community: dict | None = None            # {"sample": id} or {"site", "compartment", "trip"}
    coupling: str | None = None              # reference | gate | proportional
    state: dict | None = None                # {"<day>": {"water_mm"|"theta"|"psi_mpa"|"nacl_mol_m2": v, "version": ...}}


class VersionRelevanceRequest(BaseModel):
    model_config = ConfigDict(extra="forbid")
    target: str
    outcome: str


class ComparisonScenario(BaseModel):
    model_config = ConfigDict(extra="forbid")
    name: str = Field(min_length=1)
    forcing: list[dict[str, float]] = Field(min_length=1, max_length=366)
    disturbance_m2: float = Field(ge=0, allow_inf_nan=False)
    porosity: float | None = Field(default=None, gt=0, lt=1)


class ComparisonRequest(BaseModel):
    model_config = ConfigDict(extra="forbid")
    initial: dict[str, float]
    parameters: dict
    scenarios: list[ComparisonScenario] = Field(min_length=1, max_length=20)
    community: dict | None = None
    coupling: str | None = None


@app.exception_handler(ValueError)
async def invalid_value(request, exc):
    return JSONResponse(status_code=422, content={"detail": str(exc)})


def parameters(values: dict) -> Parameters:
    try:
        return Parameters(**values)
    except TypeError as exc:
        raise ValueError(f"Invalid/missing process parameter: {exc}") from exc


@app.get("/health")
def health():
    return {"status": "ok", "version": VERSION}


@app.get("/", response_class=HTMLResponse)
@app.get("/twin", response_class=HTMLResponse)
def console():
    return (PACKAGE/"console.html").read_text(encoding="utf-8")


@app.get("/v3/contract")
def model_contract():
    return contract()


@app.get("/v3/evidence")
def evidence(limit: int = Query(default=551, ge=1, le=551)):
    panel = pd.read_csv(DATA/"sample_evidence.tsv", sep="\t")
    return {"manifest": json.loads((DATA/"manifest.json").read_text()),
            "observations": json.loads(panel.head(limit).to_json(orient="records"))}


MAG_COLUMNS = ["rep", "domain", "phylum", "order", "family", "genus",
               "completeness", "contamination", "quality_tier", "quality_90_5",
               "cbb_steps_detected", "cbb_all_steps_detected", "h2_1l", "h2_1h", "h2_2a",
               "co_form1", "h2_carrier_union", "h2_cbb_all_steps", "h2_cbb_quality_90_5",
               "dual_h2_co_cbb_all_steps", "cbb_kos", "missing_steps"]


@lru_cache(maxsize=1)
def _mag_table() -> pd.DataFrame:
    return pd.read_csv(DATA/"mag_evidence.tsv", sep="\t", keep_default_na=False)


@app.get("/v3/mags")
def mags(candidates_only: bool = Query(default=False), limit: int = Query(default=100, ge=1, le=5438),
         q: str | None = Query(default=None, max_length=80, description="case-insensitive substring of id or taxon")):
    mag = _mag_table()
    total, candidates = len(mag), int(mag.h2_cbb_all_steps.sum())
    if candidates_only:
        mag = mag.loc[mag.h2_cbb_all_steps]
    else:
        mag = mag.sort_values(["h2_cbb_all_steps", "cbb_all_steps_detected", "completeness"], ascending=[False, False, False])
    if q:
        needle = q.strip().lower()
        hay = mag[["rep", "phylum", "order", "family", "genus"]].astype(str).agg(" ".join, axis=1).str.lower()
        mag = mag.loc[hay.str.contains(needle, regex=False)]
    present = [c for c in MAG_COLUMNS if c in mag.columns]
    return {"total_reference_mags": total,
            "h2_cbb_candidate_mags": candidates,
            "matching_mags": len(mag),
            "mags": json.loads(mag[present].head(limit).to_json(orient="records"))}


@lru_cache(maxsize=1)
def _justifications() -> dict:
    path = DATA/"justifications.json"
    return json.loads(path.read_text()) if path.exists() else {"classes": {}, "tbox_axioms": []}


@app.get("/v3/mags/{mag_id}")
def mag_detail(mag_id: str):
    """Reasoner-entailed classes, the 11 CBB steps, hydrogenase groups and OWL justifications of one MAG."""
    mag = _mag_table()
    match = mag.loc[mag.rep.eq(mag_id)]
    if match.empty:
        raise HTTPException(status_code=404, detail="Unknown reference MAG")
    row = match.iloc[0]
    classified = load_classification(DATA)
    entailed = list(classified["classes_by_mag"][mag_id])
    detected = set(str(row.cbb_kos).split(";")) if str(row.cbb_kos) else set()
    steps = [{"step": step, "present": bool(row[step]), "alternatives": list(alts),
              "detected_kos": sorted(detected & set(alts))} for step, alts in CBB_STEPS.items()]
    notes = class_annotations()
    just = _justifications()
    tbox = just.get("tbox_axioms", [])
    justifications = {}
    for cls, members in just.get("classes", {}).items():
        if mag_id in members:
            entry = members[mag_id]
            justifications[cls] = {"mag_axioms": entry["mag_axioms"],
                                   "tbox_axioms": [tbox[i] for i in entry["tbox_axiom_ids"]],
                                   "enumeration_complete": entry["enumeration_complete"],
                                   "justifications_returned": entry["justifications_returned"],
                                   "subset_minimal": bool(just.get("all_verified_subset_minimal"))}
    taxon = {k: str(row[k]) for k in ("domain", "phylum", "order", "family", "genus")}
    return {"mag_id": mag_id, "taxonomy": taxon,
            "entailed_classes": [{"class": "rak:"+c, "label": notes.get(c, (c, ""))[0], "definition": notes.get(c, (c, ""))[1]}
                                 for c in entailed],
            "cbb_steps": steps, "cbb_steps_detected": int(row.cbb_steps_detected),
            "hydrogenase_groups": {"group_1l": bool(row.h2_1l), "group_1h": bool(row.h2_1h),
                                   "group_2a": bool(row.h2_2a), "co_dehydrogenase_form_I": bool(row.co_form1)},
            "completeness": float(row.completeness), "contamination": float(row.contamination),
            "justifications": justifications,
            "ontology_sha256": classified["ontology_sha256"],
            "note": "Classes are the reasoner's entailments, not computed from table flags. Absence of an entailment is not asserted biological absence."}


@app.get("/v3/transfer")
def transfer():
    """External transfer tests (Namib, Sonoran, temperate): summary tables copied from reviews/v05."""
    path = DATA/"transfer_summary.json"
    if not path.exists():
        raise HTTPException(status_code=404, detail="transfer_summary.json is not bundled; run scripts/web/export_transfer_summary.py")
    return json.loads(path.read_text())


class AdmissibilityRequest(BaseModel):
    model_config = ConfigDict(extra="forbid")
    node: str
    value: float | None = None
    version: str | None = None            # matric | osmotic, for psi


@app.post("/v3/admissibility")
def admissibility(request: AdmissibilityRequest):
    """Verdict for clamping one graph node: zone message, and (given a value) feasibility of the regime."""
    node = request.node
    try:
        resolved = theory.resolve(node)
    except ValueError as exc:
        return {"node": node, "admissible": False, "stage": "zone", "message": str(exc)}
    if resolved == "psi" and request.version not in ("matric", "osmotic"):
        return {"node": resolved, "admissible": False, "stage": "zone",
                "message": "psi is constitutive (psi = psi_m + psi_o), so do(psi) is an ambiguous manipulation: name a version, "
                           "'matric' (handle W) or 'osmotic' (handle N), and it is translated to a state intervention."}
    if resolved in ("theta", "psi"):
        handle = "W" if resolved == "theta" or request.version == "matric" else "N"
        out = {"node": resolved, "admissible": True, "stage": "translated",
               "message": f"{resolved} is constitutive; do({resolved}) is translated into a state intervention on {handle}."}
    else:
        try:
            theory.check_admissible([node])
        except ValueError as exc:
            return {"node": resolved, "admissible": False, "stage": "zone", "message": str(exc)}
        out = {"node": resolved, "admissible": True, "stage": "zone",
               "message": f"{resolved} is a manipulable node of the causal graph."}
    if request.value is None:
        return out
    forcing = {"I": "water_mm", "c": "nacl_mol_l", "Evap": "evaporation_mm", "T": "temp_c"}
    params = {"phi": "porosity", "psi_e": "psi_entry_mpa", "b": "retention_b"}
    states = {"W": "water_mm", "N": "nacl_mol_m2"}
    if resolved in forcing:
        spec = {"forcing": {"0": {forcing[resolved]: request.value}}}
    elif resolved in params:
        spec = {"parameters": {params[resolved]: request.value}}
    else:
        if resolved in states:
            body = {states[resolved]: request.value}
        elif resolved == "theta":
            body = {"theta": request.value}
        else:
            body = {"psi_mpa": request.value, "version": request.version}
        spec = {"state": {"0": body}}
    days = [{"water_mm": 0.0, "evaporation_mm": 1.0, "temp_c": 25.0} for _ in range(14)]
    try:
        result = intervene(State(water_mm=1.8), days, Parameters(k_h2_s=0.001, k_co_s=0.001), spec, DATA)
    except ValueError as exc:
        return {"node": resolved, "admissible": False, "stage": "feasibility", "message": str(exc)}
    out["stage"] = "feasible"
    out["message"] = f"do({resolved} = {request.value:g}) is admissible and feasible on the 14-day reference scenario (initial water 1.8 mm)."
    out["h2_uptake_difference_mol_m2"] = result["gas_uptake"]["h2_mol_m2"]["difference"]
    out["trace"] = result["trace"]["targets_replaced"]
    return out


@app.get("/v3/ontology")
def ontology_reasoning():
    load_classification(DATA)
    return json.loads((DATA / "ontology_report.json").read_text())


@app.post("/v3/simulate")
def scenario(request: ScenarioRequest):
    try:
        initial = State(**request.initial)
    except TypeError as exc:
        raise ValueError(str(exc)) from exc
    return simulate(initial, request.forcing, parameters(request.parameters))


@app.post("/v3/intervene")
def intervention(request: InterventionRequest):
    try:
        initial = State(**request.initial)
    except TypeError as exc:
        raise ValueError(str(exc)) from exc
    spec = dict(request.intervention)
    for key in ("community", "coupling", "state"):
        if getattr(request, key) is not None:
            if key in spec:
                raise ValueError(f"{key} given both in the request and in the intervention")
            spec[key] = getattr(request, key)
    return intervene(initial, request.forcing, parameters(request.parameters), spec, DATA)


@app.get("/v3/graph")
def causal_graph():
    return theory.graph()


@app.get("/v3/community/{sample}")
def community_factors(sample: str):
    return comm.describe(sample)


@app.post("/v3/version-relevance")
def version_relevance(request: VersionRelevanceRequest):
    return theory.version_relevance(request.target, request.outcome)


@app.post("/v3/experiment")
def design(values: dict):
    return experiment(parameters(values))


@app.post("/v3/compare")
def comparison(values: ComparisonRequest):
    try:
        initial = State(**values.initial)
    except TypeError as exc:
        raise ValueError(str(exc)) from exc
    return compare(initial, [s.model_dump(exclude_none=True) for s in values.scenarios], parameter_ensemble(parameters(values.parameters)),
                   values.community, values.coupling)


@app.post("/v3/resource-ceiling")
def resource(values: dict):
    values = dict(values)
    mag_id = values.pop("mag_id", None)
    if not mag_id or {"pathway_supported", "co_oxidation_supported"} & values.keys():
        raise ValueError("Supply mag_id; pathway eligibility is resolved from the metagenomic evidence")
    evidence = pd.read_csv(DATA/"mag_evidence.tsv", sep="\t")
    match = evidence.loc[evidence.rep.eq(mag_id)]
    if match.empty:
        raise ValueError("Unknown reference MAG")
    row = match.iloc[0]
    classified = load_classification(DATA)
    classes = classified["classes_by_mag"][mag_id]
    try:
        supported = "H2AutotrophicCandidateMAG" in classes
        result = resource_ceiling(**values, pathway_supported=supported,
                                  co_oxidation_supported="COCarrierMAG" in classes)
        result["genomic_evidence"] = {"mag_id": mag_id, "h2_and_all_cbb_steps": supported,
                                      "owl_classes": ["rak:ReferenceMAG"] + ["rak:"+c for c in classes],
                                      "ontology_sha256": classified["ontology_sha256"],
                                      "completeness": float(row.completeness), "contamination": float(row.contamination)}
        return result
    except TypeError as exc:
        raise ValueError(str(exc)) from exc


@app.post("/v3/carbon-account")
def carbon(values: dict):
    try:
        return carbon_account(**values)
    except TypeError as exc:
        raise ValueError(str(exc)) from exc


@app.api_route("/v1/{path:path}", methods=["GET", "POST"])
def retired(path: str):
    raise HTTPException(status_code=410, detail="v0.2 heuristic endpoints are retired. Use /v3/contract and /docs for current operations; historical code remains available for reproduction.")
