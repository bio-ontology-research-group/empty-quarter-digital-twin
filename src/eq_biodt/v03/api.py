"""Current public API. Historical heuristic interventions are explicitly retired."""
from pathlib import Path
import json

from fastapi import FastAPI, HTTPException, Query
from fastapi.responses import HTMLResponse, JSONResponse
import pandas as pd
from pydantic import BaseModel, ConfigDict, Field

from .contract import contract
from .interventions import intervene
from .ontology import load_classification
from .planning import compare, experiment, parameter_ensemble
from .process import Parameters, State, simulate
from .resources import resource_ceiling, carbon_account

PACKAGE = Path(__file__).parent
DATA = PACKAGE/"data" if (PACKAGE/"data").exists() else PACKAGE.parents[2]/"data"/"v03"
app = FastAPI(title="Empty Quarter Digital Twin", version="0.4.0",
              description="Metagenomic evidence and conditional, resource-conserving scenarios. See /v3/contract.")


class ScenarioRequest(BaseModel):
    model_config = ConfigDict(extra="forbid")
    parameters: dict
    initial: dict[str, float]
    forcing: list[dict[str, float]] = Field(min_length=1, max_length=366)


class InterventionRequest(ScenarioRequest):
    intervention: dict


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
    return {"status": "ok", "version": "0.4.0"}


@app.get("/", response_class=HTMLResponse)
@app.get("/twin", response_class=HTMLResponse)
def console():
    return (PACKAGE/"console.html").read_text()


@app.get("/v3/contract")
def model_contract():
    return contract()


@app.get("/v3/evidence")
def evidence(limit: int = Query(default=551, ge=1, le=551)):
    panel = pd.read_csv(DATA/"sample_evidence.tsv", sep="\t")
    return {"manifest": json.loads((DATA/"manifest.json").read_text()),
            "observations": json.loads(panel.head(limit).to_json(orient="records"))}


@app.get("/v3/mags")
def mags(candidates_only: bool = Query(default=False), limit: int = Query(default=100, ge=1, le=5438)):
    mag = pd.read_csv(DATA/"mag_evidence.tsv", sep="\t", keep_default_na=False)
    total, candidates = len(mag), int(mag.h2_cbb_all_steps.sum())
    if candidates_only:
        mag = mag.loc[mag.h2_cbb_all_steps]
    else:
        mag = mag.sort_values(["h2_cbb_all_steps", "cbb_all_steps_detected", "completeness"], ascending=[False, False, False])
    cols = ["rep", "domain", "phylum", "order", "family", "genus",
            "completeness", "contamination", "quality_tier", "quality_90_5",
            "cbb_steps_detected", "cbb_all_steps_detected", "h2_1l", "h2_1h", "h2_2a",
            "co_form1", "h2_carrier_union", "h2_cbb_all_steps", "h2_cbb_quality_90_5",
            "dual_h2_co_cbb_all_steps", "cbb_kos", "missing_steps"]
    present = [c for c in cols if c in mag.columns]
    return {"total_reference_mags": total,
            "h2_cbb_candidate_mags": candidates,
            "mags": json.loads(mag[present].head(limit).to_json(orient="records"))}


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
    return intervene(initial, request.forcing, parameters(request.parameters), request.intervention, DATA)


@app.post("/v3/experiment")
def design(values: dict):
    return experiment(parameters(values))


@app.post("/v3/compare")
def comparison(values: ComparisonRequest):
    try:
        initial = State(**values.initial)
    except TypeError as exc:
        raise ValueError(str(exc)) from exc
    return compare(initial, [s.model_dump(exclude_none=True) for s in values.scenarios], parameter_ensemble(parameters(values.parameters)))


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
