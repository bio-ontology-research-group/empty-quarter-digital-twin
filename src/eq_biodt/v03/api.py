"""Current public API. Historical heuristic interventions are explicitly retired."""
from pathlib import Path
import json

from fastapi import FastAPI, HTTPException, Query
from fastapi.responses import HTMLResponse, JSONResponse
import pandas as pd
from pydantic import BaseModel, ConfigDict, Field

from .contract import contract
from .ontology import evaluate_ontology_reasoning
from .planning import compare, experiment, parameter_ensemble
from .process import Parameters, State, simulate
from .resources import resource_ceiling, carbon_account

PACKAGE = Path(__file__).parent
DATA = PACKAGE/"data" if (PACKAGE/"data").exists() else PACKAGE.parents[2]/"data"/"v03"
app = FastAPI(title="Empty Quarter Digital Twin", version="0.3.0",
              description="Metagenomic evidence and conditional, resource-conserving scenarios. See /v3/contract.")


class ScenarioRequest(BaseModel):
    model_config = ConfigDict(extra="forbid")
    parameters: dict
    initial: dict[str, float]
    forcing: list[dict[str, float]] = Field(min_length=1, max_length=366)


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
    return {"status": "ok", "version": "0.3.0"}


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
    return {"total_reference_mags": 5438,
            "h2_cbb_candidate_mags": int(mag.h2_cbb_all_steps.sum()) if candidates_only else 23,
            "mags": json.loads(mag[present].head(limit).to_json(orient="records"))}


@app.get("/v3/ontology")
def ontology_reasoning():
    mag = pd.read_csv(DATA/"mag_evidence.tsv", sep="\t", keep_default_na=False)
    panel = pd.read_csv(DATA/"sample_evidence.tsv", sep="\t")
    ofn_path = DATA.parent/"ontology"/"eq_biodt_v03.ofn"
    return evaluate_ontology_reasoning(mag, panel, ofn_path)


@app.post("/v3/simulate")
def scenario(request: ScenarioRequest):
    try:
        initial = State(**request.initial)
    except TypeError as exc:
        raise ValueError(str(exc)) from exc
    return simulate(initial, request.forcing, parameters(request.parameters))


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
    if not mag_id or "pathway_supported" in values:
        raise ValueError("Supply mag_id; pathway eligibility is resolved from the metagenomic evidence")
    evidence = pd.read_csv(DATA/"mag_evidence.tsv", sep="\t")
    match = evidence.loc[evidence.rep.eq(mag_id)]
    if match.empty:
        raise ValueError("Unknown reference MAG")
    row = match.iloc[0]
    try:
        result = resource_ceiling(**values, pathway_supported=bool(row.h2_cbb_all_steps))
        owl_classes = ["rak:ReferenceMAG"]
        if bool(row.h2_carrier_union):
            owl_classes.append("rak:H2CarrierMAG")
        if bool(row.cbb_all_steps_detected):
            owl_classes.append("rak:CBBCompleteMAG")
        if bool(row.h2_cbb_all_steps):
            owl_classes.append("rak:H2AutotrophicCandidateMAG")
        if bool(row.h2_cbb_quality_90_5):
            owl_classes.append("rak:HighQualityH2AutotrophicCandidateMAG")
        result["genomic_evidence"] = {"mag_id": mag_id, "h2_and_all_cbb_steps": bool(row.h2_cbb_all_steps),
                                      "owl_classes": owl_classes,
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
    raise HTTPException(status_code=410, detail="v0.2 heuristic endpoints are retired. Use /v3/contract and /docs for v0.3 operations; historical code remains available for reproduction.")
