"""Current public API. Historical heuristic interventions are explicitly retired."""
from pathlib import Path
import json

from fastapi import FastAPI, HTTPException, Query
from fastapi.responses import HTMLResponse, JSONResponse
import pandas as pd
from pydantic import BaseModel, ConfigDict, Field

from .contract import contract
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
        result["genomic_evidence"] = {"mag_id": mag_id, "h2_and_all_cbb_steps": bool(row.h2_cbb_all_steps),
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
