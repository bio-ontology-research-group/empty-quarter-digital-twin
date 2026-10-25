"""FastAPI Decision-Support Backend for the Empty Quarter Biodiversity Digital Twin (EQ-BioDT).

Provides REST endpoints for:
- Spatiotemporal site, compartment, H3, and 0.25-deg ERA5 cubes
- B-Cubed (`b3verse`) biodiversity indicators and RO-Crate v1.3 / DwC-A standards
- All 6 counterfactual management interventions (`cloud_seeding`, `revegetation`,
  `desalination_brine`, `refugium_thermal`, `megafauna_grazing`, `seismic_hydrocarbon`)
- 4-objective DESDEO-style Pareto front portfolio optimization
- Interactive Web Console served at `/` and `/twin` (compatible with rubalkhali.science).
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Dict, Optional

from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import HTMLResponse, JSONResponse
from pydantic import BaseModel, Field

from eq_biodt.data.build_twin_state_cube import TwinStateCubeBuilder
from eq_biodt.data.manifest_registry import CanonicalDataRegistry
from eq_biodt.data.sparql_bridge import RubAlKhaliSPARQLBridge
from eq_biodt.interventions import (
    CloudSeedingIntervention,
    DesalinationBrineIntervention,
    MegafaunaGrazingIntervention,
    ParetoPortfolioOptimizer,
    RefugiumThermalIntervention,
    RevegetationIntervention,
    SeismicHydrocarbonIntervention,
)
from eq_biodt.standards import B3OccurrenceCubeBuilder
from eq_biodt.validation import TwinValidationSuite

REPO_ROOT = Path(__file__).resolve().parents[3]
WEB_INDEX = REPO_ROOT / "web" / "index.html"
RO_CRATE_FILE = REPO_ROOT / "ro-crate-metadata.json"

app = FastAPI(
    title="EQ-BioDT: Empty Quarter Biodiversity Digital Twin API",
    version="0.2.0",
    description=(
        "Causal-Ontological & Mechanistic Biodiversity Digital Twin for the Rub' al Khali "
        "(Empty Quarter), conforming to EU BioDT, B-Cubed, and LTER-LIFE standards."
    ),
)

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)


class TwinRuntime:
    """Singleton state cache for fast API responses."""

    def __init__(self) -> None:
        self.registry = CanonicalDataRegistry()
        self.cube_builder = TwinStateCubeBuilder(self.registry)
        self.b3_builder = B3OccurrenceCubeBuilder(self.cube_builder)
        self.sparql = RubAlKhaliSPARQLBridge(registry=self.registry)

        self.sites_df = self.cube_builder.build_sites_cube()
        self.comp_df = self.cube_builder.build_compartment_state_cube()
        self.grid_df = self.cube_builder.build_era5_h3_grid()
        self.b3_ind_df = self.b3_builder.compute_h3_cell_indicators()

        self.cloud_sim = CloudSeedingIntervention()
        self.reveg_sim = RevegetationIntervention()
        self.desal_sim = DesalinationBrineIntervention()
        self.refugium_sim = RefugiumThermalIntervention()
        self.megafauna_sim = MegafaunaGrazingIntervention()
        self.seismic_sim = SeismicHydrocarbonIntervention()
        self.pareto_opt = ParetoPortfolioOptimizer()

    def resolve_site_compartment_row(
        self, site_id: str, compartment: str = "Surface"
    ) -> Dict[str, Any]:
        """Resolve a site (`1..60` or `RAK_1000061..RAK_1000070`) + compartment row dict."""
        comp_norm = self.cube_builder.COMPARTMENT_NORMALIZATION.get(compartment, "Surface")
        sid = str(site_id).strip()
        if sid.isdigit():
            s_num = int(sid)
            sub = self.comp_df[
                (self.comp_df["site_num"] == s_num) & (self.comp_df["compartment_norm"] == comp_norm)
            ]
            if not sub.empty:
                return sub.iloc[0].to_dict()

        # Check auxiliary site or match by site_id / site_label
        site_match = self.sites_df[
            (self.sites_df["site_id"] == sid)
            | (self.sites_df["site_label"].str.lower() == sid.lower())
            | (self.sites_df["site_iri"].str.endswith(sid))
        ]
        if site_match.empty:
            raise KeyError(f"Unknown site_id '{site_id}'")

        s_row = site_match.iloc[0].to_dict()
        if s_row.get("site_num") is not None and not (
            isinstance(s_row["site_num"], float) and s_row["site_num"] != s_row["site_num"]
        ):
            s_num = int(s_row["site_num"])
            sub = self.comp_df[
                (self.comp_df["site_num"] == s_num) & (self.comp_df["compartment_norm"] == comp_norm)
            ]
            if not sub.empty:
                return sub.iloc[0].to_dict()

        # For auxiliary sites (RAK_1000061..RAK_1000070), interpolate baseline from nearest primary site
        primary_comp = self.comp_df[self.comp_df["compartment_norm"] == comp_norm].copy()
        dists = (primary_comp["lat"] - float(s_row["lat"])) ** 2 + (
            primary_comp["lon"] - float(s_row["lon"])
        ) ** 2
        nearest = primary_comp.loc[dists.idxmin()].to_dict()
        nearest.update(
            {
                "site_id": s_row["site_id"],
                "site_label": s_row["site_label"],
                "site_iri": s_row["site_iri"],
                "habitat_regime": s_row["habitat_regime"],
                "lat": s_row["lat"],
                "lon": s_row["lon"],
                "elevation_m": s_row["elevation_m"],
                "h3_5": s_row["h3_5"],
                "h3_6": s_row["h3_6"],
            }
        )
        return nearest


_RUNTIME: Optional[TwinRuntime] = None


def get_runtime() -> TwinRuntime:
    global _RUNTIME
    if _RUNTIME is None:
        _RUNTIME = TwinRuntime()
    return _RUNTIME


class InterventionRequest(BaseModel):
    site_id: str = Field(default="35", description="Site number (1..60) or auxiliary ID (e.g., RAK_1000061)")
    compartment: str = Field(default="Surface", description="Surface, ShallowSubsurface, or RootAdjacent")
    parameters: Dict[str, Any] = Field(default_factory=dict)


class OptimizeRequest(BaseModel):
    site_id: str = Field(default="35", description="Site number (1..60) or auxiliary ID")
    compartment: str = Field(default="Surface", description="Surface, ShallowSubsurface, or RootAdjacent")
    weights: Optional[Dict[str, float]] = Field(
        default=None,
        description="Preference weights for biodiversity, chemosynthesis_cn, salinity_compaction_safety, cost_and_arg_biosafety",
    )


@app.get("/api/v1/health")
def api_health() -> Dict[str, Any]:
    rt = get_runtime()
    return {
        "status": "ok",
        "service": "EQ-BioDT",
        "version": "0.2.0",
        "sites_count": len(rt.sites_df),
        "compartment_states_count": len(rt.comp_df),
        "era5_h3_grid_cells": len(rt.grid_df),
        "virtuoso_sparql_endpoint": rt.sparql.endpoint_url,
        "virtuoso_live": rt.sparql.is_endpoint_alive(),
    }


@app.get("/api/v1/sites")
def api_list_sites() -> Dict[str, Any]:
    rt = get_runtime()
    return {
        "count": len(rt.sites_df),
        "sites": rt.sites_df.where(rt.sites_df.notna(), None).to_dict(orient="records"),
    }


@app.get("/api/v1/sites/{site_id}")
def api_site_detail(site_id: str) -> Dict[str, Any]:
    rt = get_runtime()
    try:
        comp_rows = []
        for comp in ("Surface", "ShallowSubsurface", "RootAdjacent"):
            row = rt.resolve_site_compartment_row(site_id, comp)
            comp_rows.append(row)
        return {
            "site_id": site_id,
            "site_iri": comp_rows[0]["site_iri"],
            "site_label": comp_rows[0]["site_label"],
            "compartments": comp_rows,
        }
    except KeyError as exc:
        raise HTTPException(status_code=404, detail=str(exc)) from exc


@app.get("/api/v1/grid")
def api_corridor_grid() -> Dict[str, Any]:
    rt = get_runtime()
    return {
        "count": len(rt.grid_df),
        "grid": rt.grid_df.to_dict(orient="records"),
    }


@app.get("/api/v1/cubes/b3-indicators")
def api_b3_indicators() -> Dict[str, Any]:
    rt = get_runtime()
    return {
        "standard": "B-Cubed (b3verse) H3 Biodiversity & Functional Guild Cube",
        "count": len(rt.b3_ind_df),
        "indicators": rt.b3_ind_df.to_dict(orient="records"),
    }


@app.get("/api/v1/standards/ro-crate")
def api_ro_crate() -> JSONResponse:
    if not RO_CRATE_FILE.exists():
        raise HTTPException(status_code=404, detail="ro-crate-metadata.json not found")
    data = json.loads(RO_CRATE_FILE.read_text(encoding="utf-8"))
    return JSONResponse(content=data)


@app.post("/api/v1/intervene/{intervention_type}")
def api_run_intervention(intervention_type: str, req: InterventionRequest) -> Dict[str, Any]:
    rt = get_runtime()
    try:
        site_row = rt.resolve_site_compartment_row(req.site_id, req.compartment)
    except KeyError as exc:
        raise HTTPException(status_code=404, detail=str(exc)) from exc

    itype = intervention_type.strip().lower()
    params = req.parameters or {}

    if itype == "cloud_seeding":
        return rt.cloud_sim.simulate(site_row, **params)
    if itype == "revegetation":
        return rt.reveg_sim.simulate(site_row, **params)
    if itype == "desalination_brine":
        return rt.desal_sim.simulate(site_row, **params)
    if itype == "refugium_thermal":
        return rt.refugium_sim.simulate(site_row, **params)
    if itype == "megafauna_grazing":
        return rt.megafauna_sim.simulate(site_row, **params)
    if itype == "seismic_hydrocarbon":
        return rt.seismic_sim.simulate(site_row, **params)

    raise HTTPException(
        status_code=400,
        detail=(
            f"Unknown intervention_type '{intervention_type}'. Supported: "
            "cloud_seeding, revegetation, desalination_brine, refugium_thermal, "
            "megafauna_grazing, seismic_hydrocarbon."
        ),
    )


@app.post("/api/v1/optimize")
def api_optimize_portfolio(req: OptimizeRequest) -> Dict[str, Any]:
    rt = get_runtime()
    try:
        site_row = rt.resolve_site_compartment_row(req.site_id, req.compartment)
    except KeyError as exc:
        raise HTTPException(status_code=404, detail=str(exc)) from exc
    return rt.pareto_opt.optimize_pareto_front(site_row, weights=req.weights)


@app.get("/api/v1/carbon-budget")
def api_carbon_budget() -> Dict[str, Any]:
    rt = get_runtime()
    return rt.pareto_opt.evaluate_basin_carbon_budgets(rt.comp_df, rt.grid_df)


@app.get("/api/v1/validation")
def api_validation_metrics() -> Dict[str, Any]:
    rt = get_runtime()
    suite = TwinValidationSuite(rt.cube_builder)
    return suite.run_all_validations()


@app.get("/", response_class=HTMLResponse)
@app.get("/twin", response_class=HTMLResponse)
def serve_web_console() -> HTMLResponse:
    if WEB_INDEX.exists():
        return HTMLResponse(content=WEB_INDEX.read_text(encoding="utf-8"))
    return HTMLResponse(content="<h1>EQ-BioDT API Running</h1>")
