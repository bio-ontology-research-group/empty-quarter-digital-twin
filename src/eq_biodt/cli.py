"""Command-line interface (`eq-biodt`) for the Empty Quarter Biodiversity Digital Twin."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any, Dict

import uvicorn

from eq_biodt.api.main import get_runtime
from eq_biodt.data.build_twin_state_cube import TwinStateCubeBuilder
from eq_biodt.data.manifest_registry import CanonicalDataRegistry
from eq_biodt.standards import (
    B3OccurrenceCubeBuilder,
    DwCEMoFExporter,
    FDOProvenanceRecorder,
)
from eq_biodt.validation import TwinValidationSuite


def cmd_verify_data() -> int:
    reg = CanonicalDataRegistry()
    report = reg.verify_all_sources()
    print(
        json.dumps(
            {
                "status": "VERIFIED",
                "files_verified": len(report),
                "xrf_analytes_mapped": len(reg.ontology.xrf_mappings),
                "sites_registered": len(reg.ontology.site_iris),
                "sources": report,
            },
            indent=2,
        )
    )
    return 0


def cmd_build_cubes(output_root: Path) -> int:
    reg = CanonicalDataRegistry()
    cube_builder = TwinStateCubeBuilder(reg)
    b3_builder = B3OccurrenceCubeBuilder(cube_builder)
    dwc_exporter = DwCEMoFExporter(cube_builder)

    proc_out = cube_builder.export_artifacts(output_root / "processed")
    b3_out = b3_builder.export_b3_datapackage(output_root / "standards" / "b3_cube")
    dwc_out = dwc_exporter.export_archive(output_root / "standards" / "dwca")

    print(
        json.dumps(
            {
                "status": "EXPORTED",
                "processed_cubes": {k: str(v) for k, v in proc_out.items()},
                "b3_cube": {k: str(v) for k, v in b3_out.items()},
                "dwca_emof": {k: str(v) for k, v in dwc_out.items()},
            },
            indent=2,
        )
    )
    return 0


def cmd_simulate(
    intervention: str,
    site_id: str,
    compartment: str,
    params_json: str,
    save_crate: str | None = None,
) -> int:
    rt = get_runtime()
    site_row = rt.resolve_site_compartment_row(site_id, compartment)
    params: Dict[str, Any] = json.loads(params_json) if params_json else {}

    itype = intervention.strip().lower()
    if itype == "cloud_seeding":
        res = rt.cloud_sim.simulate(site_row, **params)
    elif itype == "revegetation":
        res = rt.reveg_sim.simulate(site_row, **params)
    elif itype == "desalination_brine":
        res = rt.desal_sim.simulate(site_row, **params)
    elif itype == "refugium_thermal":
        res = rt.refugium_sim.simulate(site_row, **params)
    elif itype == "megafauna_grazing":
        res = rt.megafauna_sim.simulate(site_row, **params)
    elif itype == "seismic_hydrocarbon":
        res = rt.seismic_sim.simulate(site_row, **params)
    else:
        raise ValueError(f"Unsupported intervention '{intervention}'")

    if save_crate:
        crate = FDOProvenanceRecorder.build_run_crate(
            intervention_name=itype,
            site_iri=site_row["site_iri"],
            parameters=params,
            outputs=res["summary"],
            ocs_score=res["summary"]["ocs_score"],
        )
        FDOProvenanceRecorder.save_run_crate(crate, Path(save_crate))

    print(json.dumps(res, indent=2, default=str))
    return 0


def cmd_optimize(site_id: str, compartment: str) -> int:
    rt = get_runtime()
    site_row = rt.resolve_site_compartment_row(site_id, compartment)
    res = rt.pareto_opt.optimize_pareto_front(site_row)
    print(json.dumps(res, indent=2, default=str))
    return 0


def cmd_validate(output_dir: Path) -> int:
    suite = TwinValidationSuite()
    report = suite.run_all_validations(output_dir=output_dir)
    print(json.dumps(report, indent=2, default=str))
    return 0


def main() -> None:
    parser = argparse.ArgumentParser(
        prog="eq-biodt",
        description="Empty Quarter Biodiversity Digital Twin (EQ-BioDT) CLI",
    )
    sub = parser.add_subparsers(dest="command", required=True)

    sub.add_parser("verify-data", help="Verify canonical tables, SHA-256 digests, and ontology IRIs")

    p_val = sub.add_parser("validate", help="Run 5-tier quantitative validation suite and export metrics")
    p_val.add_argument("--output-dir", default="data/validation", help="Output directory for validation report")

    p_build = sub.add_parser("build-cubes", help="Build state cubes, B-Cubed datapackage, and DwC-A/eMoF archive")
    p_build.add_argument("--output-dir", default="data", help="Output directory (default: data)")

    p_sim = sub.add_parser("simulate", help="Run a counterfactual environmental management intervention")
    p_sim.add_argument(
        "intervention",
        choices=[
            "cloud_seeding",
            "revegetation",
            "desalination_brine",
            "refugium_thermal",
            "megafauna_grazing",
            "seismic_hydrocarbon",
        ],
    )
    p_sim.add_argument("--site", default="35", help="Site number (1..60) or auxiliary ID (e.g. RAK_1000061)")
    p_sim.add_argument("--compartment", default="Surface", help="Surface, ShallowSubsurface, or RootAdjacent")
    p_sim.add_argument("--params", default="{}", help="JSON object of intervention parameters")
    p_sim.add_argument("--save-crate", default=None, help="Optional path to save Workflow Run RO-Crate JSON-LD")

    p_opt = sub.add_parser("optimize", help="Run 4-objective DESDEO Pareto front optimization at a site")
    p_opt.add_argument("--site", default="35", help="Site number (1..60) or auxiliary ID")
    p_opt.add_argument("--compartment", default="Surface", help="Surface, ShallowSubsurface, or RootAdjacent")

    p_serve = sub.add_parser("serve", help="Start the EQ-BioDT FastAPI + Web Console server")
    p_serve.add_argument("--host", default="127.0.0.1")
    p_serve.add_argument("--port", type=int, default=8095)

    args = parser.parse_args()
    if args.command == "verify-data":
        raise SystemExit(cmd_verify_data())
    if args.command == "validate":
        raise SystemExit(cmd_validate(Path(args.output_dir)))
    if args.command == "build-cubes":
        raise SystemExit(cmd_build_cubes(Path(args.output_dir)))
    if args.command == "simulate":
        raise SystemExit(
            cmd_simulate(args.intervention, args.site, args.compartment, args.params, args.save_crate)
        )
    if args.command == "optimize":
        raise SystemExit(cmd_optimize(args.site, args.compartment))
    if args.command == "serve":
        uvicorn.run("eq_biodt.api.main:app", host=args.host, port=args.port, reload=False)


if __name__ == "__main__":
    main()
