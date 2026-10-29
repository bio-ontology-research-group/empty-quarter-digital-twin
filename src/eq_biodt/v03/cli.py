"""Reproducible public commands for v0.5."""
import argparse
import json
from pathlib import Path

import pandas as pd

from .api import DATA, resource
from .contract import contract
from .interventions import intervene
from .evidence import build_evidence, build_bundled_evidence
from .forecast import backtest, freeze, score
from .ontology import evaluate_ontology_reasoning
from .planning import compare, experiment, parameter_ensemble
from .process import Parameters, State, simulate


def main():
    parser = argparse.ArgumentParser(description="EQ-BioDT v0.5: evidence, conditional scenarios and prospective tests")
    sub = parser.add_subparsers(dest="command", required=True)
    sub.add_parser("contract")
    p = sub.add_parser("reason")
    p.add_argument("--ofn", "--output-ofn", type=Path, default=Path("eq_biodt.ofn"))
    p = sub.add_parser("build-evidence")
    for field in ("full-ko", "quality", "metadata-panel", "read-mapping", "marker-classification", "marker-ge", "output"):
        p.add_argument("--"+field, type=Path, required=True)
    p = sub.add_parser("rebuild-evidence")
    p.add_argument("--output", type=Path, required=True)
    p = sub.add_parser("backtest")
    p.add_argument("--panel", type=Path, default=DATA/"sample_evidence.tsv")
    p = sub.add_parser("freeze")
    p.add_argument("--panel", type=Path, default=DATA/"sample_evidence.tsv")
    p.add_argument("--target-date", required=True)
    p.add_argument("--output", required=True, type=Path)
    p = sub.add_parser("score")
    p.add_argument("--forecast", required=True, type=Path)
    p.add_argument("--observations", required=True, type=Path)
    for name in ("simulate", "compare", "experiment", "intervene", "resource-ceiling"):
        p = sub.add_parser(name)
        p.add_argument("--input", required=True, type=Path)
    p = sub.add_parser("serve")
    p.add_argument("--host", default="127.0.0.1")
    p.add_argument("--port", type=int, default=8000)
    args = parser.parse_args()
    if args.command == "serve":
        import uvicorn
        uvicorn.run("eq_biodt.v03.api:app", host=args.host, port=args.port)
        return 0
    if args.command == "contract":
        result = contract()
    elif args.command == "reason":
        mag = pd.read_csv(DATA/"mag_evidence.tsv", sep="\t", keep_default_na=False)
        panel = pd.read_csv(DATA/"sample_evidence.tsv", sep="\t")
        result = evaluate_ontology_reasoning(mag, panel, args.ofn)
    elif args.command == "build-evidence":
        result = build_evidence(full_ko=args.full_ko, quality=args.quality, metadata_panel=args.metadata_panel,
                                read_mapping=args.read_mapping, marker_classification=args.marker_classification,
                                marker_ge=args.marker_ge, output=args.output)
    elif args.command == "rebuild-evidence":
        result = build_bundled_evidence(args.output)
    elif args.command == "backtest":
        result = backtest(pd.read_csv(args.panel, sep="\t"))
    elif args.command == "freeze":
        result = freeze(args.panel, args.target_date, args.output)
    elif args.command == "score":
        result = score(json.loads(args.forecast.read_text()), pd.read_csv(args.observations, sep="\t"))
    else:
        values = json.loads(args.input.read_text())
        if args.command == "resource-ceiling":
            print(json.dumps(resource(values), indent=2, allow_nan=False))
            return 0
        params = Parameters(**values["parameters"])
        if args.command == "experiment":
            result = experiment(params)
        elif args.command == "intervene":
            result = intervene(State(**values["initial"]), values["forcing"], params, values["intervention"], DATA)
        elif args.command == "compare":
            result = compare(State(**values["initial"]), values["scenarios"], parameter_ensemble(params))
        else:
            result = simulate(State(**values["initial"]), values["forcing"], params)
    print(json.dumps(result, indent=2, allow_nan=False))
    return 0


if __name__ == "__main__":
    main()
