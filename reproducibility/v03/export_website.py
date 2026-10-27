"""Export computed v0.3 evidence and scenarios without browser heuristics."""
import argparse
import json
from pathlib import Path
import pandas as pd
from eq_biodt.v03.contract import contract
from eq_biodt.v03.evidence import sha256
from eq_biodt.v03.process import Parameters, State, simulate


def export(destination: Path):
    root=Path(__file__).resolve().parents[2]
    data=root/"data/v03"; report=root/"reviews/v03"
    settings=json.loads((data/"example_comparison.json").read_text())
    bundle={"version":"0.3.0", "results":json.loads((report/"results.json").read_text()),
            "manifest":json.loads((data/"manifest.json").read_text()), "contract":contract(),
            "observations":json.loads(pd.read_csv(data/"sample_evidence.tsv",sep="\t").to_json(orient="records")),
            "backtest":json.loads((report/"backtest.json").read_text())["metrics"],
            "comparison":json.loads((report/"comparison.json").read_text()),
            "examples":[{"name":s["name"],"result":simulate(State(**settings["initial"]),s["forcing"],Parameters(**settings["parameters"]))} for s in settings["scenarios"]],
            "sources":{name:sha256(data/name) for name in ["sample_evidence.tsv","mag_evidence.tsv","model_contract.json"]}}
    destination.parent.mkdir(parents=True,exist_ok=True)
    destination.write_text(json.dumps(bundle,separators=(",",":"),allow_nan=False)+"\n")
    print(destination)


if __name__=="__main__":
    parser=argparse.ArgumentParser(); parser.add_argument("destination",type=Path)
    export(parser.parse_args().destination)
