"""Reproduce the v0.7 composition experiments: digitisation, experiments, figure and table; record digests.

Usage: uv run python reproducibility/v07/reproduce.py

Result files (reviews/v07/composition/*.json, paper/tables/v07_composition.tex, paper/figures/v07_composition.{pdf,png}) contain no
timestamps and are byte-identical across runs; timings are written only to reviews/v07/reproducibility.json.
"""
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[2]


def run(relative):
    t = time.time()
    subprocess.run([sys.executable, str(ROOT / relative)], cwd=ROOT, check=True)
    return round(time.time() - t, 1)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    timings = {step: run(f"reproducibility/v07/{step}.py") for step in ("digitise_fig1", "composition", "figures")}
    src = ROOT / "src/eq_biodt/v03"
    inputs = [ROOT / "data/external/composition_constants.json", ROOT / "data/external/de_zelicourt_2018/fig1.png",
              ROOT / "data/external/de_zelicourt_2018/fig1_digitised.json", ROOT / "data/ontology/eq_biodt_v03.ofn",
              src / "process.py", src / "composition.py", src / "external_models.py", src / "community.py", src / "interventions.py",
              ROOT / "reproducibility/v06/analyse.py"] + sorted((ROOT / "reproducibility/v07").glob("*.py"))
    out = ROOT / "reviews/v07/composition"
    outputs = sorted(out.glob("*.json")) + [ROOT / "paper/tables/v07_composition.tex", ROOT / "paper/figures/v07_composition.pdf",
                                            ROOT / "paper/figures/v07_composition.png"]
    manifest = {"release": "0.7.0-dev", "timings_seconds": timings,
                "inputs": {str(p.relative_to(ROOT)): digest(p) for p in inputs},
                "outputs": {str(p.relative_to(ROOT)): digest(p) for p in outputs},
                "note": "All outputs are deterministic and byte-identical across runs (timings only here). reproduce.py itself is listed under inputs."}
    (ROOT / "reviews/v07/reproducibility.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    print("v0.7 reproduction complete", timings)


if __name__ == "__main__":
    main()
