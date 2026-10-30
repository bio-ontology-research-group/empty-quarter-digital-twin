"""Reproduce the v0.6 analyses, figures and tables, and record input/output digests.

Usage: uv run python reproducibility/v06/reproduce.py [--skip-v05]

The v0.5 reproduction (evidence rebuild, OWL reasoning with the pinned km reasoner, external-transfer analyses,
v0.5 figures) is run first unless --skip-v05 is given; it regenerates data/v03 and reviews/v05, which the v0.6
analyses read. Everything in reproducibility/v06 itself is deterministic (no random numbers except the
site bootstrap, seed 42, 2,000 draws) and takes about 15 s.
"""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[2]


def run(relative, *args):
    t = time.time()
    subprocess.run([sys.executable, str(ROOT / relative), *args], cwd=ROOT, check=True)
    return round(time.time() - t, 1)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--skip-v05", action="store_true", help="do not rerun the v0.5 reproduction")
    args = ap.parse_args()
    timings = {}
    if not args.skip_v05:
        timings["v05_reproduce"] = run("reproducibility/v05/reproduce.py")
    for step in ("analyse", "figures", "tables", "identification"):
        timings[step] = run(f"reproducibility/v06/{step}.py")
    out = ROOT / "reviews/v06"
    inputs = sorted((ROOT / "data/v03").glob("*")) + sorted((ROOT / "src/eq_biodt/v03").glob("*.py")) \
        + [ROOT / "data/ontology/eq_biodt_v03.ofn", ROOT / "reviews/v05/primary_sites.tsv", ROOT / "reproducibility/v06/data/arabia_natural_earth_50m.json",
           ROOT / "reproducibility/v03/audit.py", ROOT / "data/sources/manifest.json"] \
        + sorted((ROOT / "reproducibility/v06").glob("*.py"))
    outputs = sorted(p for p in out.glob("*") if p.name != "manifest.json" and p.suffix in (".json", ".tsv", ".png")) \
        + [ROOT / "paper/figures" / f for f in ("v06_graph.pdf", "v06_scenarios.pdf", "v06_primary_twin.pdf")] \
        + [ROOT / "paper/tables" / f for f in ("v06_numbers.tex", "v06_admissibility.tex", "v06_identification.tex")]
    manifest = {"release": "0.6.0-dev", "v05_rerun": not args.skip_v05, "timings_seconds": timings,
                "bootstrap": {"draws": 2000, "seed": 42},
                "inputs": {str(p.relative_to(ROOT)): digest(p) for p in inputs if p.exists() and p.is_file()},
                "outputs": {str(p.relative_to(ROOT)): digest(p) for p in outputs if p.exists()},
                "note": "PDF digests may vary with matplotlib version; JSON/TSV/TeX outputs are deterministic."}
    (out / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print("v0.6 reproduction complete", timings)


if __name__ == "__main__":
    main()
