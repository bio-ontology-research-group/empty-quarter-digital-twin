"""Inventory real ambiguities for a proposed experiment; no method evaluation.

Run from the repository root with ``uv run python
reproducibility/v03/ontology_design_audit.py``. Only bundled evidence is needed.
"""
from __future__ import annotations

import csv
import hashlib
import json
import math
from pathlib import Path
from statistics import median


ROOT = Path(__file__).resolve().parents[2]


def main() -> None:
    paths = [ROOT / "data/v03" / name for name in
             ("manifest.json", "mag_evidence.tsv", "sample_evidence.tsv")]
    digests = {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
    manifest = json.loads(paths[0].read_text())
    for name, expected in manifest["artifacts"].items():
        if name in digests and digests[name] != expected:
            raise ValueError(f"Evidence differs from frozen manifest: {name}")
    with paths[1].open() as handle:
        mags = list(csv.DictReader(handle, delimiter="\t"))
    with paths[2].open() as handle:
        samples = list(csv.DictReader(handle, delimiter="\t"))
    if len({m["rep"] for m in mags}) != len(mags):
        raise ValueError("Duplicate MAG identity")
    yes = lambda row, key: row[key] == "True"
    strata = {"H_and_C": [], "H_without_C": [], "C_without_H": [], "neither": []}
    overlap = 0
    marker_only = 0
    for mag in mags:
        h, c = yes(mag, "h2_carrier_union"), yes(mag, "cbb_all_steps_detected")
        key = "H_and_C" if h and c else "H_without_C" if h else "C_without_H" if c else "neither"
        strata[key].append(mag["rep"])
        overlap += sum(yes(mag, k) for k in ("h2_1l", "h2_1h", "h2_2a")) > 1
        marker_only += yes(mag, "rbcl_detected") and not c
    widths, mapped, missing = [], [], []
    for sample in samples:
        fields = [sample[k + "_share_mapped"] for k in
                  ("h2_carrier_union", "cbb_all_steps_detected", "h2_cbb_all_steps")]
        if any(x in ("NA", "") for x in fields):
            missing.append(sample["sample"])
            continue
        h, c, intersection = map(float, fields)
        if not all(math.isfinite(v) and 0 <= v <= 1 for v in (h, c, intersection)):
            raise ValueError("Invalid observed mapped-read share")
        lower, upper = max(0., h + c - 1.), min(h, c)
        if not lower - 1e-12 <= intersection <= upper + 1e-12:
            raise ValueError("Observed intersection violates same-measure bounds")
        widths.append(upper - lower)
        mapped.append({"sample": sample["sample"], "H_share": h, "C_share": c,
                       "intersection_share": intersection,
                       "marginal_only_lower": lower, "marginal_only_upper": upper})
    result = {
        "status": "design_feasibility_only_not_experiment_results",
        "inputs_sha256": digests,
        "reference_MAGs": len(mags),
        "screen_strata": {k: len(v) for k, v in strata.items()},
        "first_MAG_per_stratum_lexicographic": {k: min(v) if v else None for k, v in strata.items()},
        "multiple_hydrogenase_category_MAGs": overlap,
        "rbcl_detected_without_all_CBB_steps_MAGs": marker_only,
        "samples": len(samples), "samples_with_mapped_queries": len(mapped),
        "missing_mapped_queries": missing,
        "marginal_only_interval_width": {
            "positive_samples": sum(w > 1e-12 for w in widths),
            "minimum": min(widths), "median": median(widths), "maximum": max(widths),
            "unit": "fraction_of_mapped_read_weight",
            "interpretation": "Frechet outer bounds given H and C shares only; no biological rate bounds",
        },
        "first_sample_lexicographic": min(mapped, key=lambda s: s["sample"]),
        "evidence_scope": "Existing annotation screens and mapping weights; no independent biological truth or method comparison",
    }
    target = ROOT / "reviews/ontology-experiment/design-audit.json"
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
