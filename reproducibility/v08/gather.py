"""Gather raw UQ chunk outputs into per-question uncertainty summaries (`reviews/v08/uq/<slug>.json`)."""
from __future__ import annotations

import argparse
from dataclasses import asdict
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))

from eq_biodt.v03.composition import constants_sha256  # noqa: E402
from eq_biodt.v03.questions import QUESTIONS, evaluate  # noqa: E402
from eq_biodt.v03.uq import analyse, inputs_for_spec, task_seed  # noqa: E402
from eq_biodt.v03.uq_inputs import parameter_ranges_sha256  # noqa: E402

DEFAULT_OUT_DIR = ROOT / "reviews" / "v08" / "uq"


RELATIVE = (("h2_uptake_umol_m2", "h2_uptake_contrast_umol_m2", "h2_uptake_change_pct"),
            ("co_uptake_umol_m2", "co_uptake_contrast_umol_m2", "co_uptake_change_pct"))


def _add_relative_changes(rows: list[dict]) -> None:
    """Per-sample percentage change vs the baseline arm of the same draw: 100 * contrast / (uptake - contrast).
    Paired per row, so the rate-constant scale largely cancels. None when either value is missing or the baseline is 0."""
    for r in rows:
        for up, con, out in RELATIVE:
            if up in r or con in r:
                u, c = r.get(up), r.get(con)
                if u is None or c is None:
                    if con in r:
                        r[out] = None
                    continue
                base = u - c
                r[out] = 100.0 * c / base if base > 0 else (0.0 if c == 0 else None)


def _gather_docs(raw_dir: Path | str, boot: int = 500) -> dict[str, dict]:
    raw_path = Path(raw_dir)
    files = sorted(raw_path.glob("*.json"))
    if not files:
        raise FileNotFoundError(f"No raw UQ JSON files found in {raw_path}")

    by_level: dict[tuple[str, str, int], list[dict]] = {}
    for fp in files:
        data = json.loads(fp.read_text())
        if not isinstance(data, dict) or "slug" not in data or "outputs" not in data:
            continue
        k = (str(data["slug"]), str(data["series"]), int(data["level"]))
        by_level.setdefault(k, []).append(data)

    by_slug: dict[str, dict[tuple[str, int], dict]] = {}
    for (slug, series_id, level_idx), chunks in sorted(by_level.items()):
        chunks_sorted = sorted(chunks, key=lambda c: c["rows"][0])
        n = int(chunks_sorted[0]["n"])
        q = QUESTIONS[slug]
        s = next(sr for sr in q.series if sr.id == series_id)
        spec = s.levels[level_idx]
        inputs = inputs_for_spec(spec)
        d = len(inputs)
        expected_total = n * (d + 5)

        pos = 0
        full_outputs: list[dict] = []
        for c in chunks_sorted:
            start, stop = int(c["rows"][0]), int(c["rows"][1])
            if start != pos:
                raise ValueError(
                    f"Contiguity error for {slug}|{series_id}|{level_idx}: expected row {pos}, got {start}"
                )
            full_outputs.extend(c["outputs"])
            pos = stop
        if pos != expected_total:
            raise ValueError(
                f"Incomplete rows for {slug}|{series_id}|{level_idx}: got {pos}, expected {expected_total}"
            )

        _add_relative_changes(full_outputs)
        seed = task_seed(f"{slug}|{series_id}|{level_idx}")
        all_keys = [k for k in full_outputs[0] if not k.startswith("__")] if full_outputs else list(q.headline)
        headline = list(dict.fromkeys([*q.headline, *all_keys]))
        analysed = analyse(inputs, n, full_outputs, headline, boot=boot, seed=seed)
        central = evaluate(spec, keep_twin_result=False)

        by_slug.setdefault(slug, {})[(series_id, level_idx)] = {
            "level": level_idx,
            "label": spec.label,
            "x": spec.x,
            "n": n,
            "d": d,
            "level_seed": seed,
            "inputs": [asdict(i) for i in inputs],
            "central": {k: v for k, v in central.items() if not k.startswith("__")},
            "central_refused": central.get("__refused__"),
            "central_deduced": central.get("__deduced__", []),
            "outputs": analysed,
        }

    docs: dict[str, dict] = {}
    for slug, level_map in sorted(by_slug.items()):
        q = QUESTIONS[slug]
        series_out: list[dict] = []
        for s in q.series:
            lvls = [
                level_map[(s.id, l_idx)]
                for l_idx in range(len(s.levels))
                if (s.id, l_idx) in level_map
            ]
            if lvls:
                series_out.append(
                    {
                        "id": s.id,
                        "label_plain": s.label_plain,
                        "label_science": s.label_science,
                        "levels": lvls,
                    }
                )
        docs[slug] = {
            "slug": slug,
            "theme": q.theme,
            "headline": list(q.headline),
            "constants_sha256": constants_sha256(),
            "parameter_ranges_sha256": parameter_ranges_sha256(),
            "series": series_out,
        }
    return docs


SIG_DIGITS = 6


def _round_floats(x):
    """Round every float to SIG_DIGITS significant digits. Raw outputs differ between machines in the last bits
    (different libm; measured max relative difference 2e-11 between ws and unimatrix), so published summaries are
    rounded to make them identical across machines."""
    if isinstance(x, float):
        return float(f"{x:.{SIG_DIGITS}g}")
    if isinstance(x, dict):
        return {k: _round_floats(v) for k, v in x.items()}
    if isinstance(x, list):
        return [_round_floats(v) for v in x]
    return x


def gather_bytes(raw_dir: Path | str, boot: int = 500) -> dict[str, bytes]:
    """Return deterministic UTF-8 JSON bytes per question slug from `raw_dir` (floats at SIG_DIGITS significant digits)."""
    docs = _round_floats(_gather_docs(raw_dir, boot=boot))
    return {
        slug: (json.dumps(doc, indent=2, sort_keys=True, allow_nan=False) + "\n").encode("utf-8")
        for slug, doc in docs.items()
    }


def gather(
    raw_dir: Path | str,
    out_dir: Path | str = DEFAULT_OUT_DIR,
    boot: int = 500,
) -> dict[str, Path]:
    """Gather raw task outputs in `raw_dir` and write `<out_dir>/<slug>.json`."""
    out_path = Path(out_dir)
    out_path.mkdir(parents=True, exist_ok=True)
    rendered = gather_bytes(raw_dir, boot=boot)
    written: dict[str, Path] = {}
    for slug, content in rendered.items():
        target = out_path / f"{slug}.json"
        target.write_bytes(content)
        written[slug] = target
    return written


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Gather raw UQ task chunks into per-question JSON files.")
    parser.add_argument("--raw", type=Path, required=True, help="Directory containing raw chunk JSON files.")
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT_DIR, help="Output directory for <slug>.json.")
    parser.add_argument("--boot", type=int, default=500, help="Bootstrap resamples for Sobol CIs.")
    args = parser.parse_args(argv)
    written = gather(args.raw, args.out, boot=args.boot)
    print(f"Gathered {len(written)} question UQ files into {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
