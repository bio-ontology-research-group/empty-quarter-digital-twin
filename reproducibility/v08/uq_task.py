"""Execute a single UQ task (row slice of one question level's Saltelli design)."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))

from eq_biodt.v03.composition import constants_sha256  # noqa: E402
from eq_biodt.v03.questions import QUESTIONS  # noqa: E402
from eq_biodt.v03.uq import design, evaluate_rows, inputs_for_spec, task_seed  # noqa: E402
from eq_biodt.v03.uq_inputs import parameter_ranges_sha256  # noqa: E402


def level_seed_key(slug: str, series: str, level: int) -> str:
    return f"{slug}|{series}|{int(level)}"


def run_task(tasks: list[dict] | Path | str, task_index: int, out_dir: Path | str) -> Path:
    """Run task `task_index` from `tasks` and write `<out_dir>/<key>.json`."""
    if isinstance(tasks, (Path, str)):
        task_list = json.loads(Path(tasks).read_text())
    else:
        task_list = tasks
    t = task_list[int(task_index)]

    slug = t["slug"]
    series_id = t["series"]
    level_idx = int(t["level"])
    n = int(t.get("n", 1024))
    start, stop = int(t["rows"][0]), int(t["rows"][1])

    q = QUESTIONS[slug]
    s = next(sr for sr in q.series if sr.id == series_id)
    spec = s.levels[level_idx]
    inputs = inputs_for_spec(spec)
    seed = task_seed(level_seed_key(slug, series_id, level_idx))

    all_rows = design(inputs, n, seed)
    chunk_rows = all_rows[start:stop]
    outputs = evaluate_rows(spec, inputs, chunk_rows)

    payload = {
        "key": t["key"],
        "slug": slug,
        "series": series_id,
        "level": level_idx,
        "n": n,
        "d": len(inputs),
        "rows": [start, stop],
        "level_seed": seed,
        "constants_sha256": constants_sha256(),
        "parameter_ranges_sha256": parameter_ranges_sha256(),
        "outputs": outputs,
    }
    out_path = Path(out_dir) / f"{t['key']}.json"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(payload, indent=2, sort_keys=True, allow_nan=False) + "\n")
    return out_path


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Run one EQ-BioDT v0.8 UQ array task.")
    parser.add_argument("--task-index", type=int, required=True, help="Zero-based index in tasks.json.")
    parser.add_argument(
        "--tasks",
        type=Path,
        default=ROOT / "reviews" / "v08" / "uq" / "tasks.json",
        help="Path to tasks.json.",
    )
    parser.add_argument("--out", type=Path, required=True, help="Directory for raw task JSON outputs.")
    args = parser.parse_args(argv)
    out_path = run_task(args.tasks, args.task_index, args.out)
    print(f"Wrote {out_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
