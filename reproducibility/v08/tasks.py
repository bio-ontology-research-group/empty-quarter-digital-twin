"""Enumerate and chunk Saltelli/Jansen UQ tasks for EQ-BioDT v0.8."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))

from eq_biodt.v03.questions import QUESTIONS  # noqa: E402
from eq_biodt.v03.uq import inputs_for_spec  # noqa: E402

DEFAULT_TASKS_PATH = ROOT / "reviews" / "v08" / "uq" / "tasks.json"


def enumerate_tasks(n: int = 1024, chunk_rows: int = 2000) -> list[dict]:
    """Enumerate UQ row-chunk tasks across all 14 questions, series, and levels.

    Each task dictionary has:
    `{key: "slug|series|level|chunk", slug, series, series_index, level, rows: [start, stop], n, d, total_rows}`.
    """
    if n <= 0 or chunk_rows <= 0:
        raise ValueError("n and chunk_rows must be positive integers")
    tasks: list[dict] = []
    for slug, q in QUESTIONS.items():
        for s_idx, s in enumerate(q.series):
            for l_idx, spec in enumerate(s.levels):
                inps = inputs_for_spec(spec)
                d = len(inps)
                total_rows = n * (d + 5)
                chunk = 0
                for start in range(0, total_rows, chunk_rows):
                    stop = min(start + chunk_rows, total_rows)
                    tasks.append(
                        {
                            "key": f"{slug}|{s.id}|{l_idx}|{chunk}",
                            "slug": slug,
                            "series": s.id,
                            "series_index": s_idx,
                            "level": l_idx,
                            "rows": [start, stop],
                            "n": n,
                            "d": d,
                            "total_rows": total_rows,
                        }
                    )
                    chunk += 1
    return tasks


def merge_chunks(tasks: list[dict]) -> list[dict]:
    """Merge row-chunk tasks sharing `(slug, series, level)` into single full-level tasks."""
    grouped: dict[tuple[str, str, int], list[dict]] = {}
    for t in tasks:
        k = (t["slug"], t["series"], int(t["level"]))
        grouped.setdefault(k, []).append(t)

    merged: list[dict] = []
    for (slug, series, level), group in grouped.items():
        group_sorted = sorted(group, key=lambda x: x["rows"][0])
        start = group_sorted[0]["rows"][0]
        stop = group_sorted[-1]["rows"][1]
        base = dict(group_sorted[0])
        base["key"] = f"{slug}|{series}|{level}|0"
        base["rows"] = [start, stop]
        merged.append(base)
    return merged


def write_tasks(
    path: Path = DEFAULT_TASKS_PATH,
    n: int = 1024,
    chunk_rows: int = 2000,
) -> list[dict]:
    tasks = enumerate_tasks(n=n, chunk_rows=chunk_rows)
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(tasks, indent=2, sort_keys=True) + "\n")
    return tasks


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Enumerate EQ-BioDT v0.8 UQ tasks.")
    parser.add_argument("--n", type=int, default=1024, help="Base Sobol sample size N.")
    parser.add_argument("--chunk-rows", type=int, default=2000, help="Maximum rows per Slurm task.")
    parser.add_argument("--out", type=Path, default=DEFAULT_TASKS_PATH, help="Output tasks.json path.")
    args = parser.parse_args(argv)
    tasks = write_tasks(path=args.out, n=args.n, chunk_rows=args.chunk_rows)
    print(f"Wrote {len(tasks)} tasks to {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
