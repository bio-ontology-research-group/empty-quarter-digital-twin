"""Reproducibility entrypoint for EQ-BioDT v0.8 UQ tasks and local execution mode."""
from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))

from reproducibility.v08.gather import DEFAULT_OUT_DIR, gather  # noqa: E402
from reproducibility.v08.tasks import DEFAULT_TASKS_PATH, enumerate_tasks, write_tasks  # noqa: E402
from reproducibility.v08.uq_task import run_task  # noqa: E402


def _worker(args: tuple[list[dict], int, Path]) -> str:
    tasks, idx, raw_dir = args
    p = run_task(tasks, idx, raw_dir)
    return str(p)


def run_local(
    *,
    jobs: int = 1,
    n: int = 1024,
    chunk_rows: int = 2000,
    slugs: list[str] | None = None,
    raw_dir: Path,
    out_dir: Path = DEFAULT_OUT_DIR,
    tasks_path: Path = DEFAULT_TASKS_PATH,
    boot: int = 500,
) -> dict[str, Path]:
    """Enumerate tasks, run all (or filtered) tasks locally across `jobs` processes, then gather."""
    tasks = enumerate_tasks(n=n, chunk_rows=chunk_rows)
    if slugs:
        allowed = set(slugs)
        tasks = [t for t in tasks if t["slug"] in allowed]
    tasks_path.parent.mkdir(parents=True, exist_ok=True)
    if not slugs:
        write_tasks(path=tasks_path, n=n, chunk_rows=chunk_rows)

    raw_dir.mkdir(parents=True, exist_ok=True)
    work = [(tasks, i, raw_dir) for i in range(len(tasks))]
    if jobs <= 1:
        for item in work:
            _worker(item)
    else:
        with ProcessPoolExecutor(max_workers=jobs) as pool:
            list(pool.map(_worker, work))
    return gather(raw_dir, out_dir, boot=boot)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Reproduce EQ-BioDT v0.8 UQ outputs.")
    parser.add_argument("--local", action="store_true", help="Run UQ tasks locally with multiprocessing.")
    parser.add_argument("-j", "--jobs", type=int, default=4, help="Number of worker processes in --local mode.")
    parser.add_argument("--n", type=int, default=1024, help="Base Sobol sample size N.")
    parser.add_argument("--chunk-rows", type=int, default=2000, help="Maximum rows per task chunk.")
    parser.add_argument("--slug", action="append", default=None, help="Optional question slug(s) to restrict to.")
    parser.add_argument(
        "--raw",
        type=Path,
        default=ROOT / "reviews" / "v08" / "uq" / "raw",
        help="Directory for raw chunk JSON files.",
    )
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT_DIR, help="Output directory for <slug>.json.")
    parser.add_argument("--tasks", type=Path, default=DEFAULT_TASKS_PATH, help="Path for tasks.json.")
    parser.add_argument("--boot", type=int, default=500, help="Bootstrap resamples for Sobol CIs.")
    args = parser.parse_args(argv)

    if args.local:
        written = run_local(
            jobs=args.jobs,
            n=args.n,
            chunk_rows=args.chunk_rows,
            slugs=args.slug,
            raw_dir=args.raw,
            out_dir=args.out,
            tasks_path=args.tasks,
            boot=args.boot,
        )
        print(f"Completed local UQ run and gathered {len(written)} question files in {args.out}")
        return 0

    write_tasks(path=args.tasks, n=args.n, chunk_rows=args.chunk_rows)
    if args.raw.exists() and any(args.raw.glob("*.json")):
        written = gather(args.raw, args.out, boot=args.boot)
        print(f"Gathered {len(written)} question files in {args.out}")
    else:
        print(f"Wrote task manifest to {args.tasks} (pass --local or --raw DIR to gather)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
