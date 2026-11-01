"""Reproduce the v0.8 uncertainty quantification outputs, question files, sensitivity figure and LaTeX tables; record digests.

Usage:
  uv run python reproducibility/v08/reproduce.py
  uv run python reproducibility/v08/reproduce.py --local -j 4 [--slug SLUG] [--series SERIES] [--level LEVEL] [--n N]
"""
from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor
import hashlib
import json
from pathlib import Path
import sys
import tempfile
import time
from typing import Any

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))

from reproducibility.v08.figures import main as build_figures  # noqa: E402
from reproducibility.v08.gather import DEFAULT_OUT_DIR, _round_floats, gather, gather_bytes  # noqa: E402
from reproducibility.v08.paper_numbers import write_numbers_tex  # noqa: E402
from reproducibility.v08.questions import DEFAULT_OUT_DIR as DEFAULT_QUESTIONS_DIR, build_questions  # noqa: E402
from reproducibility.v08.tasks import enumerate_tasks  # noqa: E402
from reproducibility.v08.uq_task import run_task  # noqa: E402


def _worker(args: tuple[list[dict], int, Path]) -> str:
    tasks, idx, raw_dir = args
    p = run_task(tasks, idx, raw_dir)
    return str(p)


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _compare_summaries(
    recomputed: dict[str, bytes],
    stored_dir: Path,
    slugs: list[str] | None = None,
    series_filter: str | None = None,
    level_filter: int | None = None,
) -> list[str]:
    diffs: list[str] = []
    for slug, raw_bytes in sorted(recomputed.items()):
        if slugs and slug not in slugs:
            continue
        stored_path = stored_dir / f"{slug}.json"
        if not stored_path.exists():
            diffs.append(f"{slug}: stored summary {stored_path} does not exist")
            continue
        new_doc = _round_floats(json.loads(raw_bytes.decode("utf-8")))
        old_doc = _round_floats(json.loads(stored_path.read_text(encoding="utf-8")))
        old_series = {s["id"]: s for s in old_doc.get("series", [])}
        for s_new in new_doc.get("series", []):
            sid = s_new["id"]
            if series_filter and sid != series_filter:
                continue
            s_old = old_series.get(sid)
            if s_old is None:
                diffs.append(f"{slug}|{sid}: missing series in stored summary")
                continue
            old_levels = {int(l["level"]): l for l in s_old.get("levels", [])}
            for l_new in s_new.get("levels", []):
                lvl = int(l_new["level"])
                if level_filter is not None and lvl != level_filter:
                    continue
                l_old = old_levels.get(lvl)
                if l_old is None:
                    diffs.append(f"{slug}|{sid}|{lvl}: missing level in stored summary")
                    continue
                if l_new.get("n") != l_old.get("n"):
                    diffs.append(f"{slug}|{sid}|{lvl}: sample size n={l_new.get('n')} differs from stored n={l_old.get('n')}")
                if l_new.get("central") != l_old.get("central"):
                    diffs.append(f"{slug}|{sid}|{lvl}: central differs")
                if l_new.get("outputs") != l_old.get("outputs"):
                    for k, v_new in l_new.get("outputs", {}).items():
                        v_old = l_old.get("outputs", {}).get(k)
                        if v_new != v_old:
                            diffs.append(
                                f"{slug}|{sid}|{lvl}|{k}: median new={v_new.get('median')} vs stored={v_old.get('median') if v_old else None}, "
                                f"q05 new={v_new.get('q05')} vs stored={v_old.get('q05') if v_old else None}, "
                                f"q95 new={v_new.get('q95')} vs stored={v_old.get('q95') if v_old else None}"
                            )
    return diffs


def run_local(
    *,
    jobs: int = 1,
    n: int = 1024,
    chunk_rows: int = 2000,
    slugs: list[str] | None = None,
    series: str | None = None,
    level: int | None = None,
    raw_dir: Path | None = None,
    out_dir: Path = DEFAULT_OUT_DIR,
    boot: int = 500,
    write_outputs: bool = False,
) -> tuple[dict[str, bytes], list[str]]:
    """Enumerate tasks, run all (or filtered) tasks locally across `jobs` processes, and compare with stored summaries."""
    tasks = enumerate_tasks(n=n, chunk_rows=chunk_rows)
    if slugs:
        allowed = set(slugs)
        tasks = [t for t in tasks if t["slug"] in allowed]
    if series is not None:
        tasks = [t for t in tasks if t["series"] == series]
    if level is not None:
        tasks = [t for t in tasks if int(t["level"]) == int(level)]

    def _execute(target_raw: Path) -> dict[str, bytes]:
        target_raw.mkdir(parents=True, exist_ok=True)
        work = [(tasks, i, target_raw) for i in range(len(tasks))]
        if jobs <= 1:
            for item in work:
                _worker(item)
        else:
            with ProcessPoolExecutor(max_workers=jobs) as pool:
                list(pool.map(_worker, work))
        if write_outputs:
            gather(target_raw, out_dir, boot=boot)
        return gather_bytes(target_raw, boot=boot)

    if raw_dir is not None:
        recomputed = _execute(raw_dir)
    else:
        tmp_base = Path.home() / "tmp"
        tmp_base.mkdir(parents=True, exist_ok=True)
        with tempfile.TemporaryDirectory(dir=tmp_base, prefix="eq_biodt_v08_local_") as td:
            recomputed = _execute(Path(td))

    diffs = _compare_summaries(recomputed, out_dir, slugs=slugs, series_filter=series, level_filter=level)
    return recomputed, diffs


def write_manifest(uq_dir: Path, questions_dir: Path) -> Path:
    src = ROOT / "src/eq_biodt/v03"
    run_record_path = uq_dir / "run_record.json"
    run_record = json.loads(run_record_path.read_text(encoding="utf-8"))
    inputs = [
        ROOT / "data/external/composition_constants.json",
        ROOT / "data/external/twin_parameter_ranges.json",
        ROOT / "data/ontology/eq_biodt_v03.ofn",
        run_record_path,
        src / "community.py",
        src / "composition.py",
        src / "external_models.py",
        src / "interventions.py",
        src / "process.py",
        src / "questions.py",
        src / "uq.py",
        src / "uq_inputs.py",
    ] + sorted((ROOT / "reproducibility/v08").glob("*.py"))

    uq_outputs = [p for p in sorted(uq_dir.glob("*.json")) if p.name != "run_record.json"]
    q_outputs = sorted(questions_dir.glob("*.json"))
    paper_outputs = [
        ROOT / "paper/tables/v08_numbers.tex",
        ROOT / "paper/tables/v08_parameters.tex",
        ROOT / "paper/tables/v08_uq.tex",
        ROOT / "paper/figures/v08_sensitivity.pdf",
        ROOT / "paper/figures/v08_sensitivity.png",
    ]
    outputs = uq_outputs + q_outputs + paper_outputs

    manifest: dict[str, Any] = {
        "release": "0.8.3",
        "raw_result_digests": {
            "full_run_n1024": {
                "evaluations": run_record["full_run"]["evaluations"],
                "n": run_record["full_run"]["n"],
                "raw_sha256_of_sha256s": run_record["full_run"]["raw_sha256_of_sha256s"],
                "tasks": run_record["full_run"]["tasks"],
            },
            "convergence_rerun_n2048": {
                "n": run_record["convergence_rerun"]["n"],
                "raw_sha256_of_sha256s": run_record["convergence_rerun"]["raw_sha256_of_sha256s"],
                "still_flagged_after_rerun": run_record["convergence_rerun"]["still_flagged_after_rerun"],
                "tasks": run_record["convergence_rerun"]["tasks"],
            },
            "missing_co_genes_rerun": {
                "machine": run_record["missing_co_genes_rerun"]["machine"],
                "n": run_record["missing_co_genes_rerun"]["n"],
                "raw_sha256_of_sha256s": run_record["missing_co_genes_rerun"]["raw_sha256_of_sha256s"],
                "tasks": run_record["missing_co_genes_rerun"]["tasks"],
            },
        },
        "inputs": {str(p.relative_to(ROOT)): digest(p) for p in inputs},
        "outputs": {str(p.relative_to(ROOT)): digest(p) for p in outputs},
        "note": (
            "All generated question files, tables, figures and this manifest are deterministic and byte-identical "
            "across runs. Raw cluster evaluation digests are recorded from reviews/v08/uq/run_record.json."
        ),
    }
    manifest_path = ROOT / "reviews/v08/reproducibility.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return manifest_path


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Reproduce EQ-BioDT v0.8 UQ outputs, question files, tables and figures.")
    parser.add_argument("--local", action="store_true", help="Recompute UQ locally with multiprocessing and compare with stored summaries.")
    parser.add_argument("-j", "--jobs", type=int, default=4, help="Number of worker processes in --local mode.")
    parser.add_argument("--n", type=int, default=1024, help="Base Sobol sample size N.")
    parser.add_argument("--chunk-rows", type=int, default=2000, help="Maximum rows per task chunk.")
    parser.add_argument("--slug", action="append", default=None, help="Optional question slug(s) to restrict to in --local mode.")
    parser.add_argument("--series", type=str, default=None, help="Optional series id to restrict to in --local mode.")
    parser.add_argument("--level", type=int, default=None, help="Optional level index to restrict to in --local mode.")
    parser.add_argument("--raw", type=Path, default=None, help="Optional directory of raw chunk JSON files to gather.")
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT_DIR, help="Directory for gathered UQ <slug>.json files.")
    parser.add_argument("--questions-out", type=Path, default=DEFAULT_QUESTIONS_DIR, help="Directory for question <slug>.json files.")
    parser.add_argument("--boot", type=int, default=500, help="Bootstrap resamples for Sobol CIs.")
    args = parser.parse_args(argv)

    if args.local:
        t0 = time.time()
        recomputed, diffs = run_local(
            jobs=args.jobs,
            n=args.n,
            chunk_rows=args.chunk_rows,
            slugs=args.slug,
            series=args.series,
            level=args.level,
            raw_dir=args.raw,
            out_dir=args.out,
            boot=args.boot,
            write_outputs=False,
        )
        elapsed = round(time.time() - t0, 2)
        if not diffs:
            print(f"local UQ recomputation ({len(recomputed)} question(s), {elapsed}s): equal after 6-significant-digit rounding")
        else:
            print(f"local UQ recomputation ({len(recomputed)} question(s), {elapsed}s): {len(diffs)} difference(s) vs stored summaries:")
            for d in diffs:
                print(" -", d)
        return 0

    t0 = time.time()
    if args.raw is not None and args.raw.exists() and any(args.raw.glob("*.json")):
        gather(args.raw, args.out, boot=args.boot)
    build_questions(uq_dir=args.out, out_dir=args.questions_out)
    write_numbers_tex(args.questions_out)
    build_figures()
    manifest_path = write_manifest(args.out, args.questions_out)
    elapsed = round(time.time() - t0, 1)
    print(f"v0.8 reproduction complete in {elapsed}s -> {manifest_path.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
