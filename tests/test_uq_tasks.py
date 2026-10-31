from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from reproducibility.v08.gather import gather_bytes
from reproducibility.v08.tasks import enumerate_tasks, merge_chunks
from reproducibility.v08.uq_task import run_task


def test_chunked_equals_single(tmp_path):
    tasks = [t for t in enumerate_tasks(n=16, chunk_rows=40) if t["slug"] == "biochar" and t["level"] == 1]
    for i, _ in enumerate(tasks): run_task(tasks, i, tmp_path/"chunks")
    run_task(merge_chunks(tasks), 0, tmp_path/"single")
    assert gather_bytes(tmp_path/"chunks") == gather_bytes(tmp_path/"single")
