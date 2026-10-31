"""Parameter-uncertainty input space for joint UQ and global sensitivity (v0.8).

Combines literature constants used by a chain (triangular over [low, value, high] from
data/external/composition_constants.json) with non-intervened twin process parameters (from
data/external/twin_parameter_ranges.json).
"""
from __future__ import annotations

from dataclasses import dataclass
import hashlib
import json
from pathlib import Path
from typing import Iterable, Literal

from .composition import CompositionError, load_constants
from .process import Parameters

PARAM_REGISTRY_PATH = Path(__file__).resolve().parents[3] / "data" / "external" / "twin_parameter_ranges.json"
DISTS = {"triangular", "uniform", "loguniform"}
STATUSES = {"verified_in_source", "modeller_assumption", "source_not_checked"}
SOURCE_KEYS = {"authors", "year", "title", "doi_or_url", "page_or_table"}

_PARAM_CACHE: dict[str, dict[str, dict]] = {}


@dataclass(frozen=True)
class Input:
    name: str
    kind: Literal["constant", "parameter"]
    dist: str
    low: float
    value: float
    high: float
    status: str


def load_parameter_ranges(path: Path | None = None) -> dict[str, dict]:
    """Twin parameter range entries by id. Validates schema, bounds and reference values."""
    path = Path(path or PARAM_REGISTRY_PATH)
    key = str(path)
    if key not in _PARAM_CACHE:
        doc = json.loads(path.read_text())
        ref = Parameters(k_h2_s=.001, k_co_s=.001)
        out: dict[str, dict] = {}
        for e in doc["parameters"]:
            missing = {"id", "dist", "low", "value", "high", "unit", "meaning", "status", "source", "note"} - set(e)
            if missing:
                raise CompositionError(f"parameter range entry {e.get('id')!r} lacks {sorted(missing)}")
            pid = e["id"]
            if pid in out:
                raise CompositionError(f"duplicate parameter range id {pid!r}")
            if e["dist"] not in DISTS:
                raise CompositionError(f"parameter {pid!r}: unknown dist {e['dist']!r}")
            if e["status"] not in STATUSES:
                raise CompositionError(f"parameter {pid!r}: unknown status {e['status']!r}")
            if set(e["source"]) != SOURCE_KEYS:
                raise CompositionError(f"parameter {pid!r}: invalid source keys {sorted(e['source'])}")
            if not e["low"] <= e["value"] <= e["high"]:
                raise CompositionError(f"parameter {pid!r}: need low <= value <= high")
            if e["dist"] == "loguniform" and (e["low"] == 0 or e["high"] == 0 or (e["low"] < 0) != (e["high"] < 0)):
                raise CompositionError(f"parameter {pid!r}: loguniform bounds must have the same non-zero sign")
            if getattr(ref, pid) != e["value"]:
                raise CompositionError(f"parameter {pid!r}: value {e['value']} != reference {getattr(ref, pid)}")
            out[pid] = e
        _PARAM_CACHE[key] = out
    return _PARAM_CACHE[key]


def parameter_ranges_sha256(path: Path | None = None) -> str:
    return hashlib.sha256(Path(path or PARAM_REGISTRY_PATH).read_bytes()).hexdigest()


def inputs_for(constant_ids: Iterable[str] = (), exclude_params: set[str] | frozenset[str] = frozenset()) -> list[Input]:
    """Uncertain inputs (constants used by the chain + non-intervened twin parameters), sorted by name."""
    consts = load_constants()
    params = load_parameter_ranges()
    excl = set(exclude_params or ())
    items: list[Input] = []
    for cid in sorted(set(constant_ids)):
        if cid not in consts:
            raise CompositionError(f"unknown constant id {cid!r}")
        e = consts[cid]
        if e["low"] == e["high"]:
            continue
        items.append(Input(
            name=cid,
            kind="constant",
            dist="triangular",
            low=float(e["low"]),
            value=float(e["value"]),
            high=float(e["high"]),
            status=e.get("status", "source_not_checked"),
        ))
    for pid, e in params.items():
        if pid in excl or e["low"] == e["high"]:
            continue
        items.append(Input(
            name=pid,
            kind="parameter",
            dist=e["dist"],
            low=float(e["low"]),
            value=float(e["value"]),
            high=float(e["high"]),
            status=e["status"],
        ))
    return sorted(items, key=lambda x: x.name)
