"""Constitutive maps on observations; these maps never change causal leaves."""
from __future__ import annotations

import math
from dataclasses import dataclass


@dataclass(frozen=True)
class Observation:
    value: float | None
    unit: str
    denominator: str
    source: str
    endpoint: str
    category: str = ""

    def __post_init__(self):
        if not all((self.unit, self.denominator, self.source, self.endpoint)):
            raise ValueError("An observation requires units, denominator, source and endpoint")
        if self.value is not None and (not math.isfinite(self.value) or self.value < 0):
            raise ValueError("Observation must be nonnegative and finite, or missing")


def additive_total(parts: list[Observation], *, disjoint: bool) -> Observation:
    """Add disjoint measures with identical units, endpoint and denominator.

    Disjointness is an explicit caller assertion, never inferred from taxon labels.
    Missing is propagated, including for an empty collection.
    """
    if not parts or not disjoint:
        raise ValueError("An additive total requires identified disjoint parts")
    if any(not p.category for p in parts) or len({p.category for p in parts}) != len(parts):
        raise ValueError("Disjoint parts require distinct explicit category identifiers")
    first = parts[0]
    if any((p.unit, p.denominator, p.endpoint) !=
           (first.unit, first.denominator, first.endpoint) for p in parts):
        raise ValueError("Incompatible observation measures")
    value = None if any(p.value is None for p in parts) else sum(p.value for p in parts)
    return Observation(value, first.unit, first.denominator,
                       ";".join(sorted({p.source for p in parts})), first.endpoint,
                       "union("+",".join(sorted(p.category for p in parts))+")")


def carrier_partition(weights: dict[str, float], members: set[str]) -> dict:
    """An observed-reference partition; the complement is exactly the other MAGs.

    Overlapping carrier sets can each be measured but are never added as disjoint.
    """
    if any(not math.isfinite(v) or v < 0 for v in weights.values()):
        raise ValueError("Invalid mapped-read weight")
    total = sum(weights.values())
    selected = sum(v for k, v in weights.items() if k in members)
    return {"selected": selected, "complement": total-selected, "total": total,
            "share": selected/total if total > 0 else None}


def validate_composition(parts: dict[str, float | None], total: float = 1.) -> dict:
    """Report consistency without imposing composition priors or imputing leaves."""
    if not math.isfinite(total) or total < 0:
        raise ValueError("Invalid total")
    known = [v for v in parts.values() if v is not None]
    if any(not math.isfinite(v) or v < 0 for v in known):
        raise ValueError("Invalid composition")
    if sum(known) > total + 1e-9:
        raise ValueError("Parts exceed their shared total")
    return {"parts": dict(parts), "unassigned": total-sum(known),
            "complete": all(v is not None for v in parts.values()),
            "status": "partial" if any(v is None for v in parts.values()) else "consistent"}
