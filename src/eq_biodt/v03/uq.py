"""Saltelli / Jansen quasi-Monte Carlo uncertainty and sensitivity analysis for EQ-BioDT v0.8.

Provides:
- `task_seed(key)` — deterministic 64-bit integer seed from sha256(key).
- `design(inputs, n, seed)` — scrambled Sobol' (2d) Saltelli design with individual and status-group
  re-sampling blocks `[A, B, AB_1..AB_d, AB_G1..AB_G3]`, shape `(n * (d + 5), d)`.
- `evaluate_rows(spec, inputs, rows)` — evaluate a `LevelSpec` on each design row.
- `analyse(inputs, n, outputs, headline, boot=500, seed=0)` — joint intervals, Saltelli (2010) first-order
  and Jansen (1999) total-order Sobol' indices (individual and grouped by audit status), refusal handling
  (< 5 % variance mode, >= 5 % binary refusal-indicator mode, closed-gate exact deductive zero), and
  N/2 vs N quantile convergence diagnostics.
"""
from __future__ import annotations

from dataclasses import asdict
import hashlib
from typing import Callable

import numpy as np
from scipy.stats import qmc

from .external_models import MODELS
from .questions import QUESTIONS, LevelSpec, evaluate
from .uq_inputs import Input, inputs_for

GROUP_STATUSES: tuple[str, ...] = ("verified_in_source", "modeller_assumption", "source_not_checked")


def inputs_for_spec(spec: LevelSpec) -> list[Input]:
    """Return the sorted `Input` list for a `LevelSpec` (`inputs_for(constant_ids, spec.intervened_params)`)."""
    cids = [c for mid in (*spec.upstream, *spec.downstream) for c in MODELS[mid].constants]
    return inputs_for(cids, spec.intervened_params)


def task_seed(key: str) -> int:
    """Deterministic non-negative 64-bit seed from the first 8 bytes of sha256(key)."""
    return int.from_bytes(hashlib.sha256(key.encode("utf-8")).digest()[:8], "big")


def _icdf(u: np.ndarray, inp: Input) -> np.ndarray:
    """Map uniform u in [0, 1) through the inverse CDF of `inp`."""
    lo, mode, hi = float(inp.low), float(inp.value), float(inp.high)
    if hi == lo:
        return np.full_like(u, lo, dtype=float)
    if inp.dist == "uniform":
        return lo + u * (hi - lo)
    if inp.dist == "loguniform":
        if lo > 0.0 and hi > 0.0:
            return np.exp(np.log(lo) + u * (np.log(hi) - np.log(lo)))
        if lo < 0.0 and hi < 0.0:
            a, b = abs(hi), abs(lo)  # 0 < a <= b
            return -np.exp(np.log(a) + (1.0 - u) * (np.log(b) - np.log(a)))
        raise ValueError(f"loguniform requires bounds of the same strict sign, got [{lo}, {hi}] for {inp.name}")
    if inp.dist == "triangular":
        c = min(max(mode, lo), hi)
        span = hi - lo
        fc = (c - lo) / span
        out = np.empty_like(u, dtype=float)
        left = u <= fc
        if c > lo:
            out[left] = lo + np.sqrt(np.maximum(0.0, u[left] * span * (c - lo)))
        else:
            out[left] = lo
        if hi > c:
            out[~left] = hi - np.sqrt(np.maximum(0.0, (1.0 - u[~left]) * span * (hi - c)))
        else:
            out[~left] = hi
        return out
    raise ValueError(f"Unknown distribution {inp.dist!r} for input {inp.name!r}")


def design(inputs: list[Input], n: int, seed: int) -> np.ndarray:
    """Build the Saltelli (2010) + group-Sobol' design matrix of shape `(n * (d + 5), d)`.

    Blocks of `n` rows each:
    - Block 0: A
    - Block 1: B
    - Blocks 2 .. d+1: AB_j (A with column j from B)
    - Blocks d+2 .. d+4: AB_G for each status in GROUP_STATUSES (A with all columns of that status from B)
    """
    d = len(inputs)
    if d == 0:
        return np.empty((5 * n, 0), dtype=float)
    sampler = qmc.Sobol(d=2 * d, scramble=True, seed=seed)
    u = sampler.random(n)
    u_a, u_b = u[:, :d], u[:, d:]
    a = np.column_stack([_icdf(u_a[:, j], inp) for j, inp in enumerate(inputs)])
    b = np.column_stack([_icdf(u_b[:, j], inp) for j, inp in enumerate(inputs)])

    blocks = [a, b]
    for j in range(d):
        ab_j = a.copy()
        ab_j[:, j] = b[:, j]
        blocks.append(ab_j)
    for status in GROUP_STATUSES:
        ab_g = a.copy()
        for j, inp in enumerate(inputs):
            if inp.status == status:
                ab_g[:, j] = b[:, j]
        blocks.append(ab_g)
    return np.vstack(blocks)


def evaluate_rows(spec: LevelSpec, inputs: list[Input], rows: np.ndarray) -> list[dict]:
    """Evaluate `spec` on each row of `rows` (`shape = (m, len(inputs))`)."""
    const_idx = [(j, inp.name) for j, inp in enumerate(inputs) if inp.kind == "constant"]
    param_idx = [(j, inp.name) for j, inp in enumerate(inputs) if inp.kind == "parameter"]
    results: list[dict] = []
    for row in rows:
        consts = {name: float(row[j]) for j, name in const_idx}
        params = {name: float(row[j]) for j, name in param_idx}
        results.append(evaluate(spec, constants=consts, params=params, keep_twin_result=False))
    return results


def _sobol_block_estimates(
    mat: np.ndarray,
    d: int,
    names: list[str],
    boot: int,
    seed: int,
) -> tuple[dict[str, list[float]], dict[str, list[float]], dict[str, list[float]], dict[str, list[float]]]:
    """Compute Saltelli 2010 S1 and Jansen 1999 ST with bootstrap 95% CIs on `mat` of shape `(d + 5, m)`."""
    zero_s1 = {name: [0.0, 0.0, 0.0] for name in names}
    zero_st = {name: [0.0, 0.0, 0.0] for name in names}
    zero_gs1 = {g: [0.0, 0.0, 0.0] for g in GROUP_STATUSES}
    zero_gst = {g: [0.0, 0.0, 0.0] for g in GROUP_STATUSES}

    m = mat.shape[1]
    if m < 2:
        return zero_s1, zero_st, zero_gs1, zero_gst

    ya = mat[0]
    yb = mat[1]
    y_ab = mat[2:]  # shape (d + 3, m)
    var_y = float(np.var(np.concatenate([ya, yb]), ddof=0))
    if not np.isfinite(var_y) or var_y <= 1e-30:
        return zero_s1, zero_st, zero_gs1, zero_gst

    diff = y_ab - ya[None, :]  # shape (d + 3, m)
    term_s1 = yb[None, :] * diff
    term_st = 0.5 * (diff * diff)

    s1_est = np.mean(term_s1, axis=1) / var_y
    st_est = np.mean(term_st, axis=1) / var_y

    if boot > 0:
        rng = np.random.default_rng(seed)
        idx = rng.integers(0, m, size=(boot, m))
        counts = np.zeros((boot, m), dtype=float)
        rows_grid = np.broadcast_to(np.arange(boot)[:, None], (boot, m))
        np.add.at(counts, (rows_grid, idx), 1.0 / m)

        mean_ab = 0.5 * (counts @ (ya + yb))
        var_b = 0.5 * (counts @ (ya * ya + yb * yb)) - mean_ab * mean_ab
        valid_b = var_b > 1e-30
        if np.any(valid_b):
            counts_v = counts[valid_b]
            var_v = var_b[valid_b][:, None]
            s1_boot = (counts_v @ term_s1.T) / var_v
            st_boot = (counts_v @ term_st.T) / var_v
            s1_lo, s1_hi = np.quantile(s1_boot, [0.025, 0.975], axis=0)
            st_lo, st_hi = np.quantile(st_boot, [0.025, 0.975], axis=0)
        else:
            s1_lo = s1_hi = s1_est
            st_lo = st_hi = st_est
    else:
        s1_lo = s1_hi = s1_est
        st_lo = st_hi = st_est

    s1 = {names[j]: [float(s1_est[j]), float(s1_lo[j]), float(s1_hi[j])] for j in range(d)}
    st = {names[j]: [float(st_est[j]), float(st_lo[j]), float(st_hi[j])] for j in range(d)}
    gs1 = {
        g: [float(s1_est[d + k]), float(s1_lo[d + k]), float(s1_hi[d + k])]
        for k, g in enumerate(GROUP_STATUSES)
    }
    gst = {
        g: [float(st_est[d + k]), float(st_lo[d + k]), float(st_hi[d + k])]
        for k, g in enumerate(GROUP_STATUSES)
    }
    return s1, st, gs1, gst


def analyse(
    inputs: list[Input],
    n: int,
    outputs: list[dict],
    headline: list[str],
    boot: int = 500,
    seed: int = 0,
) -> dict:
    """Analyse the `n * (d + 5)` output rows produced on `design(inputs, n, ...)`.

    Returns `{output_name: {median, q05, q25, q75, q95, refusal_fraction, mode, S1, ST, group_S1, group_ST, convergence}}`.
    """
    d = len(inputs)
    names = [inp.name for inp in inputs]
    expected = n * (d + 5)
    if len(outputs) != expected:
        raise ValueError(f"Expected {expected} output rows for n={n}, d={d}, got {len(outputs)}")

    result: dict[str, dict] = {}
    for key in headline:
        deduced_rows = sum(1 for r in outputs if key in r.get("__deduced__", ()))
        if deduced_rows == expected or (
            deduced_rows > 0 and all(r.get(key) in (0.0, None) for r in outputs)
        ):
            result[key] = {
                "median": 0.0,
                "q05": 0.0,
                "q25": 0.0,
                "q75": 0.0,
                "q95": 0.0,
                "frac_neg": 0.0,
                "frac_zero": 1.0,
                "frac_pos": 0.0,
                "refusal_fraction": 0.0,
                "mode": "deduced",
                "S1": {name: [0.0, 0.0, 0.0] for name in names},
                "ST": {name: [0.0, 0.0, 0.0] for name in names},
                "group_S1": {g: [0.0, 0.0, 0.0] for g in GROUP_STATUSES},
                "group_ST": {g: [0.0, 0.0, 0.0] for g in GROUP_STATUSES},
                "convergence": {"half_n_shift_frac": 0.0, "flagged": False},
            }
            continue

        raw_blocks = [
            [outputs[b * n + i].get(key) for i in range(n)]
            for b in range(d + 5)
        ]
        base_vals = raw_blocks[0] + raw_blocks[1]
        refused_base = sum(1 for v in base_vals if v is None)
        refusal_fraction = float(refused_base / len(base_vals)) if base_vals else 0.0

        admitted_full = np.array([float(v) for v in base_vals if v is not None], dtype=float)
        half = max(1, n // 2)
        base_half = raw_blocks[0][:half] + raw_blocks[1][:half]
        admitted_half = np.array([float(v) for v in base_half if v is not None], dtype=float)

        if admitted_full.size > 0:
            q05, q25, median, q75, q95 = [
                float(x) for x in np.quantile(admitted_full, [0.05, 0.25, 0.50, 0.75, 0.95])
            ]
            if admitted_half.size > 0:
                q05_h, q95_h = [float(x) for x in np.quantile(admitted_half, [0.05, 0.95])]
                span = q95 - q05
                shift = float(max(abs(q05_h - q05), abs(q95_h - q95)) / span) if span > 1e-15 else 0.0
            else:
                shift = 0.0
            discrete = bool(np.all(admitted_full == np.round(admitted_full)))
            # integer-valued outputs (e.g. active days) have quantiles that jump by whole units; the 5 % rule does not apply
            convergence = {"half_n_shift_frac": shift, "discrete": discrete, "flagged": bool(shift > 0.05 and not discrete)}
        else:
            median = q05 = q25 = q75 = q95 = None
            convergence = {"half_n_shift_frac": 0.0, "flagged": False}

        if refusal_fraction >= 0.05:
            mode = "refusal_indicator"
            ind_mat = np.array(
                [[1.0 if v is None else 0.0 for v in block] for block in raw_blocks],
                dtype=float,
            )
            s1, st, gs1, gst = _sobol_block_estimates(ind_mat, d, names, boot, seed)
        elif admitted_full.size > 0 and float(np.max(admitted_full) - np.min(admitted_full)) <= 1e-9 * max(
            1.0, abs(float(np.mean(admitted_full)))
        ):
            mode = "constant"
            s1 = {}
            st = {}
            gs1 = None
            gst = None
        else:
            mode = "variance"
            valid_cols = [
                i for i in range(n) if all(raw_blocks[b][i] is not None for b in range(d + 5))
            ]
            if len(valid_cols) >= 2:
                val_mat = np.array(
                    [[float(raw_blocks[b][i]) for i in valid_cols] for b in range(d + 5)],
                    dtype=float,
                )
                s1, st, gs1, gst = _sobol_block_estimates(val_mat, d, names, boot, seed)
            else:
                s1 = {name: [0.0, 0.0, 0.0] for name in names}
                st = {name: [0.0, 0.0, 0.0] for name in names}
                gs1 = {g: [0.0, 0.0, 0.0] for g in GROUP_STATUSES}
                gst = {g: [0.0, 0.0, 0.0] for g in GROUP_STATUSES}

        if admitted_full.size == 0:
            frac_neg = frac_zero = frac_pos = None
        elif mode == "constant":
            c_val = float(median)
            frac_neg = 1.0 if c_val < 0.0 else 0.0
            frac_zero = 1.0 if c_val == 0.0 else 0.0
            frac_pos = 1.0 if c_val > 0.0 else 0.0
        else:
            frac_neg = float(f"{float(np.mean(admitted_full < 0.0)):.6g}")
            frac_zero = float(f"{float(np.mean(admitted_full == 0.0)):.6g}")
            frac_pos = float(f"{float(np.mean(admitted_full > 0.0)):.6g}")

        result[key] = {
            "median": median,
            "q05": q05,
            "q25": q25,
            "q75": q75,
            "q95": q95,
            "frac_neg": frac_neg,
            "frac_zero": frac_zero,
            "frac_pos": frac_pos,
            "refusal_fraction": refusal_fraction,
            "mode": mode,
            "S1": s1,
            "ST": st,
            "group_S1": gs1,
            "group_ST": gst,
            "convergence": convergence,
        }
    return result


# ---------------------------------------------------------------- test & convenience helpers


def uniform_inputs(d: int, low: float, high: float) -> list[Input]:
    """Create `d` synthetic uniform inputs named `x1 .. xd` for benchmark tests."""
    mid = 0.5 * (low + high)
    return [
        Input(
            name=f"x{j + 1}",
            kind="constant",
            dist="uniform",
            low=float(low),
            value=float(mid),
            high=float(high),
            status=GROUP_STATUSES[j % len(GROUP_STATUSES)],
        )
        for j in range(d)
    ]


def ishigami(x: np.ndarray, a: float = 7.0, b: float = 0.1) -> float:
    """Standard Ishigami (1990) test function on `[-pi, pi]^3`."""
    return float(np.sin(x[0]) + a * (np.sin(x[1]) ** 2) + b * (x[2] ** 4) * np.sin(x[0]))


def analyse_function(
    fn: Callable[[np.ndarray], float],
    inputs: list[Input],
    n: int = 1024,
    boot: int = 500,
    seed: int = 0,
) -> dict:
    """Run `design` and `analyse` on a pure scalar function `fn(row) -> float`."""
    rows = design(inputs, n, seed)
    outputs = [{"y": float(fn(row)), "__refused__": None, "__deduced__": []} for row in rows]
    return analyse(inputs, n, outputs, ["y"], boot=boot, seed=seed)["y"]


def level_of(
    slug: str,
    *,
    series: str | int = 0,
    x: float | str | None = None,
    sample: str | None = None,
) -> LevelSpec:
    """Locate a `LevelSpec` in `QUESTIONS[slug]` by series index/id and `x` or `sample`."""
    q = QUESTIONS[slug]
    s = q.series[series] if isinstance(series, int) else next(sr for sr in q.series if sr.id == series)
    if sample is not None:
        return next(l for l in s.levels if l.sample == sample or l.x == sample)
    if x is not None:
        return next(l for l in s.levels if l.x == x)
    return s.levels[0]


def run_level(
    spec: LevelSpec,
    n: int = 32,
    headline: list[str] | None = None,
    boot: int = 200,
    seed: int = 0,
) -> dict:
    """Design, evaluate, and analyse a single `LevelSpec` at sample size `n`."""
    inps = inputs_for_spec(spec)
    rows = design(inps, n, seed)
    raw = evaluate_rows(spec, inps, rows)
    hl = headline or [k for k in raw[0] if not k.startswith("__")]
    return {
        "n": n,
        "inputs": [asdict(i) for i in inps],
        "outputs": analyse(inps, n, raw, hl, boot=boot, seed=seed),
    }
