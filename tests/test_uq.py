from __future__ import annotations

import numpy as np

from eq_biodt.v03.uq import (
    analyse,
    analyse_function,
    design,
    ishigami,
    level_of,
    run_level,
    task_seed,
    uniform_inputs,
)
from eq_biodt.v03.uq_inputs import inputs_for

INPUTS = inputs_for()


def test_ishigami_indices():
    # analytic: S1 = [0.3139, 0.4424, 0.0], ST = [0.5576, 0.4424, 0.2437]
    res = analyse_function(ishigami, uniform_inputs(3, -np.pi, np.pi), n=4096, seed=1)
    assert np.allclose([res["S1"][k][0] for k in "x1 x2 x3".split()], [0.3139, 0.4424, 0.0], atol=0.03)
    assert np.allclose([res["ST"][k][0] for k in "x1 x2 x3".split()], [0.5576, 0.4424, 0.2437], atol=0.03)


def test_design_deterministic():
    a = design(INPUTS, 64, task_seed("q|s|0")); b = design(INPUTS, 64, task_seed("q|s|0"))
    assert np.array_equal(a, b)


def test_closed_gate_is_exact_zero():
    spec = level_of("carbon-fixation", sample="14PRr2")
    res = run_level(spec, n=32)
    out = res["outputs"]["c_fixed_upper_bound_g_m2_yr"]
    assert out["mode"] == "deduced" and out["q05"] == out["q95"] == 0.0


def test_refusal_mode_switch():
    res = analyse(INPUTS, 32, [{"y": None, "__refused__": "x"}]*10 + [{"y": 1.0}]*(32*(len(INPUTS)+5)-10), ["y"])
    assert res["y"]["refusal_fraction"] < 0.05 or res["y"]["mode"] == "refusal_indicator"
    assert res["y"]["refusal_fraction"] >= 0.05 and res["y"]["mode"] == "refusal_indicator"
