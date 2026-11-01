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
    assert out["frac_neg"] == 0.0 and out["frac_zero"] == 1.0 and out["frac_pos"] == 0.0


def test_refusal_mode_switch():
    res = analyse(INPUTS, 32, [{"y": None, "__refused__": "x"}]*10 + [{"y": 1.0}]*(32*(len(INPUTS)+5)-10), ["y"])
    assert res["y"]["refusal_fraction"] < 0.05 or res["y"]["mode"] == "refusal_indicator"
    assert res["y"]["refusal_fraction"] >= 0.05 and res["y"]["mode"] == "refusal_indicator"


def test_constant_output_mode_has_no_indices():
    total = 32 * (len(INPUTS) + 5)
    res = analyse(INPUTS, 32, [{"y": -50.0, "z": 0.0, "p": 25.0}] * total, ["y", "z", "p"])
    for key, val in [("y", -50.0), ("z", 0.0), ("p", 25.0)]:
        assert res[key]["mode"] == "constant"
        assert res[key]["median"] == val
        assert res[key]["S1"] == {}
        assert res[key]["ST"] == {}
        assert res[key]["group_S1"] is None
        assert res[key]["group_ST"] is None
    assert (res["y"]["frac_neg"], res["y"]["frac_zero"], res["y"]["frac_pos"]) == (1.0, 0.0, 0.0)
    assert (res["z"]["frac_neg"], res["z"]["frac_zero"], res["z"]["frac_pos"]) == (0.0, 1.0, 0.0)
    assert (res["p"]["frac_neg"], res["p"]["frac_zero"], res["p"]["frac_pos"]) == (0.0, 0.0, 1.0)


def test_direction_fractions_sum_to_one_and_match_raw_recount():
    import json
    from pathlib import Path
    import sys

    root = Path(__file__).resolve().parents[1]
    if str(root) not in sys.path:
        sys.path.insert(0, str(root))
    from reproducibility.v08.gather import _add_relative_changes

    uq_dir = root / "reviews" / "v08" / "uq"
    q_dir = root / "reviews" / "v08" / "questions"

    for p in sorted(uq_dir.glob("*.json")):
        if p.name == "run_record.json":
            continue
        doc = json.loads(p.read_text(encoding="utf-8"))
        for s in doc["series"]:
            for lvl in s["levels"]:
                for out_key, out in lvl["outputs"].items():
                    fn, fz, fp = out["frac_neg"], out["frac_zero"], out["frac_pos"]
                    if out["median"] is None:
                        assert fn is None and fz is None and fp is None
                    else:
                        assert 0.0 <= fn <= 1.0 and 0.0 <= fz <= 1.0 and 0.0 <= fp <= 1.0
                        assert abs((fn + fz + fp) - 1.0) < 1e-5, (p.name, s["id"], lvl["level"], out_key, fn, fz, fp)

    for p in sorted(q_dir.glob("*.json")):
        if p.name == "index.json":
            continue
        doc = json.loads(p.read_text(encoding="utf-8"))
        for s in doc["series"]:
            for lvl in s["levels"]:
                for out_key, out in lvl["outcomes"].items():
                    fn, fz, fp = out["frac_neg"], out["frac_zero"], out["frac_pos"]
                    if out["median"] is None:
                        assert fn is None and fz is None and fp is None
                    else:
                        assert abs((fn + fz + fp) - 1.0) < 1e-5, (p.name, s["id"], lvl["label"], out_key)

    raw_dir = Path("/mnt/data1/agents/projects/eq-biodt-v07/uq-raw-merged2")
    if raw_dir.exists():
        chunks = [
            json.loads(fp.read_text(encoding="utf-8"))
            for fp in sorted(raw_dir.glob("rain-pattern|rain-schedule|2|*.json"))
        ]
        assert chunks, "Expected raw chunks for rain-pattern|rain-schedule|2"
        chunks.sort(key=lambda c: c["rows"][0])
        n = int(chunks[0]["n"])
        rows: list[dict] = []
        for c in chunks:
            rows.extend(c["outputs"])
        _add_relative_changes(rows)
        base_rows = rows[: 2 * n]

        uq_rp = json.loads((uq_dir / "rain-pattern.json").read_text(encoding="utf-8"))
        lvl_uq = uq_rp["series"][0]["levels"][2]
        q_rp = json.loads((q_dir / "rain-pattern.json").read_text(encoding="utf-8"))
        lvl_q = q_rp["series"][0]["levels"][2]

        for out_key, out_uq in lvl_uq["outputs"].items():
            vals = np.array([float(r[out_key]) for r in base_rows if r.get(out_key) is not None], dtype=float)
            if out_uq["mode"] == "constant":
                c_val = float(out_uq["median"])
                exp_neg = 1.0 if c_val < 0.0 else 0.0
                exp_zero = 1.0 if c_val == 0.0 else 0.0
                exp_pos = 1.0 if c_val > 0.0 else 0.0
            else:
                exp_neg = float(f"{float(np.mean(vals < 0.0)):.6g}")
                exp_zero = float(f"{float(np.mean(vals == 0.0)):.6g}")
                exp_pos = float(f"{float(np.mean(vals > 0.0)):.6g}")
            assert out_uq["frac_neg"] == exp_neg, out_key
            assert out_uq["frac_zero"] == exp_zero, out_key
            assert out_uq["frac_pos"] == exp_pos, out_key
            if out_key in lvl_q["outcomes"]:
                assert lvl_q["outcomes"][out_key]["frac_neg"] == exp_neg
                assert lvl_q["outcomes"][out_key]["frac_zero"] == exp_zero
                assert lvl_q["outcomes"][out_key]["frac_pos"] == exp_pos
        assert lvl_uq["outputs"]["h2_uptake_change_pct"]["frac_neg"] == 0.963379
        assert lvl_uq["outputs"]["h2_uptake_change_pct"]["frac_pos"] == 0.0366211


