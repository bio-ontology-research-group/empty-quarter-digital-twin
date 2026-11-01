from __future__ import annotations

import copy
import json
from pathlib import Path
import re
import sys

import pytest

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from reproducibility.v08.questions import _augment_central_pct
from reproducibility.v08.science_texts import (
    READABLE_INPUT_NAMES,
    SCIENCE_TEMPLATES,
    _SLOT_RE,
    _fmt_pct_mag,
    _fmt_plain,
    _fmt_sig,
    _fmt_signed_pct,
    _fmt_signed_pct_endpoint,
    _load_registry_constants,
    render_science_texts,
    resolve_slot,
)
from reproducibility.v08.texts import (
    HEADLINE_SPEC,
    plain_number,
    render_question_texts,
    round_plain,
    science_central,
)

UQ_DIR = ROOT / "reviews" / "v08" / "uq"
MATH_SUBSCRIPTS = ("S_T", "γ_CO", "k_H2", "k_CO", "ψ_e")


def _load_uq(slug: str) -> dict:
    doc = json.loads((UQ_DIR / f"{slug}.json").read_text(encoding="utf-8"))
    for s in doc.get("series", []):
        for l in s.get("levels", []):
            if "central" in l and isinstance(l["central"], dict):
                _augment_central_pct(l["central"])
    return doc


def test_render_all_14_questions_succeeds():
    assert len(SCIENCE_TEMPLATES) == 14
    for slug in SCIENCE_TEMPLATES:
        uq_doc = _load_uq(slug)
        rendered = render_science_texts(slug, uq_doc)
        for field in ("title_science", "answer_science", "limits_science", "headline_science"):
            text = rendered[field]
            assert text and isinstance(text, str), f"{slug}.{field} is empty"
            assert "{" not in text and "}" not in text, f"Unresolved braces in {slug}.{field}: {text!r}"


def test_no_snake_case_and_no_raw_underscores():
    reg = _load_registry_constants()
    for key in reg:
        assert key in READABLE_INPUT_NAMES, f"Missing readable input name for {key!r}"
    for key, readable in READABLE_INPUT_NAMES.items():
        assert not re.search(r"\b[a-z][a-z0-9]*_[a-z0-9_]+\b", readable), (
            f"snake_case token in READABLE_INPUT_NAMES[{key!r}] = {readable!r}"
        )
        cleaned = readable
        for sym in MATH_SUBSCRIPTS:
            cleaned = cleaned.replace(sym, "")
        assert "_" not in cleaned, f"Unexpected '_' in READABLE_INPUT_NAMES[{key!r}] = {readable!r}"

    for slug in SCIENCE_TEMPLATES:
        uq_doc = _load_uq(slug)
        rendered = render_science_texts(slug, uq_doc)
        for field in ("title_science", "answer_science", "limits_science", "headline_science"):
            text = rendered[field]
            assert not re.search(r"\b[a-z][a-z0-9]*_[a-z0-9_]+\b", text), (
                f"snake_case token in {slug}.{field}: {text!r}"
            )
            cleaned = text
            for sym in MATH_SUBSCRIPTS:
                cleaned = cleaned.replace(sym, "")
            assert "_" not in cleaned, f"Unexpected '_' in {slug}.{field}: {text!r}"


def test_constant_outputs_have_no_indices_and_are_not_reported_with_indices():
    for slug in SCIENCE_TEMPLATES:
        uq_doc = _load_uq(slug)
        constant_targets: set[tuple[str, int, str]] = set()
        for s in uq_doc["series"]:
            for l in s["levels"]:
                l_idx = int(l["level"])
                for out_key, out_doc in l["outputs"].items():
                    if out_doc.get("mode") == "constant":
                        assert out_doc["S1"] == {}, f"{slug} {s['id']}/{l_idx}:{out_key} has non-empty S1"
                        assert out_doc["ST"] == {}, f"{slug} {s['id']}/{l_idx}:{out_key} has non-empty ST"
                        assert out_doc["group_S1"] is None, f"{slug} {s['id']}/{l_idx}:{out_key} group_S1 not None"
                        assert out_doc["group_ST"] is None, f"{slug} {s['id']}/{l_idx}:{out_key} group_ST not None"
                        constant_targets.add((s["id"], l_idx, out_key))

        tmpl = SCIENCE_TEMPLATES[slug]
        for field in ("title_science", "answer_science", "limits_science"):
            for expr in _SLOT_RE.findall(tmpl[field]):
                core = expr.split("|", 1)[0]
                if core.startswith(("const:", "special:")):
                    continue
                sl, out_key, stat = core.split(":")
                s_id, l_s = sl.split("/", 1)
                target = (s_id, int(l_s), out_key)
                if target in constant_targets:
                    assert stat not in {"st1", "st2", "st1name", "st2name", "assum"}, (
                        f"Constant output {target} queried for index stat {stat!r} in {slug}.{field}"
                    )

    ah_uq = _load_uq("air-hydrogen")
    with pytest.raises(ValueError, match="constant output"):
        resolve_slot("h2-multiplier/2:h2_uptake_change_pct:st1", ah_uq)
    with pytest.raises(ValueError, match="constant output"):
        resolve_slot("h2-multiplier/2:h2_uptake_change_pct:st1name", ah_uq)


def test_absneg_and_abspos_sign_assertions():
    for slug, tmpl in SCIENCE_TEMPLATES.items():
        uq_doc = _load_uq(slug)
        for field in ("title_science", "answer_science", "limits_science"):
            for expr in _SLOT_RE.findall(tmpl[field]):
                if "|absneg" not in expr and "|abspos" not in expr:
                    continue
                core, mod = expr.split("|", 1)
                sl, out_key, _ = core.split(":")
                s_id, l_s = sl.split("/", 1)
                l_idx = int(l_s)
                # Mutate each of the 4 guarded statistics to violate sign and verify ValueError
                for key in ("central", "median", "q05", "q95"):
                    bad_uq = copy.deepcopy(uq_doc)
                    s_doc = next(s for s in bad_uq["series"] if s["id"] == s_id)
                    l_doc = next(l for l in s_doc["levels"] if int(l["level"]) == l_idx)
                    violation = 1.0 if mod == "absneg" else -1.0
                    if key == "central":
                        l_doc["central"][out_key] = violation
                    else:
                        l_doc["outputs"][out_key][key] = violation
                    with pytest.raises(ValueError, match=f"\\|{mod} sign check failed"):
                        resolve_slot(expr, bad_uq)
                    with pytest.raises(ValueError, match=f"\\|{mod} sign check failed"):
                        render_science_texts(slug, bad_uq)


def test_word_limits_and_forbidden_phrases():
    for slug in SCIENCE_TEMPLATES:
        uq_doc = _load_uq(slug)
        rendered = render_science_texts(slug, uq_doc)
        ans_words = rendered["answer_science"].split()
        lim_words = rendered["limits_science"].split()
        assert len(ans_words) <= 170, f"{slug}.answer_science has {len(ans_words)} words (> 170)"
        assert len(lim_words) <= 90, f"{slug}.limits_science has {len(lim_words)} words (> 90)"

        for field in ("title_science", "answer_science", "limits_science"):
            text = rendered[field]
            assert "—" not in text, f"Em dash found in {slug}.{field}"
            assert not re.search(r"\bthus\b", text, flags=re.IGNORECASE), f"'thus' found in {slug}.{field}"
            assert "IQR" not in text, f"'IQR' found in {slug}.{field}"
            assert "µmol m⁻² (median" not in text, f"Old dense phrase found in {slug}.{field}"


def _independent_slot_value(expr: str, uq_doc: dict, reg: dict[str, float]) -> str:
    """Independently compute the expected string for `expr` directly from `uq_doc` and `reg`."""
    core, mod = expr.split("|", 1) if "|" in expr else (expr, None)
    if core.startswith("const:"):
        cid = core.split(":", 1)[1]
        assert mod == "plain"
        return _fmt_plain(reg[cid])
    if core.startswith("special:"):
        sname = core.split(":", 1)[1]
        levels = uq_doc["series"][0]["levels"]
        if sname == "deduced_samples":
            return ", ".join(
                str(l["x"]) for l in levels if l["outputs"]["co_uptake_umol_m2"]["mode"] == "deduced"
            )
        if sname == "open_samples":
            return ", ".join(
                str(l["x"]) for l in levels if l["outputs"]["co_uptake_umol_m2"]["mode"] == "variance"
            )
        if sname == "refused_samples":
            samps = [
                str(l["x"])
                for l in levels
                if float(l["outputs"]["co_uptake_umol_m2"]["refusal_fraction"]) == 1.0
            ]
            return f"{samps[0]} and {samps[1]}"
        raise AssertionError(f"Unexpected special slot {sname}")

    sl, out_key, stat = core.split(":")
    s_id, l_s = sl.split("/", 1)
    l_idx = int(l_s)
    s_doc = next(s for s in uq_doc["series"] if s["id"] == s_id)
    l_doc = next(l for l in s_doc["levels"] if int(l["level"]) == l_idx)
    out_doc = l_doc["outputs"][out_key]

    if stat in {"central", "median", "q05", "q95"}:
        val = float(l_doc["central"][out_key]) if stat == "central" else float(out_doc[stat])
        if mod is None:
            return _fmt_signed_pct(val) if out_key.endswith("_pct") else _fmt_plain(val)
        if mod == "absneg":
            return f"{_fmt_pct_mag(abs(val))} %"
        if mod == "abspos":
            return f"{_fmt_pct_mag(val)} %"
        if mod == "frac2pct":
            return f"{_fmt_pct_mag(100.0 * val)} %"
        if mod == "um":
            return f"{_fmt_sig(val, 3)} µmol m⁻²"
        if mod == "ug":
            return f"{_fmt_sig(val, 3)} µg C m⁻²"
        if mod == "g":
            return f"{_fmt_sig(val, 2)} g C m⁻² yr⁻¹"
        if mod == "plain":
            return _fmt_plain(val)
    if stat == "ci90":
        lo = float(out_doc["q05"])
        hi = float(out_doc["q95"])
        if mod is None:
            span = lo < 0.0 < hi
            return f"{_fmt_signed_pct_endpoint(lo, span)} to {_fmt_signed_pct_endpoint(hi, span)} %"
        if mod == "absneg":
            mags = sorted([abs(lo), abs(hi)])
            return f"{_fmt_pct_mag(mags[0])} to {_fmt_pct_mag(mags[1])} %"
        if mod == "abspos":
            vals = sorted([lo, hi])
            return f"{_fmt_pct_mag(vals[0])} to {_fmt_pct_mag(vals[1])} %"
        if mod == "frac2pct":
            vals = sorted([100.0 * lo, 100.0 * hi])
            return f"{_fmt_pct_mag(vals[0])} to {_fmt_pct_mag(vals[1])} %"
        if mod == "um":
            return f"{_fmt_sig(lo, 3)} to {_fmt_sig(hi, 3)} µmol m⁻²"
        if mod == "g":
            return f"{_fmt_sig(lo, 2)} to {_fmt_sig(hi, 2)} g C m⁻² yr⁻¹"
    if stat in {"st1", "st2", "st1name", "st2name"}:
        ranked = sorted(
            ((k, float(v[0])) for k, v in out_doc["ST"].items()),
            key=lambda kv: (-kv[1], kv[0]),
        )
        idx = 0 if stat in {"st1", "st1name"} else 1
        if stat in {"st1", "st2"}:
            return f"{ranked[idx][1]:.2f}"
        return READABLE_INPUT_NAMES[ranked[idx][0]]
    if stat == "refusal":
        return f"{_fmt_pct_mag(100.0 * float(out_doc['refusal_fraction']))} %"
    raise AssertionError(f"Unhandled slot {expr}")


def test_science_texts_contain_only_slot_resolved_numbers_matching_uq():
    reg = _load_registry_constants()
    # Allowed numeric literals in the static template prose (design facts, citations, chemical/model identifiers)
    allowed_static_numbers = {
        "0", "0.1", "0.40", "0.5", "1", "1.5", "2", "3", "4", "5", "8",
        "10", "12", "14", "20", "25", "29", "50", "51", "53", "56", "90", "187",
        "1970", "1983", "1999", "2014", "2016", "2018", "2020",
    }

    for slug, tmpl in SCIENCE_TEMPLATES.items():
        uq_doc = _load_uq(slug)
        rendered = render_science_texts(slug, uq_doc)

        for field in ("title_science", "answer_science", "limits_science"):
            raw_t = tmpl[field]
            # 1. Verify every slot resolves to the exact value from the UQ file / registry
            parts = _SLOT_RE.split(raw_t)
            reconstructed = []
            for i, part in enumerate(parts):
                if i % 2 == 0:
                    # Static template segment: check that any numbers are in allowed_static_numbers
                    for num_match in re.findall(r"\d+(?:\.\d+)?", part):
                        assert num_match in allowed_static_numbers, (
                            f"Untracked static number {num_match!r} in {slug}.{field}: {part!r}"
                        )
                    reconstructed.append(part)
                else:
                    expected_val = _independent_slot_value(part, uq_doc, reg)
                    actual_val = resolve_slot(part, uq_doc)
                    assert actual_val == expected_val, (
                        f"Slot {part!r} in {slug}.{field}: got {actual_val!r}, expected {expected_val!r}"
                    )
                    reconstructed.append(expected_val)

            assert "".join(reconstructed) == rendered[field], (
                f"Rendered {slug}.{field} does not match exact slot-resolved reconstruction"
            )


def test_plain_headline_number_and_mechanism_corrections():
    for slug, (s_id, l_idx, out_key) in HEADLINE_SPEC.items():
        uq_doc = _load_uq(slug)
        s_doc = next(s for s in uq_doc["series"] if s["id"] == s_id)
        l_doc = next(l for l in s_doc["levels"] if int(l["level"]) == l_idx)
        c_val = float(l_doc["central"][out_key])
        out_uq = l_doc["outputs"][out_key]
        txts = render_question_texts(slug, c_val, out_uq, 0.0, uq_doc)

        q_obj = {
            "slug": slug,
            "headline_value": c_val,
            "answer_plain": txts["answer_plain"],
            "answer_science": txts["answer_science"],
        }
        sci_c = science_central(q_obj)
        expected_plain = round_plain(abs(sci_c))
        assert plain_number(q_obj) == expected_plain, slug

        m_ref = re.search(
            r"reference parameterisation (?:gives )?[+−-]?(\d+(?:\.\d+)?)\s*%",
            txts["answer_science"],
        )
        if m_ref:
            assert m_ref.group(1) == expected_plain, slug
        else:
            q05 = float(out_uq["q05"])
            q95 = float(out_uq["q95"])
            if out_key.startswith("relative_yield_") and abs(q05) <= 2.0:
                q05 *= 100.0
                q95 *= 100.0
            tol = max(1e-3, 1e-5 * max(abs(q05), abs(q95)))
            assert q05 - tol <= sci_c <= q95 + tol, slug

        if slug == "rain-pattern":
            assert (
                "One soaking keeps the sand wet enough for the microbes for several days, "
                "while small daily amounts mostly evaporate the same day."
            ) in txts["answer_plain"]
            assert "pushes water deeper into the sand" not in txts["answer_plain"]
        if slug == "cloud-seeding":
            assert (
                "When extra rain does fall, the sand stays wet for longer and the microbes "
                "stay active for longer."
            ) in txts["answer_plain"]
            assert "pushes moisture deeper into the sand" not in txts["answer_plain"]
