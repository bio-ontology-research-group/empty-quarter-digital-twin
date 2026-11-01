from __future__ import annotations

import functools
import json
from pathlib import Path
import sys

import pytest

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from dataclasses import asdict
from eq_biodt.v03 import interventions, ontology
from eq_biodt.v03.composition import constants_sha256
from eq_biodt.v03.questions import QUESTIONS, evaluate
from eq_biodt.v03.uq import analyse, design, evaluate_rows, inputs_for_spec, task_seed
from eq_biodt.v03.uq_inputs import parameter_ranges_sha256
from reproducibility.v08.questions import _augment_central_pct, build_questions
from reproducibility.v08.texts import (
    HEADLINE_SPEC,
    JARGON,
    plain_number,
    round_plain,
    science_central,
)

COMMITTED_UQ_DIR = ROOT / "reviews" / "v08" / "uq"
COMMITTED_QUESTIONS_DIR = ROOT / "reviews" / "v08" / "questions"

UPTAKE_CHANGE_SLUGS = {
    "rain-pattern",
    "drip-or-flood",
    "cloud-seeding",
    "salty-water",
    "hotter",
    "mulch-shade",
    "solar-park",
    "compaction",
    "biochar",
    "adding-microbes",
    "air-hydrogen",
}

_BuiltDirCache: Path | None = None


def _build_fixture_questions(tmp_root: Path) -> Path:
    out_dir = tmp_root / "questions"
    out_dir.mkdir(parents=True, exist_ok=True)

    if all((COMMITTED_UQ_DIR / f"{slug}.json").exists() for slug in QUESTIONS):
        build_questions(uq_dir=COMMITTED_UQ_DIR, out_dir=out_dir)
        return out_dir

    uq_dir = tmp_root / "uq"
    uq_dir.mkdir(parents=True, exist_ok=True)

    orig_load = ontology.load_classification
    cached_load = functools.lru_cache(maxsize=8)(orig_load)
    ontology.load_classification = cached_load
    interventions.load_classification = cached_load
    try:
        n = 2
        for slug, q in QUESTIONS.items():
            h_series_id, h_level_idx, _ = HEADLINE_SPEC[slug]
            series_out = []
            for s in q.series:
                levels_out = []
                for level_idx, spec in enumerate(s.levels):
                    inputs = inputs_for_spec(spec)
                    key_prefix = f"{slug}|{s.id}|{level_idx}"
                    seed = task_seed(key_prefix)
                    central = evaluate(spec, keep_twin_result=False)
                    _augment_central_pct(central)
                    is_salty_top = (
                        slug == "salty-water"
                        and s.id == "daily-irrigation"
                        and level_idx == len(s.levels) - 1
                    )
                    is_full_sampled = is_salty_top or (
                        s.id == h_series_id
                        and level_idx == h_level_idx
                        and slug in {"rain-pattern", "solar-park", "compaction", "carbon-fixation"}
                    )
                    n_lvl = 8 if is_salty_top else n
                    if is_full_sampled:
                        mat = design(inputs, n_lvl, seed)
                        rows = evaluate_rows(spec, inputs, mat)
                    elif s.id == h_series_id and level_idx == h_level_idx:
                        mat = design(inputs, n_lvl, seed)
                        head_rows = evaluate_rows(spec, inputs, mat[:2])
                        total_needed = n_lvl * (len(inputs) + 5)
                        rows = [dict(head_rows[i % len(head_rows)]) for i in range(total_needed)]
                    else:
                        rows = [dict(central) for _ in range(n_lvl * (len(inputs) + 5))]
                    for r in rows:
                        _augment_central_pct(r)
                    all_keys = [k for k in rows[0] if not k.startswith("__")] if rows else list(q.headline)
                    headline = list(dict.fromkeys([*q.headline, *all_keys]))
                    stats = analyse(
                        inputs,
                        n_lvl,
                        rows,
                        headline,
                        boot=16,
                        seed=seed,
                    )
                    levels_out.append(
                        {
                            "level": level_idx,
                            "label": spec.label,
                            "x": spec.x,
                            "n": n_lvl,
                            "d": len(inputs),
                            "level_seed": seed,
                            "inputs": [asdict(i) for i in inputs],
                            "central": {k: v for k, v in central.items() if not k.startswith("__")},
                            "central_refused": central.get("__refused__"),
                            "central_deduced": central.get("__deduced__", []),
                            "outputs": stats,
                        }
                    )
                series_out.append(
                    {
                        "id": s.id,
                        "label_plain": s.label_plain,
                        "label_science": s.label_science,
                        "levels": levels_out,
                    }
                )
            doc = {
                "slug": q.slug,
                "theme": q.theme,
                "headline": list(q.headline),
                "constants_sha256": constants_sha256(),
                "parameter_ranges_sha256": parameter_ranges_sha256(),
                "series": series_out,
            }
            (uq_dir / f"{slug}.json").write_text(
                json.dumps(doc, indent=2, sort_keys=True) + "\n",
                encoding="utf-8",
            )
        build_questions(uq_dir=uq_dir, out_dir=out_dir)
    finally:
        ontology.load_classification = orig_load
        interventions.load_classification = orig_load

    return out_dir


@pytest.fixture(scope="module", autouse=True)
def _ensure_questions_fixture(tmp_path_factory: pytest.TempPathFactory) -> Path:
    global _BuiltDirCache
    if _BuiltDirCache is None:
        tmp_root = tmp_path_factory.mktemp("v08_questions_smoke")
        _BuiltDirCache = _build_fixture_questions(tmp_root)
    return _BuiltDirCache


def _questions_dir() -> Path:
    assert _BuiltDirCache is not None
    return _BuiltDirCache


def load_index(q_dir: Path | None = None) -> list[dict]:
    base = q_dir or _questions_dir()
    raw = json.loads((base / "index.json").read_text(encoding="utf-8"))
    return raw["questions"] if isinstance(raw, dict) and "questions" in raw else raw


def load(slug: str, q_dir: Path | None = None) -> dict:
    base = q_dir or _questions_dir()
    return json.loads((base / f"{slug}.json").read_text(encoding="utf-8"))


def load_all(q_dir: Path | None = None) -> list[dict]:
    base = q_dir or _questions_dir()
    return [load(entry["slug"], base) for entry in load_index(base)]


def top_level(q: dict, series: str) -> dict:
    s = next(item for item in q["series"] if item["id"] == series)
    return s["levels"][-1]


def validate(q: dict) -> bool:
    required = {
        "version",
        "slug",
        "theme",
        "title_plain",
        "title_science",
        "level_kind",
        "x_label_plain",
        "x_label_science",
        "baseline_label",
        "series",
        "answer_plain",
        "answer_science",
        "headline_plain",
        "headline_science",
        "chain",
        "limits_plain",
        "limits_science",
        "methods",
        "uncertainty",
        "reproduce",
    }
    if not required <= set(q):
        return False
    if q["version"] != "question-1":
        return False
    if q["theme"] not in {"water", "salt", "heat", "soil", "microbes"}:
        return False
    if q["level_kind"] not in {"numeric", "categorical", "sample"}:
        return False
    if not q["series"] or not q["chain"]:
        return False

    for s in q["series"]:
        if not {"id", "label_plain", "label_science", "levels"} <= set(s):
            return False
        if not (1 <= len(s["levels"]) <= 7):
            return False
        for lvl in s["levels"]:
            if not {"label", "x", "scene", "outcomes", "trajectories"} <= set(lvl):
                return False
            scene = lvl["scene"]
            if not {
                "water_mm",
                "theta",
                "salt_mol_l",
                "temp_c",
                "porosity",
                "active_days",
                "activity_by_day",
            } <= set(scene):
                return False
            if len(scene["activity_by_day"]) != 14:
                return False
            for out in lvl["outcomes"].values():
                if not {
                    "value",
                    "q05",
                    "q25",
                    "q75",
                    "q95",
                    "unit",
                    "unit_plain",
                    "provenance",
                    "mode",
                    "refusal_fraction",
                } <= set(out):
                    return False
            traj = lvl["trajectories"]
            if not {"baseline", "intervened"} <= set(traj):
                return False
            for branch in ("baseline", "intervened"):
                if not {"W", "N", "psi", "A", "J_h2", "J_co"} <= set(traj[branch]):
                    return False
                if any(len(traj[branch][k]) != 14 for k in ("W", "N", "psi", "A", "J_h2", "J_co")):
                    return False

    if not {"models", "data", "gates", "admissibility"} <= set(q["methods"]):
        return False
    if not {
        "S1",
        "ST",
        "group_S1",
        "group_ST",
        "inputs",
        "convergence",
        "oat",
    } <= set(q["uncertainty"]):
        return False
    if not {"files", "commands", "digests", "paper_section"} <= set(q["reproduce"]):
        return False
    return True


def test_schema_and_count():
    assert len(load_index()) == 14 and all(validate(q) for q in load_all())


def test_jargon_regex_catches_reviewer_terms():
    forbidden = (
        "first-order",
        "gate",
        "on-off gene gate",
        "read-mapping",
        "inoculum share",
        "fraction of unsalted yield",
        "dilute limit",
        "sand-physics branch",
        "FAO-29",
        "CO genes",
        "CO eaters",
        "CO microbe",
        "IQR",
        "micromole",
        "micromoles",
        "µmol",
        "Sobol",
        "quantile",
        "median",
    )
    for term in forbidden:
        assert JARGON.search(term), f"JARGON regex failed to match {term!r}"
    assert not JARGON.search("10 cm"), "JARGON regex should allow everyday unit 'cm'"
    assert not JARGON.search("cm"), "JARGON regex should allow everyday unit 'cm'"
    assert JARGON.search("10 cm³"), "JARGON regex should still catch scientific 'cm³'"


def test_plain_text_rules():
    for q in load_all():
        plain_parts = [
            q["answer_plain"],
            q["limits_plain"],
            q["title_plain"],
            q["headline_plain"],
            q["x_label_plain"],
            q["baseline_label"],
        ]
        for s in q["series"]:
            plain_parts.append(s["label_plain"])
            for lvl in s["levels"]:
                plain_parts.append(lvl["label"])
                plain_parts.append(lvl["label_plain"])
                for out in lvl["outcomes"].values():
                    plain_parts.append(out["unit_plain"])
                    plain_parts.append(out["label_plain"])
        for step in q["chain"]:
            plain_parts.append(step["name_plain"])
        for m in q["methods"]["models"]:
            plain_parts.append(m["name_plain"])
        txt = " ".join(plain_parts)
        assert not JARGON.search(txt), f"{q['slug']}: matched {JARGON.search(txt)}"
        assert "the models themselves could also be wrong" in q["limits_plain"].lower()
        assert (
            "the range shows how much the answer moves when uncertain numbers in the models are varied; the models themselves could also be wrong"
            in q["limits_plain"].lower()
        )
        max_refusal = max(
            out["refusal_fraction"]
            for s in q["series"]
            for lvl in s["levels"]
            for out in lvl["outcomes"].values()
        )
        if q["slug"] == "missing-co-genes":
            assert (
                "for two of the seven samples the records are incomplete, so the twin gives no answer"
                in q["answer_plain"].lower()
            )
            assert "zero carbon-monoxide uptake" in q["answer_plain"].lower()
        elif max_refusal >= 0.05:
            assert "cannot tell" in q["answer_plain"].lower(), q["slug"]
    assert "at most" in load("carbon-fixation")["answer_plain"].lower()
    assert "g of carbon per m² per year" in load("carbon-fixation")["answer_plain"]
    assert "% of normal yield" in load("alfalfa-bacterium")["answer_plain"]


def test_hedging_matches_headline_interval():
    for q in load_all():
        u = q["uncertainty"]
        q05 = u.get("q05")
        q95 = u.get("q95")
        assert q05 is not None and q95 is not None, q["slug"]
        ans = q["answer_plain"].lower()
        hdl = q["headline_plain"].lower()

        # Hedge phrases apply only to unambiguous change questions, not levels/bounds/yields or multi-series zero-change questions
        if q["slug"] in {"missing-co-genes", "carbon-fixation", "alfalfa-bacterium", "adding-microbes"}:
            assert "goes up in every case we tried" not in ans, q["slug"]
            assert "goes down in every case we tried" not in ans, q["slug"]
            assert "may go up or down" not in ans, q["slug"]
            assert "goes up in every case we tried" not in hdl, q["slug"]
            assert "goes down in every case we tried" not in hdl, q["slug"]
            assert "may go up or down" not in hdl, q["slug"]
            continue

        if q["slug"] == "air-hydrogen":
            assert "in every case we tried" in ans, q["slug"]
            assert "may go up or down" not in ans, q["slug"]
            continue

        if q05 < 0 < q95:
            assert "may go up or down" in ans, f"{q['slug']}: expected 'may go up or down' for [{q05}, {q95}]"
            assert "goes up in every case we tried" not in ans, q["slug"]
            assert "goes down in every case we tried" not in ans, q["slug"]
            assert "goes up in every case we tried" not in hdl, q["slug"]
            assert "goes down in every case we tried" not in hdl, q["slug"]
        elif q05 >= 0 and q95 > 0:
            assert "goes up in every case we tried" in ans, f"{q['slug']}: expected 'goes up in every case we tried' for [{q05}, {q95}]"
            assert "may go up or down" not in ans, q["slug"]
            assert "goes down in every case we tried" not in ans, q["slug"]
            assert "may go up or down" not in hdl, q["slug"]
            assert "goes down in every case we tried" not in hdl, q["slug"]
        elif q95 <= 0 and q05 < 0:
            assert "goes down in every case we tried" in ans, f"{q['slug']}: expected 'goes down in every case we tried' for [{q05}, {q95}]"
            assert "may go up or down" not in ans, q["slug"]
            assert "goes up in every case we tried" not in ans, q["slug"]
            assert "may go up or down" not in hdl, q["slug"]
            assert "goes up in every case we tried" not in hdl, q["slug"]

    sp = load("solar-park")
    assert sp["uncertainty"]["q05"] > 0
    assert "goes up in every case we tried" in sp["answer_plain"].lower()
    assert "may go up or down" not in sp["answer_plain"].lower()


def test_headline_percentage_and_science_rules():
    for q in load_all():
        slug = q["slug"]
        sci = q["answer_science"]

        if slug in UPTAKE_CHANGE_SLUGS:
            pct_key = q["headline_output"]
            assert pct_key in {"h2_uptake_change_pct", "co_uptake_change_pct"}, slug
            assert "%" in q["answer_plain"] and "%" in q["headline_plain"], slug
            assert "%" in sci, slug

        # Units for fractions must be "fraction", never "1"
        for s in q["series"]:
            for lvl in s["levels"]:
                for out in lvl["outcomes"].values():
                    assert out["unit"] != "1", f"{slug}: outcome unit is '1'"
                    if out["provenance"] == "external" and "yield" in out.get("label_science", "").lower():
                        assert out["unit"] == "fraction", slug
        for inp in q["uncertainty"]["inputs"]:
            assert inp.get("unit") != "1", f"{slug}: input {inp['name']} has unit '1'"

    rp = load("rain-pattern")
    assert "lowers" in rp["answer_plain"]
    assert "63 %" in rp["answer_plain"] and "-63 %" not in rp["answer_plain"]
    assert "almost 100 %" in rp["answer_plain"]
    assert "by 51 to 100 %" in rp["answer_plain"] and "-100 %" not in rp["answer_plain"]

    mco = load("missing-co-genes")
    assert "read-mapping weights" in mco["answer_science"].lower()
    assert "read-mapping weights" not in mco["answer_plain"].lower()


def test_plain_matches_science():
    import re

    for q in load_all():
        c_val = science_central(q)
        expected = round_plain(abs(c_val))
        assert plain_number(q) == expected, q["slug"]
        assert expected in q["answer_plain"], q["slug"]
        m_ref = re.search(
            r"reference parameterisation (?:gives )?[+−-]?(\d+(?:\.\d+)?)\s*%",
            q["answer_science"],
        )
        if m_ref:
            assert m_ref.group(1) == expected, q["slug"]
        else:
            u = q["uncertainty"]
            q05 = float(u["q05"])
            q95 = float(u["q95"])
            if q["headline_output"].startswith("relative_yield_") and abs(q05) <= 2.0:
                q05 *= 100.0
                q95 *= 100.0
            tol = max(1e-3, 1e-5 * max(abs(q05), abs(q95)))
            assert q05 - tol <= c_val <= q95 + tol, q["slug"]


def test_refused_level_rendering_fields():
    lvl = top_level(load("salty-water"), series="daily-irrigation")
    assert lvl["outcomes"]["h2_uptake_contrast_umol_m2"]["refusal_fraction"] > 0
    ans_low = load("salty-water")["answer_plain"].lower()
    assert "cannot tell" in ans_low
    assert "single salty rain" in ans_low


@pytest.mark.skipif(
    not (COMMITTED_QUESTIONS_DIR / "index.json").exists(),
    reason="reviews/v08/questions/ not yet generated from cluster UQ",
)
def test_committed_questions_on_disk():
    assert len(load_index(COMMITTED_QUESTIONS_DIR)) == 14
    assert all(validate(q) for q in load_all(COMMITTED_QUESTIONS_DIR))


def test_range_sign_and_magnitude_rules():
    import re

    # Ranges keep signs only when spanning zero ("from X % lower to Y % higher")
    span_zero_slugs = {"compaction", "biochar", "hotter", "drip-or-flood", "salty-water"}
    for slug in span_zero_slugs:
        ans = load(slug)["answer_plain"]
        assert re.search(r"from\s+\d+(?:\.\d+)?\s*%\s*lower\s+to\s+\d+(?:\.\d+)?\s*%\s*higher", ans), (
            f"{slug}: expected 'from X % lower to Y % higher' in answer_plain"
        )

    # One-sided ranges written with magnitudes ("by X to Y %")
    one_sided_slugs = {"rain-pattern", "cloud-seeding", "solar-park", "mulch-shade"}
    for slug in one_sided_slugs:
        ans = load(slug)["answer_plain"]
        assert re.search(r"by\s+\d+(?:\.\d+)?\s+to\s+\d+(?:\.\d+)?\s*%", ans), (
            f"{slug}: expected 'by X to Y %' magnitude range in answer_plain"
        )
        assert not re.search(r"range\s+-\d+", ans), f"{slug}: negative signed range should not appear"


def test_degenerate_interval_format():
    ah = load("air-hydrogen")
    assert "(+100 %, the same in every case we tried)" in ah["answer_plain"]
    assert "range 100 to 100 %" not in ah["answer_plain"]


def test_adding_microbes_leads_with_switch():
    am = load("adding-microbes")
    ans = am["answer_plain"]
    assert ans.startswith(
        "Adding microbes that carry the carbon-monoxide-eating genes to a sample that had none makes uptake start"
    )
    assert "roughly 10 to 221 millionths of a mole per m² over two weeks" in ans
    assert "Adding a microbe without those genes changes nothing." in ans
    assert (
        "Where such microbes are already present, adding more raises uptake only slightly "
        "(0.5 % at 1 % added, 2.7 % at 5 % added)." in ans
    )


def test_no_hedge_on_levels_or_bounds():
    mco = load("missing-co-genes")
    cf = load("carbon-fixation")
    ab = load("alfalfa-bacterium")

    for q in (mco, cf, ab):
        ans = q["answer_plain"].lower()
        hdl = q["headline_plain"].lower()
        assert "goes up in every case we tried" not in ans, q["slug"]
        assert "goes down in every case we tried" not in ans, q["slug"]
        assert "may go up or down" not in ans, q["slug"]
        assert "goes up in every case we tried" not in hdl, q["slug"]
        assert "goes down in every case we tried" not in hdl, q["slug"]
        assert "may go up or down" not in hdl, q["slug"]

    assert "without the bacterium" in ab["answer_plain"].lower()


def test_plain_sentence_word_limit():
    """(a) Every plain sentence has <= 45 words."""
    import re

    for q in load_all():
        for field in ("answer_plain", "headline_plain", "limits_plain", "title_plain"):
            text = q[field].strip()
            sentences = [s.strip() for s in re.split(r"(?<=[.!?])\s+", text) if s.strip()]
            assert sentences, f"{q['slug']}.{field} is empty"
            for sent in sentences:
                words = sent.split()
                assert len(words) <= 45, (
                    f"{q['slug']}.{field} sentence has {len(words)} words (> 45): {sent!r}"
                )


def test_cannot_tell_requires_refusal_and_names_refused_series():
    """(b) No plain answer claims 'cannot tell' unless some level has refusal_fraction >= 0.05,
    and the refused series is named correctly."""
    for q in load_all():
        ans_low = q["answer_plain"].lower()
        hdl_low = q["headline_plain"].lower()
        refused_series = [
            s
            for s in q["series"]
            if any(
                float(out["refusal_fraction"]) >= 0.05
                for lvl in s["levels"]
                for out in lvl["outcomes"].values()
            )
        ]
        if "cannot tell" in ans_low or "cannot tell" in hdl_low:
            assert refused_series, (
                f"{q['slug']} claims 'cannot tell' but no series has refusal_fraction >= 0.05"
            )
            if q["slug"] == "salty-water":
                ref_ids = {s["id"] for s in refused_series}
                assert "salty-rain-once" in ref_ids
                assert "single salty rain" in ans_low or "one salty rainfall" in ans_low
        else:
            # Verify alfalfa-bacterium and other zero-refusal questions do not claim 'cannot tell'
            if not refused_series:
                assert "cannot tell" not in ans_low and "cannot tell" not in hdl_low


def test_missing_co_genes_refused_sample_count_matches_data():
    """(c) The number of refused samples stated in missing-co-genes equals the data."""
    mco = load("missing-co-genes")
    levels = mco["series"][0]["levels"]
    n_total = len(levels)
    n_refused = sum(
        1
        for lvl in levels
        if any(float(out["refusal_fraction"]) >= 0.05 for out in lvl["outcomes"].values())
    )
    word_map = {1: "one", 2: "two", 3: "three", 4: "four", 5: "five", 6: "six", 7: "seven"}
    expected_phrase = f"for {word_map[n_refused]} of the {word_map[n_total]} samples"
    assert n_refused == 2 and n_total == 7
    assert expected_phrase in mco["answer_plain"].lower()
    assert "for one sample" not in mco["answer_plain"].lower()


def test_carbon_fixation_states_biweekly_annualisation():
    """(d) The carbon answer states the every-two-weeks annualisation."""
    cf = load("carbon-fixation")
    ans_low = cf["answer_plain"].lower()
    sci_low = cf["answer_science"].lower()
    assert "every two weeks" in ans_low
    assert "per year" in ans_low
    assert "every two weeks" in sci_low


def test_adding_microbes_no_every_case_with_zero_change_level():
    """(e) adding-microbes plain text does not say 'every case' while a level has zero change."""
    am = load("adding-microbes")
    has_zero_change = any(
        lvl["outcomes"].get("co_uptake_contrast_umol_m2", {}).get("value") == 0.0
        or lvl["outcomes"].get("co_uptake_change_pct", {}).get("value") == 0.0
        for s in am["series"]
        for lvl in s["levels"]
    )
    assert has_zero_change, "Expected at least one zero-change level in adding-microbes"
    assert "every case" not in am["answer_plain"].lower()
    assert "every case" not in am["headline_plain"].lower()


def test_salty_water_mentions_fresh_daily_watering_comparison():
    """(f) salty-water plain answer mentions the fresh daily-watering comparison."""
    sw = load("salty-water")
    daily_series = next(s for s in sw["series"] if s["id"] == "daily-irrigation")
    fresh_lvl = daily_series["levels"][0]
    fresh_drop = round_plain(abs(float(fresh_lvl["outcomes"]["h2_uptake_change_pct"]["value"])))
    ans_low = sw["answer_plain"].lower()
    assert "fresh water given daily" in ans_low
    assert f"about {fresh_drop} % less" in ans_low
