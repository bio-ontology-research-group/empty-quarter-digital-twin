from __future__ import annotations

import functools
import json
from pathlib import Path
import re
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
from reproducibility.v08.questions import build_questions
from reproducibility.v08.texts import (
    HEADLINE_SPEC,
    JARGON,
    plain_number,
    round_plain,
    science_central,
)

COMMITTED_QUESTIONS_DIR = ROOT / "reviews" / "v08" / "questions"

_BuiltDirCache: Path | None = None


def _build_fixture_questions(tmp_root: Path) -> Path:
    uq_dir = tmp_root / "uq"
    out_dir = tmp_root / "questions"
    uq_dir.mkdir(parents=True, exist_ok=True)
    out_dir.mkdir(parents=True, exist_ok=True)

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
                for out in lvl["outcomes"].values():
                    plain_parts.append(out["unit_plain"])
        for step in q["chain"]:
            plain_parts.append(step["name_plain"])
        txt = " ".join(plain_parts)
        assert not JARGON.search(txt), f"{q['slug']}: matched {JARGON.search(txt)}"
        assert "the models themselves could also be wrong" in q["limits_plain"].lower()
        assert (
            "the range shows how much the answer moves when uncertain numbers in the models are varied; the models themselves could also be wrong"
            in q["limits_plain"].lower()
        )
        u = q["uncertainty"]
        if u.get("q05") is not None and u.get("q95") is not None and u["q05"] < 0 < u["q95"]:
            assert "may go up or down" in q["answer_plain"].lower(), q["slug"]
        max_refusal = max(
            out["refusal_fraction"]
            for s in q["series"]
            for lvl in s["levels"]
            for out in lvl["outcomes"].values()
        )
        if max_refusal >= 0.05:
            assert "the model cannot tell for such" in q["answer_plain"].lower(), q["slug"]
    assert "at most" in load("carbon-fixation")["answer_plain"].lower()


def test_plain_matches_science():
    for q in load_all():
        assert plain_number(q) == round_plain(science_central(q)), q["slug"]
        assert round_plain(science_central(q)) in q["answer_plain"], q["slug"]


def test_refused_level_rendering_fields():
    lvl = top_level(load("salty-water"), series="daily-irrigation")
    assert lvl["outcomes"]["h2_uptake_contrast_umol_m2"]["refusal_fraction"] > 0
    assert "the model cannot tell for such salty water" in load("salty-water")["answer_plain"].lower()


@pytest.mark.skipif(
    not (COMMITTED_QUESTIONS_DIR / "index.json").exists(),
    reason="reviews/v08/questions/ not yet generated from cluster UQ",
)
def test_committed_questions_on_disk():
    assert len(load_index(COMMITTED_QUESTIONS_DIR)) == 14
    assert all(validate(q) for q in load_all(COMMITTED_QUESTIONS_DIR))
