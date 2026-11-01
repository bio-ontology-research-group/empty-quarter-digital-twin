"""Build per-question plain/science JSON files (`reviews/v08/questions/<slug>.json` and `index.json`)."""
from __future__ import annotations

import argparse
import functools
import hashlib
import json
from pathlib import Path
import sys
from typing import Any

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
if str(ROOT / "src") not in sys.path:
    sys.path.insert(0, str(ROOT / "src"))

from dataclasses import asdict, is_dataclass
from eq_biodt.v03 import interventions, ontology, process  # noqa: E402
from eq_biodt.v03.composition import constants_sha256, load_constants  # noqa: E402
from eq_biodt.v03.external_models import MODELS  # noqa: E402
from eq_biodt.v03.questions import QUESTIONS, LevelSpec, Question, evaluate  # noqa: E402
from eq_biodt.v03.uq_inputs import (  # noqa: E402
    load_parameter_ranges,
    parameter_ranges_sha256,
)
from reproducibility.v08.texts import (  # noqa: E402
    HEADLINE_SPEC,
    MODEL_META,
    OUTCOME_META,
    PLAIN_BASELINE_LABELS,
    PLAIN_TITLES,
    plain_level_label,
    render_question_texts,
)

DEFAULT_UQ_DIR = ROOT / "reviews" / "v08" / "uq"
DEFAULT_OUT_DIR = ROOT / "reviews" / "v08" / "questions"

PAPER_SECTIONS: dict[str, str] = {
    "rain-pattern": "§3.2 Paired hydrological regimes & §3.6 Joint parameter uncertainty",
    "drip-or-flood": "§3.5 External model composition (Experiment 2) & §3.6",
    "cloud-seeding": "§3.5 External model composition (Experiment 5) & §3.6",
    "salty-water": "§3.2 Paired regimes, §3.5 (Experiment 1: Brackish irrigation) & §3.6",
    "alfalfa-bacterium": "§3.5 External model composition (Experiment 10: SA187 inoculation) & §3.6",
    "hotter": "§3.2 Paired regimes, §3.5 (Experiment 3: Climate warming) & §3.6",
    "mulch-shade": "§3.5 External model composition (Experiment 4: Mulch and shade) & §3.6",
    "solar-park": "§3.5 External model composition (Experiment 8: Agrivoltaics) & §3.6",
    "compaction": "§3.2 Paired regimes, §3.5 (Experiment 6: Off-road compaction) & §3.6",
    "biochar": "§3.5 External model composition (Experiment 7: Biochar amendment) & §3.6",
    "missing-co-genes": "§3.1 Ontological capacity gates & §3.2 Closed-gate communities",
    "adding-microbes": "§3.5 External model composition (Experiment 11: Constitutive inoculation) & §3.6",
    "air-hydrogen": "§3.5 External model composition (Experiment 9: Atmospheric H₂) & §3.6",
    "carbon-fixation": "§3.5 External model composition (Experiment 12: Knallgas C-fixation bound) & §3.6",
}


def _extract_branch_trajectory(sim: dict[str, Any] | None) -> dict[str, list[float]]:
    if not sim or "trajectory" not in sim or not sim["trajectory"]:
        zeros = [0.0] * 14
        return {
            "W": list(zeros),
            "N": list(zeros),
            "psi": list(zeros),
            "A": list(zeros),
            "J_h2": list(zeros),
            "J_co": list(zeros),
        }
    traj = sim["trajectory"]
    w_list: list[float] = []
    n_list: list[float] = []
    psi_list: list[float] = []
    a_list: list[float] = []
    jh2_list: list[float] = []
    jco_list: list[float] = []
    for day_entry in traj[:14]:
        w_list.append(float(day_entry.get("water_mm") or 0.0))
        n_list.append(float(day_entry.get("salt_mol_m2") or 0.0))
        gas = day_entry.get("gas") or {}
        psi_val = gas.get("total_water_potential_mpa")
        if psi_val is None:
            psi_val = gas.get("matric_potential_mpa")
        psi_list.append(float(psi_val if psi_val is not None else -85.0))
        a_list.append(float(gas.get("activation") or 0.0))
        h2_s = gas.get("h2_mol_m2_s")
        co_s = gas.get("co_mol_m2_s")
        jh2_list.append(float(h2_s or 0.0) * 86400.0 * 1e6)
        jco_list.append(float(co_s or 0.0) * 86400.0 * 1e6)
    return {
        "W": w_list,
        "N": n_list,
        "psi": psi_list,
        "A": a_list,
        "J_h2": jh2_list,
        "J_co": jco_list,
    }


def _extract_scene(central_res: dict[str, Any], intervened_traj: dict[str, list[float]]) -> dict[str, Any]:
    tw = central_res.get("__twin__")
    sim = tw["result"]["intervened"] if (tw and tw.get("result")) else None
    params = sim.get("parameters", {}) if sim else {}
    forcing = sim.get("forcing", []) if sim else []
    summary = sim.get("summary", {}) if sim else {}

    porosity = float(central_res.get("porosity") or params.get("porosity") or 0.4)
    depth_mm = float(params.get("depth_m", 0.1)) * 1000.0
    if summary.get("water_input_mm") is not None:
        water_mm = float(summary["water_input_mm"])
    elif forcing:
        water_mm = float(sum(float(f.get("water_mm", 0.0)) for f in forcing))
    else:
        water_mm = 12.0

    w_arr = intervened_traj["W"]
    mean_w = sum(w_arr) / len(w_arr) if w_arr else 0.0
    theta = max(0.0, min(porosity, mean_w / depth_mm if depth_mm > 0 else 0.0))
    salt_mol_l = float(central_res.get("mean_nacl_mol_l") or 0.0)
    temp_c = (
        float(sum(float(f.get("temp_c", 25.0)) for f in forcing) / len(forcing))
        if forcing
        else 25.0
    )
    act_by_day = [float(a) for a in intervened_traj["A"]]
    active_days = int(
        central_res.get("active_days")
        if central_res.get("active_days") is not None
        else sum(1 for a in act_by_day if a > 0.01)
    )
    return {
        "water_mm": water_mm,
        "theta": theta,
        "salt_mol_l": salt_mol_l,
        "temp_c": temp_c,
        "porosity": porosity,
        "active_days": active_days,
        "activity_by_day": act_by_day,
    }


def _build_chain_and_methods(q: Question) -> tuple[list[dict[str, str]], dict[str, Any]]:
    upstream_ids: list[str] = []
    downstream_ids: list[str] = []
    for s in q.series:
        for lvl in s.levels:
            for m_id in lvl.upstream:
                if m_id not in upstream_ids:
                    upstream_ids.append(m_id)
            for m_id in lvl.downstream:
                if m_id not in downstream_ids:
                    downstream_ids.append(m_id)

    ordered: list[tuple[str, str]] = (
        [(m, "upstream") for m in upstream_ids]
        + [("eq_biodt_v03", "twin")]
        + [(m, "downstream") for m in downstream_ids]
    )

    chain: list[dict[str, str]] = []
    models_out: list[dict[str, Any]] = []
    const_reg = load_constants()

    for m_id, kind in ordered:
        meta = MODEL_META[m_id]
        chain.append(
            {
                "id": m_id,
                "name_plain": meta["name_plain"],
                "name_science": meta["name_science"],
                "kind": kind,
            }
        )
        if m_id == "eq_biodt_v03":
            models_out.append(
                {
                    "id": m_id,
                    "name": meta["name_science"],
                    "name_plain": meta["name_plain"],
                    "kind": "twin",
                    "equation_tex": meta["equation_tex"],
                    "citation": "EQ-BioDT v0.3 process-ontological soil model (Campbell 1974; Moldrup 2000; van 't Hoff)",
                    "validity": {
                        "osmotic_domain_mol_l": [0.0, 0.1],
                        "porosity": [0.05, 0.90],
                        "water_mm": [0.0, 200.0],
                    },
                    "constants": [],
                }
            )
        else:
            ext = MODELS[m_id]
            c_list = []
            for cid in ext.constants:
                if cid in const_reg:
                    c_entry = dict(const_reg[cid])
                    c_entry["id"] = cid
                    c_list.append(c_entry)
            models_out.append(
                {
                    "id": m_id,
                    "name": ext.name,
                    "name_plain": meta["name_plain"],
                    "kind": ext.kind,
                    "equation_tex": meta["equation_tex"],
                    "validity": {k: list(v) for k, v in ext.validity.items()},
                    "constants": c_list,
                }
            )

    gates = [
        {
            "id": "Group1hNiFeHydrogenase",
            "pathway": "Atmospheric H₂ oxidation",
            "rule": "Closed when sample has 0 MAGs entailed as carrying Group 1h [NiFe]-hydrogenase; fixes J_H2 = 0 exactly.",
        },
        {
            "id": "FormICODehydrogenaseSubunit",
            "pathway": "Atmospheric CO oxidation",
            "rule": "Closed when sample has 0 MAGs entailed as carrying Form I CO dehydrogenase; fixes J_CO = 0 exactly.",
        },
        {
            "id": "CBBAutotrophicCarbonFixationPathway",
            "pathway": "Knallgas autotrophic CO₂ fixation",
            "rule": "Closed when sample has 0 MAGs carrying both Group 1h [NiFe]-hydrogenase and CBB RuBisCO; fixes C_fixed = 0 exactly.",
        },
    ]

    methods = {
        "models": models_out,
        "data": [
            "data/v03/entailed_memberships.json",
            "data/v03/mag_evidence.tsv",
            "data/external/composition_constants.json",
            "data/external/twin_parameter_ranges.json",
        ],
        "gates": gates,
        "admissibility": {
            "dilute_nacl_limit_mol_l": 0.1,
            "refusal_threshold_fraction": 0.05,
            "scope": "Joint parameter and literature-constant uncertainty only; excludes model-structure error.",
        },
    }
    return chain, methods


def _compute_oat(
    spec: LevelSpec,
    out_key: str,
    raw_inputs: list[dict[str, Any]],
) -> dict[str, dict[str, float | None]]:
    oat: dict[str, dict[str, float | None]] = {}
    base_consts = {i["name"]: float(i["value"]) for i in raw_inputs if i["kind"] == "constant"}
    base_params = {i["name"]: float(i["value"]) for i in raw_inputs if i["kind"] == "parameter"}

    for inp in raw_inputs:
        name = inp["name"]
        kind = inp["kind"]
        vals: list[float | None] = []
        for bound in (float(inp["low"]), float(inp["high"])):
            c_override = dict(base_consts)
            p_override = dict(base_params)
            if kind == "constant":
                c_override[name] = bound
            else:
                p_override[name] = bound
            res = evaluate(spec, constants=c_override, params=p_override, keep_twin_result=False)
            v = res.get(out_key) if res.get("__refused__") is None else None
            vals.append(float(v) if isinstance(v, (int, float)) else None)
        lo_v, hi_v = vals[0], vals[1]
        span = abs(hi_v - lo_v) if (lo_v is not None and hi_v is not None) else None
        oat[name] = {"low": lo_v, "high": hi_v, "span": span}
    return oat


def _enrich_inputs(raw_inputs: list[dict[str, Any]]) -> list[dict[str, Any]]:
    const_reg = load_constants()
    param_reg = load_parameter_ranges()
    enriched: list[dict[str, Any]] = []
    for inp in raw_inputs:
        name = inp["name"]
        entry = dict(inp)
        src_meta = const_reg.get(name) if inp["kind"] == "constant" else param_reg.get(name)
        if src_meta:
            entry["unit"] = src_meta.get("unit", "")
            entry["meaning"] = src_meta.get("meaning", "")
            entry["source"] = src_meta.get("source", {})
            if "note" in src_meta:
                entry["note"] = src_meta["note"]
        enriched.append(entry)
    return enriched


def build_questions(
    uq_dir: Path | str = DEFAULT_UQ_DIR,
    out_dir: Path | str = DEFAULT_OUT_DIR,
) -> dict[str, Path]:
    """Build all 14 `<out_dir>/<slug>.json` files and `<out_dir>/index.json`."""
    uq_path = Path(uq_dir)
    out_path = Path(out_dir)
    out_path.mkdir(parents=True, exist_ok=True)

    orig_load = ontology.load_classification
    cached_load = functools.lru_cache(maxsize=8)(orig_load)
    ontology.load_classification = cached_load
    interventions.load_classification = cached_load

    orig_sim = process.simulate
    sim_cache: dict[Any, Any] = {}

    def _cached_sim(initial: Any, forcing: list[dict[str, Any]], p: Any, state_overrides: Any = None) -> Any:
        if state_overrides:
            return orig_sim(initial, forcing, p, state_overrides)
        i_k = tuple(sorted((asdict(initial) if is_dataclass(initial) else dict(initial)).items()))
        f_k = tuple(tuple(sorted((asdict(f) if is_dataclass(f) else dict(f)).items())) for f in forcing)
        p_k = tuple(sorted((asdict(p) if is_dataclass(p) else dict(p)).items()))
        key = (i_k, f_k, p_k)
        if key not in sim_cache:
            sim_cache[key] = orig_sim(initial, forcing, p, None)
        return sim_cache[key]

    process.simulate = _cached_sim
    interventions.simulate = _cached_sim

    written: dict[str, Path] = {}
    index_entries: list[dict[str, Any]] = []

    try:
        for slug, q in QUESTIONS.items():
            uq_file = uq_path / f"{slug}.json"
            if not uq_file.exists():
                raise FileNotFoundError(f"Missing UQ file for {slug}: {uq_file}")
            uq_raw_bytes = uq_file.read_bytes()
            uq_sha = hashlib.sha256(uq_raw_bytes).hexdigest()
            uq_doc = json.loads(uq_raw_bytes.decode("utf-8"))

            uq_series_map = {s["id"]: s for s in uq_doc["series"]}
            series_out: list[dict[str, Any]] = []
            max_refusal_fraction = 0.0

            for s in q.series:
                uq_s = uq_series_map[s.id]
                uq_levels_map = {int(lvl["level"]): lvl for lvl in uq_s["levels"]}
                levels_out: list[dict[str, Any]] = []

                for l_idx, spec in enumerate(s.levels):
                    uq_lvl = uq_levels_map[l_idx]
                    central_res = evaluate(spec, keep_twin_result=True)
                    tw = central_res.get("__twin__")
                    base_sim = tw["result"]["baseline"] if (tw and tw.get("result")) else None
                    intv_sim = tw["result"]["intervened"] if (tw and tw.get("result")) else None

                    base_traj = _extract_branch_trajectory(base_sim)
                    intv_traj = _extract_branch_trajectory(intv_sim)
                    scene = _extract_scene(central_res, intv_traj)

                    central_vals = uq_lvl.get("central", {})
                    central_deduced = set(uq_lvl.get("central_deduced") or [])
                    central_refused = uq_lvl.get("central_refused")

                    outcomes_out: dict[str, Any] = {}
                    out_keys = [
                        k for k in q.headline if k in OUTCOME_META and k in uq_lvl["outputs"]
                    ]
                    for out_key in out_keys:
                        uq_o = uq_lvl["outputs"][out_key]
                        meta = OUTCOME_META[out_key]
                        c_val = central_vals.get(out_key, central_res.get(out_key))
                        if c_val is not None:
                            c_val = float(c_val)
                        ref_f = float(uq_o.get("refusal_fraction", 0.0))
                        if out_key in q.headline:
                            max_refusal_fraction = max(max_refusal_fraction, ref_f)
                        prov = (
                            "deduced"
                            if (out_key in central_deduced or uq_o.get("mode") == "deduced")
                            else meta["provenance"]
                        )
                        outcomes_out[out_key] = {
                            "value": c_val,
                            "median": uq_o.get("median"),
                            "q05": uq_o.get("q05"),
                            "q25": uq_o.get("q25"),
                            "q75": uq_o.get("q75"),
                            "q95": uq_o.get("q95"),
                            "unit": meta["unit"],
                            "unit_plain": meta["unit_plain"],
                            "label_plain": meta["label_plain"],
                            "label_science": meta["label_science"],
                            "provenance": prov,
                            "mode": uq_o.get("mode", "variance"),
                            "refusal_fraction": ref_f,
                            "refusal_reason": central_refused,
                        }

                    lbl_plain = plain_level_label(slug, s.id, spec.label, spec.x)
                    levels_out.append(
                        {
                            "label": lbl_plain,
                            "label_plain": lbl_plain,
                            "label_science": spec.label,
                            "x": spec.x,
                            "sample": spec.sample,
                            "scene": scene,
                            "outcomes": outcomes_out,
                            "trajectories": {
                                "baseline": base_traj,
                                "intervened": intv_traj,
                            },
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

            h_series_id, h_level_idx, h_out_key = HEADLINE_SPEC[slug]
            h_series_spec = next(sr for sr in q.series if sr.id == h_series_id)
            h_level_spec = h_series_spec.levels[h_level_idx]
            h_uq_series = uq_series_map[h_series_id]
            h_uq_lvl = next(l for l in h_uq_series["levels"] if int(l["level"]) == h_level_idx)
            h_uq_out = h_uq_lvl["outputs"][h_out_key]
            h_central_val = float(h_uq_lvl["central"][h_out_key])

            texts = render_question_texts(
                slug=slug,
                central_val=h_central_val,
                outcome_uq=h_uq_out,
                max_refusal_fraction=max_refusal_fraction,
            )

            chain, methods = _build_chain_and_methods(q)
            enriched_inputs = _enrich_inputs(h_uq_lvl["inputs"])
            oat = _compute_oat(h_level_spec, h_out_key, h_uq_lvl["inputs"])

            uncertainty = {
                "headline_output": h_out_key,
                "headline_series": h_series_id,
                "headline_level": h_level_idx,
                "headline_label": h_level_spec.label,
                "n": int(h_uq_lvl["n"]),
                "d": int(h_uq_lvl["d"]),
                "median": h_uq_out.get("median"),
                "q05": h_uq_out.get("q05"),
                "q25": h_uq_out.get("q25"),
                "q75": h_uq_out.get("q75"),
                "q95": h_uq_out.get("q95"),
                "mode": h_uq_out.get("mode", "variance"),
                "refusal_fraction": float(h_uq_out.get("refusal_fraction", 0.0)),
                "S1": h_uq_out.get("S1", {}),
                "ST": h_uq_out.get("ST", {}),
                "group_S1": h_uq_out.get("group_S1", {}),
                "group_ST": h_uq_out.get("group_ST", {}),
                "inputs": enriched_inputs,
                "convergence": h_uq_out.get("convergence", {}),
                "oat": oat,
            }

            reproduce_block = {
                "files": [
                    "src/eq_biodt/v03/questions.py",
                    "src/eq_biodt/v03/uq.py",
                    "reproducibility/v08/questions.py",
                    f"reviews/v08/uq/{slug}.json",
                ],
                "commands": [
                    "uv run python reproducibility/v08/reproduce.py",
                    "uv run python reproducibility/v08/questions.py --uq-dir reviews/v08/uq --out reviews/v08/questions",
                ],
                "digests": {
                    "constants_sha256": constants_sha256(),
                    "parameter_ranges_sha256": parameter_ranges_sha256(),
                    "uq_sha256": uq_sha,
                },
                "paper_section": PAPER_SECTIONS[slug],
            }

            title_plain = PLAIN_TITLES[slug]
            baseline_plain = PLAIN_BASELINE_LABELS[slug]

            q_doc = {
                "version": "question-1",
                "twin_version": "0.8.0",
                "slug": slug,
                "theme": q.theme,
                "title_plain": title_plain,
                "title_science": q.title_science,
                "level_kind": q.level_kind,
                "x_label_plain": q.x_label_plain,
                "x_label_science": q.x_label_science,
                "baseline_label": baseline_plain,
                "headline": list(q.headline),
                "headline_output": h_out_key,
                "headline_level": {
                    "series": h_series_id,
                    "level": h_level_idx,
                    "label": h_level_spec.label,
                    "x": h_level_spec.x,
                },
                "headline_value": h_central_val,
                "headline_plain": texts["headline_plain"],
                "headline_science": texts["headline_science"],
                "answer_plain": texts["answer_plain"],
                "answer_science": texts["answer_science"],
                "series": series_out,
                "chain": chain,
                "limits_plain": texts["limits_plain"],
                "limits_science": texts["limits_science"],
                "methods": methods,
                "uncertainty": uncertainty,
                "reproduce": reproduce_block,
            }

            q_out_file = out_path / f"{slug}.json"
            q_out_file.write_text(
                json.dumps(q_doc, indent=2, sort_keys=True, allow_nan=False) + "\n",
                encoding="utf-8",
            )
            written[slug] = q_out_file

            index_entries.append(
                {
                    "slug": slug,
                    "theme": q.theme,
                    "title_plain": title_plain,
                    "title_science": q.title_science,
                    "headline_plain": texts["headline_plain"],
                    "headline_science": texts["headline_science"],
                    "answer_plain": texts["answer_plain"],
                    "answer_science": texts["answer_science"],
                    "thumbnail": None,
                }
            )

        index_doc = {
            "version": "question-1",
            "twin_version": "0.8.0",
            "questions": index_entries,
        }
        index_file = out_path / "index.json"
        index_file.write_text(
            json.dumps(index_doc, indent=2, sort_keys=True, allow_nan=False) + "\n",
            encoding="utf-8",
        )
        written["index"] = index_file
    finally:
        ontology.load_classification = orig_load
        interventions.load_classification = orig_load
        process.simulate = orig_sim
        interventions.simulate = orig_sim

    return written


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--uq-dir",
        type=Path,
        default=DEFAULT_UQ_DIR,
        help="Directory containing gathered UQ JSON files (<slug>.json)",
    )
    parser.add_argument(
        "--out",
        type=Path,
        default=DEFAULT_OUT_DIR,
        help="Output directory for question JSON files and index.json",
    )
    args = parser.parse_args(argv)
    build_questions(uq_dir=args.uq_dir, out_dir=args.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
