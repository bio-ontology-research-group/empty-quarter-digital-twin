#!/usr/bin/env python3
"""Export the static web bundle of the public Twin page (rubalkhali.science/model).

Writes ``<out>/{bundle.json, manifest.json, py/...}``. The browser runs the released model code
(``eq_biodt.v03.process`` and ``eq_biodt.v03.theory``) unchanged under Pyodide; every number shown on
the page is read from the committed ``reviews/`` and ``data/v03`` outputs, never recomputed from
scratch by a different implementation. Nothing here modifies the model, the manuscript or the data.

    python scripts/web/export_web_bundle.py --out ../empty-quarter/website/frontend/public/twin-data/eq_biodt_v06
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src"
DATA = ROOT / "data" / "v03"
V05 = ROOT / "reviews" / "v05"
V06 = ROOT / "reviews" / "v06"
PAPER_TABLES = ROOT / "paper" / "tables"
sys.path.insert(0, str(SRC))

from eq_biodt.v03 import theory  # noqa: E402

DEFAULT_OUT = ROOT.parent / "empty-quarter" / "website" / "frontend" / "public" / "twin-data" / "eq_biodt_v06"
VERBATIM = ("v03/process.py", "v03/theory.py")
SCHEDULES = ("intact", "single_pulse", "four_pulses", "daily_pulses")
MODES = ("gate", "proportional")
CBB_STEPS = {  # eq_biodt.v03.evidence.CBB_STEPS (evidence.py imports pandas, so it is read by value and tested)
    "carboxylation": ("K01601",), "phosphoglycerate_kinase": ("K00927",),
    "gap_dehydrogenase": ("K05298", "K00150", "K00134"), "triose_isomerase": ("K01803",),
    "aldolase": ("K01623", "K01624"), "fructose_bisphosphatase": ("K03841", "K02446", "K11532", "K01086"),
    "transketolase": ("K00615",),
    "sedoheptulose_regeneration": ("K01100", "K11532", "K01086", "K00616", "K13810"),
    "ribose_isomerase": ("K01807", "K01808"), "ribulose_epimerase": ("K01783",),
    "phosphoribulokinase": ("K00855", "K28462"),
}
SBPASE = {"K01100", "K11532", "K01086"}
TRANSALDOLASE = {"K00616", "K13810"}
STEP_LABEL = {
    "carboxylation": "RuBisCO large subunit", "phosphoglycerate_kinase": "Phosphoglycerate kinase",
    "gap_dehydrogenase": "GAPDH", "triose_isomerase": "Triose-phosphate isomerase", "aldolase": "Fructose-bisphosphate aldolase",
    "fructose_bisphosphatase": "Fructose-1,6-bisphosphatase", "transketolase": "Transketolase",
    "sedoheptulose_regeneration": "Sedoheptulose-bisphosphatase or transaldolase route",
    "ribose_isomerase": "Ribose-5-phosphate isomerase", "ribulose_epimerase": "Ribulose-phosphate 3-epimerase",
    "phosphoribulokinase": "Phosphoribulokinase",
}
CLASSES = ("H2CarrierMAG", "COCarrierMAG", "H2AutotrophicCandidateMAG", "CBBCompleteMAG")


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def rnd(x, sig=6):
    """Round to significant digits (JSON-safe; None for NaN)."""
    if x is None or (isinstance(x, float) and x != x):
        return None
    if isinstance(x, (int, bool)):
        return x
    return float(f"{float(x):.{sig}g}")


def git(*args) -> str:
    try:
        return subprocess.run(["git", "-C", str(ROOT), *args], capture_output=True, text=True, check=True).stdout.strip()
    except Exception:
        return ""


def tex_numbers() -> dict:
    out = {}
    for f in sorted(PAPER_TABLES.glob("v0*_numbers.tex")):
        for name, value in re.findall(r"\\newcommand\{\\([A-Za-z]+)\}\{([^}]*)\}", f.read_text()):
            out[name] = value
    return out


def project_version() -> str:
    m = re.search(r'^version\s*=\s*"([^"]+)"', (ROOT / "pyproject.toml").read_text(), re.M)
    return m.group(1) if m else "0.6.0"


# --------------------------------------------------------------------------- readable axioms
KO = re.compile(r"kegg:(K\d{5})")


def readable(axiom: str) -> str:
    """Render an OWL functional-syntax axiom in plain language (lossless for the axiom shapes used here)."""
    a = axiom.replace("rak:", "").replace("kegg:", "")
    a = a.replace("RAK_GUILD_", "guild ")
    m = re.fullmatch(r"SubClassOf\(\?mag ObjectSomeValuesFrom\(hasKO (K\d{5})\)\)", a)
    if m:
        return f"this MAG has KO {m[1]}"
    m = re.fullmatch(r"SubClassOf\(\?mag ObjectSomeValuesFrom\(hasMarkerGuild guild (\w+)\)\)", a)
    if m:
        return f"this MAG has marker guild {m[1]}"
    m = re.fullmatch(r"SubClassOf\((K\d{5}) StepKO_(\w+)\)", a)
    if m:
        return f"{m[1]} is a KO of CBB step {m[2]}"
    m = re.fullmatch(r"EquivalentClasses\(CBBStep_(\w+) ObjectSomeValuesFrom\(hasKO StepKO_\w+\)\)", a)
    if m:
        return f"CBB step {m[1]} = has some KO of that step"
    m = re.fullmatch(r"EquivalentClasses\(H2CarrierMAG ObjectSomeValuesFrom\(hasMarkerGuild guild TOTAL_H2ASE\)\)", a)
    if m:
        return "H2CarrierMAG = has a hydrogenase marker guild (groups 1l, 1h or 2a)"
    m = re.fullmatch(r"EquivalentClasses\(COCarrierMAG ObjectSomeValuesFrom\(hasMarkerGuild guild COXL_FORMI\)\)", a)
    if m:
        return "COCarrierMAG = has the form I CO-dehydrogenase marker guild"
    if a.startswith("EquivalentClasses(CBBCompleteMAG"):
        return "CBBCompleteMAG = has all 11 CBB steps"
    if a.startswith("EquivalentClasses(H2AutotrophicCandidateMAG"):
        return "H2AutotrophicCandidateMAG = H2CarrierMAG and CBBCompleteMAG"
    if a.startswith("DisjointUnion(guild TOTAL_H2ASE"):
        return "hydrogenase guild = exactly one of groups 1l, 1h, 2a"
    return a


# --------------------------------------------------------------------------- sections
def build_samples() -> tuple[list, dict, dict]:
    ev = pd.read_csv(DATA / "sample_evidence.tsv", sep="\t")
    paired = pd.read_csv(V06 / "paired_regimes_by_sample.tsv", sep="\t").set_index("sample")
    # the median that defines gamma is that of q over the 550 mapped samples
    medians = {"H2CarrierMAG": float(paired.q_h2.median()), "COCarrierMAG": float(paired.q_co.median())}
    samples = []
    for r in ev.itertuples(index=False):
        d = {"id": str(r.sample), "site": int(r.site), "compartment": r.compartment, "trip": int(r.trip_num),
             "date": str(r.collection_date), "lat": rnd(r.lat, 8), "lon": rnd(r.lon, 8), "mapped": r.sample in paired.index}
        if d["mapped"]:
            p = paired.loc[r.sample]
            d["q"] = {"H2CarrierMAG": rnd(p.q_h2, 8), "COCarrierMAG": rnd(p.q_co, 8),
                      "H2AutotrophicCandidateMAG": rnd(r.h2_cbb_all_steps_share_mapped, 8),
                      "CBBCompleteMAG": rnd(r.cbb_all_steps_detected_share_mapped, 8)}
            d["gamma"] = {m: {"h2": rnd(p[f"gamma_h2_{m}"]), "co": rnd(p[f"gamma_co_{m}"])} for m in MODES}
            d["uptake_umol_m2"] = {m: {g: [rnd(p[f"{g}_{m}_{s}"] * 1e6) for s in SCHEDULES] for g in ("h2", "co")} for m in MODES}
            d["effect_umol_m2"] = {m: {g: [rnd(p[f"{g}eff_{m}_{s}"] * 1e6) for s in SCHEDULES[1:]] for g in ("h2", "co")} for m in MODES}
        samples.append(d)
    return samples, medians, {"schedules": list(SCHEDULES), "effect_schedules": list(SCHEDULES[1:]), "unit": "umol m-2 over 14 days"}


def build_sites(samples: list) -> dict:
    sites = pd.read_csv(V05 / "primary_sites.tsv", sep="\t")
    camp = pd.read_csv(V05 / "primary_campaigns.tsv", sep="\t")
    conn = json.loads((V05 / "primary_connection.json").read_text())
    aux = set(conn["auxiliary_sites"])
    by_site = {}
    for s in samples:
        by_site.setdefault(s["site"], set()).add(s["trip"])
    return {"sites": [{"site": int(r.site), "lat": rnd(r.latitude, 8), "lon": rnd(r.longitude, 8), "samples": int(r.samples),
                       "campaigns": int(r.campaigns), "trips": sorted(by_site[int(r.site)]), "auxiliary": int(r.site) in aux}
                      for r in sites.itertuples()],
            "campaigns": [{"trip": int(r.trip_num), "samples": int(r.samples), "sites": int(r.sites),
                           "first_date": r.first_date, "last_date": r.last_date} for r in camp.itertuples()],
            "primary_connection": {k: conn[k] for k in ("samples", "usable_mapping_samples", "sites", "core_sites", "auxiliary_sites",
                                                        "site_campaigns", "compartments", "latitude_range", "longitude_range",
                                                        "ambiguous_marginals", "median_overlap_interval_percentage_points",
                                                        "co_admission_ceiling_ratio", "data_boundary")}}


def build_graph() -> dict:
    g = theory.graph()
    adm = []
    for n in g["nodes"]:
        try:
            theory.check_admissible([n["id"]])
            v, why = "admissible", n["physical_action"] or ""
        except ValueError as e:
            v, why = "refused", str(e)
        adm.append({"node": n["id"], "verdict": v, "reason": why})
    layout = {  # hand-tuned layout of the manuscript figure (reproducibility/v06/figures.py); fixed constants are omitted from the drawing
        "xpos": {"K": 36, "E": 68, "gamma_h2": 20, "q_h2": 36, "D_K": 52, "q_co": 68, "gamma_co": 84,
                 "Evap": 8, "I": 19, "c": 30, "T": 41, "psi_e": 53, "b": 64, "phi": 76,
                 "W_max": 10, "theta": 20.5, "psi_o": 31, "psi": 42, "psi_m": 53, "P_cap": 64, "P": 75, "D_h2": 86, "D_co": 98,
                 "W": 13, "N": 29, "A": 42, "k_h2": 58, "k_co": 71, "J_h2": 86, "J_co": 98},
        "row_of": {"ontological": 0, "evidence": 0, "deduced": 1, "exogenous_manipulable": 2, "constitutive": 3, "causal": 4},
        "band_titles": ["Ontology K and evidence E", "Deduced D_K(E)", "Exogenous inputs", "Constitutive maps \u03a6", "Causal mechanisms"],
        "symbols": {"K": "K", "E": "E", "D_K": "D_K(E)", "q_h2": "q_H\u2082", "q_co": "q_CO", "gamma_h2": "\u03b3_H\u2082", "gamma_co": "\u03b3_CO",
                    "I": "I_t", "c": "c_t", "Evap": "Ev_t", "T": "T_t", "phi": "\u03c6", "psi_e": "\u03c8_e", "b": "b", "W_max": "W_max",
                    "theta": "\u03b8", "psi_m": "\u03c8_m", "psi_o": "\u03c8_o", "psi": "\u03c8", "D_h2": "D_H\u2082", "D_co": "D_CO",
                    "A": "A_t", "k_h2": "k_H\u2082", "k_co": "k_CO", "J_h2": "J_H\u2082", "J_co": "J_CO", "P": "P", "P_cap": "P_cap",
                    "W": "W_t", "N": "N_t"}}
    assert set(layout["xpos"]) == {n["id"] for n in g["nodes"] if n["zone"] != "exogenous_fixed"}
    return {**g, "layout": layout, "admissibility": adm, "version_relevance_psi": {
        o: {k: r[k] for k in ("handles", "relevant", "bypass_paths", "reading")}
        for o in ("A", "D_h2", "J_h2", "J_co") for r in [theory.version_relevance("psi", o)]}}


def build_mags() -> dict:
    mag = pd.read_csv(DATA / "mag_evidence.tsv", sep="\t")
    ent = json.loads((DATA / "entailed_memberships.json").read_text())["classes_by_mag"]
    just = json.loads((DATA / "justifications.json").read_text())
    tbox = [readable(a) for a in just["tbox_axioms"]]
    sel = mag[mag.h2_cbb_all_steps | mag.co_form1]
    out = []
    for r in sel.itertuples(index=False):
        kos = set(filter(None, (r.cbb_kos or "").split(";"))) if isinstance(r.cbb_kos, str) else set()
        steps = {}
        for step, alts in CBB_STEPS.items():
            hit = [k for k in alts if k in kos]
            steps[step] = {"present": bool(getattr(r, step)), "found": hit, "accepted": list(alts)}
        route = ("SBPase (K01100/K11532/K01086)" if kos & SBPASE else None,
                 "transaldolase (K00616/K13810)" if kos & TRANSALDOLASE else None)
        classes = ent.get(r.rep, [])
        d = {"id": r.rep, "domain": r.domain, "phylum": r.phylum if isinstance(r.phylum, str) else "",
             "order": r.order if isinstance(r.order, str) else "", "family": r.family if isinstance(r.family, str) else "",
             "genus": r.genus if isinstance(r.genus, str) else "", "completeness": rnd(r.completeness, 4),
             "contamination": rnd(r.contamination, 4), "quality": r.quality_tier, "steps": steps,
             "steps_present": int(sum(s["present"] for s in steps.values())),
             "regeneration_routes": [x for x in route if x],
             "hydrogenase": [g for g, f in (("1l", r.h2_1l), ("1h", r.h2_1h), ("2a", r.h2_2a)) if bool(f)],
             "co_form1": bool(r.co_form1), "classes": classes, "justification": {}}
        for cls in ("H2AutotrophicCandidateMAG", "COCarrierMAG"):
            j = just["classes"].get(cls, {}).get(r.rep)
            if j:
                d["justification"][cls] = {"tbox": j["tbox_axiom_ids"], "abox": [readable(a) for a in j["mag_axioms"]],
                                           "minimal": bool(just["all_verified_subset_minimal"]), "enumeration_complete": j["enumeration_complete"]}
        out.append(d)
    out.sort(key=lambda d: (not ("H2AutotrophicCandidateMAG" in d["classes"]), d["id"]))
    rep = json.loads((DATA / "ontology_report.json").read_text())
    return {"mags": out, "tbox": tbox, "step_labels": STEP_LABEL, "step_order": list(CBB_STEPS),
            "counts": rep["entailed_class_counts"],
            "reasoner": {k: rep[k] for k in ("reasoner", "reasoner_version", "selected_route", "expressivity", "logical_axioms",
                                             "declared_classes", "consistent", "unsatisfiable_classes", "ontology_sha256")},
            "justification_note": "Each justification is a minimal subset of the TBox and of the MAG's own assertions that entails the class "
                                  "membership (all_verified_subset_minimal=%s); other justifications may exist." % just["all_verified_subset_minimal"],
            "all_reference_mags": int(len(mag))}


def build_ontology() -> dict:
    maps = pd.read_csv(DATA / "ontology_mappings.tsv", sep="\t")
    obs = json.loads((DATA / "observation_maps.json").read_text())
    pairs = {}
    for key, p in obs["pairs"].items():
        pairs[key] = {"a": p["a"], "b": p["b"], "relation": p["relation"],
                      "maps": sorted({(m["map"], tuple(m["assumptions"])) for m in p["declared_maps"]})}
        pairs[key]["maps"] = [{"map": m, "assumptions": list(a)} for m, a in pairs[key]["maps"]]
    return {"mappings": maps.fillna("").to_dict("records"),
            "endpoints": {k: {"label": v["label"], "definition": v["definition"], "unit_ucum": v["unit_ucum"],
                              "denominator": v["denominator"]} for k, v in obs["endpoints"].items()},
            "declared_maps": obs["declared_maps"], "undeclared_maps": obs["undeclared_maps"], "pairs": pairs,
            "probe_count": obs["probe_count"], "disjoint_pairs": obs["disjoint_pairs"],
            "ontology_sha256": obs["ontology_sha256"]}


def build_transfer() -> dict:
    ext = json.loads((V05 / "external_transfer.json").read_text())
    namib = pd.read_csv(V05 / "namib_hydration_contrasts.tsv", sep="\t")
    pred = pd.read_csv(V05 / "namib_conditional_predictions.tsv", sep="\t")
    son = pd.read_csv(V05 / "sonoran_hydration_contrasts.tsv", sep="\t")
    lab = pd.read_csv(V05 / "sonoran_lab_observations.tsv", sep="\t")
    fp = pd.read_csv(V05 / "sonoran_field_predictions.tsv", sep="\t")
    names = {"frozen_reference": "Frozen reference", "calibrated_kinetics": "Kinetic coefficient calibrated",
             "calibrated_diffusion_only": "Diffusion-only coefficient calibrated", "training_mean": "Training mean",
             "linear_moisture_temperature": "Linear moisture and temperature", "calibrated_retention": "Soil-response parameters calibrated"}
    order = list(names)
    field = [{"model": names[m["model"]], "key": m["model"], "rmse": rnd(m["rmse"], 4), "mae": rnd(m["mae"], 4)}
             for m in sorted(ext["sonoran"]["field_metrics"], key=lambda m: order.index(m["model"]))]
    per_campaign = []
    for (camp, model), g in fp.groupby(["campaign", "model"]):
        per_campaign.append({"campaign": camp, "model": names.get(model, model),
                             "rmse": rnd(float(((g.observed - g.predicted) ** 2).mean() ** .5), 4)})
    temperate = {}
    for r in ext["temperate_leave_soil_out"]:
        temperate.setdefault(r["soil"].replace("PB_sand", "Sand"), {})[r["model"]] = rnd(r["rmse"], 4)
    lab_means = [{"plot": plot, "gwc": rnd(gwc, 4), "mean": rnd(g.uptake_nmol_m2_s.mean(), 4), "n": int(len(g))}
                 for (plot, gwc), g in lab.groupby(["Treatment", "gwc"])]
    return {
        "namib": {"contrasts": [{k: (None if (isinstance(v, float) and v != v) else rnd(v, 5) if isinstance(v, float) else v)
                                 for k, v in r.items()} for r in namib.to_dict("records")],
                  "conditional": [{"site": r["site"], "bulk_density": rnd(r["assumed_bulk_density_g_cm3"], 3),
                                   "dry_flux": rnd(r["dry_reference_flux"], 4), "wet_flux": rnd(r["wet_reference_flux"], 4),
                                   "direction": int(r["predicted_direction"])} for r in pred.to_dict("records")],
                  "sites": ext["namib"]["sites"], "dry_quantified": ext["namib"]["dry_sites_with_quantified_atmospheric_rates"],
                  "wet_quantified": ext["namib"]["wet_sites_with_quantified_atmospheric_rates"],
                  "reference_failure": ext["namib"]["reference_failure"], "scope": ext["namib"]["prediction_scope"]},
        "sonoran_lab": {"contrasts": [{k: (rnd(v, 4) if isinstance(v, float) else v) for k, v in r.items()} for r in son.to_dict("records")],
                        "model_default": {k: rnd(v, 4) for k, v in ext["sonoran"]["lab_default_predictions"].items()},
                        "assumed_density_g_cm3": ext["sonoran"]["lab_default_assumed_density_g_cm3"],
                        "observed_means": lab_means, "scope": ext["sonoran"]["lab_scope"]},
        "sonoran_field": {"table": field, "per_campaign": per_campaign, "collars": ext["sonoran"]["field_collars"],
                          "plot_campaign_means": ext["sonoran"]["field_plot_campaign_means"],
                          "validation": ext["sonoran"]["field_validation"], "unit": "nmol H2 m-2 s-1"},
        "temperate": {"rmse": temperate, "unit": "bottle rate (min-1)", "models": ["activation_diffusion", "diffusion_only", "training_mean"]},
        "protocol": "reviews/v05/external-evaluation-protocol.md", "retrospective": True}


def build_scenarios() -> dict:
    summ = json.loads((V06 / "paired_regimes_summary.json").read_text())
    tv = json.loads((V06 / "treatment_versions.json").read_text())
    vr = json.loads((V06 / "version_relevance.json").read_text())
    res = json.loads((V06 / "results.json").read_text())
    members = [{k: (rnd(v, 6) if isinstance(v, float) else v) for k, v in m.items()} for m in tv["members"]]
    comp = {}
    for k, v in summ["by_compartment"].items():
        gas, mode, sched, c = k.split("|")
        comp.setdefault(f"{gas}|{mode}|{sched}", {})[c] = {"n": v["n"], **{q: rnd(v[q] * 1e6, 5) for q in ("median", "q1", "q3")}}
    return {"paired_summary": {"samples": summ["samples"], "closed_gates": summ["closed_gates"],
                               "by_compartment_umol_m2": comp, "contrasts": summ["contrasts"],
                               "gamma_proportional_quantiles": summ["gamma_proportional_quantiles"],
                               "reference_uptake_mol_m2": summ["reference_uptake_mol_m2"]},
            "treatment_versions": {"summary": tv["summary"], "members": members},
            "version_relevance": vr, "results": res}


def build_key_numbers(paper: dict, ont: dict, conn: dict, ref: dict) -> dict:
    return {"metagenomes": 551, "mapped_samples": int(paper["MappedSamples"]), "campaigns": 5, "site_campaign_visits": 217,
            "sites": 61, "core_sites": 59, "auxiliary_sites": 2, "reference_mags": int(paper["ReferenceMAGs"]),
            "cbb_complete": int(paper["CBBFull"]), "h2_cbb_candidates": int(paper["HtwoCBB"]),
            "h2_cbb_hq": int(paper["HtwoCBBStrict"]), "dual_h2_co": int(paper["DualHtwoCO"]),
            "h2_carriers": int(paper["HtwoCarriers"]), "co_carriers": ont["counts"]["COCarrierMAG"],
            "ambiguous_samples": int(paper["AmbiguousSamples"]), "median_interval_width_pp": rnd(conn["median_overlap_interval_percentage_points"], 3),
            "graph_nodes": int(paper["GraphNodes"]), "admissible_nodes": int(paper["AdmissibleNodes"]),
            "admissible_inputs": int(paper["AdmissibleInputs"]), "admissible_states": int(paper["AdmissibleStates"]),
            "refused_nodes": int(paper["RefusedNodes"]), "closed_co_gates": int(paper["ClosedCOGates"]),
            "closed_h2_gates": int(paper["ClosedHtwoGates"]), "closed_co_zero_effects": int(paper["ClosedCOZeroEffects"]),
            "ref_h2_single_umol": float(paper["RefHtwoSingle"]), "ref_h2_four_umol": float(paper["RefHtwoFour"]),
            "ref_h2_daily_umol": rnd(ref["daily_pulses"]["h2"] * 1e6, 4), "ref_co_single_umol": float(paper["RefCOSingle"]),
            "ref_co_four_umol": float(paper["RefCOFour"]), "version_members": int(paper["VersionMembers"]),
            "single_first": int(paper["VersionSingleFirst"]), "four_first": int(paper["VersionFourFirst"]),
            "daily_lowest": int(paper["VersionDailyLowest"]), "co_admission_ceiling_ratio": rnd(conn["co_admission_ceiling_ratio"], 3)}


def build_identification() -> dict:
    """Compact view of reviews/v06/identification.json for the page: claim, plain and closure tests, derivation chain, adjustment sets."""
    r = json.loads((V06 / "identification.json").read_text())

    def steps(closure):
        return [{"rule": d["rule"], "premises": d["premises"], "conclusion": d["node"], "value": d["value"]}
                for d in closure["derivations"] if d["premises"]]

    def constants(closure):
        return sorted(d["node"] for d in closure["derivations"] if not d["premises"])

    dn, ox, fb = r["deductive_null"], r["oxygen_context"], r["field_backdoor"]
    ec = ox["executable_check"]
    return {"note": r["note"], "queries": [
        {"query": "deductive_null", "claim": dn["claim"], "rule": dn["rule"],
         "plain_dsep_holds": dn["result"]["graph_test"], "closure_dsep_holds": dn["result"]["graph_test_Z_star_csi"],
         "open_path": dn["result"]["open_path_plain"], "closure_steps": steps(dn["closure"]), "closure_constants": constants(dn["closure"]),
         "adjustment_sets": None, "numerical_check": None, "verdict": dn["verdict"]},
        {"query": "oxygen_context", "claim": ox["claim"], "rule": "Rule 3",
         "plain_dsep_holds": ox["without_stratum"]["graph_test"], "closure_dsep_holds": ox["stratum"]["graph_test_Z_star_csi"],
         "open_path": ox["without_stratum"]["open_path_plain"], "closure_steps": steps(ox["closure"]), "closure_constants": constants(ox["closure"]),
         "adjustment_sets": None,
         "numerical_check": {"h2_invariant_gamma_co_0": ec["gamma_co=0"]["h2_invariant"], "h2_invariant_gamma_co_1": ec["gamma_co=1"]["h2_invariant"],
                             "do_water_still_changes_h2": not ox["control_water_still_affects_h2"]["holds"]},
         "verdict": ox["verdict"]},
        {"query": "field_backdoor", "claim": fb["estimand"], "rule": "back-door",
         "plain_dsep_holds": False, "closure_dsep_holds": None,
         "open_path": fb["open_backdoor_path_without_adjustment"], "closure_steps": [], "closure_constants": [],
         "adjustment_sets": {k: fb[k]["minimal_sets"] for k in ("all_measured", "T_unmeasured", "T_and_Ev_unmeasured", "phi_unmeasured")},
         "positivity_caveat": fb["positivity_caveat"], "limits": fb["limits"], "numerical_check": None, "verdict": fb["verdict"]}]}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", type=Path, default=DEFAULT_OUT)
    out = ap.parse_args().out.resolve()
    if out.exists():
        shutil.rmtree(out)
    (out / "py" / "eq_biodt" / "v03").mkdir(parents=True)

    version = project_version()
    for rel in VERBATIM:
        shutil.copyfile(SRC / "eq_biodt" / rel, out / "py" / "eq_biodt" / rel)
    shutil.copyfile(Path(__file__).with_name("live_scenario.py"), out / "py" / "eq_biodt" / "live.py")
    (out / "py" / "eq_biodt" / "__init__.py").write_text(f'"""EQ-BioDT web runtime (minimal package for the browser)."""\n__version__ = "{version}"\n')
    (out / "py" / "eq_biodt" / "v03" / "__init__.py").write_text(
        f'"""EQ-BioDT causal-ontological twin (web runtime; released process and theory modules are copied verbatim)."""\nVERSION = "{version}"\n')

    paper = tex_numbers()
    samples, medians, columns = build_samples()
    sites = build_sites(samples)
    mags = build_mags()
    ontology = build_ontology()
    scenarios = build_scenarios()
    commit = git("rev-parse", "HEAD")
    dirty = bool(git("status", "--porcelain", "--untracked-files=no", "--", "src", "data/v03", "reviews", "paper/tables"))
    bundle = {
        "schema": 1, "version": version, "title": "EQ-BioDT v0.6",
        "provenance": {"source_repository": "https://github.com/bio-ontology-research-group/empty-quarter-digital-twin",
                       "commit": commit, "branch": git("rev-parse", "--abbrev-ref", "HEAD"), "tree_dirty": dirty,
                       "generated_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
                       "ontology_sha256": ontology["ontology_sha256"],
                       "inputs_sha256": {str(p.relative_to(ROOT)): sha256(p) for p in sorted(
                           [*DATA.glob("*.tsv"), *DATA.glob("*.json"), V06 / "paired_regimes_by_sample.tsv", V06 / "treatment_versions.json",
                            V05 / "external_transfer.json"])}},
        "notice": ("Conditional process model. It reports what follows from the evidence, the ontology and stated mechanisms; "
                   "reaction rate coefficients are supplied scenario parameters, not fitted to Rub' al-Khali flux data, and no net carbon removal is estimated."),
        "key_numbers": build_key_numbers(paper, mags, sites["primary_connection"], scenarios["paired_summary"]["reference_uptake_mol_m2"]),
        "paper_numbers": paper,
        "graph": build_graph(), "identification": build_identification(),
        "medians": medians, "sample_columns": columns, "samples": samples, **{k: sites[k] for k in ("sites", "campaigns")},
        "primary_connection": sites["primary_connection"],
        "scenarios": scenarios, "candidates": mags, "ontology": ontology, "transfer": build_transfer(),
        "reference_schedule": {"horizon_days": 14, "total_water_mm": 12, "initial_water_mm": 1.8, "evaporation_mm": 1.0,
                               "temperature_c": 25.0, "k_h2_s": 0.001, "k_co_s": 0.001,
                               "schedules": {"intact": [], "single_pulse": [0], "four_pulses": [0, 3, 7, 10],
                                             "daily_pulses": list(range(14))}},
    }
    text = json.dumps(bundle, separators=(",", ":"), allow_nan=False, ensure_ascii=False)
    (out / "bundle.json").write_text(text)

    files = sorted(p for p in out.rglob("*") if p.is_file() and p.name != "manifest.json")
    manifest = {"version": version, "commit": commit, "tree_dirty": dirty, "generated_utc": bundle["provenance"]["generated_utc"],
                "files": {str(p.relative_to(out)): {"sha256": sha256(p), "bytes": p.stat().st_size} for p in files},
                "verbatim_from_src": {rel: sha256(SRC / "eq_biodt" / rel) for rel in VERBATIM},
                "pyodide": "see /pyodide/MANIFEST.json"}
    (out / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(f"wrote {out}: bundle.json {len(text) / 1e6:.2f} MB, {len(files)} files, commit {commit[:8]}{' (dirty tree)' if dirty else ''}")


if __name__ == "__main__":
    main()
