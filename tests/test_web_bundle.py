"""The web bundle of the public Twin page: faithful copies, numbers equal to reviews/ outputs, schema."""
import hashlib
import importlib.util
import json
import re
import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
EXPORTER = ROOT / "scripts" / "web" / "export_web_bundle.py"
LIVE = ROOT / "scripts" / "web" / "live_scenario.py"


@pytest.fixture(scope="module")
def out(tmp_path_factory):
    d = tmp_path_factory.mktemp("web") / "eq_biodt_v06"
    subprocess.run([sys.executable, str(EXPORTER), "--out", str(d)], check=True, cwd=ROOT)
    return d


@pytest.fixture(scope="module")
def bundle(out):
    return json.loads((out / "bundle.json").read_text())


@pytest.fixture(scope="module")
def live():
    spec = importlib.util.spec_from_file_location("eq_live", LIVE)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def test_python_files_are_verbatim_copies(out):
    for rel in ("v03/process.py", "v03/theory.py"):
        assert (out / "py" / "eq_biodt" / rel).read_bytes() == (ROOT / "src" / "eq_biodt" / rel).read_bytes()
    assert (out / "py" / "eq_biodt" / "live.py").read_bytes() == LIVE.read_bytes()
    init = (out / "py" / "eq_biodt" / "v03" / "__init__.py").read_text()
    assert re.search(r'^VERSION = "0\.6\.\d+', init, re.M)


def test_manifest_digests(out):
    manifest = json.loads((out / "manifest.json").read_text())
    assert {"py/eq_biodt/__init__.py", "py/eq_biodt/v03/process.py", "py/eq_biodt/v03/theory.py", "bundle.json"} <= set(manifest["files"])
    for rel, meta in manifest["files"].items():
        data = (out / rel).read_bytes()
        assert hashlib.sha256(data).hexdigest() == meta["sha256"] and len(data) == meta["bytes"]
    assert len(manifest["commit"]) == 40


def test_schema_and_size(out, bundle):
    for key in ("version", "provenance", "notice", "key_numbers", "graph", "samples", "sites", "campaigns", "medians", "scenarios",
                "candidates", "ontology", "transfer", "reference_schedule"):
        assert key in bundle
    assert (out / "bundle.json").stat().st_size < 3_000_000
    assert {"nodes", "edges", "zones", "admissibility"} <= set(bundle["graph"])
    assert set(bundle["scenarios"]) >= {"paired_summary", "treatment_versions", "version_relevance"}
    assert bundle["key_numbers"]["metagenomes"] == 551 and len(bundle["samples"]) == 551
    assert len(bundle["sites"]) == 61 and sum(s["auxiliary"] for s in bundle["sites"]) == 2


def test_key_numbers_equal_manuscript_macros(bundle):
    macros = {}
    for f in (ROOT / "paper" / "tables").glob("v0*_numbers.tex"):
        macros.update(re.findall(r"\\newcommand\{\\([A-Za-z]+)\}\{([^}]*)\}", f.read_text()))
    k = bundle["key_numbers"]
    for name, key in (("ReferenceMAGs", "reference_mags"), ("CBBFull", "cbb_complete"), ("HtwoCBB", "h2_cbb_candidates"),
                      ("HtwoCBBStrict", "h2_cbb_hq"), ("DualHtwoCO", "dual_h2_co"), ("AmbiguousSamples", "ambiguous_samples"),
                      ("ClosedCOGates", "closed_co_gates"), ("GraphNodes", "graph_nodes"), ("AdmissibleNodes", "admissible_nodes"),
                      ("VersionSingleFirst", "single_first"), ("VersionFourFirst", "four_first"), ("MappedSamples", "mapped_samples")):
        assert k[key] == int(macros[name]), name
    assert (k["ref_h2_single_umol"], k["ref_h2_four_umol"]) == (float(macros["RefHtwoSingle"]), float(macros["RefHtwoFour"]))
    assert round(k["ref_h2_daily_umol"], 2) == 2.87            # "Daily additions totalling 12 mm yield 2.87 umol H2 m-2"
    assert k["co_carriers"] == 198 and k["h2_carriers"] == 1094
    assert bundle["candidates"]["counts"]["H2AutotrophicCandidateMAG"] == 136


def test_samples_equal_paired_regimes(bundle):
    t = pd.read_csv(ROOT / "reviews/v06/paired_regimes_by_sample.tsv", sep="\t").set_index("sample")
    mapped = [s for s in bundle["samples"] if s["mapped"]]
    assert len(mapped) == len(t) == 550
    for s in mapped:
        r = t.loc[s["id"]]
        for mode in ("gate", "proportional"):
            assert s["gamma"][mode]["h2"] == pytest.approx(r[f"gamma_h2_{mode}"], rel=1e-5, abs=1e-12)
            assert s["gamma"][mode]["co"] == pytest.approx(r[f"gamma_co_{mode}"], rel=1e-5, abs=1e-12)
            for gas in ("h2", "co"):
                for i, sched in enumerate(("intact", "single_pulse", "four_pulses", "daily_pulses")):
                    assert s["uptake_umol_m2"][mode][gas][i] == pytest.approx(r[f"{gas}_{mode}_{sched}"] * 1e6, rel=1e-5, abs=1e-12)
                for i, sched in enumerate(("single_pulse", "four_pulses", "daily_pulses")):
                    assert s["effect_umol_m2"][mode][gas][i] == pytest.approx(r[f"{gas}eff_{mode}_{sched}"] * 1e6, rel=1e-5, abs=1e-12)
    closed = json.loads((ROOT / "reviews/v06/results.json").read_text())["paired"]["closed_gates"]["co_closed"]
    assert sorted(s["id"] for s in mapped if s["gamma"]["gate"]["co"] == 0) == sorted(closed)
    assert all(s["effect_umol_m2"][m]["co"] == [0, 0, 0] for s in mapped if s["id"] in closed for m in ("gate", "proportional"))
    med = bundle["medians"]
    for s in mapped:       # gamma_proportional = q / median(q)
        assert s["gamma"]["proportional"]["h2"] == pytest.approx(s["q"]["H2CarrierMAG"] / med["H2CarrierMAG"], rel=1e-4, abs=1e-9)


def test_overlap_interval_positive_in_549(bundle):
    wide = 0
    for s in bundle["samples"]:
        if s["mapped"]:
            a, b = s["q"]["H2CarrierMAG"], s["q"]["CBBCompleteMAG"]
            lo, hi = max(0, a + b - 1), min(a, b)
            wide += hi - lo > 1e-12
    assert wide == bundle["key_numbers"]["ambiguous_samples"] == 549


def test_treatment_versions_equal_review(bundle):
    t = pd.read_csv(ROOT / "reviews/v06/treatment_versions.tsv", sep="\t")
    m = bundle["scenarios"]["treatment_versions"]["members"]
    assert len(m) == 24
    for row, r in zip(m, t.itertuples()):
        assert row["uptake_single_pulse"] == pytest.approx(r.uptake_single_pulse, rel=1e-5)
        assert row["ranking"] == r.ranking
    s = bundle["scenarios"]["treatment_versions"]["summary"]
    assert (s["members_where_single_beats_four"], s["members_where_four_beats_single"], s["members_where_daily_is_lowest"]) == (22, 2, 24)


def test_candidates(bundle):
    c = bundle["candidates"]["mags"]
    h2 = [m for m in c if "H2AutotrophicCandidateMAG" in m["classes"]]
    co = [m for m in c if "COCarrierMAG" in m["classes"]]
    assert len(h2) == 136 and len(co) == 198
    assert all(m["steps_present"] == 11 for m in h2)
    assert all("H2AutotrophicCandidateMAG" in m["justification"] for m in h2)
    assert all("COCarrierMAG" in m["justification"] for m in co)
    assert len(c) == 136 + 198 - 34 and len({m["id"] for m in c}) == len(c)   # 34 dual candidates are in both groups


def test_transfer_equals_review_and_paper(bundle):
    tr = bundle["transfer"]
    n = pd.read_csv(ROOT / "reviews/v05/namib_hydration_contrasts.tsv", sep="\t").set_index("site")
    for r in tr["namib"]["contrasts"]:
        if r["fold"] is not None:
            assert r["fold"] == pytest.approx(n.loc[r["site"], "fold"], rel=1e-4)
            assert r["fold_lcl"] == pytest.approx(n.loc[r["site"], "fold_lcl"], rel=1e-4)
    c14 = {r["site"]: r for r in tr["namib"]["contrasts"]}
    assert (round(c14["C14-2"]["fold"], 1), round(c14["C14-2"]["fold_lcl"], 1), round(c14["C14-2"]["fold_ucl"], 1)) == (41.3, 31.1, 56.0)
    assert (round(c14["C14-4"]["fold"], 1), round(c14["C14-4"]["fold_lcl"], 1), round(c14["C14-4"]["fold_ucl"], 1)) == (19.4, 15.6, 24.0)
    table = (ROOT / "paper/tables/v05_external.tex").read_text()
    for row in tr["sonoran_field"]["table"]:
        assert f"{row['model']} & {row['rmse']:.3f} & {row['mae']:.3f}" in table.replace("0.6921", "0.692").replace("0.5291", "0.529") or \
            f"{row['model']} & {round(row['rmse'], 3):.3f} & {round(row['mae'], 3):.3f}" in table
    assert tr["sonoran_field"]["collars"] == 64 and tr["sonoran_field"]["plot_campaign_means"] == 16
    ext = json.loads((ROOT / "reviews/v05/external_transfer.json").read_text())["temperate_leave_soil_out"]
    for r in ext:
        soil = r["soil"].replace("PB_sand", "Sand")
        assert tr["temperate"]["rmse"][soil][r["model"]] == pytest.approx(r["rmse"], rel=1e-3)


def test_graph_matches_theory(bundle):
    from eq_biodt.v03 import theory
    g = theory.graph()
    assert [n["id"] for n in bundle["graph"]["nodes"]] == [n["id"] for n in g["nodes"]]
    adm = {a["node"]: a["verdict"] for a in bundle["graph"]["admissibility"]}
    assert sum(v == "admissible" for v in adm.values()) == 11 and len(adm) == 42   # v0.7 adds C_h2, C_co


def test_identification_equals_review(bundle):
    r = json.loads((ROOT / "reviews" / "v06" / "identification.json").read_text())
    q = {x["query"]: x for x in bundle["identification"]["queries"]}
    assert list(q) == ["deductive_null", "oxygen_context", "field_backdoor"]
    dn, ox, fb = r["deductive_null"], r["oxygen_context"], r["field_backdoor"]
    assert q["deductive_null"]["claim"] == dn["claim"] and q["deductive_null"]["rule"] == dn["rule"]
    assert q["deductive_null"]["verdict"] == dn["verdict"]
    assert q["deductive_null"]["open_path"] == dn["result"]["open_path_plain"]
    assert q["deductive_null"]["plain_dsep_holds"] is False and q["deductive_null"]["closure_dsep_holds"] is True
    chain = [(s["conclusion"], s["premises"]) for s in q["deductive_null"]["closure_steps"]]
    assert chain == [("gamma_co", ["q_co"]), ("k_co@0", ["gamma_co"]), ("J_co@0", ["k_co@0"])]
    assert all(s["rule"] for s in q["deductive_null"]["closure_steps"])
    assert q["oxygen_context"]["claim"] == ox["claim"] and q["oxygen_context"]["verdict"] == ox["verdict"]
    assert q["oxygen_context"]["open_path"] == ox["without_stratum"]["open_path_plain"]
    assert q["oxygen_context"]["plain_dsep_holds"] is False and q["oxygen_context"]["closure_dsep_holds"] is True
    assert q["oxygen_context"]["numerical_check"] == {"h2_invariant_gamma_co_0": True, "h2_invariant_gamma_co_1": False, "do_water_still_changes_h2": True}
    assert q["field_backdoor"]["claim"] == fb["estimand"] and q["field_backdoor"]["verdict"] == fb["verdict"]
    assert q["field_backdoor"]["open_path"] == fb["open_backdoor_path_without_adjustment"]
    assert q["field_backdoor"]["adjustment_sets"] == {k: fb[k]["minimal_sets"] for k in ("all_measured", "T_unmeasured", "T_and_Ev_unmeasured", "phi_unmeasured")}
    assert q["field_backdoor"]["adjustment_sets"]["T_and_Ev_unmeasured"] == []
    assert q["field_backdoor"]["positivity_caveat"] == fb["positivity_caveat"]


# ---- the browser driver, run natively against the same released modules ----
def test_live_reference_daily_matches_paper(live):
    r = json.loads(live.run(json.dumps({"mode": "reference", "members": False,
                                        "intervention": {"water": {"total_mm": 12, "version": "daily_pulses"}}})))
    assert round(r["effect"]["h2_umol_m2"], 2) == 2.87


def test_live_closed_gate_is_zero_by_deduction(live, bundle):
    sample = next(s for s in bundle["samples"] if s["id"] == "51PRr3")
    r = json.loads(live.run(json.dumps({"mode": "gate", "shares": sample["q"], "medians": bundle["medians"], "members": False,
                                        "intervention": {"water": {"total_mm": 12, "version": "single_pulse"}}})))
    assert r["nulls"]["co"]["zero_by_deduction"] and r["nulls"]["co"]["numerically_confirmed"]
    assert r["effect"]["co_umol_m2"] == 0 and r["effect"]["h2_umol_m2"] == pytest.approx(1102.3, abs=0.05)


def test_live_equals_precomputed_sample_regimes(live, bundle):
    s = next(s for s in bundle["samples"] if s["mapped"] and s["gamma"]["proportional"]["h2"] not in (0, 1))
    for i, ver in enumerate(("single_pulse", "four_pulses", "daily_pulses")):
        r = json.loads(live.run(json.dumps({"mode": "proportional", "shares": s["q"], "medians": bundle["medians"], "members": False,
                                            "intervention": {"water": {"total_mm": 12, "version": ver}}})))
        assert r["effect"]["h2_umol_m2"] == pytest.approx(s["effect_umol_m2"]["proportional"]["h2"][i], rel=1e-5)
        assert r["effect"]["co_umol_m2"] == pytest.approx(s["effect_umol_m2"]["proportional"]["co"][i], rel=1e-5)


def test_live_member_grid_equals_planning_and_versions(live, bundle):
    from eq_biodt.v03.planning import parameter_ensemble
    base = live.base_parameters()
    assert live.parameter_members(base) == parameter_ensemble(base)
    r = json.loads(live.run(json.dumps({"mode": "reference", "intervention": {"water": {"total_mm": 12, "version": "single_pulse"}}})))
    ref = {m["member"]: m for m in bundle["scenarios"]["treatment_versions"]["members"]}
    for m in r["members"]:
        assert m["h2_intervened"] * 1e-6 == pytest.approx(ref[m["member"]]["uptake_single_pulse"], rel=1e-5)


def test_live_refusals_and_psi_versions(live):
    v = {x["target"]: x for x in json.loads(live.admissible(json.dumps(["J_h2", "q_h2", "D_K", "theta", "psi", "W", "I"])))}
    assert [v[k]["verdict"] for k in ("J_h2", "q_h2", "D_K", "theta", "psi")] == ["refused"] * 5
    assert "deduced zone" in v["q_h2"]["message"] and "causal outcome" in v["J_h2"]["message"]
    assert v["W"]["verdict"] == v["I"]["verdict"] == "admissible"
    with pytest.raises(ValueError, match="specify version"):
        live.run(json.dumps({"mode": "reference", "baseline": {"water_mm": 12, "version": "single_pulse"},
                             "intervention": {"state": {"target": "psi_mpa", "value": -1.0, "day": 1}}}))
    r = json.loads(live.run(json.dumps({"mode": "reference", "members": False, "baseline": {"water_mm": 12, "version": "single_pulse"},
                                        "intervention": {"state": {"target": "psi_mpa", "value": -1.0, "version": "matric", "day": 1}}})))
    assert round(r["psi_versions"]["matric"]["state_after"]["water_mm"], 1) == 7.1       # paper: 7.1 mm
    assert round(r["psi_versions"]["osmotic"]["nacl_mol_l_after"], 3) == 0.182            # paper: 0.182 mol/L
    assert r["relevance"]["J_h2"]["relevant"] and not r["relevance"]["A"]["relevant"]
