"""v0.6 analyses of the stratified causal-ontological model (deterministic; bootstrap seed 42).

Writes reviews/v06/{paired_regimes_by_sample.tsv, paired_regimes_summary.json,
treatment_versions.json, treatment_versions.tsv, version_relevance.json, admissibility.tsv, results.json}.
"""
from dataclasses import replace
import importlib.util
import json
from pathlib import Path
import sys
import time

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "reviews/v06"
sys.path.insert(0, str(ROOT / "src"))
from eq_biodt.v03 import community as comm, theory                       # noqa: E402
from eq_biodt.v03.interventions import community_factors, intervene       # noqa: E402
from eq_biodt.v03.planning import compare, parameter_ensemble             # noqa: E402
from eq_biodt.v03.process import Parameters, State, gas_exchange, simulate  # noqa: E402

DATA = ROOT / "data/v03"
SCHEDULES = ("intact", "single_pulse", "four_pulses", "daily_pulses")
COMPARTMENTS = ("Surface", "Deep", "PlantRoot")
MODES = ("gate", "proportional")


def reference_schedules():
    """The four reference schedules of reproducibility/v03/reproduce.py (14 d, 12 mm total)."""
    out = []
    for name, n in zip(SCHEDULES, (0, 1, 4, 14)):
        days = {int(i * 14 / n) for i in range(n)} if n else set()
        forcing = [{"water_mm": 12 / n if day in days else 0, "nacl_mol_l": 0., "evaporation_mm": 1., "temp_c": 25.}
                   for day in range(14)]
        out.append({"name": name, "forcing": forcing, "disturbance_m2": 0. if n == 0 else 1.})
    return out


def load_cluster_ci():
    spec = importlib.util.spec_from_file_location("v03_audit", ROOT / "reproducibility/v03/audit.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod.cluster_ci


def gas_uptake(result):
    rows = result["trajectory"]
    return (sum(r["gas"]["h2_mol_m2_s"] * 86400 for r in rows), sum(r["gas"]["co_mol_m2_s"] * 86400 for r in rows))


def paired_regimes(schedules, base):
    table = comm._load()["table"]
    panel = comm._load()["panel"].set_index("sample")
    cache = {}
    rows = []
    for sample in sorted(table.index):
        shares = comm.carrier_shares(sample)["shares"]
        row = {"sample": sample, "site": str(panel.loc[sample, "site"]), "trip": int(panel.loc[sample, "trip_num"]),
               "compartment": panel.loc[sample, "compartment"],
               "q_h2": shares["H2CarrierMAG"], "q_co": shares["COCarrierMAG"]}
        for mode in MODES:
            g = comm.gamma(shares, mode)
            row[f"gamma_h2_{mode}"], row[f"gamma_co_{mode}"] = g["h2"], g["co"]
            key = (g["h2"], g["co"])
            if key not in cache:
                p = replace(base, gamma_h2=g["h2"], gamma_co=g["co"])
                cache[key] = {s["name"]: gas_uptake(simulate(State(1.8), s["forcing"], p)) for s in schedules}
            up = cache[key]
            for s in SCHEDULES:
                row[f"h2_{mode}_{s}"], row[f"co_{mode}_{s}"] = up[s]
                if s != "intact":
                    row[f"h2eff_{mode}_{s}"] = up[s][0] - up["intact"][0]
                    row[f"coeff_{mode}_{s}"] = up[s][1] - up["intact"][1]
        rows.append(row)
    return pd.DataFrame(rows), len(cache)


def iqr(x):
    return [float(v) for v in np.quantile(x, [.25, .5, .75])]


def summarise(df, panel_checks):
    out = {"samples": int(len(df)), "by_compartment": {}}
    for mode in MODES:
        for s in SCHEDULES[1:]:
            for gas in ("h2", "co"):
                col = f"{gas}eff_{mode}_{s}"
                for comp in COMPARTMENTS:
                    x = df.loc[df.compartment.eq(comp), col].to_numpy()
                    q1, med, q3 = iqr(x)
                    out["by_compartment"][f"{gas}|{mode}|{s}|{comp}"] = {"n": int(len(x)), "median": med, "q1": q1, "q3": q3}
    return out


def contrasts(df, cluster_ci):
    rng = np.random.default_rng(42)
    res = []
    for mode, col in (("proportional", "h2eff_proportional_single_pulse"), ("gate", "h2eff_gate_single_pulse")):
        grouped = df.groupby(["site", "trip", "compartment"])[col].mean().unstack()
        rng = np.random.default_rng(42)
        for a, b in (("Deep", "PlantRoot"), ("Deep", "Surface"), ("PlantRoot", "Surface")):
            paired = grouped[[a, b]].dropna().reset_index()
            paired["difference"] = paired[a] - paired[b]
            lo, hi = cluster_ci(paired, ["difference"], lambda x: x[:, 0], rng)
            res.append({"coupling": mode, "contrast": f"{a} minus {b}", "matched_site_campaigns": int(len(paired)),
                        "sites": int(paired.site.nunique()), "mean_difference_mol_m2": float(paired.difference.mean()),
                        "site_bootstrap_95_percentile": [lo, hi], "bootstrap_draws": 2000, "seed": 42})
    return res


def treatment_versions_export(schedules, base):
    comp = compare(State(1.8), schedules, parameter_ensemble(base))
    tv = comp["treatment_versions"]
    names = [s["name"] for s in schedules]
    up = {(r["member"], r["scenario"]): r["h2_uptake_mol_m2"] for r in comp["outcomes"]}
    n = len(comp["parameters"])
    members = []
    group = next(g for g in tv["groups"] if abs(g["total_water_mm"] - 12) < 1e-9)
    for m in group["members"]:
        i = m["member"]
        par = comp["parameters"][i]
        members.append({"member": i, "retention_b": par["retention_b"], "dry_threshold_mpa": par["dry_threshold_mpa"],
                        "activation_shape": par["activation_shape"], "k_h2_s": par["k_h2_s"],
                        "control": up[(i, tv["control"])],
                        **{f"uptake_{v}": m["uptake_mol_m2"][v] for v in group["versions"]},
                        "ranking": ">".join(m["ranking"]), "top": m["ranking"][0],
                        **{f"effect_{v}": m["uptake_mol_m2"][v] - up[(i, tv["control"])] for v in group["versions"]}})
    mem = pd.DataFrame(members)
    effects = [mem[f"effect_{v}"] for v in group["versions"]]
    summary = {"members": n, "total_water_mm": group["total_water_mm"], "versions": group["versions"], "control": tv["control"],
               "ranking_changes_across_members": group["ranking_changes_across_members"],
               "version_relevant": tv["version_relevant"],
               "distinct_rankings": sorted(mem.ranking.unique()), "ranking_counts": mem.ranking.value_counts().to_dict(),
               "members_where_single_beats_four": int((mem.uptake_single_pulse > mem.uptake_four_pulses).sum()),
               "members_where_four_beats_single": int((mem.uptake_four_pulses > mem.uptake_single_pulse).sum()),
               "members_where_daily_is_lowest": int(((mem.uptake_daily_pulses < mem.uptake_single_pulse) & (mem.uptake_daily_pulses < mem.uptake_four_pulses)).sum()),
               "effect_vs_control_range_mol_m2": tv["effect_vs_control_range_mol_m2"],
               "effect_range_by_version_mol_m2": {v: [float(mem[f"effect_{v}"].min()), float(mem[f"effect_{v}"].max())] for v in group["versions"]},
               "pairwise": [{k: r[k] for k in ("a", "b", "uptake_difference_range_mol_m2", "rank_reversal")} for r in comp["pairwise_uptake_comparisons"]],
               "pareto_in_every_member": comp["pareto_in_every_member"]}
    return {"summary": summary, "members": members}, mem


def version_relevance_demo(base):
    """do(psi) at the start of day 1, after a 12 mm pulse on day 0, under the matric and osmotic versions.

    Instantaneous A, D_h2 and J_h2 are evaluated at the intervened state; the 3-day uptake uses intervene() and is
    reported only where the process model stays inside its ideal-dilute NaCl domain on all three days."""
    forcing = [{"water_mm": 12. if d == 0 else 0., "nacl_mol_l": 0., "evaporation_mm": 1., "temp_c": 25.} for d in range(14)]
    traj = simulate(State(1.8), forcing, base)["trajectory"]
    state = State(traj[0]["water_mm"], traj[0]["salt_mol_m2"])
    out = {"state_entering_day_1": {"water_mm": state.water_mm, "salt_mol_m2": state.salt_mol_m2},
           "relevance": {o: theory.version_relevance("psi", o) for o in ("A", "D_h2", "J_h2", "J_co")},
           "targets": []}
    reference = None
    for psi in (-1.0, -0.4):
        entry = {"psi_mpa": psi, "versions": {}}
        for version in ("matric", "osmotic"):
            t = theory.translate_constituted("psi", psi, state, base, version, 25.)
            w, n = t["water_mm"], t["salt_mol_m2"]
            theta = w / (1000 * base.depth_m)
            psi_o = -theory._osmotic_coefficient(base, 25.) * n / w
            g = gas_exchange(theta, 25., base, osmotic_mpa=psi_o)
            v = {"handle": t["handle"], "water_mm": w, "salt_mol_m2": n, "theta": theta, "nacl_mol_l": n / w,
                 "psi_total_mpa": g["total_water_potential_mpa"], "A": g["activation"],
                 "D_h2_m2_s": base.d_h2_m2_s * g["gas_diffusivity_ratio"], "J_h2_mol_m2_s": g["h2_mol_m2_s"],
                 "within_ideal_dilute_domain": (n / w) <= base.osmotic_domain_mol_l}
            if v["within_ideal_dilute_domain"]:
                r = intervene(state, forcing[1:4], base, {"state": {"0": {"psi_mpa": psi, "version": version}}}, DATA)
                v["h2_uptake_3d_mol_m2"] = r["gas_uptake"]["h2_mol_m2"]["intervened"]
                v["h2_uptake_3d_baseline_mol_m2"] = r["gas_uptake"]["h2_mol_m2"]["baseline"]
                v["intervene_checks_passed"] = all(r["checks"].values())
                v["intervene_reports_version_relevance"] = {k: x["relevant"] for k, x in r["version_relevance"].items()}
            else:
                v["h2_uptake_3d_mol_m2"] = None
                v["note"] = "NaCl concentration exceeds the ideal-dilute domain (0.1 mol/L): the process model abstains"
            entry["versions"][version] = v
        a = [entry["versions"][k]["A"] for k in ("matric", "osmotic")]
        j = [entry["versions"][k]["J_h2_mol_m2_s"] for k in ("matric", "osmotic")]
        entry["A_equal"] = abs(a[0] - a[1]) < 1e-9
        entry["J_h2_ratio_matric_over_osmotic"] = j[0] / j[1]
        out["targets"].append(entry)
    return out


def admissibility():
    rows = []
    for n in theory.graph()["nodes"]:
        try:
            theory.check_admissible([n["id"]])
            verdict, reason = "admissible", ""
        except ValueError as e:
            verdict, reason = "refused", str(e)
        rows.append({"node": n["id"], "label": n["label"], "zone": n["zone"], "manipulable": n["manipulable"],
                     "verdict": verdict, "reason": reason, "physical_action": n["physical_action"] or ""})
    return pd.DataFrame(rows)


def main():
    t0 = time.time()
    OUT.mkdir(parents=True, exist_ok=True)
    base = Parameters(k_h2_s=.001, k_co_s=.001)
    schedules = reference_schedules()
    cluster_ci = load_cluster_ci()

    df, distinct = paired_regimes(schedules, base)
    df.to_csv(OUT / "paired_regimes_by_sample.tsv", sep="\t", index=False, float_format="%.12g")

    # cross-check one sample against the public intervene() interface
    cross = []
    for sample in (df["sample"].iloc[0], df.loc[df.q_co.idxmax(), "sample"]):
        for mode in MODES:
            r = intervene(State(1.8), schedules[1]["forcing"], base, {"community": {"sample": sample}, "coupling": mode,
                                                                       "forcing": {"0": {"water_mm": 0.}}}, DATA)
            cross.append(r["checks"]["gamma_equal_across_arms"])
    closed = comm.closed_gates()
    co_closed = set(closed["co_closed"])
    h2_closed = set(closed["h2_closed"])
    closed_rows = df[df["sample"].isin(co_closed)]
    co_effect_cols = [c for c in df if c.startswith("coeff_")]
    co_zero = bool((closed_rows[co_effect_cols] == 0).all().all()) if len(closed_rows) else True
    co_uptake_cols = [c for c in df if c.startswith("co_gate_")]
    co_uptake_zero = bool((closed_rows[co_uptake_cols] == 0).all().all()) if len(closed_rows) else True
    h2_closed_rows = df[df["sample"].isin(h2_closed)]
    h2_zero = bool((h2_closed_rows[[c for c in df if c.startswith("h2eff_")]] == 0).all().all()) if len(h2_closed_rows) else True
    summary = summarise(df, None)
    summary["closed_gates"] = {"mapped_samples": closed["samples"], "co_closed": sorted(co_closed), "h2_closed": sorted(h2_closed),
                               "co_closed_count": len(co_closed), "h2_closed_count": len(h2_closed),
                               "co_effect_exactly_zero_all_schedules_gate_and_proportional": co_zero,
                               "co_uptake_exactly_zero_gate": co_uptake_zero, "h2_effect_exactly_zero_when_h2_closed": h2_zero}
    summary["contrasts"] = contrasts(df, cluster_ci)
    summary["gamma_proportional_quantiles"] = {g: [float(v) for v in df[f"gamma_{g}_proportional"].quantile([0, .25, .5, .75, 1])] for g in ("h2", "co")}
    summary["distinct_gamma_pairs_simulated"] = distinct
    summary["intervene_cross_check_gamma_equal"] = all(cross)
    ref = {s["name"]: gas_uptake(simulate(State(1.8), s["forcing"], base)) for s in schedules}
    summary["reference_uptake_mol_m2"] = {k: {"h2": v[0], "co": v[1]} for k, v in ref.items()}
    (OUT / "paired_regimes_summary.json").write_text(json.dumps(summary, indent=2, allow_nan=False) + "\n")

    tv, mem = treatment_versions_export(schedules, base)
    (OUT / "treatment_versions.json").write_text(json.dumps(tv, indent=2, allow_nan=False) + "\n")
    mem.to_csv(OUT / "treatment_versions.tsv", sep="\t", index=False, float_format="%.12g")

    vr = version_relevance_demo(base)
    (OUT / "version_relevance.json").write_text(json.dumps(vr, indent=2, allow_nan=False) + "\n")

    adm = admissibility()
    adm.to_csv(OUT / "admissibility.tsv", sep="\t", index=False)

    results = {"paired": {k: summary[k] for k in ("samples", "closed_gates", "contrasts", "distinct_gamma_pairs_simulated")},
               "treatment_versions": tv["summary"], "admissibility_counts": adm.groupby(["zone", "verdict"]).size().reset_index().to_dict("records"),
               "runtime_seconds": round(time.time() - t0, 1)}
    (OUT / "results.json").write_text(json.dumps(results, indent=2, allow_nan=False) + "\n")
    print(json.dumps({"samples": len(df), "co_closed": len(co_closed), "h2_closed": len(h2_closed), "co_zero": co_zero,
                      "runtime_s": results["runtime_seconds"]}))


if __name__ == "__main__":
    main()
