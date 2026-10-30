"""LaTeX macros and the supplement admissibility table from reviews/v06 (no new computation)."""
from pathlib import Path
import json

import sys

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "reviews/v06"
TAB = ROOT / "paper/tables"
sys.path.insert(0, str(ROOT / "src"))
UM = 1e6  # mol m-2 -> umol m-2
COMP = {"Surface": "Surface", "Deep": "Deep", "PlantRoot": "Root"}
SCHED = {"single_pulse": "Single", "four_pulses": "Four", "daily_pulses": "Daily"}
SYM = {"K": "K", "E": "E", "D_K": "D_K(E)", "q_h2": "q_{\\mathrm{H}_2}", "q_co": "q_{\\mathrm{CO}}",
       "gamma_h2": "\\gamma_{\\mathrm{H}_2}", "gamma_co": "\\gamma_{\\mathrm{CO}}", "I": "I_t", "c": "c_t", "Evap": "\\mathrm{Ev}_t",
       "T": "T_t", "phi": "\\phi", "psi_e": "\\psi_e", "b": "b", "L": "L", "C_a": "C_a", "D0": "D_0", "k_ref": "k_{\\mathrm{ref}}",
       "Q10": "Q_{10}", "act_par": "\\psi_d,\\psi_w,a", "drain_par": "f,d", "osm_par": "\\varphi,R", "O2": "S_{\\mathrm{O}_2}",
       "capture": "\\eta", "W": "W_t", "N": "N_t", "W_max": "W_{\\max}", "theta": "\\theta", "psi_m": "\\psi_m", "psi_o": "\\psi_o",
       "psi": "\\psi", "D_h2": "D_{\\mathrm{H}_2}", "D_co": "D_{\\mathrm{CO}}", "A": "A_t", "k_h2": "k_{\\mathrm{H}_2}",
       "k_co": "k_{\\mathrm{CO}}", "J_h2": "J_{\\mathrm{H}_2}", "J_co": "J_{\\mathrm{CO}}", "P": "P", "P_cap": "P_{\\mathrm{cap}}"}
ZONE = {"ontological": "Ontology $K$", "evidence": "Evidence $E$", "deduced": "Deduced", "exogenous_manipulable": "Input",
        "exogenous_fixed": "Constant", "constitutive": "Constitutive", "causal": "Causal"}
NOTE = {"ontological": "model revision", "evidence": "observation update", "deduced": "deduced from $K$, $E$; observation update",
        "exogenous_fixed": "fixed in this release", "causal": "outcome mechanism; act on its inputs"}


def f(x, d=1):
    return f"{x:.{d}f}"


def sci(x, d=2):
    m, e = f"{x:.{d}e}".split("e")
    return f"{m}\\times10^{{{int(e)}}}"


def main():
    summ = json.loads((OUT / "paired_regimes_summary.json").read_text())
    tv = json.loads((OUT / "treatment_versions.json").read_text())["summary"]
    vr = json.loads((OUT / "version_relevance.json").read_text())
    adm = pd.read_csv(OUT / "admissibility.tsv", sep="\t", keep_default_na=False)
    m = {}
    cg = summ["closed_gates"]
    m["MappedSamples"] = summ["samples"]
    m["ClosedCOGates"] = cg["co_closed_count"]
    m["ClosedHtwoGates"] = cg["h2_closed_count"]
    m["ClosedCOZeroEffects"] = cg["co_closed_count"] * 3 * 2          # samples x schedules x couplings, verified exactly 0
    m["ClosedCOAllZero"] = "yes" if cg["co_effect_exactly_zero_all_schedules_gate_and_proportional"] else "NO"
    for key, v in summ["by_compartment"].items():
        gas, mode, sched, comp = key.split("|")
        if sched not in SCHED:
            continue
        tag = {"h2": "Htwo", "co": "CO"}[gas] + {"gate": "Gate", "proportional": "Prop"}[mode] + SCHED[sched]
        for stat, name in (("median", "Median"), ("q1", "QOne"), ("q3", "QThree")):
            m[f"{tag}{name}{COMP[comp]}"] = f(v[stat] * UM, 1)
        m[f"{tag}N{COMP[comp]}"] = v["n"]
    ref = summ["reference_uptake_mol_m2"]
    for s, name in SCHED.items():
        m[f"RefHtwo{name}"] = f(ref[s]["h2"] * UM, 1)
        m[f"RefCO{name}"] = f(ref[s]["co"] * UM, 1)
    names = {"Deep minus PlantRoot": "DeepRoot", "Deep minus Surface": "DeepSurface", "PlantRoot minus Surface": "RootSurface"}
    for c in summ["contrasts"]:
        if c["coupling"] != "proportional":
            continue
        n = names[c["contrast"]]
        m[f"Contrast{n}Est"] = f(c["mean_difference_mol_m2"] * UM, 0)
        m[f"Contrast{n}Low"] = f(c["site_bootstrap_95_percentile"][0] * UM, 0)
        m[f"Contrast{n}High"] = f(c["site_bootstrap_95_percentile"][1] * UM, 0)
        m[f"Contrast{n}Visits"] = c["matched_site_campaigns"]
        m[f"Contrast{n}Sites"] = c["sites"]
    for c in summ["contrasts"]:
        if c["coupling"] == "gate":
            m[f"GateContrast{names[c['contrast']]}Est"] = f(c["mean_difference_mol_m2"] * UM, 0)
    m["GammaPropMinCO"] = f(summ["gamma_proportional_quantiles"]["co"][0], 2)
    m["GammaPropMaxCO"] = f(summ["gamma_proportional_quantiles"]["co"][4], 1)
    m["GammaPropMinHtwo"] = f(summ["gamma_proportional_quantiles"]["h2"][0], 3)
    m["GammaPropMaxHtwo"] = f(summ["gamma_proportional_quantiles"]["h2"][4], 2)
    # treatment versions
    m["VersionMembers"] = tv["members"]
    m["VersionSingleFirst"] = tv["members_where_single_beats_four"]
    m["VersionFourFirst"] = tv["members_where_four_beats_single"]
    m["VersionDailyLowest"] = tv["members_where_daily_is_lowest"]
    m["VersionRankings"] = len(tv["distinct_rankings"])
    lo, hi = tv["effect_vs_control_range_mol_m2"]
    m["VersionEffectLow"], m["VersionEffectHigh"] = f(lo * UM, 1), f(hi * UM, 0)
    for v, tag in (("single_pulse", "Single"), ("four_pulses", "Four"), ("daily_pulses", "Daily")):
        a, b = tv["effect_range_by_version_mol_m2"][v]
        m[f"VersionEffect{tag}Low"], m[f"VersionEffect{tag}High"] = f(a * UM, 1), f(b * UM, 0)
    # version relevance demo
    for t in vr["targets"]:
        tag = "Strong" if t["psi_mpa"] == -1.0 else "Mild"
        for ver, vtag in (("matric", "Matric"), ("osmotic", "Osmotic")):
            d = t["versions"][ver]
            m[f"Psi{tag}{vtag}A"] = f(d["A"], 3)
            m[f"Psi{tag}{vtag}DHtwo"] = f(d["D_h2_m2_s"] * 1e6, 2)
            m[f"Psi{tag}{vtag}JHtwo"] = f(d["J_h2_mol_m2_s"] * 1e9, 3)
            m[f"Psi{tag}{vtag}Water"] = f(d["water_mm"], 1)
            m[f"Psi{tag}{vtag}Salt"] = f(d["nacl_mol_l"], 3)
        m[f"Psi{tag}JRatio"] = f(t["J_h2_ratio_matric_over_osmotic"], 3)
        m[f"Psi{tag}AEqual"] = "yes" if t["A_equal"] else "no"
        d = t["versions"]
        if all(d[v]["h2_uptake_3d_mol_m2"] is not None for v in d):
            m[f"Psi{tag}UptakeMatric"] = f(d["matric"]["h2_uptake_3d_mol_m2"] * UM, 1)
            m[f"Psi{tag}UptakeOsmotic"] = f(d["osmotic"]["h2_uptake_3d_mol_m2"] * UM, 1)
    m["PsiStateWater"] = f(vr["state_entering_day_1"]["water_mm"], 1)
    m["PsiRelevantA"] = "yes" if vr["relevance"]["A"]["relevant"] else "no"
    m["PsiRelevantJHtwo"] = "yes" if vr["relevance"]["J_h2"]["relevant"] else "no"
    m["PsiRelevantDHtwo"] = "yes" if vr["relevance"]["D_h2"]["relevant"] else "no"
    # admissibility
    m["GraphNodes"] = len(adm)
    m["AdmissibleNodes"] = int((adm.verdict == "admissible").sum())
    m["RefusedNodes"] = int((adm.verdict == "refused").sum())
    m["AdmissibleInputs"] = int(((adm.verdict == "admissible") & adm.zone.eq("exogenous_manipulable")).sum())
    m["AdmissibleStates"] = int(((adm.verdict == "admissible") & adm.zone.eq("causal")).sum())
    TAB.mkdir(exist_ok=True)
    (TAB / "v06_numbers.tex").write_text("% Generated by reproducibility/v06/tables.py; do not edit by hand.\n" +
                                         "\n".join(f"\\newcommand{{\\{k}}}{{{v}}}" for k, v in m.items()) + "\n")
    rows = []
    handles = {}
    from eq_biodt.v03.theory import ZONES
    adm = adm.assign(_z=adm.zone.map(ZONES.index)).sort_values("_z", kind="stable")
    for r in adm.itertuples():
        sym = f"${SYM[r.node]}$"
        if r.verdict == "admissible":
            note = r.physical_action
        elif r.zone == "constitutive":
            note = "defined by its arguments; " + ("translated to $W$ or $N$" if r.node in ("theta", "psi") else "act on its constituents")
        else:
            note = NOTE[r.zone]
        rows.append(f"{sym} & {ZONE[r.zone]} & {'yes' if r.manipulable else 'no'} & {'admissible' if r.verdict == 'admissible' else 'refused'} & {note} \\\\")
    (TAB / "v06_admissibility.tex").write_text("\n".join(rows) + "\n")
    print(f"{len(m)} macros, {len(rows)} admissibility rows")


if __name__ == "__main__":
    main()
