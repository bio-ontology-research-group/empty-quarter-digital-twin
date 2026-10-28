"""Generate current evidence, OWL reasoning, backtests, planning products and manuscript artifacts.

Run from the repository root: uv run python reproducibility/v03/reproduce.py
An optional --paper directory receives generated tables, assumptions and figures.
"""
from pathlib import Path
import argparse
import json
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from eq_biodt.v03.contract import contract
from eq_biodt.v03.evidence import sha256
from eq_biodt.v03.forecast import backtest, visits
from eq_biodt.v03.ontology import evaluate_ontology_reasoning
from eq_biodt.v03.planning import compare, experiment, parameter_ensemble
from eq_biodt.v03.process import Parameters, State, gas_exchange

sys.path.insert(0, str(Path(__file__).parent))
from external_moisture import analyse


def write(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")


def tex(text):
    mapping = {"\\": r"\textbackslash{}", "_": r"\_", "%": r"\%", "&": r"\&", "#": r"\#",
               "^": r"\textasciicircum{}", "{": r"\{", "}": r"\}", "$": r"\$"}
    return "".join(mapping.get(c, c) for c in str(text)).replace("->", r"$\rightarrow$")


def tex_statement(text: str) -> str:
    """Render contract clause statements with proper LaTeX math and chemical formulas."""
    replacements = [
        ("matric potential = entry potential times (theta/porosity)^(-b)",
         r"$\psi_m = \psi_e(\theta/\phi)^{-b}$"),
        ("-70 to -100 MPa", r"$-70$ to $-100$~MPa"),
        ("D = D0 (porosity-theta)^(10/3)/porosity^2",
         r"$D = D_0(\phi-\theta)^{10/3}/\phi^2$"),
        ("D C'' - k C = 0", r"$D C'' - k C = 0$"),
        ("C_atm sqrt(D k) tanh(L sqrt(k/D))",
         r"$C_{\mathrm{atm}}\sqrt{D k}\tanh(L\sqrt{k/D})$"),
        ("Q10^((T-Tref)/10)", r"$Q_{10}^{(T-T_{\mathrm{ref}})/10}$"),
        ("0-50 C", r"$0$--$50~^\circ\mathrm{C}$"),
        ("1000 L porosity mm", r"$1000\,L\phi$~mm"),
        ("-2 phi c R T", r"$-2\varphi c R T$"),
        ("H2+0.5O2 -> H2O; CO+0.5O2 -> CO2; CO2+2H2 -> CH2O+H2O",
         r"$\mathrm{H_2}+0.5\mathrm{O_2}\rightarrow\mathrm{H_2O}$; "
         r"$\mathrm{CO}+0.5\mathrm{O_2}\rightarrow\mathrm{CO_2}$; "
         r"$\mathrm{CO_2}+2\mathrm{H_2}\rightarrow\mathrm{CH_2O}+\mathrm{H_2O}$"),
        ("non-CO2", r"non-$\mathrm{CO_2}$"),
        ("Log1p", r"$\log(1+\mathrm{GE})$"),
        ("alpha=10", r"$\alpha=10$"),
        ("[-1,1]", r"$[-1,1]$"),
        ("+/-14 days", r"$\pm 14$~days"),
        ("(-150,-85,-30,-3,-0.3,-0.03 MPa)",
         r"($-150$, $-85$, $-30$, $-3$, $-0.3$, $-0.03$~MPa)"),
        ("(0.3,0.4,0.5)", r"($0.3$, $0.4$, $0.5$)"),
        ("-3 MPa", r"$-3$~MPa"),
        ("25 C", r"$25~^\circ\mathrm{C}$"),
        ("1e-9 mol/m2/s", r"$10^{-9}\text{ mol m}^{-2}\text{ s}^{-1}$"),
        ("b in {3,5}", r"$b\in\{3,5\}$"),
        ("{-70,-100} MPa", r"$\{-70,-100\}$~MPa"),
        ("{0.5,2}", r"$\{0.5,2\}$"),
        ("{0.1,1,10}", r"$\{0.1,1,10\}$"),
        ("CO2", r"$\mathrm{CO_2}$"),
        ("H2", r"$\mathrm{H_2}$"),
        ("O2", r"$\mathrm{O_2}$"),
    ]
    out = str(text)
    tokens = {}
    for idx, (src, dst) in enumerate(replacements):
        if src in out:
            tok = f"@@TOK{chr(65 + idx // 26)}{chr(65 + idx % 26)}@@"
            out = out.replace(src, tok)
            tokens[tok] = dst
    out = tex(out)
    for tok, dst in tokens.items():
        out = out.replace(tok, dst)
    return out


def tex_default(value) -> str:
    if value is None:
        return "Unavailable"
    if value == "REQUIRED":
        return "User supplied"
    if value == "activation_diffusion":
        return "Activation and diffusion"
    if isinstance(value, float):
        if value < 0:
            return f"$-{abs(value):g}$"
        if 0 < value < 1e-3:
            mant, exp = f"{value:.2e}".split("e")
            mant = mant.rstrip("0").rstrip(".")
            return rf"${mant}\times 10^{{{int(exp)}}}$"
        return f"${value:g}$"
    return tex(value)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--paper", type=Path)
    ap.add_argument("--reji-source", type=Path)
    args = ap.parse_args()
    root = Path(__file__).resolve().parents[2]
    out = root / "reviews/v03"
    out.mkdir(exist_ok=True)
    panel = pd.read_csv(root / "data/v03/sample_evidence.tsv", sep="\t")
    manifest = json.loads((root / "data/v03/manifest.json").read_text())
    for name, digest in manifest["artifacts"].items():
        if sha256(root / "data/v03" / name) != digest:
            raise ValueError("Evidence artifact digest mismatch")
    mag = pd.read_csv(root / "data/v03/mag_evidence.tsv", sep="\t", keep_default_na=False)
    ont = evaluate_ontology_reasoning(mag, panel, root / "data/ontology/eq_biodt_v03.ofn")
    write(out / "ontology_reasoning.json", ont)
    bt = backtest(panel)
    write(out / "backtest.json", bt)
    external = analyse(out, args.reji_source)
    p = Parameters(k_h2_s=.001, k_co_s=.001)
    forces = []
    for name, n in [("intact", 0), ("single_pulse", 1), ("four_pulses", 4), ("daily_pulses", 14)]:
        days = {int(i * 14 / n) for i in range(n)}
        forcing = [{"water_mm": 12 / n if day in days else 0, "nacl_mol_l": 0., "evaporation_mm": 1., "temp_c": 25.} for day in range(14)]
        forces.append({"name": name, "forcing": forcing, "disturbance_m2": 0. if n == 0 else 1.})
    settings = {"parameters": {"k_h2_s": .001, "k_co_s": .001}, "initial": {"water_mm": 1.8}, "scenarios": forces}
    write(root / "data/v03/example_comparison.json", settings)
    write(root / "data/v03/example_scenario.json", {k: v for k, v in settings.items() if k != "scenarios"} | {"forcing": forces[2]["forcing"]})
    comparison = compare(State(1.8), forces, parameter_ensemble(p))
    write(out / "comparison.json", comparison)
    design = experiment(p)
    write(out / "experiment.json", design)
    pd.DataFrame(design["schedule"]).to_csv(out / "experiment_schedule.tsv", sep="\t", index=False)
    model = contract()
    write(root / "data/v03/model_contract.json", model)

    pred_df = pd.DataFrame(bt["predictions"])
    ridge_cov = float(pred_df.loc[pred_df.model.eq("geography_season_ridge"), "covered_90"].mean())
    summary = {
        "version": "0.3.0",
        "samples": len(panel),
        "reference_MAGs": len(mag),
        "mapped_samples": int(panel.mapped_read_pct.notna().sum()),
        "forecast_visits": len(visits(panel)),
        "test_visits": sum(r["n_visits"] for r in bt["metrics"] if r["model"] == "persistence"),
        "conformal_90_coverage_ridge": round(ridge_cov, 4),
        "CBB_screens": {
            k: int(mag[k].sum())
            for k in ["rbcl_detected", "rbcl_prk_detected", "three_gene_screen", "cbb_all_steps_detected", "h2_cbb_all_steps", "h2_cbb_quality_90_5"]
        },
        "owl_reasoning": {
            "reasoner": ont["reasoner"],
            "selected_route": ont["selected_route"],
            "logical_axioms": ont["logical_axioms"],
            "classification_wall_seconds": ont["classification_wall_seconds"],
            "justification_axiom_count": ont["example_candidate_justification"]["axiom_count"],
        },
        "external_bottles": external["bottles"],
        "ensemble_members": len(comparison["parameters"]),
        "uptake_rank_reversals": sum(r["rank_reversal"] for r in comparison["pairwise_uptake_comparisons"]),
        "source_code_sha256": {str(f.relative_to(root)): sha256(f) for f in sorted((root / "src/eq_biodt/v03").glob("*.py"))},
    }
    write(out / "results.json", summary)

    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9, "axes.spines.top": False, "axes.spines.right": False, "svg.fonttype": "none"})
    fig, axes = plt.subplots(1, 2, figsize=(9, 3.3), layout="constrained")
    theta = np.linspace(0, .4, 401)
    for b, color in [(3, "#a75d20"), (4, "#007b83"), (5, "#59669f")]:
        pp = Parameters(k_h2_s=.001, k_co_s=.001, retention_b=b)
        axes[0].plot(theta, [gas_exchange(x, 25, pp)["h2_mol_m2_s"] * 1e9 for x in theta], label=f"Retention b = {b}", color=color)
    axes[0].set(xlabel="Volumetric water content", ylabel="Conditional H₂ uptake (nmol m⁻² s⁻¹)", title="A   Dry activation and pore transport")
    axes[0].legend(frameon=False, fontsize=8)
    ef = pd.DataFrame(external["fits"]).groupby(["soil", "model"]).rmse_per_min.mean().unstack()
    x = np.arange(len(ef))
    axes[1].bar(x - .18, ef.diffusion_only, .36, color="#a75d20", label="Diffusion envelope")
    axes[1].bar(x + .18, ef.activation_diffusion, .36, color="#007b83", label="Activation × diffusion envelope")
    soil_labels = {"Forest": "Forest", "Meadow": "Meadow", "PB_sand": "Pine Barrens"}
    axes[1].set(xticks=x, xticklabels=[soil_labels[s] for s in ef.index], ylabel="Held-replicate bottle-rate RMSE (min⁻¹)", title="B   Published bottle observations")
    axes[1].legend(frameon=False, fontsize=8)
    fig.savefig(out / "process_validation.pdf")
    fig.savefig(out / "process_validation.png", dpi=180)
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(9, 3.3), layout="constrained")
    scores = pd.DataFrame(bt["metrics"])
    for forecast_model, color, label in [("geography_season_ridge", "#007b83", "Geography + season"), ("persistence", "#a75d20", "Persistence"), ("compartment_mean", "#59669f", "Compartment mean")]:
        g = scores.loc[scores.model.eq(forecast_model)]
        axes[0].plot(g.campaign, g.rmse_log1p_ge, "o-", color=color, label=label)
    axes[0].set(xlabel="Held-out campaign", ylabel="RMSE of log(1 + marker GE)", xticks=[2, 3, 4, 5], title="A   Future-campaign prediction")
    axes[0].legend(frameon=False, fontsize=8)
    outcomes = pd.DataFrame(comparison["outcomes"])
    for i, name in enumerate([s["name"] for s in forces]):
        v = outcomes.loc[outcomes.scenario.eq(name), "h2_uptake_mol_m2"] * 1e6
        axes[1].plot([v.min(), v.max()], [i, i], lw=5, color="#007b83")
    axes[1].set(yticks=range(4), yticklabels=["Intact control", "One 12 mm pulse", "Four 3 mm pulses", "Daily 12/14 mm pulses"], xlabel="14-day H₂ uptake range (µmol m⁻²)", title="B   Explicit parameter sensitivity")
    fig.savefig(out / "forecast_scenarios.pdf")
    fig.savefig(out / "forecast_scenarios.png", dpi=180)
    plt.close(fig)

    if args.paper:
        dest = args.paper
        (dest / "figures").mkdir(exist_ok=True)
        (dest / "tables").mkdir(exist_ok=True)
        for filename in ["process_validation.pdf", "forecast_scenarios.pdf"]:
            (dest / "figures" / filename).write_bytes((out / filename).read_bytes())
        body = [
            "% Generated from the executable v0.3 contract; do not edit by hand.",
            r"\section{Complete axioms, observation rules and assumptions}",
            "The following stable identifiers enumerate the implemented model contract. Ontological statements specify quantity types and constitutive relations. Causal statements specify conditional mechanisms; predictive and decision assumptions have separate roles.",
        ]
        body += [r"\subsection*{" + tex(c["id"] + ": " + c["kind"]) + "}\n" + tex_statement(c["statement"]) for c in model["clauses"]]
        parameter_labels = {
            "k_h2_s": r"$\mathrm{H_2}$ reference reaction coefficient ($\mathrm{s}^{-1}$)",
            "k_co_s": r"$\mathrm{CO}$ reference reaction coefficient ($\mathrm{s}^{-1}$)",
            "porosity": r"Porosity $\phi$ (fraction)",
            "depth_m": r"Soil layer depth $L$ (m)",
            "psi_entry_mpa": r"Entry water potential $\psi_e$ (MPa)",
            "retention_b": r"Retention exponent $b$",
            "dry_threshold_mpa": r"Dry activation threshold $\psi_d$ (MPa)",
            "wet_activation_mpa": r"Full-activation endpoint $\psi_w$ (MPa)",
            "activation_shape": r"Activation shape exponent $a$",
            "q10": r"Temperature response $Q_{10}$",
            "reference_temp_c": r"Reference temperature $T_{\mathrm{ref}}$ ($^\circ\mathrm{C}$)",
            "d_h2_m2_s": r"$\mathrm{H_2}$ free-air diffusivity $D_{0,\mathrm{H_2}}$ ($\mathrm{m}^2\,\mathrm{s}^{-1}$)",
            "d_co_m2_s": r"$\mathrm{CO}$ free-air diffusivity $D_{0,\mathrm{CO}}$ ($\mathrm{m}^2\,\mathrm{s}^{-1}$)",
            "h2_mol_m3": r"Atmospheric $\mathrm{H_2}$ concentration $C_{a,\mathrm{H_2}}$ ($\mathrm{mol}\,\mathrm{m}^{-3}$)",
            "co_mol_m3": r"Atmospheric $\mathrm{CO}$ concentration $C_{a,\mathrm{CO}}$ ($\mathrm{mol}\,\mathrm{m}^{-3}$)",
            "oxygen_supply_mol_m2_s": r"Oxygen supply limit ($\mathrm{mol}\,\mathrm{m}^{-2}\,\mathrm{s}^{-1}$)",
            "energy_capture_fraction": r"Captured energy fraction $\eta$",
            "maintenance_w_m2": r"Independent maintenance demand $M$ ($\mathrm{W}\,\mathrm{m}^{-2}$)",
            "field_capacity_fraction": r"Field-capacity fraction $f$ of pore storage",
            "drainage_fraction_day": r"Daily drainage fraction $d$",
            "osmotic_coefficient": r"Osmotic coefficient $\varphi$",
            "osmotic_domain_mol_l": r"Dilute NaCl domain limit ($\mathrm{mol}\,\mathrm{L}^{-1}$)",
            "mechanism": "Physiological mechanism",
        }
        assert parameter_labels.keys() == model["parameter_defaults"].keys()
        body += [r"\section{Process parameters}", r"\begin{longtable}{p{0.63\textwidth}p{0.27\textwidth}}\toprule Parameter & Default \\\midrule\endhead"]
        body += [parameter_labels[k] + " & " + tex_default(v) + r" \\" for k, v in model["parameter_defaults"].items()]
        body += [
            r"\bottomrule\end{longtable}",
            "Required kinetic coefficients have no fitted default. The example workflows supply $0.001~\\mathrm{s}^{-1}$ for each gas. All remaining defaults are scenario assumptions or physical constants with the scope stated in C1--C12.",
            r"\section{Carbon-fixation annotation rules}",
            r"\begin{longtable}{p{0.43\textwidth}p{0.49\textwidth}}\toprule Reaction step & Alternative KOs (logical OR) \\\midrule\endhead",
        ]
        body += [tex(k.replace("_", " ").capitalize()) + " & " + tex(", ".join(v)) + r" \\" for k, v in model["CBB_step_alternatives"].items()]
        body += [r"\bottomrule\end{longtable}"]
        (dest / "assumptions.tex").write_text("\n\n".join(body) + "\n")
        rows = []
        for campaign, g in scores.groupby("campaign"):
            vals = g.set_index("model").rmse_log1p_ge
            cov = g.set_index("model").loc["geography_season_ridge", "coverage_90"] * 100
            rows.append(
                f"{int(campaign)} & {int(g.iloc[0].n_visits)} & {vals['geography_season_ridge']:.3f} & {vals['persistence']:.3f} & {vals['compartment_mean']:.3f} & {vals['global_mean']:.3f} & {cov:.1f}\\%" + r" \\"
            )
        (dest / "tables/v03_backtest.tex").write_text("\n".join(rows) + "\n")
        macros = {
            "ReferenceMAGs": len(mag),
            "CBBFull": summary["CBB_screens"]["cbb_all_steps_detected"],
            "HtwoCBB": summary["CBB_screens"]["h2_cbb_all_steps"],
            "HtwoCBBStrict": summary["CBB_screens"]["h2_cbb_quality_90_5"],
            "ArchaealRuBisCO": ont["entailed_class_counts"]["ArchaealRuBisCOMAG"],
            "DualHtwoCO": ont["entailed_class_counts"]["DualH2COAutotrophicCandidateMAG"],
            "HtwoCarriers": ont["entailed_class_counts"]["H2CarrierMAG"],
            "MultiHtwoMAGs": ont["ablations"]["additive_hydrogenase_guilds_without_carrier_union"]["multi_hydrogenase_mags_double_counted"],
            "JustAxioms": ont["example_candidate_justification"]["axiom_count"],
            "AmbiguousSamples": ont["ablations"]["sample_marginals_without_genome_linkage"]["nonzero_ambiguity_samples"],
            "ForecastVisits": summary["forecast_visits"],
            "TestVisits": summary["test_visits"],
            "ConformalCoverage": f"{ridge_cov * 100:.1f}\\%",
            "RankReversals": summary["uptake_rank_reversals"],
        }
        (dest / "tables/v03_numbers.tex").write_text("\n".join("\\newcommand{\\" + k + "}{" + str(v) + "}" for k, v in macros.items()) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
