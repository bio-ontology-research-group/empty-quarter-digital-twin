"""Generate current evidence, backtests, planning products and manuscript artifacts.

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
from eq_biodt.v03.planning import compare, experiment, parameter_ensemble
from eq_biodt.v03.process import Parameters, State, gas_exchange

sys.path.insert(0, str(Path(__file__).parent))
from external_moisture import analyse


def write(path, value):
    path.write_text(json.dumps(value,indent=2,allow_nan=False)+"\n")


def tex(text):
    mapping = {"\\":r"\textbackslash{}", "_":r"\_", "%":r"\%", "&":r"\&", "#":r"\#",
               "^":r"\textasciicircum{}", "{":r"\{", "}":r"\}", "$":r"\$"}
    return "".join(mapping.get(c,c) for c in str(text)).replace("->", r"$\rightarrow$")


def main():
    ap=argparse.ArgumentParser();ap.add_argument("--paper",type=Path);ap.add_argument("--reji-source",type=Path)
    args=ap.parse_args()
    root=Path(__file__).resolve().parents[2];out=root/"reviews/v03";out.mkdir(exist_ok=True)
    panel=pd.read_csv(root/"data/v03/sample_evidence.tsv",sep="\t")
    manifest=json.loads((root/"data/v03/manifest.json").read_text())
    for name,digest in manifest["artifacts"].items():
        if sha256(root/"data/v03"/name) != digest:
            raise ValueError("Evidence artifact digest mismatch")
    bt=backtest(panel);write(out/"backtest.json",bt)
    external=analyse(out,args.reji_source)
    p=Parameters(k_h2_s=.001,k_co_s=.001)
    forces=[]
    for name,n in [("intact",0),("single_pulse",1),("four_pulses",4),("daily_pulses",14)]:
        days={int(i*14/n) for i in range(n)}
        forcing=[{"water_mm":12/n if day in days else 0,"nacl_mol_l":0.,"evaporation_mm":1.,"temp_c":25.} for day in range(14)]
        forces.append({"name":name,"forcing":forcing,"disturbance_m2":0. if n==0 else 1.})
    settings={"parameters":{"k_h2_s":.001,"k_co_s":.001},"initial":{"water_mm":1.8},"scenarios":forces}
    write(root/"data/v03/example_comparison.json",settings)
    write(root/"data/v03/example_scenario.json",{k:v for k,v in settings.items() if k!="scenarios"}|{"forcing":forces[2]["forcing"]})
    comparison=compare(State(1.8),forces,parameter_ensemble(p));write(out/"comparison.json",comparison)
    design=experiment(p);write(out/"experiment.json",design)
    pd.DataFrame(design["schedule"]).to_csv(out/"experiment_schedule.tsv",sep="\t",index=False)
    model=contract();write(root/"data/v03/model_contract.json",model)
    mag=pd.read_csv(root/"data/v03/mag_evidence.tsv",sep="\t")
    summary={"version":"0.3.0","samples":len(panel),"reference_MAGs":len(mag),
             "mapped_samples":int(panel.mapped_read_pct.notna().sum()),
             "forecast_visits":len(visits(panel)),
             "test_visits":sum(r["n_visits"] for r in bt["metrics"] if r["model"]=="persistence"),
             "CBB_screens":{k:int(mag[k].sum()) for k in ["rbcl_detected","rbcl_prk_detected","three_gene_screen","cbb_all_steps_detected","h2_cbb_all_steps","h2_cbb_quality_90_5"]},
             "external_bottles":external["bottles"],"ensemble_members":len(comparison["parameters"]),
             "uptake_rank_reversals":sum(r["rank_reversal"] for r in comparison["pairwise_uptake_comparisons"]),
             "unissued_forecast_reason":"No next sampling date has been scheduled.",
             "source_code_sha256":{str(f.relative_to(root)):sha256(f) for f in sorted((root/"src/eq_biodt/v03").glob("*.py"))}}
    write(out/"results.json",summary)
    plt.rcParams.update({"font.family":"DejaVu Sans","font.size":9,"axes.spines.top":False,"axes.spines.right":False,"svg.fonttype":"none"})
    fig,axes=plt.subplots(1,2,figsize=(9,3.3),layout="constrained")
    theta=np.linspace(0,.4,401)
    for b,color in [(3,"#a75d20"),(4,"#007b83"),(5,"#59669f")]:
        pp=Parameters(k_h2_s=.001,k_co_s=.001,retention_b=b)
        axes[0].plot(theta,[gas_exchange(x,25,pp)["h2_mol_m2_s"]*1e9 for x in theta],label=f"Retention b = {b}",color=color)
    axes[0].set(xlabel="Volumetric water content",ylabel="Conditional H₂ uptake (nmol m⁻² s⁻¹)",title="A   Dry activation and pore transport")
    axes[0].legend(frameon=False,fontsize=8)
    ef=pd.DataFrame(external["fits"]).groupby(["soil","model"]).rmse_per_min.mean().unstack()
    x=np.arange(len(ef))
    axes[1].bar(x-.18,ef.diffusion_only,.36,color="#a75d20",label="Diffusion envelope")
    axes[1].bar(x+.18,ef.activation_diffusion,.36,color="#007b83",label="Activation × diffusion envelope")
    axes[1].set(xticks=x,xticklabels=["Forest","Meadow","Desert sand"],ylabel="Held-replicate bottle-rate RMSE (min⁻¹)",title="B   Published bottle observations")
    axes[1].legend(frameon=False,fontsize=8)
    fig.savefig(out/"process_validation.pdf");fig.savefig(out/"process_validation.png",dpi=180);plt.close(fig)
    fig,axes=plt.subplots(1,2,figsize=(9,3.3),layout="constrained")
    scores=pd.DataFrame(bt["metrics"])
    for forecast_model,color,label in [("geography_season_ridge","#007b83","Geography + season"),("persistence","#a75d20","Persistence"),("compartment_mean","#59669f","Compartment mean")]:
        g=scores.loc[scores.model.eq(forecast_model)]
        axes[0].plot(g.campaign,g.rmse_log1p_ge,"o-",color=color,label=label)
    axes[0].set(xlabel="Held-out campaign",ylabel="RMSE of log(1 + marker GE)",xticks=[2,3,4,5],title="A   Future-campaign prediction")
    axes[0].legend(frameon=False,fontsize=8)
    outcomes=pd.DataFrame(comparison["outcomes"])
    for i,name in enumerate([s["name"] for s in forces]):
        v=outcomes.loc[outcomes.scenario.eq(name),"h2_uptake_mol_m2"]*1e6
        axes[1].plot([v.min(),v.max()],[i,i],lw=5,color="#007b83")
    axes[1].set(yticks=range(4),yticklabels=["Intact control","One 12 mm pulse","Four 3 mm pulses","Daily 12/14 mm pulses"],xlabel="14-day H₂ uptake range (µmol m⁻²)",title="B   Explicit parameter sensitivity")
    fig.savefig(out/"forecast_scenarios.pdf");fig.savefig(out/"forecast_scenarios.png",dpi=180);plt.close(fig)
    if args.paper:
        dest=args.paper; (dest/"figures").mkdir(exist_ok=True);(dest/"tables").mkdir(exist_ok=True)
        for filename in ["process_validation.pdf","forecast_scenarios.pdf"]:
            (dest/"figures"/filename).write_bytes((out/filename).read_bytes())
        body=["% Generated from the executable v0.3 contract; do not edit by hand.",r"\section{Complete axioms, observation rules and assumptions}",
              "The following stable identifiers enumerate the implemented model contract. Ontological statements specify quantity types and constitutive relations. Causal statements specify conditional mechanisms; predictive and decision assumptions have separate roles."]
        body += [r"\subsection*{"+tex(c["id"]+": "+c["kind"])+"}\n"+tex(c["statement"]) for c in model["clauses"]]
        parameter_labels = {
            "k_h2_s":"H2 reference reaction coefficient (1/s)", "k_co_s":"CO reference reaction coefficient (1/s)",
            "porosity":"Porosity (fraction)", "depth_m":"Soil layer depth (m)",
            "psi_entry_mpa":"Entry water potential (MPa)", "retention_b":"Retention exponent b",
            "dry_threshold_mpa":"Dry activation threshold (MPa)", "wet_activation_mpa":"Full-activation endpoint (MPa)",
            "activation_shape":"Activation shape exponent", "q10":"Temperature response Q10",
            "reference_temp_c":"Reference temperature (degrees C)", "d_h2_m2_s":"H2 free-air diffusivity (m2/s)",
            "d_co_m2_s":"CO free-air diffusivity (m2/s)", "h2_mol_m3":"Atmospheric H2 concentration (mol/m3)",
            "co_mol_m3":"Atmospheric CO concentration (mol/m3)", "oxygen_supply_mol_m2_s":"Oxygen supply limit (mol/m2/s)",
            "energy_capture_fraction":"Captured energy fraction", "maintenance_w_m2":"Independent maintenance demand (W/m2)",
            "field_capacity_fraction":"Field-capacity fraction of pore storage", "drainage_fraction_day":"Daily drainage fraction",
            "osmotic_coefficient":"Osmotic coefficient", "osmotic_domain_mol_l":"Dilute NaCl domain limit (mol/L)",
            "mechanism":"Physiological mechanism"}
        assert parameter_labels.keys() == model["parameter_defaults"].keys()
        body += [r"\section{Process parameters}",r"\begin{longtable}{p{0.63\textwidth}p{0.27\textwidth}}\toprule Parameter & Default \\\midrule\endhead"]
        def default_label(value):
            if value is None: return "Unavailable"
            if value == "REQUIRED": return "User supplied"
            if value == "activation_diffusion": return "Activation and diffusion"
            return str(value)
        body += [tex(parameter_labels[k])+" & "+tex(default_label(v))+r" \\" for k,v in model["parameter_defaults"].items()]
        body += [r"\bottomrule\end{longtable}","Required kinetic coefficients have no fitted default. The example workflows supply 0.001 per second for each gas. All remaining defaults are scenario assumptions or physical constants with the scope stated in C1--C12.",r"\section{Carbon-fixation annotation rules}",r"\begin{longtable}{p{0.43\textwidth}p{0.49\textwidth}}\toprule Reaction step & Alternative KOs (logical OR) \\\midrule\endhead"]
        body += [tex(k.replace("_", " ").capitalize())+" & "+tex(", ".join(v))+r" \\" for k,v in model["CBB_step_alternatives"].items()]
        body += [r"\bottomrule\end{longtable}"]
        (dest/"assumptions.tex").write_text("\n\n".join(body)+"\n")
        rows=[]
        for campaign,g in scores.groupby("campaign"):
            vals=g.set_index("model").rmse_log1p_ge
            rows.append(f"{int(campaign)} & {int(g.iloc[0].n_visits)} & {vals['geography_season_ridge']:.3f} & {vals['persistence']:.3f} & {vals['compartment_mean']:.3f} & {vals['global_mean']:.3f}"+r" \\")
        (dest/"tables/v03_backtest.tex").write_text("\n".join(rows)+"\n")
        macros={"ReferenceMAGs":len(mag),"CBBFull":summary["CBB_screens"]["cbb_all_steps_detected"],"HtwoCBB":summary["CBB_screens"]["h2_cbb_all_steps"],"HtwoCBBStrict":summary["CBB_screens"]["h2_cbb_quality_90_5"],"ForecastVisits":summary["forecast_visits"],"TestVisits":summary["test_visits"],"RankReversals":summary["uptake_rank_reversals"]}
        (dest/"tables/v03_numbers.tex").write_text("\n".join("\\newcommand{\\"+k+"}{"+str(v)+"}" for k,v in macros.items())+"\n")
    print(json.dumps(summary,indent=2))


if __name__=="__main__":
    main()
