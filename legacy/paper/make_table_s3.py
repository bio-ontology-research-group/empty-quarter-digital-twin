"""Write the objective-score rows of Supplementary Table S3 (six portfolios at Site 30).

Run from the repository root: `uv run python paper/make_table_s3.py`.
Writes `paper/tables/table_s3_scores.tex`, which `supplementary.tex` inputs.
"""

from pathlib import Path

from eq_biodt.api.main import get_runtime

# Portfolio labels used in Supplementary Table S3, keyed by the optimizer's policy names.
PORTFOLIO_LABELS = {
    "Baseline Hyperarid Do-Nothing": ("P1", "Do-nothing baseline"),
    "Integrated Eco-Restoration & Oryx Rewilding (Balanced Optimum)": (
        "P2",
        "Integrated eco-restoration and oryx rewilding",
    ),
    "Low-Water Chemosynthetic Crust & Refugium Conservation": (
        "P3",
        "Low-water chemosynthetic crust and refugium conservation",
    ),
    "High-Canopy Oasis Revegetation + Site 60 Ca-SO4 Irrigation": (
        "P4",
        "High-canopy oasis revegetation with Ca-SO$_4$ irrigation",
    ),
    "Active Energy-Corridor Nodal Seismic + Oil-Spill Bioremediation": (
        "P5",
        "Energy corridor with nodal seismic and oil-spill bioremediation",
    ),
    "Unmitigated Heavy Seeding + Brine Disposal + Camel Overgrazing + Vibroseis": (
        "P6",
        "Unmitigated seeding, brine, camel overgrazing and Vibroseis",
    ),
}


def main():
    rt = get_runtime()
    row30 = rt.resolve_site_compartment_row("30", "Surface")
    opt = rt.pareto_opt.optimize_pareto_front(row30)

    # The optimizer returns portfolios ranked by Pareto membership, then utility.
    ranked = opt["all_evaluated_portfolios"]
    rank_of = {ev["policy"]["name"]: i + 1 for i, ev in enumerate(ranked)}

    rows = []
    for ev in sorted(ranked, key=lambda e: PORTFOLIO_LABELS[e["policy"]["name"]][0]):
        name = ev["policy"]["name"]
        code, label = PORTFOLIO_LABELS[name]
        o = ev["objectives"]
        rob = ev["probabilistic_robustness"]
        rows.append(
            f"{code} & {label} & {'yes' if ev['is_pareto_optimal'] else 'no'} & "
            f"{o['J1_biodiversity_refugium_max']:.4f} & "
            f"{o['J2_chemosynthesis_cn_flux_max']:.4f} & "
            f"{o['J3_salinity_compaction_risk_min']:.4f} & "
            f"{o['J4_water_energy_arg_risk_min']:.4f} & "
            f"{ev['desdeo_utility_score']:.4f} & "
            f"{rob['pareto_dominance_prob']:.3f} & "
            f"{rob['rank1_prob']:.3f} & "
            f"{rank_of[name]} & {ev['ocs_score']:.1f} \\\\"
        )

    out = Path(__file__).resolve().parent / "tables" / "table_s3_scores.tex"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text("\n".join(rows) + "\n", encoding="utf-8")
    print(f"wrote {out}")
    print("\n".join(rows))


if __name__ == "__main__":
    main()
