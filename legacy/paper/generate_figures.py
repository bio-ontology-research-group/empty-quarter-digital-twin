"""Generate publication-quality multi-panel vector PDF figures for EQ-BioDT."""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from eq_biodt.interventions.cloud_seeding import CloudSeedingIntervention
from eq_biodt.interventions.desalination_brine import DesalinationBrineIntervention
from eq_biodt.interventions.pareto_optimizer import ParetoPortfolioOptimizer
from eq_biodt.interventions.revegetation import RevegetationIntervention

mpl.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.size": 8.5,
        "axes.labelsize": 9.0,
        "axes.titlesize": 9.5,
        "legend.fontsize": 7.5,
        "xtick.labelsize": 8.0,
        "ytick.labelsize": 8.0,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    }
)


def generate_fig1(root: Path, out_dir: Path) -> None:
    sites = pd.read_parquet(root / "data/processed/eq_sites_cube.parquet")
    comp = pd.read_parquet(root / "data/processed/eq_compartment_state_cube.parquet")
    grid = pd.read_parquet(root / "data/processed/eq_era5_h3_corridor_grid.parquet")

    fig, axes = plt.subplots(2, 2, figsize=(12.2, 8.2), constrained_layout=True)

    # Panel A: 1,043 km Transect spatial map + empirical domain-of-support grid mask
    ax = axes[0, 0]
    supp_grid = grid[grid["in_empirical_support"]]
    abst_grid = grid[~grid["in_empirical_support"]]
    ax.scatter(
        abst_grid["lon"],
        abst_grid["lat"],
        c="#d9d9d9",
        marker="s",
        s=14,
        alpha=0.55,
        label=f"Abstained grid (n={len(abst_grid)})",
    )
    ax.scatter(
        supp_grid["lon"],
        supp_grid["lat"],
        c="#bdd7e7",
        marker="s",
        s=14,
        alpha=0.65,
        label=f"Supported grid (n={len(supp_grid)})",
    )
    evap_by_site = comp.groupby("site_num")["xrf_evaporite_index"].mean()
    evap = sites["site_num"].map(evap_by_site).fillna(0.0)
    primary = sites[sites["is_primary_transect"]]
    aux = sites[~sites["is_primary_transect"]]
    sc = ax.scatter(
        primary["lon"],
        primary["lat"],
        c=evap.loc[primary.index],
        cmap="viridis",
        s=42,
        edgecolors="black",
        linewidths=0.45,
        zorder=4,
        label="Primary transect (n=60)",
    )
    ax.scatter(
        aux["lon"],
        aux["lat"],
        c=evap.loc[aux.index],
        cmap="viridis",
        marker="^",
        s=52,
        edgecolors="crimson",
        linewidths=0.75,
        zorder=5,
        label="Auxiliary/Sabkha (n=10)",
    )
    cbar = fig.colorbar(sc, ax=ax, shrink=0.84, pad=0.02)
    cbar.set_label("XRF Evaporite Index (CLR Ca+Mg+Na+S+Cl vs Si)")
    ax.set_xlabel("Longitude (°E)")
    ax.set_ylabel("Latitude (°N)")
    ax.set_title("a  Rub' al-Khali 1,043 km Transect & Support Mask")
    ax.legend(loc="lower right", frameon=True, ncol=2)
    ax.grid(alpha=0.25, linestyle="--")

    # Panel B: Compartment distribution of Group 1l [NiFe]-hydrogenase (hydB_1l GE at 60% floor)
    ax = axes[0, 1]
    order = ["Surface", "ShallowSubsurface", "RootAdjacent"]
    labels = [
        "Surface\n(0-2 cm, n=45)",
        "Subsurface\n(5-10 cm, n=51)",
        "Root-adjacent\n(Rhizosphere, n=59)",
    ]
    colors = ["#d95f02", "#1b9e77", "#7570b3"]
    obs_comp = comp[comp["wgs_observed"]].copy()
    data_hyd = [
        np.log2(obs_comp.loc[obs_comp["compartment_norm"] == c, "hydB_1l_ge_raw"].dropna().values + 0.01)
        for c in order
    ]
    bp = ax.boxplot(
        data_hyd,
        tick_labels=labels,
        patch_artist=True,
        widths=0.50,
        medianprops={"color": "black", "linewidth": 1.5},
    )
    for patch, color in zip(bp["boxes"], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.74)
    ax.set_ylabel(r"$\log_2(\mathrm{hydB\_1l}\ \mathrm{GE} + 0.01)$ (60% floor)")
    ax.set_title("b  Compartment Trace-Gas Capacity (155 Observed States)")
    ax.grid(axis="y", alpha=0.25, linestyle="--")

    # Panel C: 16S Shannon diversity vs Group 1l hydrogenase across compartments
    ax = axes[1, 0]
    for c, lbl, col, mk in zip(
        order,
        ["Surface (n=45)", "Subsurface 5-10 cm (n=51)", "Root-adjacent (n=59)"],
        colors,
        ["o", "s", "^"],
    ):
        sub = obs_comp[obs_comp["compartment_norm"] == c].dropna(
            subset=["hydB_1l_ge_raw", "shannon_16s"]
        )
        ax.scatter(
            np.log2(sub["hydB_1l_ge_raw"] + 0.01),
            sub["shannon_16s"],
            c=col,
            marker=mk,
            s=30,
            alpha=0.80,
            label=lbl,
        )
    ax.set_xlabel(r"$\log_2(\mathrm{hydB\_1l}\ \mathrm{GE} + 0.01)$")
    ax.set_ylabel("16S rRNA Shannon Diversity (H')")
    ax.set_title("c  Chemolithoautotrophy vs Community Diversity")
    ax.legend(loc="lower left", frameon=True)
    ax.grid(alpha=0.25, linestyle="--")

    # Panel D: Empirical Phylum Composition (16S ASV n=180 vs CoverM MAG n=155) & Unmapped Read Fraction
    ax = axes[1, 1]
    x_pos = np.arange(len(order))
    w_bar = 0.34
    phyla_pairs = [
        ("actino_16s_rel_raw", "mag_actinomycetota_raw", "Actinomycetota", "#2b83ba"),
        ("thermoplasmatota_16s_rel_raw", "mag_thermoplasmatota_raw", "Thermoplasmatota", "#f46d43"),
        ("halobacteriota_16s_rel_raw", "mag_halobacteriota_raw", "Halobacteriota", "#d53e4f"),
    ]
    bot_16s = np.zeros(len(order))
    bot_mag = np.zeros(len(order))
    for col_16s, col_mag, p_label, p_col in phyla_pairs:
        v16 = np.array(
            [
                float(comp.loc[comp["compartment_norm"] == c, col_16s].dropna().mean()) * 100.0
                for c in order
            ]
        )
        vmag = np.array(
            [
                float(obs_comp.loc[obs_comp["compartment_norm"] == c, col_mag].dropna().mean()) * 100.0
                for c in order
            ]
        )
        ax.bar(
            x_pos - w_bar / 2,
            v16,
            bottom=bot_16s,
            width=w_bar,
            color=p_col,
            edgecolor="black",
            linewidth=0.35,
            label=f"{p_label} (16S / MAG)",
        )
        ax.bar(
            x_pos + w_bar / 2,
            vmag,
            bottom=bot_mag,
            width=w_bar,
            color=p_col,
            edgecolor="black",
            linewidth=0.35,
            hatch="//",
            alpha=0.82,
        )
        bot_16s += v16
        bot_mag += vmag

    ax.bar(
        x_pos - w_bar / 2,
        np.maximum(0.0, 100.0 - bot_16s),
        bottom=bot_16s,
        width=w_bar,
        color="#e0e0e0",
        edgecolor="black",
        linewidth=0.3,
        label="Other phyla (solid=16S, hatch=MAG)",
    )
    ax.bar(
        x_pos + w_bar / 2,
        np.maximum(0.0, 100.0 - bot_mag),
        bottom=bot_mag,
        width=w_bar,
        color="#e0e0e0",
        edgecolor="black",
        linewidth=0.3,
        hatch="//",
        alpha=0.82,
    )
    unmapped_means = [
        float(obs_comp.loc[obs_comp["compartment_norm"] == c, "unmapped_read_pct_raw"].dropna().mean())
        for c in order
    ]
    ax2 = ax.twinx()
    ax2.plot(
        x_pos,
        unmapped_means,
        "ko--",
        linewidth=1.5,
        markersize=5.5,
        label="Unmapped WGS reads (%)",
    )
    ax2.set_ylim(0, 100)
    ax2.set_ylabel("Unmapped Shotgun Reads (%)")
    ax.set_xticks(x_pos)
    ax.set_xticklabels(["Surface", "Subsurface (5-10 cm)", "Root-adjacent"])
    ax.set_ylim(0, 100)
    ax.set_ylabel("Phylum Share (% of 16S or Mapped MAG Reads)")
    ax.set_title("d  16S ASV (Solid) vs CoverM MAG (Hatched) & Unmapped WGS")
    ax.legend(loc="upper left", frameon=True, fontsize=7.0)
    ax2.legend(loc="upper right", frameon=True, fontsize=7.0)
    ax.grid(axis="y", alpha=0.25, linestyle="--")

    out_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_dir / "fig1_transect_cube_architecture.pdf", dpi=300)
    plt.close(fig)


def generate_fig2(root: Path, out_dir: Path) -> None:
    rep = json.loads((root / "data/validation/validation_report.json").read_text())

    fig, axes = plt.subplots(2, 2, figsize=(12.4, 8.2), constrained_layout=True)

    # Panel A: Tier 1 Spatial CV, Tier 2 Campaign LOCO CV & Prospective Trip 5 Forecast
    ax = axes[0, 0]
    t1 = rep["tier1_spatial_block_cv"]["targets"]
    t2 = rep["tier2_temporal_loco_cv"]["pooled_metrics"]["log2_hydB_1l"]
    t_prosp = rep["tier2_temporal_loco_cv"]["trip5_prospective_forecast_benchmark"]
    p_lvl = t_prosp["level_forecasts_log2_hydB_1l"]["longitudinal_twin"]
    p_del = t_prosp["temporal_anomaly_forecasts_delta_trip5_minus_history"]["longitudinal_twin"]

    names = [
        "hydB_1l\n(Spatial CV)",
        "16S Shannon\n(Spatial CV)",
        "vOTU Richness\n(Spatial CV)",
        "hydB_1l\n(Campaign LOCO)",
        "Trip 5 Level\n(Prospective)",
        "Trip 5 Delta\n(Within-Site)",
    ]
    rhos = [
        t1["log2_hydB_1l"]["pearson_r"],
        t1["shannon_16s_site_compartment"]["pearson_r"],
        t1["log1p_viral_richness"]["pearson_r"],
        t2["pearson_r"],
        p_lvl["pearson_r"],
        p_del["pearson_r"],
    ]
    spearmans = [
        t1["log2_hydB_1l"]["spearman_rho"],
        t1["shannon_16s_site_compartment"]["spearman_rho"],
        t1["log1p_viral_richness"]["spearman_rho"],
        t2["spearman_rho"],
        p_lvl["spearman_rho"],
        p_del["spearman_rho"],
    ]
    r2s = [
        t1["log2_hydB_1l"]["r2"],
        t1["shannon_16s_site_compartment"]["r2"],
        t1["log1p_viral_richness"]["r2"],
        t2["r2"],
        p_lvl["r2"],
        p_del["r2"],
    ]
    x = np.arange(len(names))
    w = 0.25
    ax.bar(x - w, rhos, width=w, color="#2b83ba", label="Pearson r")
    ax.bar(
        x,
        spearmans,
        width=w,
        color="#abdda4",
        edgecolor="black",
        linewidth=0.4,
        label=r"Spearman $\rho$",
    )
    ax.bar(x + w, r2s, width=w, color="#fdae61", label=r"Out-of-sample $R^2$")
    ax.set_xticks(x)
    ax.set_xticklabels(names, rotation=16, ha="right")
    ax.set_ylim(0, 0.72)
    ax.set_ylabel("Out-of-Sample Metric")
    ax.set_title("a  Spatial CV, Campaign LOCO & Trip 5 Holdout Forecast")
    ax.legend(loc="upper right", frameon=True)
    ax.grid(axis="y", alpha=0.25, linestyle="--")

    # Panel B: Tier 3 Ontological Polytope Projection & MinT Reconciliation (n=550)
    ax = axes[0, 1]
    t3 = rep["tier3_causal_ontological_projection"]
    cats = ["Mean CVM\n(Violation, GE)", "All-Node\nRMSE (GE)", "Leaf-Guild\nRMSE (GE)"]
    unconstrained = [
        t3["mean_cvm_unconstrained"],
        t3["rmse_unconstrained"],
        t3["leaf_rmse_unconstrained"],
    ]
    projected_euclid = [
        t3["mean_cvm_projected_two_layer"],
        t3["rmse_projected_two_layer"],
        t3["leaf_rmse_projected_two_layer"],
    ]
    projected_mint = [
        0.0,
        t3["rmse_nonnegative_mint"],
        t3["leaf_rmse_nonnegative_mint"],
    ]
    x3 = np.arange(len(cats))
    w3 = 0.25
    ax.bar(
        x3 - w3,
        unconstrained,
        width=w3,
        color="#d7191c",
        alpha=0.82,
        label=r"Unconstrained $\hat{y}$ (OCS=0.00)",
    )
    ax.bar(
        x3,
        projected_euclid,
        width=w3,
        color="#1a9641",
        alpha=0.85,
        label=r"Euclidean $\Pi_{\mathcal{C}(\mathcal{T}_{\mathrm{ont}})}$ (OCS=1.00)",
    )
    ax.bar(
        x3 + w3,
        projected_mint,
        width=w3,
        color="#2b83ba",
        alpha=0.85,
        label=r"Nonnegative Shrinkage-MinT (OCS=1.00)",
    )
    ax.set_xticks(x3)
    ax.set_xticklabels(cats)
    ax.set_ylim(0, 0.56)
    ax.set_ylabel("Genome Equivalents (GE)")
    ax.set_title("b  OWL 2 Polytope Projection & MinT Reconciliation (n=550)")
    ax.legend(loc="upper left", frameon=True)
    ax.grid(axis="y", alpha=0.25, linestyle="--")

    # Panel C: Tier 4 Paired Calibration + 97.4% Lineage-Turnover & CBB Co-Carriage Decomposition
    ax = axes[1, 0]
    r_paradox = rep["tier4_natural_experiment_calibration"]["root_vs_bulk_rhizosphere_paradox"]
    s_refugium = rep["tier4_natural_experiment_calibration"]["surface_vs_subsurface_refugium"]
    cbb_sub = {
        r["quantity"]: r
        for r in rep["tier4_natural_experiment_calibration"]["h2_oxidation_vs_cbb_autotrophy_decoupling"][
            "cbb_subsets_position_and_aridity"
        ]
        if r["contrast"] == "root/bulk"
    }

    y_pos = np.array([0, 1, 2, 3])
    emp_med = np.array(
        [
            r_paradox["empirical_median_log2fc_hydB_1l"],
            cbb_sub["l1_cbb"]["median_log2"],
            cbb_sub["l1_nocbb"]["median_log2"],
            s_refugium["empirical_median_log2fc_hydB_1l"],
        ]
    )
    ci_low = np.array(
        [
            r_paradox["empirical_95ci"][0],
            cbb_sub["l1_cbb"]["ci_lo"],
            cbb_sub["l1_nocbb"]["ci_lo"],
            s_refugium["empirical_95ci"][0],
        ]
    )
    ci_high = np.array(
        [
            r_paradox["empirical_95ci"][1],
            cbb_sub["l1_cbb"]["ci_hi"],
            cbb_sub["l1_nocbb"]["ci_hi"],
            s_refugium["empirical_95ci"][1],
        ]
    )

    ax.errorbar(
        emp_med,
        y_pos + 0.08,
        xerr=[emp_med - ci_low, ci_high - emp_med],
        fmt="o",
        color="black",
        capsize=3.5,
        label="Empirical paired median (95% CI)",
    )
    ax.scatter(
        [
            r_paradox["twin_predicted_median_log2fc_hydB_1l"],
            s_refugium["twin_predicted_median_log2fc_hydB_1l"],
        ],
        [0 - 0.12, 3 - 0.12],
        color="#d95f02",
        marker="D",
        s=52,
        zorder=4,
        label="EQ-BioDT scenario prediction",
    )
    ax.axvline(0.0, color="gray", linestyle=":", linewidth=1.0)
    ax.set_yticks(y_pos)
    ax.set_yticklabels(
        [
            "Root vs Bulk All 1l\n(97.4% genus turnover)",
            "Root vs Bulk 1l+CBB\n(101/558 MAGs, 18.1%)",
            "Root vs Bulk 1l no-CBB\n(457/558 MAGs, 81.9%)",
            "Surface vs 5-10 cm\n(n=68 paired visits)",
        ]
    )
    ax.set_xlabel(r"Paired Contrast ($\log_2$ Fold Change)")
    ax.set_title("c  Paired Contrasts, Turnover & CBB Co-Carriage")
    ax.legend(loc="lower left", frameon=True)
    ax.grid(axis="x", alpha=0.25, linestyle="--")

    # Panel D: Space-for-Time Transport Test (Between-site spatial vs Within-site temporal slopes on P90_era5)
    ax = axes[1, 1]
    sft = rep["tier2_temporal_loco_cv"]["space_for_time_transport"]["h123_within_vs_between_hydB_1l"]
    p90_rows = {r["compartment"]: r for r in sft if r["exposure"] == "P90_era5"}
    comps_sft = ["Surface", "PlantRoot", "Deep"]
    lbls_sft = [
        "Surface (0-2 cm)\n(n=79 visits, q=0.013)",
        "Root-adjacent\n(n=193 visits, q=0.020)",
        "Subsurface (5-10 cm)\n(n=103 visits, q=0.900)",
    ]
    between_b = [p90_rows[c]["between_slope"] for c in comps_sft]
    within_b = [p90_rows[c]["within_slope"] for c in comps_sft]
    within_lo = [p90_rows[c]["within_boot_lo"] for c in comps_sft]
    within_hi = [p90_rows[c]["within_boot_hi"] for c in comps_sft]
    y_sft = np.arange(len(comps_sft))
    ax.barh(
        y_sft - 0.16,
        between_b,
        height=0.28,
        color="#fc8d59",
        edgecolor="black",
        linewidth=0.4,
        label="Between-site spatial slope",
    )
    ax.barh(
        y_sft + 0.16,
        within_b,
        height=0.28,
        color="#2b83ba",
        edgecolor="black",
        linewidth=0.4,
        xerr=[
            np.array(within_b) - np.array(within_lo),
            np.array(within_hi) - np.array(within_b),
        ],
        capsize=3,
        label="Within-site temporal slope (95% CI)",
    )
    ax.axvline(0.0, color="black", linestyle=":", linewidth=0.9)
    ax.set_yticks(y_sft)
    ax.set_yticklabels(lbls_sft)
    ax.set_xlabel(r"Slope of $\log_2(\mathrm{hydB\_1l})$ on $\log(1 + P_{90}^{\mathrm{ERA5}})$")
    ax.set_title("d  Space-for-Time Transport Test Across Revisits")
    ax.legend(loc="lower left", frameon=True)
    ax.grid(axis="x", alpha=0.25, linestyle="--")

    fig.savefig(out_dir / "fig2_validation_causal_polytope.pdf", dpi=300)
    plt.close(fig)


def generate_fig3(root: Path, out_dir: Path) -> None:
    comp = pd.read_parquet(root / "data/processed/eq_compartment_state_cube.parquet")
    s30_row = (
        comp[(comp["site_num"] == 30) & (comp["compartment_norm"] == "Surface")].iloc[0].to_dict()
    )
    s60_row = (
        comp[(comp["site_num"] == 60) & (comp["compartment_norm"] == "Surface")].iloc[0].to_dict()
    )

    fig, axes = plt.subplots(2, 2, figsize=(12.4, 8.2), constrained_layout=True)

    # Panel A: Cloud Seeding 14-day Impulse Response across Precipitation Pulses
    ax = axes[0, 0]
    cloud = CloudSeedingIntervention()
    for mm, col, ls in [(6.0, "#1b9e77", "-"), (15.0, "#d95f02", "--"), (25.0, "#7570b3", "-.")]:
        res = cloud.simulate(s30_row, rain_pulse_mm=mm, horizon_days=14)
        t = [p["day"] for p in res["trajectory"]]
        h2_pw = [p["h2_oxidation_power_uw_m3"] for p in res["trajectory"]]
        ax.plot(t, h2_pw, color=col, linestyle=ls, linewidth=1.8, label=f"+{int(mm)} mm pulse")
    ax.set_xlabel("Time Since Prescribed Rainfall Pulse (days)")
    ax.set_ylabel(r"Atmospheric $\mathrm{H}_2$ Oxidation Power ($\mu\mathrm{W\ m^{-3}}$)")
    ax.set_title("a  Millington-Quirk Moisture Pulse Response (S30)")
    ax.legend(loc="lower right", frameon=True)
    ax.grid(alpha=0.25, linestyle="--")

    # Panel B: Revegetation Rhizosphere Paradox (Area-Weighted, Biomass-Weighted & Absolute Copy Density)
    ax = axes[0, 1]
    reveg = RevegetationIntervention()
    fracs = np.linspace(0.05, 0.80, 12)
    area_wt = []
    biom_wt = []
    abs_mult = []
    for f in fracs:
        s = reveg.simulate(
            s30_row,
            shrub_cover_fraction=float(f),
            rak_pr_pgpr_dose=0.85,
            site60_rbcl_hyab_plasmid_dose=0.75,
        )["summary"]
        area_wt.append(s["landscape_hydB_1l_ge"])
        biom_wt.append(s["biomass_weighted_hydB_1l_ge"])
        abs_mult.append(s["absolute_hydB_1l_density_multiplier"])
    ax.plot(
        fracs * 100,
        area_wt,
        "o-",
        color="#1b9e77",
        linewidth=1.6,
        label="Area-weighted hydB_1l (GE)",
    )
    ax.plot(
        fracs * 100,
        biom_wt,
        "s--",
        color="#d95f02",
        linewidth=1.6,
        label="Biomass-weighted hydB_1l (GE)",
    )
    ax.plot(
        fracs * 100,
        abs_mult,
        "^:",
        color="#7570b3",
        linewidth=1.7,
        label="Absolute copy density (x baseline)",
    )
    ax.set_xlabel("Shrub Canopy Cover (%)")
    ax.set_ylabel("Capacity (GE) or Density Multiplier (x)")
    ax.set_title("b  Rhizosphere Paradox: Relative vs Absolute Capacity")
    ax.legend(loc="upper right", frameon=True)
    ax.grid(alpha=0.25, linestyle="--")

    # Panel C: Common-Baseline 4x2x2 Factorial Brine Analysis (S60 Baseline)
    ax = axes[1, 0]
    desal = DesalinationBrineIntervention()
    f_grid = desal.evaluate_factorial_brine_grid(s60_row)
    ec_levels = [10.65, 35.0, 85.0, 216.2]
    for q_val, inoc_val, lbl, col, mk, ls in [
        (1.4, 0.6, "Ca-SO4 (q=1.4) + ect/rbcL inoc", "#1a9641", "o", "-"),
        (1.4, 0.0, "Ca-SO4 (q=1.4), uninoculated", "#66bd63", "s", "--"),
        (0.1, 0.6, "Na-Cl (q=0.1) + ect/rbcL inoc", "#fdae61", "^", "-."),
        (0.1, 0.0, "Na-Cl (q=0.1), uninoculated", "#d7191c", "v", ":"),
    ]:
        sub_rows = [
            r
            for r in f_grid
            if abs(r["ca_so4_to_nacl_ratio"] - q_val) < 1e-6
            and abs(r["sphingomonas_ect_kdp_plasmid_dose"] - inoc_val) < 1e-6
        ]
        sub_rows.sort(key=lambda r: r["discharge_salinity_ms_cm"])
        c_flux = [r["autotrophic_carbon_fixation_flux"] for r in sub_rows]
        ax.plot(
            range(len(ec_levels)),
            c_flux,
            marker=mk,
            linestyle=ls,
            color=col,
            linewidth=1.6,
            label=lbl,
        )
    ax.set_xticks(range(len(ec_levels)))
    ax.set_xticklabels(["10.65\n(Site 60)", "35.0\n(Brackish)", "85.0\n(Saline)", "216.2\n(Site 59)"])
    ax.set_xlabel(r"Discharge Electrical Conductivity ($\mathrm{mS\ cm^{-1}}$)")
    ax.set_ylabel("Autotrophic C Fixation Flux (model units)")
    ax.set_title("c  Common-Baseline 4x2x2 Factorial Brine Audit")
    ax.legend(loc="upper right", frameon=True)
    ax.grid(alpha=0.25, linestyle="--")

    # Panel D: 4-Objective Portfolio Frontier & Monte Carlo Dominance Probabilities (B=200)
    ax = axes[1, 1]
    opt = ParetoPortfolioOptimizer()
    pareto = opt.optimize_pareto_front(s30_row)
    all_p = pareto["all_evaluated_portfolios"]
    j1 = [p["objectives"]["J1_biodiversity_refugium_max"] for p in all_p]
    j2 = [p["objectives"]["J2_chemosynthesis_cn_flux_max"] for p in all_p]
    p_dom = [p["probabilistic_robustness"]["pareto_dominance_prob"] for p in all_p]
    sc = ax.scatter(
        j2,
        j1,
        c=p_dom,
        cmap="viridis",
        vmin=0.0,
        vmax=1.0,
        s=68,
        alpha=0.92,
        edgecolors="black",
        linewidths=0.6,
    )
    short_codes = {
        "Baseline Hyperarid Do-Nothing": "P1",
        "Integrated Eco-Restoration & Oryx Rewilding (Balanced Optimum)": "P2",
        "Low-Water Chemosynthetic Crust & Refugium Conservation": "P3",
        "High-Canopy Oasis Revegetation + Site 60 Ca-SO4 Irrigation": "P4",
        "Active Energy-Corridor Nodal Seismic + Oil-Spill Bioremediation": "P5",
        "Unmitigated Heavy Seeding + Brine Disposal + Camel Overgrazing + Vibroseis": "P6",
    }
    for p in all_p:
        code = short_codes.get(p["policy"]["name"], "")
        dom = p["probabilistic_robustness"]["pareto_dominance_prob"]
        ax.annotate(
            f"{code} ({dom:.2f})",
            (
                p["objectives"]["J2_chemosynthesis_cn_flux_max"],
                p["objectives"]["J1_biodiversity_refugium_max"],
            ),
            textcoords="offset points",
            xytext=(4, 4),
            fontsize=7.5,
        )
    rec = pareto["recommended_policy"]
    ax.scatter(
        [rec["objectives"]["J2_chemosynthesis_cn_flux_max"]],
        [rec["objectives"]["J1_biodiversity_refugium_max"]],
        color="cyan",
        marker="*",
        s=185,
        edgecolors="black",
        linewidths=0.9,
        zorder=5,
        label="Top utility & dominance (P3)",
    )
    cbar = fig.colorbar(sc, ax=ax, shrink=0.84, pad=0.02)
    cbar.set_label("Monte Carlo Pareto Dominance Prob. (B=200)")
    ax.set_xlabel("Chemosynthetic & C/N Flux Score (J2)")
    ax.set_ylabel("Biodiversity & Refugium Score (J1)")
    ax.set_title("d  Multi-Objective Frontier & Probabilistic Robustness")
    ax.legend(loc="lower right", frameon=True)
    ax.grid(alpha=0.25, linestyle="--")

    fig.savefig(out_dir / "fig3_interventions_pareto.pdf", dpi=300)
    plt.close(fig)


if __name__ == "__main__":
    repo_root = Path(__file__).resolve().parents[1]
    figures_dir = repo_root / "paper" / "figures"
    generate_fig1(repo_root, figures_dir)
    generate_fig2(repo_root, figures_dir)
    generate_fig3(repo_root, figures_dir)
    print("Generated all 3 4-panel vector PDF figures in", figures_dir)
