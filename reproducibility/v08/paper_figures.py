"""Main-text figures and answer table for the restructured EQ-BioDT v0.8 manuscript.

Generates:
  - paper/figures/v08_composition_main.{pdf,png}  (Figure A: How composition extends the twin)
  - paper/figures/v08_answers.{pdf,png}           (Figure B: Management answers under joint uncertainty)
  - paper/tables/v08_answers.tex                  (Main-text Table: 14 management questions)

Also exports READABLE_NAMES, QUESTION_TITLES_SHORT, and helper utilities shared with
reproducibility/v08/figures.py for supplementary v08_sensitivity.{pdf,png} and v08_uq.tex.
"""
from __future__ import annotations

import json
from pathlib import Path
import sys
from typing import Any

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.lines import Line2D
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from reproducibility.v08.paper_numbers import (  # noqa: E402
    _fmt_sig_latex,
    fmt_frac_pct,
    fmt_pct_value,
    fmt_st,
)
QDIR = ROOT / "reviews" / "v08" / "questions"
UQDIR = ROOT / "reviews" / "v08" / "uq"
V07_COMP_DIR = ROOT / "reviews" / "v07" / "composition"
DE_ZELICOURT_PATH = ROOT / "data" / "external" / "de_zelicourt_2018" / "fig1_digitised.json"
FIG_DIR = ROOT / "paper" / "figures"
TAB_DIR = ROOT / "paper" / "tables"

MM_TO_INCH = 1.0 / 25.4
FIG_WIDTH_IN = 180.0 * MM_TO_INCH  # 180 mm Nature Communications full-width figure

plt.rcParams.update(
    {
        "font.family": "DejaVu Sans",
        "font.size": 6.8,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "savefig.dpi": 300,
        "mathtext.fontset": "dejavusans",
        "svg.hashsalt": "eq-biodt-v08-paper",
    }
)

INK = "#26323A"
MUTED = "#5A6872"
GRID_GREY = "#D9E0E6"
LIGHT_BG = "#F6F8FA"

#Consistent publication palette across v0.6-v0.8
COL_TWIN = "#D2693C"       # terracotta / rust
COL_EXTERNAL = "#4E8FA6"   # steel blue
COL_DEDUCED = "#4F9A82"    # sage green
COL_ONTOLOGY = "#8E7CB3"   # muted purple
COL_ACCENT = "#D9A441"     # ochre / gold

# Direction colours in Figure B
COL_DECREASE = "#D2693C"      # decrease in >=95% of draws
COL_INCREASE = "#4F9A82"      # increase in >=95% of draws
COL_UNDETERMINED = "#4E8FA6"  # direction undetermined

# Single centralized mapping from code parameter/constant identifiers to readable scientific names
READABLE_NAMES: dict[str, str] = {
    # Twin process-model parameters
    "k_h2_s": r"H$_2$ rate constant $k_{\mathrm{H_2}}$",
    "k_co_s": r"CO rate constant $k_{\mathrm{CO}}$",
    "retention_b": r"retention exponent $b$",
    "psi_entry_mpa": r"air-entry potential $\psi_e$",
    "q10": r"temperature coefficient $Q_{10}$",
    "dry_threshold_mpa": r"activation threshold $\psi_d$",
    "psi_d": r"activation threshold $\psi_d$",
    "wet_activation_mpa": r"full-activation potential $\psi_w$",
    "psi_w": r"full-activation potential $\psi_w$",
    "activation_shape": r"activation shape $a$",
    "a": r"activation shape $a$",
    "porosity": r"porosity $\phi$",
    "drainage_fraction_day": "drainage fraction",
    "drainage": "drainage fraction",
    "field_capacity_fraction": "field-capacity fraction",
    "field_capacity": "field-capacity fraction",
    "d_h2_m2_s": "free-air diffusivity",
    "d_co_m2_s": "free-air diffusivity",
    "D_H2": "free-air diffusivity",
    "D_CO": "free-air diffusivity",
    # External composition constants (derived from registry descriptions)
    "cloud_seeding_rain_increase_frac": "rain increase per seeding",
    "solar_panel_evap_reduction_frac": "evaporation reduction under panels",
    "solar_panel_soil_temp_reduction_c": "soil cooling under panels",
    "panel_wash_l_m2_per_wash": "panel wash water per wash",
    "panel_wash_interval_days": "panel wash interval",
    "evap_scaling_per_c": "evaporation increase per °C",
    "soil_temp_per_air_temp_warming": "soil-to-air warming ratio",
    "mulch_evap_reduction_frac": "evaporation reduction under mulch",
    "shade_net_soil_temp_reduction_c": "soil cooling under shade net",
    "biochar_porosity_increase_per_pct": "porosity increase per % biochar",
    "offroad_bd_max_increase_g_cm3": "maximum bulk-density increase",
    "offroad_bd_passes_efold": "compaction e-folding passes",
    "particle_density_g_cm3": "mineral particle density",
    "knallgas_co2_per_h2": r"CO$_2$ fixed per H$_2$ oxidised",
    "molar_mass_c_g_mol": "carbon molar mass",
    "days_per_year": "days per year",
    "biocrust_net_c_fixation_g_c_m2_yr": "biocrust net carbon uptake",
    "h2_background_ppm": r"background H$_2$ mixing ratio",
    "h2_hydrogen_economy_global_ppm": r"hydrogen-economy H$_2$ mixing ratio",
    "nacl_mol_l_per_dS_m": r"NaCl concentration per dS m$^{-1}$",
    "fao29_ece_per_ecw": "FAO leaching concentration factor",
    "fao29_ecsw_per_ece": "soil-water to extract EC ratio",
    "crop_etc_mm_per_day": "crop water demand",
    "alfalfa_etc_mm_per_day": "alfalfa irrigation demand",
    "mh_date_palm_threshold_dS_m": "date-palm salinity threshold",
    "mh_date_palm_slope_pct_per_dS_m": "date-palm salinity slope",
    "mh_barley_threshold_dS_m": "barley salinity threshold",
    "mh_barley_slope_pct_per_dS_m": "barley salinity slope",
    "mh_alfalfa_threshold_dS_m": "alfalfa salinity threshold",
    "mh_alfalfa_slope_pct_per_dS_m": "alfalfa salinity slope",
    "mh_tomato_threshold_dS_m": "tomato salinity threshold",
    "mh_tomato_slope_pct_per_dS_m": "tomato salinity slope",
    "sa187_mock_ratio_fresh": "trial yield ratio (fresh)",
    "sa187_mock_ratio_dry": "trial yield ratio (dry)",
    "sa187_benefit_b": "SA187 yield benefit",
    "sa187_effect_fresh_low_s1": "SA187 low-salinity effect (fresh, S1)",
    "sa187_effect_fresh_low_s2": "SA187 low-salinity effect (fresh, S2)",
    "sa187_effect_dry_low_s1": "SA187 low-salinity effect (dry, S1)",
    "sa187_effect_dry_low_s2": "SA187 low-salinity effect (dry, S2)",
    "sa187_ecw_low_dS_m": "trial low-salinity irrigation EC",
    "sa187_ecw_high_dS_m": "trial high-salinity irrigation EC",
    "sa187_soil_ec_dS_m": "trial site soil EC",
}

# Short readable titles for each question slug (used in supplementary figure and table)
QUESTION_TITLES_SHORT: dict[str, str] = {
    "rain-pattern": "Rain pattern",
    "drip-or-flood": "Drip vs flood timing",
    "cloud-seeding": "Cloud seeding",
    "salty-water": "Brackish irrigation",
    "alfalfa-bacterium": "Alfalfa with SA187",
    "hotter": "Climate warming",
    "mulch-shade": "Mulch and shade",
    "solar-park": "Solar park",
    "compaction": "Vehicle compaction",
    "biochar": "Biochar addition",
    "missing-co-genes": "Missing CO genes",
    "adding-microbes": "Adding CO microbes",
    "air-hydrogen": r"Atmospheric H$_2$",
    "carbon-fixation": "Carbon fixation",
}

# Readable headline level labels for supplementary table (no code identifiers or sample IDs)
HEADLINE_LABELS_READABLE: dict[str, str] = {
    "rain-pattern": "Four 3 mm showers",
    "drip-or-flood": "Irrigation every 10 days",
    "cloud-seeding": "Cloud seeding +10% rain",
    "salty-water": r"Daily irrigation 8 dS m$^{-1}$",
    "alfalfa-bacterium": r"7.81 dS m$^{-1}$ (fitted leaching)",
    "hotter": "+4 °C warming, rain unchanged",
    "mulch-shade": "Full organic mulch",
    "solar-park": "Full solar-panel cover",
    "compaction": "50 vehicle passes",
    "biochar": "Biochar 5% w/w",
    "missing-co-genes": "Open-gate community",
    "adding-microbes": "+1% CO-carrier share",
    "air-hydrogen": r"2$\times$ background H$_2$",
    "carbon-fixation": "Single 12 mm storm (open gate)",
}


def readable_name(param_id: str, *, latex: bool = False) -> str:
    """Return the readable scientific name for a parameter or external constant identifier."""
    if param_id not in READABLE_NAMES:
        raise KeyError(f"Unmapped parameter/constant identifier: {param_id!r}")
    s = READABLE_NAMES[param_id]
    if latex:
        return to_latex(s)
    return s


def to_latex(s: str) -> str:
    """Convert a readable label (which may already contain $...$ math) to safe LaTeX."""
    s = s.replace("°C", r"$^\circ$C").replace("°", r"$^\circ$")
    # Escape % only when not already escaped
    out: list[str] = []
    for i, ch in enumerate(s):
        if ch == "%" and (i == 0 or s[i - 1] != "\\"):
            out.append(r"\%")
        elif ch == "&" and (i == 0 or s[i - 1] != "\\"):
            out.append(r"\&")
        else:
            out.append(ch)
    return "".join(out)


def _save_fig(fig: plt.Figure, stem: str) -> tuple[Path, Path]:
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    pdf_path = FIG_DIR / f"{stem}.pdf"
    png_path = FIG_DIR / f"{stem}.png"
    fig.savefig(
        pdf_path,
        metadata={"CreationDate": None, "ModDate": None, "Producer": None},
    )
    fig.savefig(
        png_path,
        dpi=300,
        metadata={"Software": None},
    )
    plt.close(fig)
    return pdf_path, png_path


# ==============================================================================
# Figure A: paper/figures/v08_composition_main.{pdf,png}
# ==============================================================================
def _draw_provenance_tag(
    ax: plt.Axes,
    x: float,
    y: float,
    tag: str,
    color: str,
) -> None:
    """Draw a compact coloured provenance pill tag on schematic panel a."""
    w, h = 10.5, 4.4
    pill = FancyBboxPatch(
        (x - w / 2.0, y - h / 2.0),
        w,
        h,
        boxstyle="round,pad=0.2,rounding_size=1.8",
        facecolor=color,
        edgecolor="none",
        alpha=0.92,
        zorder=6,
    )
    ax.add_patch(pill)
    ax.text(
        x,
        y,
        tag,
        ha="center",
        va="center",
        fontsize=5.3,
        fontweight="bold",
        color="white",
        zorder=7,
    )


def _draw_schematic_box(
    ax: plt.Axes,
    x: float,
    y: float,
    w: float,
    h: float,
    title: str,
    body: str,
    *,
    edge_color: str = INK,
    face_color: str = "white",
    prov_tag: tuple[str, str] | None = None,
    linestyle: str = "-",
) -> None:
    box = FancyBboxPatch(
        (x, y),
        w,
        h,
        boxstyle="round,pad=0.3,rounding_size=2.0",
        facecolor=face_color,
        edgecolor=edge_color,
        linewidth=1.0,
        linestyle=linestyle,
        zorder=3,
    )
    ax.add_patch(box)
    ax.text(
        x + w / 2.0,
        y + h - 4.3,
        title,
        ha="center",
        va="top",
        fontsize=6.1,
        fontweight="bold",
        color=INK,
        zorder=5,
    )
    ax.text(
        x + w / 2.0,
        y + (h - 5.5) / 2.0,
        body,
        ha="center",
        va="center",
        fontsize=5.6,
        color=INK,
        linespacing=1.18,
        zorder=5,
    )
    if prov_tag is not None:
        tag_text, tag_col = prov_tag
        _draw_provenance_tag(ax, x + w - 5.8, y + h + 0.4, tag_text, tag_col)


def _draw_arrow(
    ax: plt.Axes,
    p1: tuple[float, float],
    p2: tuple[float, float],
    *,
    color: str = INK,
    ls: str = "-",
    rad: float = 0.0,
) -> None:
    arr = FancyArrowPatch(
        p1,
        p2,
        arrowstyle="-|>,head_length=3.5,head_width=2.2",
        connectionstyle=f"arc3,rad={rad}",
        color=color,
        linewidth=1.0,
        linestyle=ls,
        zorder=4,
    )
    ax.add_patch(arr)


def build_figure_composition_main() -> tuple[Path, Path]:
    """Generate Figure A: paper/figures/v08_composition_main.{pdf,png}."""
    alf_doc = json.loads((V07_COMP_DIR / "inoculated_alfalfa_salinity.json").read_text(encoding="utf-8"))
    brack_doc = json.loads((V07_COMP_DIR / "brackish_irrigation.json").read_text(encoding="utf-8"))
    dig_doc = json.loads(DE_ZELICOURT_PATH.read_text(encoding="utf-8"))

    # Port validity interval from brackish_irrigation.json caveats ("0.1 mol/L")
    import re as _re
    m_dom = _re.search(r"(\d+(?:\.\d+)?)\s*mol/L", " ".join(brack_doc.get("caveats", [])))
    valid_hi = float(m_dom.group(1)) if m_dom else 0.1

    # Panel b numbers from reviews/v07/composition/inoculated_alfalfa_salinity.json
    # (cross-checked with data/external/de_zelicourt_2018/fig1_digitised.json)
    obs_fresh = alf_doc["observed"]["fresh"]
    obs_dry = alf_doc["observed"]["dry"]
    assert abs(float(obs_fresh["value"]) - float(dig_doc["season1_mock_ratio_high_over_low"]["fresh"])) < 1e-4
    assert abs(float(obs_dry["value"]) - float(dig_doc["season1_mock_ratio_high_over_low"]["dry"])) < 1e-4

    pred_fao = alf_doc["validation_mock_ratio"]["A_fao_fixed_leaching"]["predicted_ratio"]
    pred_twin = alf_doc["validation_mock_ratio"]["B_twin_physics"]["predicted_ratio"]

    # Panel c numbers from reviews/v07/composition/inoculated_alfalfa_salinity.json
    lvl_312 = next(lv for lv in alf_doc["levels"] if abs(float(lv["ec_w_dS_m"]) - 3.12) < 1e-6)
    fao_ece_312 = lvl_312["downstream"]["fao29_leaching"]["ece_dS_m"]
    fit_fresh = alf_doc["fitted_leaching"]["fresh"]
    fit_dry = alf_doc["fitted_leaching"]["dry"]
    measured_soil_ec = float(alf_doc["fitted_leaching"]["comparison_with_measured_soil_ec"]["measured_soil_ec_dS_m"])

    fig = plt.figure(figsize=(FIG_WIDTH_IN, 4.65))
    gs_top = fig.add_gridspec(
        1,
        1,
        left=0.045,
        right=0.985,
        top=0.94,
        bottom=0.50,
    )
    gs_bot = fig.add_gridspec(
        1,
        2,
        width_ratios=[1.08, 0.92],
        left=0.195,
        right=0.985,
        top=0.41,
        bottom=0.11,
        wspace=0.30,
    )

    ax_a = fig.add_subplot(gs_top[0, 0])
    ax_b = fig.add_subplot(gs_bot[0, 0])
    ax_c = fig.add_subplot(gs_bot[0, 1])

    # ---------------- Panel a: Schematic of composed chain ----------------
    ax_a.set_xlim(0, 180)
    ax_a.set_ylim(0, 56)
    ax_a.axis("off")
    ax_a.text(
        -5.5,
        54.8,
        "a",
        fontsize=9.5,
        fontweight="bold",
        color=INK,
        ha="left",
        va="top",
    )
    ax_a.text(
        0.5,
        54.3,
        "How composition extends the twin across management inputs, typed ports, and endpoints",
        fontsize=7.3,
        fontweight="bold",
        color=INK,
        ha="left",
        va="top",
    )

    # Top row boxes (brackish irrigation chain)
    y_top = 24.5
    h_top = 22.0
    # 1. Management input
    _draw_schematic_box(
        ax_a,
        0.5,
        y_top,
        28.0,
        h_top,
        "Management input",
        "Irrigation water EC\n(dS m$^{-1}$)",
        edge_color=MUTED,
        face_color=LIGHT_BG,
    )
    # 2. Upstream translation models
    _draw_schematic_box(
        ax_a,
        33.5,
        y_top,
        32.0,
        h_top,
        "Upstream translation",
        "EC to NaCl;\nirrigation schedule",
        edge_color=COL_EXTERNAL,
        face_color="#F2F7FA",
        prov_tag=("external", COL_EXTERNAL),
    )
    # 3. Typed ports
    _draw_schematic_box(
        ax_a,
        70.5,
        y_top,
        35.5,
        h_top,
        "Typed ports",
        "NaCl in added water\n(mol L$^{-1}$);\nwater depth (mm d$^{-1}$)",
        edge_color=INK,
        face_color="white",
    )
    # 4. Twin intervention & trajectory
    _draw_schematic_box(
        ax_a,
        111.0,
        y_top,
        33.5,
        h_top,
        "Twin intervention",
        "do(water, salt) $\\rightarrow$\nwater storage, salt,\nH$_2$ uptake",
        edge_color=COL_TWIN,
        face_color="#FDF5F1",
        prov_tag=("twin", COL_TWIN),
    )
    # 5. Downstream endpoint models
    _draw_schematic_box(
        ax_a,
        149.5,
        y_top,
        30.0,
        h_top,
        "Downstream models",
        "Soil-solution to EC$_e$;\nMaas–Hoffman\ncrop yield",
        edge_color=COL_EXTERNAL,
        face_color="#F2F7FA",
        prov_tag=("external", COL_EXTERNAL),
    )

    # Horizontal arrows along top chain
    y_mid = y_top + h_top / 2.0
    _draw_arrow(ax_a, (28.5, y_mid), (33.5, y_mid))
    _draw_arrow(ax_a, (65.5, y_mid), (70.5, y_mid))
    _draw_arrow(ax_a, (106.0, y_mid), (111.0, y_mid))
    _draw_arrow(ax_a, (144.5, y_mid), (149.5, y_mid))

    # Port validity check callout below Typed ports
    val_box = FancyBboxPatch(
        (49.5, 2.0),
        60.5,
        15.0,
        boxstyle="round,pad=0.25,rounding_size=1.8",
        facecolor="#FFF8EC",
        edgecolor=COL_ACCENT,
        linewidth=0.95,
        linestyle="--",
        zorder=3,
    )
    ax_a.add_patch(val_box)
    ax_a.text(
        79.75,
        9.5,
        f"Port validity check: dilute osmotic domain [0, {valid_hi:g}] mol L$^{{-1}}$\n"
        r"Refusal when outside interval (e.g. $>0.1$ mol L$^{-1}$)",
        ha="center",
        va="center",
        fontsize=5.3,
        color=INK,
        linespacing=1.18,
        zorder=5,
    )
    _draw_arrow(ax_a, (85.0, y_top), (85.0, 17.2), color=COL_ACCENT, ls="--")

    # Second branch: Deductive gate on knallgas carbon fixation
    gate_box = FancyBboxPatch(
        (114.5, 2.0),
        65.0,
        15.0,
        boxstyle="round,pad=0.25,rounding_size=1.8",
        facecolor="#F1F8F5",
        edgecolor=COL_DEDUCED,
        linewidth=0.95,
        zorder=3,
    )
    ax_a.add_patch(gate_box)
    _draw_provenance_tag(ax_a, 173.2, 17.3, "deduced", COL_DEDUCED)
    ax_a.text(
        145.5,
        9.5,
        "Deductive gate on knallgas carbon fixation:\n"
        r"H$_2$-autotrophic candidate share $= 0 \rightarrow$ deduced $0$",
        ha="center",
        va="center",
        fontsize=5.3,
        color=INK,
        linespacing=1.18,
        zorder=5,
    )
    _draw_arrow(ax_a, (128.0, y_top), (128.0, 17.2), color=COL_DEDUCED)

    # Provenance legend stacked cleanly in bottom-left of Panel a
    ax_a.text(0.5, 18.2, "Provenance tags:", fontsize=5.6, fontweight="bold", color=INK, va="center")
    _draw_provenance_tag(ax_a, 6.2, 13.0, "twin", COL_TWIN)
    ax_a.text(12.8, 13.0, "process model", fontsize=5.3, color=MUTED, va="center")
    _draw_provenance_tag(ax_a, 6.2, 7.6, "external", COL_EXTERNAL)
    ax_a.text(12.8, 7.6, "literature map", fontsize=5.3, color=MUTED, va="center")
    _draw_provenance_tag(ax_a, 6.2, 2.2, "deduced", COL_DEDUCED)
    ax_a.text(12.8, 2.2, "ontological closure", fontsize=5.3, color=MUTED, va="center")

    # ---------------- Panel b: Hada Al-Sham alfalfa trial test ----------------
    ax_b.text(
        -0.38,
        1.13,
        "b",
        transform=ax_b.transAxes,
        fontsize=9.5,
        fontweight="bold",
        color=INK,
        ha="left",
        va="top",
    )
    ax_b.set_title(
        "Hada Al-Sham alfalfa trial (7.81 vs 3.12 dS m$^{-1}$)",
        fontsize=7.0,
        fontweight="bold",
        color=INK,
        loc="left",
        pad=5,
    )

    b_items = [
        (
            "Observed trial\n(fresh biomass)",
            float(obs_fresh["value"]),
            float(obs_fresh["low"]),
            float(obs_fresh["high"]),
            INK,
            "o",
            f"{obs_fresh['value']:.2f} ({obs_fresh['low']:.3f}–{obs_fresh['high']:.3f})",
        ),
        (
            "Observed trial\n(dry biomass)",
            float(obs_dry["value"]),
            float(obs_dry["low"]),
            float(obs_dry["high"]),
            INK,
            "s",
            f"{obs_dry['value']:.2f} ({obs_dry['low']:.3f}–{obs_dry['high']:.3f})",
        ),
        (
            "Standard fixed leaching\n(FAO, $f = 1.5$)",
            float(pred_fao["value"]),
            float(pred_fao["low"]),
            float(pred_fao["high"]),
            COL_EXTERNAL,
            "D",
            f"{pred_fao['value']:.2f} ({pred_fao['low']:.2f}–{pred_fao['high']:.3f})",
        ),
        (
            "Twin bare-sand\nsalt balance",
            float(pred_twin["value"]),
            float(pred_twin["low"]),
            float(pred_twin["high"]),
            COL_TWIN,
            "^",
            f"{pred_twin['value']:.3f} ({pred_twin['low']:.3f}–{pred_twin['high']:.3f})",
        ),
    ]

    # Shaded observed envelope band
    ax_b.axvspan(
        float(obs_fresh["low"]),
        float(obs_dry["high"]),
        color=GRID_GREY,
        alpha=0.55,
        zorder=1,
    )
    ax_b.axhline(1.5, color=GRID_GREY, lw=0.7, ls="--", zorder=1)

    y_pos_b = np.arange(len(b_items))[::-1]
    for yi, (lbl, val, lo, hi, col, marker, txt) in zip(y_pos_b, b_items):
        ax_b.errorbar(
            val,
            yi,
            xerr=[[max(0.0, val - lo)], [max(0.0, hi - val)]],
            fmt=marker,
            color=col,
            ecolor=col,
            markersize=4.8,
            elinewidth=1.3,
            capsize=2.5,
            zorder=4,
        )
        # Place numerical annotation cleanly without overlapping error bar
        text_x = hi + 0.025 if hi < 0.70 else lo - 0.025
        ha_align = "left" if hi < 0.70 else "right"
        ax_b.text(
            text_x,
            yi,
            txt,
            va="center",
            ha=ha_align,
            fontsize=5.6,
            color=col,
            fontweight="bold" if col == COL_TWIN else "normal",
        )

    ax_b.set_yticks(y_pos_b)
    ax_b.set_yticklabels([item[0] for item in b_items], fontsize=6.0, color=INK)
    ax_b.set_xlim(-0.03, 1.05)
    ax_b.set_ylim(-0.5, 3.5)
    ax_b.set_xlabel("High-to-low salinity mock yield ratio", fontsize=6.5, color=INK)
    ax_b.tick_params(axis="both", labelsize=6.0)

    # ---------------- Panel c: Implied root-zone ECe at 3.12 dS m^-1 ----------------
    ax_c.text(
        -0.18,
        1.13,
        "c",
        transform=ax_c.transAxes,
        fontsize=9.5,
        fontweight="bold",
        color=INK,
        ha="left",
        va="top",
    )
    ax_c.set_title(
        "Implied root-zone EC$_e$ at 3.12 dS m$^{-1}$",
        fontsize=6.8,
        fontweight="bold",
        color=INK,
        loc="left",
        pad=5,
    )

    c_items = [
        (
            "Standard\nleaching rule\n($f = 1.5$)",
            float(fao_ece_312["value"]),
            float(fao_ece_312["low"]),
            float(fao_ece_312["high"]),
            COL_EXTERNAL,
            f"{fao_ece_312['value']:.2f}",
        ),
        (
            "Fitted leaching\nfactor\n(fresh)",
            float(fit_fresh["ece_at_3_12_dS_m"]),
            float(fit_fresh["ece_at_3_12_range"][0]),
            float(fit_fresh["ece_at_3_12_range"][1]),
            COL_DEDUCED,
            f"{fit_fresh['ece_at_3_12_dS_m']:.2f}",
        ),
        (
            "Fitted leaching\nfactor\n(dry)",
            float(fit_dry["ece_at_3_12_dS_m"]),
            float(fit_dry["ece_at_3_12_range"][0]),
            float(fit_dry["ece_at_3_12_range"][1]),
            COL_TWIN,
            f"{fit_dry['ece_at_3_12_dS_m']:.2f}",
        ),
    ]

    # Measured site soil conductivity horizontal line + band
    ax_c.axhspan(
        measured_soil_ec - 0.12,
        measured_soil_ec + 0.12,
        color=COL_ACCENT,
        alpha=0.22,
        zorder=1,
    )
    ax_c.axhline(
        measured_soil_ec,
        color=INK,
        lw=1.0,
        ls="--",
        zorder=2,
    )
    ax_c.text(
        0.18,
        measured_soil_ec + 0.18,
        f"measured soil EC ({measured_soil_ec:.2f} dS m$^{{-1}}$)",
        ha="left",
        va="bottom",
        fontsize=5.4,
        color=INK,
        fontstyle="italic",
        zorder=5,
    )

    x_pos_c = np.arange(len(c_items))
    for xi, (lbl, val, lo, hi, col, txt) in zip(x_pos_c, c_items):
        ax_c.errorbar(
            xi,
            val,
            yerr=[[max(0.0, val - lo)], [max(0.0, hi - val)]],
            fmt="o",
            color=col,
            ecolor=col,
            markersize=5.0,
            elinewidth=1.4,
            capsize=3.0,
            zorder=4,
        )
        y_lbl = hi + 0.22 if xi == 0 else lo - 0.18
        va_lbl = "bottom" if xi == 0 else "top"
        ax_c.text(
            xi,
            y_lbl,
            f"{txt}\n[{lo:.2f}–{hi:.2f}]",
            ha="center",
            va=va_lbl,
            fontsize=5.4,
            color=col,
            fontweight="bold",
        )

    ax_c.set_xticks(x_pos_c)
    ax_c.set_xticklabels([item[0] for item in c_items], fontsize=5.8, color=INK)
    ax_c.set_xlim(-0.45, 2.5)
    ax_c.set_ylim(0.55, 7.65)
    ax_c.set_ylabel("Root-zone EC$_e$ (dS m$^{-1}$)", fontsize=6.5, color=INK)
    ax_c.tick_params(axis="both", labelsize=6.0)

    return _save_fig(fig, "v08_composition_main")


# ==============================================================================
# Figure B: paper/figures/v08_answers.{pdf,png}
# ==============================================================================
# 15 management intervention rows grouped as specified in the prompt
FIG_B_GROUPS: list[tuple[str, list[tuple[str, str, str, int]]]] = [
    (
        "Water timing",
        [
            ("four 3 mm showers", "rain-pattern", "rain-schedule", 2),
            ("daily drizzle (same 12 mm)", "rain-pattern", "rain-schedule", 3),
            ("daily irrigation (fresh)", "drip-or-flood", "irrigation-interval", 0),
            ("irrigation every 3 days", "drip-or-flood", "irrigation-interval", 1),
            ("irrigation every 10 days", "drip-or-flood", "irrigation-interval", 2),
        ],
    ),
    (
        "Water amount and quality",
        [
            ("cloud seeding +10 % rain", "cloud-seeding", "seeding-intensity", 2),
            ("daily brackish irrigation 8 dS m$^{-1}$", "salty-water", "daily-irrigation", 3),
        ],
    ),
    (
        "Surface cover and microclimate",
        [
            ("full organic mulch", "mulch-shade", "mulch-only", 2),
            ("full shade net", "mulch-shade", "shade-only", 2),
            ("full mulch + shade", "mulch-shade", "mulch-and-shade", 2),
            ("full solar-panel cover", "solar-park", "panel-cover", 2),
        ],
    ),
    (
        "Soil structure",
        [
            ("biochar 5 % w/w", "biochar", "biochar-dose", 3),
            ("50 vehicle passes", "compaction", "offroad-passes", 3),
        ],
    ),
    (
        "Warming",
        [
            ("soil warming +4 °C (fixed moisture)", "hotter", "warmer-soil-only", 3),
            ("climate warming +4 °C (rain unchanged)", "hotter", "warmer-climate", 3),
        ],
    ),
]

HEATMAP_COLS = [
    r"Rate constant $k_{\mathrm{H_2}}$",
    r"Retention exponent $b$",
    r"Air-entry potential $\psi_e$",
    r"Temperature coeff. $Q_{10}$",
    r"Activation ($\psi_d, \psi_w, a$)",
    r"Porosity & drainage",
    "External constants",
]


def _extract_level_st(
    q_doc: dict[str, Any],
    slug: str,
    series_id: str,
    level_idx: int,
    out_key: str = "h2_uptake_change_pct",
) -> tuple[dict[str, list[float]], dict[str, list[float]], str]:
    """Return (ST, group_ST, mode) for a specific level's H2 percentage-change output.

    Reads from reviews/v08/questions/<slug>.json when the level is the question's headline
    or if outcome-level ST is embedded; falls back to reviews/v08/uq/<slug>.json for
    non-headline levels.
    """
    s = next(sr for sr in q_doc["series"] if sr["id"] == series_id)
    lvl = s["levels"][level_idx]
    out = lvl["outcomes"][out_key]
    mode = str(out.get("mode", "variance"))
    if mode != "variance":
        return {}, {}, mode

    if "ST" in out and "group_ST" in out and out["ST"]:
        return out["ST"], out["group_ST"], mode

    u = q_doc.get("uncertainty", {})
    if (
        u.get("headline_series") == series_id
        and int(u.get("headline_level", -1)) == level_idx
        and u.get("headline_output") == out_key
        and u.get("ST")
    ):
        return u["ST"], u.get("group_ST", {}), mode

    uq_path = UQDIR / f"{slug}.json"
    if uq_path.exists():
        uq_doc = json.loads(uq_path.read_text(encoding="utf-8"))
        uq_s = next(sr for sr in uq_doc["series"] if sr["id"] == series_id)
        uq_l = next(lv for lv in uq_s["levels"] if int(lv["level"]) == level_idx)
        uq_o = uq_l["outputs"][out_key]
        return uq_o.get("ST", {}), uq_o.get("group_ST", {}), str(uq_o.get("mode", "variance"))

    return {}, {}, mode


def _compute_heatmap_row(
    st_map: dict[str, list[float]],
    group_st: dict[str, list[float]],
    mode: str,
) -> list[float]:
    if mode != "variance" or not st_map:
        return [np.nan] * 7

    def get_st(k: str) -> float:
        return float(st_map[k][0]) if k in st_map else 0.0

    k_val = get_st("k_h2_s")
    b_val = get_st("retention_b")
    psi_e_val = get_st("psi_entry_mpa")
    q10_val = get_st("q10")
    act_val = max(
        get_st("dry_threshold_mpa"),
        get_st("wet_activation_mpa"),
        get_st("activation_shape"),
    )
    soil_val = max(
        get_st("porosity"),
        get_st("drainage_fraction_day"),
        get_st("field_capacity_fraction"),
    )
    ext_val = float(group_st.get("verified_in_source", [0.0])[0]) + float(
        group_st.get("source_not_checked", [0.0])[0]
    )
    return [k_val, b_val, psi_e_val, q10_val, act_val, soil_val, ext_val]


def build_figure_answers(q_dir: Path = QDIR) -> tuple[Path, Path]:
    """Generate Figure B: paper/figures/v08_answers.{pdf,png}."""
    q_docs: dict[str, dict[str, Any]] = {}
    for _, items in FIG_B_GROUPS:
        for _, slug, _, _ in items:
            if slug not in q_docs:
                q_docs[slug] = json.loads((q_dir / f"{slug}.json").read_text(encoding="utf-8"))

    # Flatten rows while recording group boundaries
    rows_data: list[dict[str, Any]] = []
    group_spans: list[tuple[str, int, int]] = []
    idx_cursor = 0
    for grp_name, items in FIG_B_GROUPS:
        start_idx = idx_cursor
        for label, slug, series_id, level_idx in items:
            q_doc = q_docs[slug]
            s = next(sr for sr in q_doc["series"] if sr["id"] == series_id)
            lvl = s["levels"][level_idx]
            out = lvl["outcomes"]["h2_uptake_change_pct"]
            st_map, group_st, mode = _extract_level_st(q_doc, slug, series_id, level_idx)
            hm_vals = _compute_heatmap_row(st_map, group_st, mode)
            frac_neg = float(out["frac_neg"])
            frac_pos = float(out["frac_pos"])
            if frac_neg >= 0.95:
                cat = "decrease in ≥95 % of draws"
                col = COL_DECREASE
            elif frac_pos >= 0.95:
                cat = "increase in ≥95 % of draws"
                col = COL_INCREASE
            else:
                cat = "direction undetermined"
                col = COL_UNDETERMINED
            rows_data.append(
                {
                    "label": label,
                    "group": grp_name,
                    "median": float(out["median"]),
                    "q05": float(out["q05"]),
                    "q25": float(out["q25"]),
                    "q75": float(out["q75"]),
                    "q95": float(out["q95"]),
                    "frac_neg_pct": 100.0 * frac_neg,
                    "category": cat,
                    "color": col,
                    "heatmap": hm_vals,
                }
            )
            idx_cursor += 1
        group_spans.append((grp_name, start_idx, idx_cursor - 1))

    n_rows = len(rows_data)
    # Assign y-coordinates with a vertical gap between groups so group titles never collide with row labels
    group_gap = 0.72
    y_list: list[float] = []
    cur_y = 0.0
    for grp_idx, (_, i_start, i_end) in enumerate(group_spans):
        if grp_idx > 0:
            cur_y -= group_gap
        for _ in range(i_start, i_end + 1):
            y_list.append(cur_y)
            cur_y -= 1.0
    # Shift so bottom row is at y = 0.0
    min_y = min(y_list)
    y_coords = np.array([y - min_y for y in y_list], dtype=float)
    max_y = float(np.max(y_coords))

    fig = plt.figure(figsize=(FIG_WIDTH_IN, 5.55))
    gs = fig.add_gridspec(
        1,
        3,
        width_ratios=[1.55, 0.55, 1.35],
        left=0.272,
        right=0.985,
        top=0.855,
        bottom=0.16,
        wspace=0.10,
    )
    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1], sharey=ax_a)
    ax_c = fig.add_subplot(gs[0, 2], sharey=ax_a)

    # Panel labels a, b, c (bold lowercase)
    ax_a.text(
        -0.66,
        1.07,
        "a",
        transform=ax_a.transAxes,
        fontsize=9.5,
        fontweight="bold",
        color=INK,
        ha="left",
        va="top",
    )
    ax_b.text(
        -0.08,
        1.07,
        "b",
        transform=ax_b.transAxes,
        fontsize=9.5,
        fontweight="bold",
        color=INK,
        ha="left",
        va="top",
    )
    ax_c.text(
        -0.05,
        1.07,
        "c",
        transform=ax_c.transAxes,
        fontsize=9.5,
        fontweight="bold",
        color=INK,
        ha="left",
        va="top",
    )

    # ---------------- Panel a: Forest plot ----------------
    ax_a.axvline(0.0, color=INK, lw=0.85, ls="--", zorder=2)
    for yi, r in zip(y_coords, rows_data):
        c = r["color"]
        # 5-95% interval (thin line)
        ax_a.plot(
            [r["q05"], r["q95"]],
            [yi, yi],
            color=c,
            lw=1.1,
            solid_capstyle="round",
            zorder=3,
        )
        # Interquartile range q25-q75 (thick line)
        ax_a.plot(
            [r["q25"], r["q75"]],
            [yi, yi],
            color=c,
            lw=3.0,
            solid_capstyle="butt",
            zorder=4,
        )
        # Median point
        ax_a.plot(
            r["median"],
            yi,
            marker="o",
            color=c,
            markersize=4.5,
            markeredgecolor="white",
            markeredgewidth=0.4,
            zorder=5,
        )

    ax_a.set_xscale("symlog", linthresh=10.0, linscale=0.55)
    ax_a.set_xlim(-105.0, 220.0)
    ax_a.set_xticks([-100, -30, -10, 0, 10, 40, 200])
    ax_a.set_xticklabels(["$-$100", "$-$30", "$-$10", "0", "+10", "+40", "+200"], fontsize=5.6)
    ax_a.set_yticks(y_coords)
    ax_a.set_yticklabels([r["label"] for r in rows_data], fontsize=6.0, color=INK)
    ax_a.set_ylim(-0.65, max_y + 0.85)
    ax_a.set_xlabel("14-day H$_2$ uptake change vs 12 mm storm (%)", fontsize=6.4, color=INK)
    ax_a.set_title("Effect & 5–95 % / IQR intervals", fontsize=6.8, fontweight="bold", color=INK, loc="left", pad=6)

    # Group separators and labels in the vertical gaps above each group
    for grp_idx, (grp_name, i_start, i_end) in enumerate(group_spans):
        y_top_grp = y_coords[i_start]
        if grp_idx > 0:
            sep_y = y_top_grp + 0.82
            ax_a.axhline(sep_y, color=GRID_GREY, lw=0.8, ls="-", zorder=1)
            ax_b.axhline(sep_y, color=GRID_GREY, lw=0.8, ls="-", zorder=1)
            ax_c.axhline(sep_y, color=GRID_GREY, lw=0.8, ls="-", zorder=1)
        ax_a.text(
            -0.03,
            y_top_grp + 0.48,
            grp_name,
            transform=ax_a.get_yaxis_transform(),
            ha="right",
            va="center",
            fontsize=5.8,
            fontweight="bold",
            fontstyle="italic",
            color=MUTED,
        )

    # Direction legend centered at top of figure
    legend_handles = [
        Line2D([0], [0], color=COL_DECREASE, marker="o", lw=1.8, ms=4, label="decrease in ≥95 % of draws"),
        Line2D([0], [0], color=COL_INCREASE, marker="o", lw=1.8, ms=4, label="increase in ≥95 % of draws"),
        Line2D([0], [0], color=COL_UNDETERMINED, marker="o", lw=1.8, ms=4, label="direction undetermined"),
    ]
    fig.legend(
        handles=legend_handles,
        loc="upper center",
        bbox_to_anchor=(0.54, 0.985),
        ncol=3,
        frameon=False,
        fontsize=6.0,
        handlelength=1.5,
        columnspacing=1.4,
    )

    # ---------------- Panel b: Share of draws with a decrease ----------------
    ax_b.axvline(95.0, color=GRID_GREY, lw=0.7, ls=":", zorder=1)
    ax_b.axvline(5.0, color=GRID_GREY, lw=0.7, ls=":", zorder=1)
    for yi, r in zip(y_coords, rows_data):
        fn = r["frac_neg_pct"]
        ax_b.barh(
            yi,
            fn,
            height=0.58,
            color=r["color"],
            alpha=0.82,
            edgecolor="none",
            zorder=3,
        )
        txt_str = f"{fn:.0f}%" if abs(fn - round(fn)) < 0.05 else f"{fn:.1f}%"
        x_txt = fn + 3.0 if fn < 68.0 else fn - 3.0
        ha_txt = "left" if fn < 68.0 else "right"
        col_txt = INK if fn < 68.0 else "white"
        ax_b.text(
            x_txt,
            yi,
            txt_str,
            va="center",
            ha=ha_txt,
            fontsize=5.4,
            fontweight="bold" if fn >= 68.0 else "normal",
            color=col_txt,
            zorder=5,
        )

    ax_b.set_xlim(0.0, 108.0)
    ax_b.set_xticks([0, 50, 100])
    ax_b.set_xticklabels(["0", "50", "100"], fontsize=5.8)
    ax_b.tick_params(axis="y", left=False, labelleft=False)
    ax_b.set_xlabel("Draws < 0 (%)", fontsize=6.4, color=INK)
    ax_b.set_title("Decrease share", fontsize=6.8, fontweight="bold", color=INK, loc="left", pad=6)

    # ---------------- Panel c: Heatmap of total-order Sobol' indices ----------------
    hm_matrix = np.array([r["heatmap"] for r in rows_data], dtype=float)
    cmap = LinearSegmentedColormap.from_list(
        "eq_seq",
        ["#F7F9FB", "#B7D4DF", "#4E8FA6", "#26323A"],
    ).with_extremes(bad="#F0F2F4")

    n_cols = len(HEATMAP_COLS)
    # Plot cells aligned with y_coords
    for row_idx, yi in enumerate(y_coords):
        for col_idx in range(n_cols):
            val = hm_matrix[row_idx, col_idx]
            face = "#F0F2F4" if np.isnan(val) else cmap(np.clip(val / 0.85, 0.0, 1.0))
            rect = plt.Rectangle(
                (col_idx - 0.5, yi - 0.5),
                1.0,
                1.0,
                facecolor=face,
                edgecolor="white",
                linewidth=0.6,
                zorder=2,
            )
            ax_c.add_patch(rect)
            if not np.isnan(val) and val >= 0.1:
                txt_col = "white" if val >= 0.45 else INK
                ax_c.text(
                    col_idx,
                    yi,
                    f"{val:.2f}",
                    ha="center",
                    va="center",
                    fontsize=5.3,
                    color=txt_col,
                    fontweight="bold" if val >= 0.35 else "normal",
                    zorder=4,
                )

    ax_c.set_xlim(-0.5, n_cols - 0.5)
    ax_c.set_xticks(np.arange(n_cols))
    ax_c.set_xticklabels(HEATMAP_COLS, rotation=32, ha="right", fontsize=5.6, color=INK)
    ax_c.tick_params(axis="y", left=False, labelleft=False)
    ax_c.set_title("Total-order Sobol' sensitivity $S_T$", fontsize=6.8, fontweight="bold", color=INK, loc="left", pad=6)

    return _save_fig(fig, "v08_answers")


# ==============================================================================
# Main-text Table: paper/tables/v08_answers.tex
# ==============================================================================
# 14 questions ordered in the 6 narrative groups of the restructured paper:
#   1. Water timing and amount (rain-pattern, drip-or-flood, cloud-seeding, salty-water)
#   2. Surface and microclimate (mulch-shade, solar-park, hotter)
#   3. Soil structure (compaction, biochar)
#   4. Atmosphere (air-hydrogen)
#   5. Community and carbon (missing-co-genes, adding-microbes, carbon-fixation)
#   6. Crop (alfalfa-bacterium)
TABLE_SPECS: list[tuple[str, str, str, str, bool]] = [
    # (slug, plain_question, readable_intervention, readable_endpoint, end_of_group)
    (
        "rain-pattern",
        r"Does splitting rain into showers change H$_2$ uptake?",
        "Four 3 mm showers",
        r"14-d H$_2$ uptake change (\%)",
        False,
    ),
    (
        "drip-or-flood",
        r"Does drip or flood irrigation timing change H$_2$ uptake?",
        "Irrigation every 10 days",
        r"14-d H$_2$ uptake change (\%)",
        False,
    ),
    (
        "cloud-seeding",
        r"How much does cloud seeding change H$_2$ uptake?",
        r"Cloud seeding (+10\% rain)",
        r"14-d H$_2$ uptake change (\%)",
        False,
    ),
    (
        "salty-water",
        r"How does brackish irrigation change H$_2$ uptake?",
        r"Daily irrigation (8 dS m$^{-1}$)",
        r"14-d H$_2$ uptake change (\%)",
        True,
    ),
    (
        "mulch-shade",
        r"How do surface mulch and shade nets change H$_2$ uptake?",
        "Full organic mulch",
        r"14-d H$_2$ uptake change (\%)",
        False,
    ),
    (
        "solar-park",
        r"How do solar panels change desert soil H$_2$ uptake?",
        "Full solar-panel cover",
        r"14-d H$_2$ uptake change (\%)",
        False,
    ),
    (
        "hotter",
        r"How does climate warming change H$_2$ uptake?",
        r"Warming +4 $^\circ$C (rain unchanged)",
        r"14-d H$_2$ uptake change (\%)",
        True,
    ),
    (
        "compaction",
        r"How does off-road vehicle traffic change H$_2$ uptake?",
        "50 vehicle passes",
        r"14-d H$_2$ uptake change (\%)",
        False,
    ),
    (
        "biochar",
        r"Does adding biochar to desert sand change H$_2$ uptake?",
        r"Biochar 5\% w/w",
        r"14-d H$_2$ uptake change (\%)",
        True,
    ),
    (
        "air-hydrogen",
        r"How does higher atmospheric H$_2$ change soil H$_2$ uptake?",
        r"2$\times$ background H$_2$",
        r"14-d H$_2$ uptake change (\%)",
        True,
    ),
    (
        "missing-co-genes",
        "What happens to CO uptake when CO-oxidising genes are absent?",
        "Closed- vs open-gate communities",
        r"14-d CO uptake ($\mu$mol m$^{-2}$)",
        False,
    ),
    (
        "adding-microbes",
        "Does adding CO-oxidising microbes increase soil CO uptake?",
        r"+1\% CO-carrier share",
        r"14-d CO uptake change (\%)",
        False,
    ),
    (
        "carbon-fixation",
        r"How much carbon can H$_2$-oxidising autotrophs fix per year?",
        "Single 12 mm storm (annualised)",
        r"C-fixation ceiling (g C m$^{-2}$ yr$^{-1}$)",
        True,
    ),
    (
        "alfalfa-bacterium",
        "How much alfalfa yield does SA187 recover under brackish water?",
        r"7.81 dS m$^{-1}$ + SA187 (fitted leaching)",
        "Relative yield (fraction)",
        False,
    ),
]


def _format_signed_pct(val: float) -> str:
    s = fmt_pct_value(val)
    if s.startswith("-"):
        return "$-$" + s[1:]
    if float(val) > 0:
        return "+" + s
    return s


def _format_median_and_direction(slug: str, q_doc: dict[str, Any]) -> tuple[str, str]:
    u = q_doc["uncertainty"]
    med = float(u["median"])
    q05 = float(u["q05"])
    q95 = float(u["q95"])
    fn = float(u["frac_neg"])
    fp = float(u["frac_pos"])
    mode = str(u.get("mode", "variance"))

    if slug == "missing-co-genes":
        levels = q_doc["series"][0]["levels"]
        n_closed = sum(
            1 for lv in levels if lv["outcomes"]["co_uptake_umol_m2"].get("mode") == "deduced"
        )
        med_s = _fmt_sig_latex(med, 3)
        lo_s = _fmt_sig_latex(q05, 3)
        hi_s = _fmt_sig_latex(q95, 3)
        med_col = f"0 [0, 0] ({n_closed} closed); {med_s} [{lo_s}, {hi_s}] (open)"
        dir_txt = "all draws" if fp == 1.0 else f"{fmt_frac_pct(fp)}\\% of draws"
        dir_col = f"exact zero by deduction ({n_closed} samples); positive, {dir_txt} (open)"
        return med_col, dir_col

    if slug == "air-hydrogen" or (mode == "constant" and abs(q95 - q05) < 1e-9):
        s_val = _format_signed_pct(med)
        med_col = f"{s_val} [{s_val}, {s_val}]"
        dir_txt = "all draws" if fp == 1.0 else f"{fmt_frac_pct(fp)}\\% of draws"
        return med_col, f"constant ({s_val}\\%, {dir_txt})"

    if slug == "carbon-fixation":
        med_s = _fmt_sig_latex(med, 2)
        lo_s = _fmt_sig_latex(q05, 2)
        hi_s = _fmt_sig_latex(q95, 2)
        med_col = f"{med_s} [{lo_s}, {hi_s}]"
        dir_txt = "all draws" if fp == 1.0 else f"{fmt_frac_pct(fp)}\\% of draws"
        return med_col, f"positive bound, {dir_txt}"

    if slug == "alfalfa-bacterium":
        med_s = _fmt_sig_latex(med, 3)
        lo_s = _fmt_sig_latex(q05, 3)
        hi_s = _fmt_sig_latex(q95, 3)
        med_col = f"{med_s} [{lo_s}, {hi_s}]"
        dir_txt = "all draws" if fp == 1.0 else f"{fmt_frac_pct(fp)}\\% of draws"
        return med_col, f"positive yield, {dir_txt}"

    # Percentage-change questions
    med_s = _format_signed_pct(med)
    lo_s = _format_signed_pct(q05)
    hi_s = _format_signed_pct(q95)
    med_col = f"{med_s} [{lo_s}, {hi_s}]"

    if fn >= 0.95:
        if fn == 1.0:
            dir_col = "decrease, all draws"
        else:
            dir_col = f"decrease, {fmt_frac_pct(fn)}\\% of draws"
    elif fp >= 0.95:
        if fp == 1.0:
            dir_col = "increase, all draws"
        else:
            dir_col = f"increase, {fmt_frac_pct(fp)}\\% of draws"
    else:
        if fn >= fp:
            dir_col = f"undetermined ({fmt_frac_pct(fn)}\\% decrease)"
        else:
            dir_col = f"undetermined ({fmt_frac_pct(fp)}\\% increase)"
    return med_col, dir_col


def _format_main_drivers(q_doc: dict[str, Any]) -> str:
    st_map = q_doc["uncertainty"].get("ST") or {}
    ranked = sorted(((k, float(v[0])) for k, v in st_map.items()), key=lambda kv: (-kv[1], kv[0]))
    top2 = ranked[:2]
    parts = [f"{readable_name(name, latex=True)} ({fmt_st(val)})" for name, val in top2]
    return ", ".join(parts)


def write_answers_table(
    questions_dir: Path = QDIR,
    out_file: Path = TAB_DIR / "v08_answers.tex",
) -> Path:
    """Write paper/tables/v08_answers.tex (14 data rows, booktabs, \\small, tabularx)."""
    rows_tex: list[str] = []
    for slug, q_plain, intv_label, endpoint_label, end_of_group in TABLE_SPECS:
        q_doc = json.loads((questions_dir / f"{slug}.json").read_text(encoding="utf-8"))
        med_col, dir_col = _format_median_and_direction(slug, q_doc)
        drivers_col = _format_main_drivers(q_doc)
        row_line = (
            f"{q_plain} & {intv_label} & {endpoint_label} & "
            f"{med_col} & {dir_col} & {drivers_col} \\\\"
        )
        if end_of_group:
            row_line += "\n\\addlinespace[3pt]"
        rows_tex.append(row_line)

    table_str = (
        "% Generated by reproducibility/v08/paper_figures.py; do not edit by hand.\n"
        "{\\small\n"
        "\\begin{tabularx}{\\textwidth}{@{} >{\\raggedright\\arraybackslash}p{3.0cm} "
        ">{\\raggedright\\arraybackslash}p{2.2cm} "
        ">{\\raggedright\\arraybackslash}p{2.1cm} "
        ">{\\raggedright\\arraybackslash}p{2.1cm} "
        ">{\\raggedright\\arraybackslash}p{2.3cm} "
        ">{\\raggedright\\arraybackslash}X @{} }\n"
        "\\toprule\n"
        "Question & Intervention & Endpoint & Median [5--95\\,\\%] & Direction & Main drivers ($S_T$) \\\\\n"
        "\\midrule\n"
        + "\n".join(rows_tex)
        + "\n\\bottomrule\n"
        "\\end{tabularx}\n"
        "}\n"
    )
    out_file.parent.mkdir(parents=True, exist_ok=True)
    out_file.write_text(table_str, encoding="utf-8")
    return out_file


def main() -> None:
    build_figure_composition_main()
    build_figure_answers()
    write_answers_table()


if __name__ == "__main__":
    main()
