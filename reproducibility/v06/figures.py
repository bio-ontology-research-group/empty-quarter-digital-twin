"""v0.6 publication figures: stratified graph, community-coupled scenarios and the improved primary-twin figure."""
from pathlib import Path
import json
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Patch, Rectangle

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
from eq_biodt.v03 import theory                                    # noqa: E402
sys.path.insert(0, str(ROOT / "reproducibility/v06"))
import v06_graph                                                   # noqa: E402

OUT = ROOT / "reviews/v06"
FIG = ROOT / "paper/figures"
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8, "axes.spines.top": False, "axes.spines.right": False,
                     "pdf.fonttype": 42, "savefig.dpi": 200, "mathtext.fontset": "dejavusans"})
INK = "#26323A"
ZONE_COLOUR = {"ontological": "#8E7CB3", "evidence": "#4E8FA6", "deduced": "#4F9A82", "exogenous_manipulable": "#D9A441",
               "exogenous_fixed": "#BDB6A6", "constitutive": "#6C8EBF", "causal": "#D2693C"}
BAND = [("ontological", "evidence"), ("deduced",), ("exogenous_manipulable",), ("constitutive",), ("causal",)]
BAND_TITLE = {0: "Ontology $K$ and\nevidence $E$", 1: "Deduced\n$D_K(E)$", 2: "Exogenous\ninputs",
              3: "Constitutive\nmaps $\\Phi$", 4: "Causal\nmechanisms"}
SYMBOL = {"K": "K", "E": "E", "D_K": "D_K(E)", "q_h2": "q_{\\mathrm{H}_2}", "q_co": "q_{\\mathrm{CO}}",
          "gamma_h2": "\\gamma_{\\mathrm{H}_2}", "gamma_co": "\\gamma_{\\mathrm{CO}}", "I": "I_t", "c": "c_t", "Evap": "\\mathrm{Ev}_t",
          "T": "T_t", "phi": "\\phi", "psi_e": "\\psi_e", "b": "b", "L": "L", "C_a": "C_a", "D0": "D_0", "k_ref": "k_{\\mathrm{ref}}",
          "Q10": "Q_{10}", "act_par": "\\psi_d,\\psi_w,a", "drain_par": "f,\\,d", "osm_par": "\\varphi,R", "O2": "S_{\\mathrm{O}_2}",
          "capture": "\\eta", "W": "W_t", "N": "N_t", "W_max": "W_{\\max}", "theta": "\\theta", "psi_m": "\\psi_m", "psi_o": "\\psi_o",
          "psi": "\\psi", "D_h2": "D_{\\mathrm{H}_2}", "D_co": "D_{\\mathrm{CO}}", "A": "A_t", "k_h2": "k_{\\mathrm{H}_2}",
          "k_co": "k_{\\mathrm{CO}}", "J_h2": "J_{\\mathrm{H}_2}", "J_co": "J_{\\mathrm{CO}}", "P": "P", "P_cap": "P_{\\mathrm{cap}}"}


def save(fig, name):
    FIG.mkdir(exist_ok=True)
    OUT.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG / f"{name}.pdf", bbox_inches="tight", metadata={"CreationDate": None, "ModDate": None})
    fig.savefig(OUT / f"{name}.png", bbox_inches="tight", dpi=170)
    plt.close(fig)


# ------------------------------------------------------------------ Fig: stratified graph
ROWS = 5
ROW_OF = {"ontological": 0, "evidence": 0, "deduced": 1, "exogenous_manipulable": 2, "constitutive": 3, "causal": 4}
# Hand-tuned x positions (fixed constants are omitted from the figure; see the caption).
XPOS = {"K": 36, "E": 68,
        "gamma_h2": 20, "q_h2": 36, "D_K": 52, "q_co": 68, "gamma_co": 84,
        "Evap": 8, "I": 19, "c": 30, "T": 41, "psi_e": 53, "b": 64, "phi": 76,
        "W_max": 10, "theta": 20.5, "psi_o": 31, "psi": 42, "psi_m": 53, "P_cap": 64, "P": 75, "D_h2": 86, "D_co": 98,
        "W": 13, "N": 29, "A": 42, "k_h2": 58, "k_co": 71, "J_h2": 86, "J_co": 98}


def layout(g):
    return {n["id"]: (XPOS[n["id"]], ROW_OF[n["zone"]]) for n in g["nodes"] if n["zone"] != "exogenous_fixed"}


def draw_graph():
    from matplotlib.patches import PathPatch
    from matplotlib.path import Path as MPath
    g = v06_graph.graph(theory)
    POS = layout(g)
    nodes = {n["id"]: n for n in g["nodes"] if n["id"] in POS}
    fig, ax = plt.subplots(figsize=(7.0, 6.9))
    ax.set_xlim(-15, 103)
    ax.set_ylim(ROWS - 0.42, -0.55)
    ax.axis("off")
    for r in range(ROWS):
        ax.add_patch(Rectangle((-15, r - .5), 118, 1.0, facecolor=ZONE_COLOUR[BAND[r][0]], alpha=.11, lw=0, zorder=0))
        ax.text(-14, r, BAND_TITLE[r], ha="left", va="center", fontsize=8.4, fontweight="bold", color=INK, zorder=5)
    style = {"mechanism": dict(ls="-", lw=1.25, color="#A8431F", alpha=.8),
             "definition": dict(ls=(0, (4, 2)), lw=1.05, color="#34567F", alpha=.75),
             "deduction": dict(ls=(0, (1, 1.6)), lw=1.4, color="#2F7A63", alpha=.95)}
    hw, hh = 3.7, .22

    def anchor(src, dst):
        (x1, y1), (x2, y2) = POS[src], POS[dst]
        if y1 == y2:
            sgn = 1 if x2 > x1 else -1
            return (x1 + sgn * hw, y1), (x2 - sgn * hw, y2)
        sgn = 1 if y2 > y1 else -1
        # attach along the box edge, shifted towards the partner to spread the fan
        sx = x1 + np.clip((x2 - x1) * .08, -hw * .6, hw * .6)
        dx = x2 + np.clip((x1 - x2) * .08, -hw * .6, hw * .6)
        return (sx, y1 + sgn * hh), (dx, y2 - sgn * hh)

    gaps = {}
    for r in range(ROWS):
        xs = sorted(POS[i][0] for i in POS if POS[i][1] == r)
        gaps[r] = [xs[0] - 5] + [(a + b) / 2 for a, b in zip(xs, xs[1:])] + [xs[-1] + 5]
    used = {}
    for e in g["edges"]:
        if e["source"] not in POS:
            continue
        if e["kind"] == "lag":
            x, y = POS[e["target"]]
            path = MPath([(x - 1.6, y + hh), (x - 3.4, y + hh + .5), (x + 3.4, y + hh + .5), (x + 1.6, y + hh + .02)],
                         [MPath.MOVETO, MPath.CURVE4, MPath.CURVE4, MPath.CURVE4])
            ax.add_patch(FancyArrowPatch(path=path, arrowstyle="-|>", mutation_scale=7, lw=1.3, color="#5E3A87", zorder=4,
                                         shrinkA=0, shrinkB=0))
            continue
        s_, d_ = anchor(e["source"], e["target"])
        st = dict(style[e["kind"]])
        if nodes[e["source"]]["zone"] == "exogenous_fixed":
            st.update(alpha=.55, lw=.85)
        (x1, y1), (x2, y2) = POS[e["source"]], POS[e["target"]]
        pts = [s_]
        if abs(y2 - y1) > 1:
            lo, hi = (y1, y2) if y2 > y1 else (y2, y1)
            rows_between = list(range(int(lo) + 1, int(hi))) if y2 > y1 else list(range(int(hi) - 1, int(lo), -1))
            for r in rows_between:
                t = (r - y1) / (y2 - y1)
                gx = min(gaps[r], key=lambda g: abs(g - (s_[0] + t * (d_[0] - s_[0]))))
                used[(r, gx)] = used.get((r, gx), 0) + 1
                gx += .55 * ((used[(r, gx)] + 1) // 2) * (1 if used[(r, gx)] % 2 else -1)
                sgn = 1 if y2 > y1 else -1
                pts += [(gx, r - sgn * .5 + sgn * .02), (gx, r + sgn * .5 - sgn * .02)]
        pts.append(d_)
        path = MPath(pts, [MPath.MOVETO] + [MPath.LINETO] * (len(pts) - 1))
        ax.add_patch(FancyArrowPatch(path=path, arrowstyle="-|>", mutation_scale=6, zorder=1, shrinkA=0, shrinkB=0, **st))
    for i, n in nodes.items():
        x, y = POS[i]
        manip = n["manipulable"]
        fixed = n["zone"] == "exogenous_fixed"
        ax.add_patch(FancyBboxPatch((x - hw, y - hh), 2 * hw, 2 * hh, boxstyle="round,pad=0,rounding_size=.12",
                                    facecolor="#EEEAE0" if fixed else ZONE_COLOUR[n["zone"]],
                                    edgecolor=INK if manip else "#55606A", linewidth=2.2 if manip else .7, zorder=3))
        dark = n["zone"] in ("ontological", "evidence", "deduced", "constitutive", "causal")
        ax.text(x, y, f"${SYMBOL[i]}$", ha="center", va="center", fontsize=10,
                color="white" if dark else INK, fontweight="bold" if manip else "normal", zorder=6)
    handles = [Line2D([], [], color="#A8431F", lw=1.5, label="mechanism (do replaces it)"),
               Line2D([], [], color="#34567F", lw=1.3, ls=(0, (4, 2)), label="definition (invariant)"),
               Line2D([], [], color="#2F7A63", lw=1.6, ls=(0, (1, 1.6)), label="deduction from $K$, $E$"),
               Line2D([], [], color="#5E3A87", lw=1.4, label="lag $t\\to t{+}1$ (curved)"),
               Patch(facecolor="white", edgecolor=INK, linewidth=2.0, label="physical handle (admissible do target)")]
    ax.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, -0.005), ncol=3, frameon=False, fontsize=8,
              handlelength=2.4, columnspacing=1.4)
    save(fig, "v06_graph")


# ------------------------------------------------------------------ Fig: scenarios
COL = ["#196B75", "#D27A36", "#72618C", "#677B40"]
COMP = [("Surface", "Surface"), ("Deep", "Subsurface"), ("PlantRoot", "Root-adjacent")]


def jitter(n, width=.34):
    return ((np.arange(n) * 0.6180339887) % 1 - .5) * 2 * width


def draw_scenarios():
    tv = pd.read_csv(OUT / "treatment_versions.tsv", sep="\t")
    df = pd.read_csv(OUT / "paired_regimes_by_sample.tsv", sep="\t")
    summary = json.loads((OUT / "paired_regimes_summary.json").read_text())
    ref = summary["reference_uptake_mol_m2"]
    fig, axes = plt.subplots(1, 3, figsize=(7.2, 3.1), gridspec_kw={"width_ratios": [1.35, 1, 1.1]}, layout="constrained")
    # A: members sorted by single-pulse uptake
    ax = axes[0]
    t = tv.sort_values("uptake_single_pulse").reset_index(drop=True)
    x = np.arange(1, len(t) + 1)
    series = [("single_pulse", "One 12 mm pulse", COL[0], "o"), ("four_pulses", "Four 3 mm pulses", COL[1], "s"),
              ("daily_pulses", "Daily 12/14 mm", COL[2], "^")]
    rev = t.index[t.uptake_four_pulses > t.uptake_single_pulse]
    for i in rev:
        ax.axvspan(i + .5, i + 1.5, color="#F3D9C4", lw=0, zorder=0)
    ax.axhline(0, color=COL[3], lw=1.6, label="Intact control (0)", zorder=2)
    for key, lab, c, m in series:
        ax.plot(x, t[f"uptake_{key}"] * 1e6, marker=m, ms=3.4, lw=.8, color=c, label=lab, zorder=3)
    ax.set_yscale("symlog", linthresh=1, linscale=.35)
    ax.set_yticks([0, 1, 10, 100, 1000, 10000])
    ax.set_yticklabels(["0", "1", "10", "$10^2$", "$10^3$", "$10^4$"])
    ax.set(xlabel="Parameter member (sorted by one-pulse uptake)", ylabel="14-day H$_2$ uptake ($\\mu$mol m$^{-2}$)", xlim=(.4, 24.6))
    ax.legend(frameon=False, fontsize=6.6, loc="upper left", handlelength=1.4, borderaxespad=.1)
    if len(rev):
        ax.text(rev.max() + 1.05, ax.get_ylim()[0] if False else 0.55, "four pulses\nexceed one", fontsize=6.4, ha="right", va="bottom", color="#9B5A25")
    ax.set_title("A  Treatment versions, 24 members", loc="left", fontweight="bold", fontsize=8)
    # B: single-pulse H2 effect by compartment, proportional coupling
    ax = axes[1]
    for k, (comp, lab) in enumerate(COMP):
        v = df.loc[df.compartment.eq(comp), "h2eff_proportional_single_pulse"].sort_values().to_numpy() * 1e6
        ax.scatter(k + jitter(len(v)), v, s=3.5, color=COL[k], alpha=.45, lw=0, rasterized=True)
        q1, med, q3 = np.quantile(v, [.25, .5, .75])
        ax.plot([k - .36, k + .36], [med, med], color=INK, lw=1.6, zorder=4)
        ax.plot([k, k], [q1, q3], color=INK, lw=3.4, alpha=.55, solid_capstyle="butt", zorder=3)
    ax.axhline(ref["single_pulse"]["h2"] * 1e6, color="#8B8B8B", lw=.9, ls="--")
    ax.text(2.45, ref["single_pulse"]["h2"] * 1e6 + 40, "reference\ncommunity", fontsize=6.2, color="#666", ha="right", va="bottom")
    ax.set(xticks=range(3), xticklabels=[l.replace("-", "-\n") if "-" in l else l for _, l in COMP], xlim=(-.55, 2.55),
           ylabel="One-pulse H$_2$ effect ($\\mu$mol m$^{-2}$)")
    ax.tick_params(axis="x", labelsize=7)
    ax.set_title("B  Proportional coupling", loc="left", fontweight="bold", fontsize=8)
    # C: CO effect vs carrier share
    ax = axes[2]
    closed = df.q_co <= 0
    ax.scatter(df.q_co[~closed] * 100, df.coeff_proportional_single_pulse[~closed] * 1e6, s=4, color="#9AA7A9", lw=0, rasterized=True, label="Proportional")
    ax.scatter(df.q_co[~closed] * 100, df.coeff_gate_single_pulse[~closed] * 1e6, s=4, color=COL[0], lw=0, rasterized=True, label="Gate (open)")
    ax.scatter(df.q_co[closed] * 100, df.coeff_gate_single_pulse[closed] * 1e6, s=26, color="#C0392B", edgecolor="white", lw=.5, zorder=5,
               label=f"Gate closed (n = {int(closed.sum())}): exactly 0")
    ax.scatter(df.q_co[closed] * 100, df.coeff_proportional_single_pulse[closed] * 1e6, s=26, marker="x", color="#C0392B", lw=1, zorder=5)
    ax.set_xscale("symlog", linthresh=.1, linscale=.35)
    ax.set_xticks([0, .1, 1, 10]); ax.set_xticklabels(["0", "0.1", "1", "10"])
    ax.set(xlabel="CO-carrier share of mapped reads (%)", ylabel="One-pulse CO effect ($\\mu$mol m$^{-2}$)", xlim=(-.02, 60))
    ax.legend(frameon=False, fontsize=6.4, loc="upper left", handletextpad=.3, borderaxespad=.1)
    ax.set_title("C  CO under deductive gate", loc="left", fontweight="bold", fontsize=8)
    save(fig, "v06_scenarios")


# ------------------------------------------------------------------ Fig 1 (improved)
def _polys(geom):
    if geom["type"] == "Polygon":
        return [geom["coordinates"]]
    return geom["coordinates"]


def draw_coast(ax, names=None, lw=.5):
    fc = json.loads((ROOT / "reproducibility/v06/data/arabia_natural_earth_50m.json").read_text())["features"]
    for f in fc:
        if f["name"] in ("United States of America", "France"):
            continue
        for poly in _polys(f):
            xy = np.array(poly[0])
            ax.fill(xy[:, 0], xy[:, 1], facecolor="#F4EEDF", edgecolor="#A99F8A", lw=lw, zorder=1)


def box(ax, x, y, w, h, text, fc, ec="#9AA6A8", fs=7.2, bold=False, lw=.9, tc=INK):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.004,rounding_size=.018", facecolor=fc, edgecolor=ec, lw=lw,
                                transform=ax.transAxes, zorder=3))
    ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=fs, transform=ax.transAxes, zorder=4,
            fontweight="bold" if bold else "normal", color=tc, linespacing=1.15)


def arrow(ax, a, b, color="#55606A", style="-|>", lw=1.1, rad=0., ls="-"):
    ax.add_patch(FancyArrowPatch(a, b, transform=ax.transAxes, arrowstyle=style, mutation_scale=8, lw=lw, color=color, zorder=5,
                                 connectionstyle=f"arc3,rad={rad}", shrinkA=0, shrinkB=0, ls=ls))


def draw_primary():
    sites = pd.read_csv(ROOT / "reviews/v05/primary_sites.tsv", sep="\t")
    panel = pd.read_csv(ROOT / "data/v03/sample_evidence.tsv", sep="\t")
    fig = plt.figure(figsize=(9, 7.6), layout="constrained")
    gs = fig.add_gridspec(2, 2, height_ratios=[1.12, 1], width_ratios=[.82, 1.18])
    ax = fig.add_subplot(gs[0, :])
    lon0, lon1, lat0, lat1 = 43.4, 56.6, 17.7, 22.3
    ax.set_facecolor("#DCE8EC")
    draw_coast(ax)
    ax.grid(color="white", lw=.7, alpha=.9, zorder=0.5)
    k = 4.6
    pts = ax.scatter(sites.longitude, sites.latitude, c=sites.campaigns, s=k * sites.samples, cmap="viridis", vmin=1, vmax=5,
                     edgecolor="white", linewidth=.5, zorder=4, alpha=.95)
    aux = sites[sites.site.gt(60)]
    ax.scatter(aux.longitude, aux.latitude, s=k * aux.samples + 85, facecolor="none", edgecolor="#9B452F", linewidth=1.3, zorder=5)
    ax.set_xlim(lon0, lon1)
    ax.set_ylim(lat0, lat1)
    ax.set_aspect(1 / np.cos(np.radians(20.2)))
    ax.set_xlabel("Longitude (°E)")
    ax.set_ylabel("Latitude (°N)")
    for name, (x, y) in {"SAUDI ARABIA": (49.2, 22.05), "OMAN": (55.3, 18.4), "YEMEN": (46.2, 17.95)}.items():
        ax.text(x, y, name, fontsize=7, color="#8A826F", ha="center", va="center", style="italic", zorder=2)
    # scale bar, 100 km at 20 N
    dlon = 100 / (111.32 * np.cos(np.radians(20.2)))
    x0, y0 = 54.2, 17.95
    ax.plot([x0, x0 + dlon], [y0, y0], color=INK, lw=2.2, solid_capstyle="butt", zorder=6)
    ax.text(x0 + dlon / 2, y0 + .12, "100 km", ha="center", va="bottom", fontsize=7, zorder=6)
    # legends: size and ring
    for n, (lx, ly) in zip((1, 10, 25), ((44.2, 17.95), (44.9, 17.95), (45.7, 17.95))):
        ax.scatter([lx], [ly], s=k * n, color="#B8C2C4", edgecolor="#55606A", lw=.5, zorder=6)
        ax.text(lx, ly + .28 if n == 1 else ly + .38, str(n), ha="center", fontsize=6.6, zorder=6)
    ax.text(44.0, 18.62, "samples per site (point area)", fontsize=6.8, ha="left", va="bottom", zorder=6)
    ax.scatter([47.3], [17.95], s=k * 3 + 85, facecolor="none", edgecolor="#9B452F", lw=1.3, zorder=6)
    ax.text(47.7, 17.95, "auxiliary sites", fontsize=6.8, va="center", zorder=6)
    ax.set_title("A  Primary sampling footprint, Rub’ al-Khali", loc="left", fontweight="bold", fontsize=9, pad=8)
    ax.text(.995, .975, "551 samples · 61 sites · five campaigns", transform=ax.transAxes, ha="right", va="top", fontsize=7.4,
            bbox={"facecolor": "white", "alpha": .8, "lw": 0, "pad": 2.5}, zorder=7)
    cb = fig.colorbar(pts, ax=ax, shrink=.6, pad=.012, ticks=[1, 2, 3, 4, 5], aspect=18)
    cb.set_label("Campaigns sampled")
    ins = ax.inset_axes([.005, .60, .27, .39])
    ins.set_facecolor("#DCE8EC")
    draw_coast(ins, lw=.35)
    ins.set_xlim(33, 62)
    ins.set_ylim(11, 33)
    ins.add_patch(Rectangle((lon0, lat0), lon1 - lon0, lat1 - lat0, fill=False, edgecolor="#C0392B", lw=1.1, zorder=5))
    ins.set_xticks([]); ins.set_yticks([])
    for sp in ins.spines.values():
        sp.set_visible(True); sp.set_linewidth(.7); sp.set_color("#55606A")
    ins.set_aspect(1 / np.cos(np.radians(22)))
    ax_b = fig.add_subplot(gs[1, 0])
    counts = pd.crosstab(panel.trip_num, panel.compartment)
    bottom = np.zeros(5)
    for i, c in enumerate(["Surface", "Deep", "PlantRoot"]):
        ax_b.bar(np.arange(5), counts[c], bottom=bottom, label={"PlantRoot": "Root-adjacent", "Deep": "Subsurface", "Surface": "Surface"}[c], color=COL[i])
        bottom += counts[c].to_numpy()
    ax_b.set(xticks=np.arange(5), xticklabels=["Mar ’23", "Jul ’23", "Feb ’24", "Aug ’24", "Oct ’25"], ylabel="Metagenomes")
    ax_b.tick_params(axis="x", rotation=25)
    ax_b.legend(frameon=False, fontsize=7.5)
    ax_b.set_title("B  Repeated biological observations", loc="left", fontweight="bold", fontsize=9, pad=8)
    ax = fig.add_subplot(gs[1, 1])
    ax.axis("off")
    ax.set_title("C  Two zones: deduction conditions, do() intervenes", loc="left", fontweight="bold", fontsize=9, pad=8)
    ax.add_patch(FancyBboxPatch((0, 0), .485, 1, boxstyle="round,pad=0.003,rounding_size=.02", facecolor="#EAF2F0", edgecolor="#9DB7B1",
                                lw=.8, transform=ax.transAxes, zorder=0))
    ax.add_patch(FancyBboxPatch((.515, 0), .485, 1, boxstyle="round,pad=0.003,rounding_size=.02", facecolor="#F8EEE6", edgecolor="#D9B79F",
                                lw=.8, transform=ax.transAxes, zorder=0))
    ax.text(.2425, .965, "Invariant zone", ha="center", va="top", fontsize=7.6, fontweight="bold", color="#2F6F63", transform=ax.transAxes)
    ax.text(.7575, .965, "Causal zone", ha="center", va="top", fontsize=7.6, fontweight="bold", color="#A8431F", transform=ax.transAxes)
    box(ax, .03, .75, .20, .15, "Evidence $E$\nMAGs, KOs,\nmapped reads", "#FFFFFF", fs=6.8)
    box(ax, .255, .75, .20, .15, "Ontology $K$\nOWL TBox,\ncarrier classes", "#FFFFFF", fs=6.8)
    box(ax, .03, .54, .425, .14, "Deduced memberships $D_K(E)$\ncarrier shares $q_s(C)$", "#D8EBE5", fs=7)
    box(ax, .03, .37, .425, .11, "Community factor $\\gamma_g$", "#BFE0D6", fs=7.4, bold=True)
    box(ax, .03, .03, .425, .26, "Constitutive maps $\\Phi$\n(invariant definitions)\n$\\theta{=}W/1000L$, $\\psi$, $D_g$, $P$", "#E3E9F4", fs=6.8)
    box(ax, .545, .75, .41, .15, "Exogenous inputs $u$\nwater, salt, $T$, evaporation, soil", "#FBE5B8", fs=6.8)
    box(ax, .545, .22, .41, .42, "Mechanisms\n$W_t,\\,N_t \\rightarrow A_t \\rightarrow k_g \\rightarrow J_g$\n$\\gamma_g$ enters $k_g$", "#F1C3AB", fs=7.4, bold=True)
    box(ax, .545, .02, .41, .12, "do$(u,W,N)$ regimes: paired effects on $J_g$", "#FFFFFF", ec="#A8431F", fs=6.6, lw=1.2)
    arrow(ax, (.13, .75), (.13, .68)); arrow(ax, (.355, .75), (.30, .68))
    arrow(ax, (.2425, .54), (.2425, .48))
    arrow(ax, (.455, .425), (.545, .425), color="#2F7A63", lw=1.6)
    arrow(ax, (.75, .75), (.75, .64))
    arrow(ax, (.59, .22), (.59, .14), color="#A8431F", lw=1.4)
    ax.text(.62, .18, "do replaces $u$, $W$, $N$;\n$E$, $K$, $\\gamma_g$, $\\Phi$ unchanged", ha="left", va="center", fontsize=5.8,
            color="#6A3A22", transform=ax.transAxes)
    arrow(ax, (.455, .26), (.545, .26), style="<|-|>", color="#34567F")
    save(fig, "v06_primary_twin")


if __name__ == "__main__":
    draw_graph()
    draw_scenarios()
    draw_primary()
