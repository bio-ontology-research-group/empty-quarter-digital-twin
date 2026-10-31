"""v0.7 composition figure (small multiples, one panel per experiment) and LaTeX table, read from reviews/v07/composition.

Writes paper/figures/v07_composition.{pdf,png} and paper/tables/v07_composition.tex. No computation beyond formatting."""
from pathlib import Path
import json

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "reviews/v07/composition"
FIG = ROOT / "paper/figures"
TAB = ROOT / "paper/tables"
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 7, "axes.spines.top": False, "axes.spines.right": False,
                     "pdf.fonttype": 42, "savefig.dpi": 200, "mathtext.fontset": "dejavusans", "svg.hashsalt": "eq-biodt-v07"})
INK = "#26323A"
COLOURS = ["#D2693C", "#4E8FA6", "#4F9A82", "#8E7CB3", "#D9A441"]
ORDER = ["brackish_irrigation", "drip_vs_flood", "mulch_shade", "biochar", "inoculation", "carbon_fixation_bound", "h2_leakage",
         "solar_park", "cloud_seeding", "compaction", "warming", "inoculated_alfalfa_salinity"]


def panel(ax, n, doc):
    f = doc["figure"]
    for i, s in enumerate(f["series"]):
        pts = [(x, y, lo, hi) for x, y, lo, hi in zip(s["x"], s["y"], s["low"], s["high"]) if y is not None]
        if not pts:
            continue
        x, y, lo, hi = zip(*pts)
        c = COLOURS[i % len(COLOURS)]
        if len(x) > 1:
            ax.plot(x, y, "-o", color=c, lw=1.2, ms=3, label=s["label"])
            band = [(a, l, h) for a, l, h in zip(x, lo, hi) if l is not None and h is not None]
            if band:
                ax.fill_between([b[0] for b in band], [b[1] for b in band], [b[2] for b in band], color=c, alpha=.2, lw=0)
        else:
            ax.errorbar(x, y, yerr=[[y[0] - lo[0] if lo[0] is not None else 0], [hi[0] - y[0] if hi[0] is not None else 0]], fmt="o", color=c,
                        ms=4, capsize=2, label=s["label"])
    if f.get("baseline_y") is not None:
        ax.axhline(f["baseline_y"], color=INK, lw=.7, ls=":")
    ax.set_title(f["title"], fontsize=7.5, loc="left", color=INK)
    ax.set_xlabel(f["xlabel"], fontsize=6.5)
    ax.set_ylabel(f["ylabel"], fontsize=6.5)
    ax.tick_params(labelsize=6)
    if len(f["series"]) > 1:
        ax.legend(fontsize=5, frameon=False, loc="best")


def figure():
    docs = {i: json.loads((SRC / f"{i}.json").read_text()) for i in ORDER}
    fig, axes = plt.subplots(3, 4, figsize=(11, 7.2))
    for ax, (i, d) in zip(axes.flat, docs.items()):
        panel(ax, i, d)
    fig.tight_layout()
    FIG.mkdir(exist_ok=True)
    fig.savefig(FIG / "v07_composition.pdf", bbox_inches="tight", metadata={"CreationDate": None, "ModDate": None, "Producer": None})
    fig.savefig(FIG / "v07_composition.png", bbox_inches="tight", dpi=170, metadata={"Software": None})
    plt.close(fig)


def tex(s):
    return (s.replace("\\", "\\textbackslash{}").replace("&", "\\&").replace("%", "\\%").replace("_", "\\_").replace("#", "\\#")
             .replace("$", "\\$"))


def num(x):
    if x is None:
        return "--"
    if x < 0:
        return "$-$" + num(-x)
    if x != 0 and (abs(x) >= 1e4 or abs(x) < 1e-2):
        m, e = f"{x:.3g}".split("e") if "e" in f"{x:.3g}" else (f"{x:.3g}", None)
        return f"${m}\\times10^{{{int(e)}}}$" if e else m
    return f"{x:.3g}"


PROV = {"twin": "twin", "deduced": "deduced"}
UNIT = {"% of baseline H2 uptake": "\\% of baseline", "fraction of unsalted yield": "fraction", "umol m-2 over 14 d": "$\\mu$mol m$^{-2}$",
        "ug C m-2 yr-1 (annualised)": "$\\mu$g C m$^{-2}$ yr$^{-1}$", "fraction of unsalted, uninoculated yield": "fraction"}


def table():
    s = json.loads((SRC / "summary.json").read_text())
    rows = []
    for r in s["experiments"]:
        h = r["headline"]
        prov = h["provenance"].replace("external:", "ext.: ").replace("_", "\\_")
        rows.append(f"{r['n']} & {tex(r['id'].replace('_', ' '))} & {tex(h['label'])} & {num(h['value'])} [{num(h['low'])}, {num(h['high'])}] & "
                    f"{UNIT.get(h['unit'], tex(h['unit']))} & {prov} \\\\")
    body = "\n".join(rows)
    t = ("\\begin{tabular}{@{}rp{2.1cm}p{4.4cm}p{3.3cm}p{2cm}p{2.3cm}@{}}\n\\toprule\n\\# & Experiment & Headline & Central [low, high] & Unit & Provenance \\\\\n\\midrule\n"
         + body + "\n\\bottomrule\n\\end{tabular}\n")
    TAB.mkdir(exist_ok=True)
    (TAB / "v07_composition.tex").write_text(t)


if __name__ == "__main__":
    figure()
    table()
