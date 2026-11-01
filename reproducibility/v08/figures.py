"""v0.8 UQ sensitivity figure and LaTeX tables, read from data/external/twin_parameter_ranges.json and reviews/v08/questions.

Writes paper/figures/v08_sensitivity.{pdf,png}, paper/tables/v08_parameters.tex and paper/tables/v08_uq.tex.
No computation beyond formatting."""
from __future__ import annotations

import json
from pathlib import Path
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from reproducibility.v08.paper_figures import (  # noqa: E402
    FIG_WIDTH_IN,
    HEADLINE_LABELS_READABLE,
    QUESTION_TITLES_SHORT,
    readable_name,
    to_latex,
)

QDIR = ROOT / "reviews/v08/questions"
PARAM_PATH = ROOT / "data/external/twin_parameter_ranges.json"
FIG = ROOT / "paper/figures"
TAB = ROOT / "paper/tables"

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 6.8,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "savefig.dpi": 300,
    "mathtext.fontset": "dejavusans",
    "svg.hashsalt": "eq-biodt-v08",
})

INK = "#26323A"
STATUS_COLOURS = {
    "modeller_assumption": "#D2693C",
    "verified_in_source": "#4F9A82",
    "source_not_checked": "#4E8FA6",
}
STATUS_LABELS = {
    "modeller_assumption": "modeller assumption",
    "verified_in_source": "verified in source",
    "source_not_checked": "source not checked",
}

PARAM_UNITS = {
    "1/s": "s$^{-1}$",
    "MPa": "MPa",
    "-": "--",
    "m2/s": "m$^2$ s$^{-1}$",
    "1/day": "day$^{-1}$",
}

OUTCOME_UNITS = {
    "%": "\\%",
    "fraction": "fraction",
    "µmol m⁻²": "$\\mu$mol m$^{-2}$",
    "g C m⁻² yr⁻¹": "g C m$^{-2}$ yr$^{-1}$",
}


def tex(s: str) -> str:
    return (
        s.replace("\\", "\\textbackslash{}")
        .replace("&", "\\&")
        .replace("%", "\\%")
        .replace("_", "\\_")
        .replace("#", "\\#")
        .replace("$", "\\$")
    )


def num(x: float | None) -> str:
    if x is None:
        return "--"
    if x < 0:
        return "$-$" + num(-x)
    if x != 0 and (abs(x) >= 1e4 or abs(x) < 1e-3):
        s = f"{x:.3g}"
        if "e" in s:
            m, e = s.split("e")
            return f"${m}\\times10^{{{int(e)}}}$"
        return s
    return f"{x:.3g}"


def pnum(x: float) -> str:
    if x < 0:
        return "$-$" + pnum(-x)
    if x != 0 and abs(x) < 1e-4:
        s = f"{x:.3g}"
        if "e" in s:
            m, e = s.split("e")
            return f"${m}\\times10^{{{int(e)}}}$"
    return f"{x:g}"


def parameters_table() -> None:
    reg = json.loads(PARAM_PATH.read_text())["parameters"]
    rows = []
    for i, p in enumerate(reg, start=1):
        pname = readable_name(p["id"], latex=True)
        if p["id"] == "d_h2_m2_s":
            pname = r"free-air diffusivity (H$_2$)"
        elif p["id"] == "d_co_m2_s":
            pname = "free-air diffusivity (CO)"
        dist = tex(p["dist"])
        val_str = f"{pnum(p['value'])} [{pnum(p['low'])}, {pnum(p['high'])}]"
        unit = PARAM_UNITS.get(p["unit"], tex(p["unit"]))
        status = STATUS_LABELS.get(p["status"], tex(p["status"]))
        rows.append(f"{i} & {pname} & {dist} & {val_str} & {unit} & {status} \\\\")
    body = "\n".join(rows)
    t = (
        "\\begin{tabular}{@{}rlp{1.9cm}p{4.1cm}p{1.4cm}p{3.1cm}@{}}\n"
        "\\toprule\n"
        "\\# & Parameter & Distribution & Reference [low, high] & Unit & Status \\\\\n"
        "\\midrule\n"
        + body
        + "\n\\bottomrule\n\\end{tabular}\n"
    )
    TAB.mkdir(exist_ok=True)
    (TAB / "v08_parameters.tex").write_text(t)


def _headline_outcome(doc: dict) -> dict:
    u = doc["uncertainty"]
    ser = next(s for s in doc["series"] if s["id"] == u["headline_series"])
    lvl = ser["levels"][u["headline_level"]]
    return lvl["outcomes"][u["headline_output"]]


def _top_contributors(doc: dict, k: int = 2) -> list[dict]:
    u = doc["uncertainty"]
    status_by_name = {inp["name"]: inp["status"] for inp in u["inputs"]}
    ranked = sorted(u["ST"].items(), key=lambda kv: (-kv[1][0], kv[0]))
    return [
        {"id": name, "ST": st, "status": status_by_name[name]}
        for name, st in ranked[:k]
    ]


def uq_table() -> None:
    idx = json.loads((QDIR / "index.json").read_text())["questions"]
    rows = []
    for i, qinfo in enumerate(idx, start=1):
        slug = qinfo["slug"]
        doc = json.loads((QDIR / f"{slug}.json").read_text())
        u = doc["uncertainty"]
        out = _headline_outcome(doc)
        unit = OUTCOME_UNITS.get(out["unit"], tex(out["unit"]))
        top = _top_contributors(doc, 2)
        top_str = ", ".join(
            f"{readable_name(t['id'], latex=True)} ({num(t['ST'][0])})"
            for t in top
        )
        assum = u["group_ST"]["modeller_assumption"][0]
        q_title = QUESTION_TITLES_SHORT.get(slug, tex(slug))
        hl_label = to_latex(HEADLINE_LABELS_READABLE.get(slug, u["headline_label"]))
        label = f"{hl_label} ({unit})"
        ref_val = out.get("reference", out["value"])
        val_str = f"{num(ref_val)} & {num(u['median'])} [{num(u['q05'])}, {num(u['q95'])}]"
        rows.append(
            f"{i} & {q_title} & {label} & {val_str} & {top_str} & {num(assum)} \\\\"
        )
    body = "\n".join(rows)
    t = (
        "\\begin{tabular}{@{}rp{2.3cm}p{3.4cm}rp{2.5cm}p{4.4cm}r@{}}\n"
        "\\toprule\n"
        "\\# & Question & Headline level (unit) & Central & Median [5\\%, 95\\%] & Top two $S_T$ contributors & $S_T^{\\mathrm{assum}}$ \\\\\n"
        "\\midrule\n"
        + body
        + "\n\\bottomrule\n\\end{tabular}\n"
    )
    TAB.mkdir(exist_ok=True)
    (TAB / "v08_uq.tex").write_text(t)


def sensitivity_figure() -> None:
    idx = json.loads((QDIR / "index.json").read_text())["questions"]
    docs = [json.loads((QDIR / f"{q['slug']}.json").read_text()) for q in idx]
    n = len(docs)
    y = np.arange(n)[::-1]

    fig = plt.figure(figsize=(FIG_WIDTH_IN, 4.4))
    gs = fig.add_gridspec(
        1,
        2,
        width_ratios=[1.0, 1.25],
        left=0.20,
        right=0.985,
        top=0.91,
        bottom=0.11,
        wspace=0.08,
    )
    ax_grp = fig.add_subplot(gs[0, 0])
    ax_top = fig.add_subplot(gs[0, 1], sharey=ax_grp)

    statuses = ["modeller_assumption", "verified_in_source", "source_not_checked"]
    bar_h = 0.24
    offsets = [bar_h, 0.0, -bar_h]

    for st, off in zip(statuses, offsets):
        vals = [d["uncertainty"]["group_ST"][st][0] for d in docs]
        los = [max(0.0, v - d["uncertainty"]["group_ST"][st][1]) for v, d in zip(vals, docs)]
        his = [max(0.0, d["uncertainty"]["group_ST"][st][2] - v) for v, d in zip(vals, docs)]
        ax_grp.barh(
            y + off,
            vals,
            height=bar_h * 0.92,
            color=STATUS_COLOURS[st],
            label=STATUS_LABELS[st],
            alpha=0.9,
        )
        ax_grp.errorbar(
            vals,
            y + off,
            xerr=[los, his],
            fmt="none",
            ecolor=INK,
            elinewidth=0.6,
            capsize=1.5,
        )

    ax_grp.text(
        -0.44,
        1.06,
        "a",
        transform=ax_grp.transAxes,
        fontsize=9.5,
        fontweight="bold",
        color=INK,
        ha="left",
        va="top",
    )
    ax_grp.axvline(1.0, color=INK, lw=0.6, ls=":")
    ax_grp.set_yticks(y)
    ax_grp.set_yticklabels(
        [f"{i + 1}. {QUESTION_TITLES_SHORT.get(d['slug'], d['slug'])}" for i, d in enumerate(docs)],
        fontsize=6.2,
        color=INK,
    )
    ax_grp.set_xlabel("Group total-order Sobol' index $S_T^{\\mathrm{group}}$", fontsize=6.6, color=INK)
    ax_grp.set_title("Sensitivity by verification status", fontsize=7.0, fontweight="bold", loc="left", color=INK)
    ax_grp.set_xlim(0, 1.42)
    ax_grp.legend(fontsize=5.6, frameon=False, loc="lower right")

    top_off = [0.16, -0.16]
    for rank, off in enumerate(top_off):
        vals, los, his, cols, names = [], [], [], [], []
        for d in docs:
            tc = _top_contributors(d, 2)[rank]
            v, lo, hi = tc["ST"]
            vals.append(v)
            los.append(max(0.0, v - lo))
            his.append(max(0.0, hi - v))
            cols.append(STATUS_COLOURS[tc["status"]])
            names.append(readable_name(tc["id"], latex=False))
        ax_top.barh(y + off, vals, height=0.28, color=cols, alpha=0.88 if rank == 0 else 0.58)
        ax_top.errorbar(vals, y + off, xerr=[los, his], fmt="none", ecolor=INK, elinewidth=0.6, capsize=1.5)
        for yi, v, hi, nm in zip(y + off, vals, his, names):
            if v > 0.02:
                ax_top.text(v + hi + 0.015, yi, nm, va="center", ha="left", fontsize=5.1, color=INK)

    ax_top.text(
        -0.06,
        1.06,
        "b",
        transform=ax_top.transAxes,
        fontsize=9.5,
        fontweight="bold",
        color=INK,
        ha="left",
        va="top",
    )
    ax_top.axvline(1.0, color=INK, lw=0.6, ls=":")
    ax_top.tick_params(axis="y", left=False, labelleft=False)
    ax_top.set_xlabel("Individual total-order Sobol' index $S_T$", fontsize=6.6, color=INK)
    ax_top.set_title("Top two parameter contributors", fontsize=7.0, fontweight="bold", loc="left", color=INK)
    ax_top.set_xlim(0, 1.55)

    FIG.mkdir(exist_ok=True)
    fig.savefig(
        FIG / "v08_sensitivity.pdf",
        metadata={"CreationDate": None, "ModDate": None, "Producer": None},
    )
    fig.savefig(FIG / "v08_sensitivity.png", dpi=300, metadata={"Software": None})
    plt.close(fig)


def main() -> None:
    parameters_table()
    uq_table()
    sensitivity_figure()


if __name__ == "__main__":
    main()
