"""Supplementary tables for the composition section of the EQ-BioDT manuscript.

Writes paper/tables/v08_constants_rows.tex (table rows for the 42 external constants, read from
data/external/composition_constants.json) and paper/tables/v08_envelopes.tex (the one-at-a-time constant envelopes of the 12
composition experiments, read from reviews/v07/composition/summary.json). No computation beyond formatting.

Usage: uv run python reproducibility/v08/supplement_tables.py
"""
from __future__ import annotations

import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from reproducibility.v08.paper_figures import readable_name  # noqa: E402

CONSTANTS = ROOT / "data" / "external" / "composition_constants.json"
SUMMARY = ROOT / "reviews" / "v07" / "composition" / "summary.json"
TAB = ROOT / "paper" / "tables"

STATUS_TEXT = {
    "verified_in_source": "verified in source",
    "modeller_assumption": "modeller assumption",
    "source_not_checked": "source not checked",
}

UNIT_TEXT = {
    "mol/L per dS/m": r"mol L$^{-1}$ per dS m$^{-1}$",
    "dimensionless": "--",
    "dS/m": r"dS m$^{-1}$",
    "% yield per dS/m": r"\% per dS m$^{-1}$",
    "mm/day": r"mm day$^{-1}$",
    "fraction": "fraction",
    "degC": r"$^\circ$C",
    "porosity (absolute) per % w/w": r"per \% w/w",
    "mol CO2 per mol H2": r"mol CO$_2$ per mol H$_2$",
    "g/mol": r"g mol$^{-1}$",
    "day/year": r"day yr$^{-1}$",
    "g C m-2 yr-1": r"g C m$^{-2}$ yr$^{-1}$",
    "ppm (mole fraction)": "ppm",
    "L per m2 of panel per wash": r"L m$^{-2}$ per wash",
    "day": "day",
    "g/cm3": r"g cm$^{-3}$",
    "passes": "passes",
    "fraction per degC": r"$^\circ$C$^{-1}$",
    "degC per degC": "--",
    "-": "--",
}

ORGANISATIONS = {
    "IUPAC": "IUPAC",
    "IAU 1976": "IAU",
    "IAU": "IAU",
    "World Meteorological Organization": "WMO",
    "Assumption": "assumed, no source",
}


def num(x: float, digits: int = 4) -> str:
    """Compact number: up to `digits` significant digits, the exact repr when it is short."""
    s = f"{x:.{digits}g}"
    if digits == 4 and float(s) != x and len(repr(float(x))) <= 7:
        s = repr(float(x))
    if s.startswith("-"):
        s = "$-$" + s[1:]
    return s


def short_source(source: dict) -> str:
    authors = source["authors"].split(" (")[0].strip()
    year = source["year"]
    if authors in ORGANISATIONS:
        name = ORGANISATIONS[authors]
        return name if name == "assumed, no source" else f"{name} ({year})"
    et_al = authors.endswith("et al.")
    people = [p.strip() for p in authors.replace(" et al.", "").split(",") if p.strip()]
    surnames = [p.split()[0] for p in people]
    if surnames[0] == "de":
        surnames = ["de Z\\'elicourt"] + surnames[1:]
    if surnames[0].startswith("Hassanpour"):
        surnames[0] = "Hassanpour Adeh"
    if et_al or len(surnames) > 2:
        return f"{surnames[0]} et al. ({year})"
    if len(surnames) == 2:
        return f"{surnames[0]} and {surnames[1]} ({year})"
    return f"{surnames[0]} ({year})"


def constants_rows() -> str:
    consts = json.loads(CONSTANTS.read_text())["constants"]
    lines = []
    for i, c in enumerate(consts, start=1):
        name = readable_name(c["id"], latex=True)
        unit = UNIT_TEXT[c["unit"]]
        src = short_source(c["source"])
        lines.append(
            f"{i} & {name} & {num(c['value'])} [{num(c['low'])}, {num(c['high'])}] & {unit} & {src} & {STATUS_TEXT[c['status']]} \\\\"
        )
    return "\n".join(lines) + "\n"


EXPERIMENTS = {
    "brackish_irrigation": ("Brackish irrigation", r"Alfalfa relative yield at 8 dS m$^{-1}$, daily irrigation", "fraction"),
    "drip_vs_flood": ("Drip or flood", r"H$_2$ uptake change, irrigation every 10 days", r"\% of baseline"),
    "mulch_shade": ("Mulch and shade", r"H$_2$ uptake change, full mulch cover", r"\% of baseline"),
    "biochar": ("Biochar", r"H$_2$ uptake change, 5\% w/w biochar", r"\% of baseline"),
    "inoculation": ("Inoculation", r"CO uptake after 1\% CO-carrier inoculation of a closed-gate sample", r"$\mu$mol m$^{-2}$"),
    "carbon_fixation_bound": ("Carbon fixation bound", "Upper bound of carbon fixed by hydrogen oxidizers, sample with the highest autotroph share", r"g C m$^{-2}$ yr$^{-1}$"),
    "h2_leakage": ("Hydrogen leakage", r"H$_2$ uptake change at twice the background mixing ratio", r"\% of baseline"),
    "solar_park": ("Solar park", r"H$_2$ uptake change, full panel cover", r"\% of baseline"),
    "cloud_seeding": ("Cloud seeding", r"H$_2$ uptake change, full-intensity seeding (+10\% rain on rain days)", r"\% of baseline"),
    "compaction": ("Vehicle compaction", r"H$_2$ uptake change, 50 passes", r"\% of baseline"),
    "warming": ("Climate warming", r"H$_2$ uptake change, +4 $^\circ$C with unchanged rainfall", r"\% of baseline"),
    "inoculated_alfalfa_salinity": ("Inoculated alfalfa", r"Inoculated alfalfa relative yield at 7.81 dS m$^{-1}$, salinity fitted to the field trial (fresh weight)", "fraction"),
}

PROVENANCE = {"twin": "twin"}


def envelope_table() -> str:
    exps = json.loads(SUMMARY.read_text())["experiments"]
    lines = [
        r"\begin{tabular}{@{}rp{2.3cm}p{5.0cm}p{3.1cm}p{2.1cm}p{1.8cm}@{}}",
        r"\toprule",
        r"\# & Experiment & Headline & Central [low, high] & Unit & Provenance \\",
        r"\midrule",
    ]
    for e in exps:
        name, headline, unit = EXPERIMENTS[e["id"]]
        h = e["headline"]
        v, lo, hi = h["value"], h["low"], h["high"]
        if e["id"] == "carbon_fixation_bound":
            v, lo, hi = v / 1e6, lo / 1e6, hi / 1e6
        prov = PROVENANCE.get(h["provenance"], "external model")
        lines.append(
            f"{e['n']} & {name} & {headline} & {num(v, 3)} [{num(lo, 3)}, {num(hi, 3)}] & {unit} & {prov} \\\\"
        )
    lines += [r"\bottomrule", r"\end{tabular}"]
    return "\n".join(lines) + "\n"


def main() -> None:
    TAB.mkdir(exist_ok=True)
    (TAB / "v08_constants_rows.tex").write_text("% Generated by reproducibility/v08/supplement_tables.py; do not edit by hand.\n" + constants_rows())
    (TAB / "v08_envelopes.tex").write_text(envelope_table())
    print("wrote v08_constants_rows.tex and v08_envelopes.tex")


if __name__ == "__main__":
    main()
