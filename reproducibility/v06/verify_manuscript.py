"""Verify every hard-coded number in paper/main.tex and paper/supplementary.tex against the generated outputs.

Usage: uv run python reproducibility/v06/verify_manuscript.py [-v] [--list-literals]

Design
------
1. Both .tex files are normalised (comments, equations, citations, labels, column specs removed; `~` and `\\%` resolved).
2. Every numeric literal in the remaining prose, captions and table headers is extracted (`manuscript_text.NUM`),
   after blanking equation constants, identifiers (C14-2, K01100, M00165, S1), years, dates, DOIs and sub/superscripts.
3. REGISTRY maps a context template (text around the number, `[[x]]` marks a literal) to a computation path
   (an expression on the generated files) and the rounding rule: the literal is compared with the computed value
   rounded half-up to the literal's own number of decimals (so "6.13" must equal round(value, 2); "1,078" an integer;
   "8.61%" is read as 8.61).  A registry entry that no longer matches the text fails (stale registry); a mismatch fails.
4. ALLOW lists literals that are not outputs: fixed design inputs, definitions and arithmetic constants.  Each entry carries
   a reason, and, where the input is defined in code or in the model contract, a cross-check against that definition.
5. Any numeric literal in the text covered by neither list fails.
6. Every \\macro used in the .tex files must be defined (paper/tables/*.tex or the preamble), and every value macro
   used in the text is recomputed independently from the generated outputs (MACROS).  The generated row tables
   (contrasts, paired errors, external transfer, campaigns) are cross-checked against their sources.

Not covered: numbers written as words ("five campaigns", "two auxiliary locations", "ten sites"), the generated
paper/assumptions.tex (rendered from the executable contract by reproducibility/v03/reproduce.py), and the bibliography.
"""
import argparse
import ast
import json
import re
import sys
from decimal import ROUND_HALF_UP, Decimal
from itertools import product
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from manuscript_text import NUM, NUMP, macros_used, mask, normalise  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
PAPER = ROOT / "paper"
FILES = ("main.tex", "supplementary.tex")
UM = 1e6   # mol m-2 -> umol m-2
NM = 1e9   # mol m-2 s-1 -> nmol m-2 s-1


# ----------------------------------------------------------------------------------------------- sources

def _ledger_core_sites():
    """Numbered core sites in the March 2023 field ledger (the full transect)."""
    import csv
    rows = csv.DictReader(open(ROOT / "data/sources/trip1-field-ledger.tsv"), delimiter="\t")
    return {int(r["site"]) for r in rows if r["site"].isdigit()}


def _missing_core_sites(d):
    """Ledger core sites without a metagenome in the released cohort."""
    import csv
    have = {int(r["site"]) for r in csv.DictReader(open(ROOT / "data/v03/sample_evidence.tsv"), delimiter="\t") if r["site"].isdigit()}
    return sorted(_ledger_core_sites() - have)

def J(p):
    return json.loads((ROOT / p).read_text())


def T(p):
    return pd.read_csv(ROOT / p, sep="\t")


class D:
    """Generated outputs and design sources, loaded once."""
    def __init__(self):
        self.res = J("reviews/v03/results.json")
        self.aud = J("reviews/v03/audit.json")
        self.ont = J("data/v03/ontology_report.json")
        self.jus = J("data/v03/justifications.json")
        self.con = J("data/v03/model_contract.json")
        self.ev = T("data/v03/sample_evidence.tsv")
        self.pc = J("reviews/v05/primary_connection.json")
        self.ext = J("reviews/v05/external_transfer.json")
        self.trace = J("reviews/v05/deduction_causal_trace.json")
        self.nam_c = T("reviews/v05/namib_hydration_contrasts.tsv")
        self.nam_p = T("reviews/v05/namib_conditional_predictions.tsv")
        self.nam_o = T("reviews/v05/namib_observations.tsv")
        self.son_h = T("reviews/v05/sonoran_hydration_contrasts.tsv")
        self.son_p = T("reviews/v05/sonoran_field_predictions.tsv")
        self.camp = T("reviews/v05/primary_campaigns.tsv")
        self.v6 = J("reviews/v06/results.json")
        self.psum = J("reviews/v06/paired_regimes_summary.json")
        self.tvj = J("reviews/v06/treatment_versions.json")
        self.tvt = T("reviews/v06/treatment_versions.tsv")
        self.vr = J("reviews/v06/version_relevance.json")
        self.adm = T("reviews/v06/admissibility.tsv")
        self.ps = T("reviews/v06/paired_regimes_by_sample.tsv")
        self.defaults = self.con["parameter_defaults"]
        self.clauses = {c["id"]: c["statement"] for c in self.con["clauses"]}
        self.src = {p: (ROOT / p).read_text() for p in (
            "reproducibility/v06/analyse.py", "reproducibility/v05/external_transfer.py", "src/eq_biodt/v03/resources.py")}

    # derived helpers -------------------------------------------------------------------------------------------
    def cl(self, cid, rx):
        m = re.search(rx, self.clauses[cid])
        assert m, f"contract clause {cid} does not match {rx!r}"
        return float(m.group(1))

    def nc(self, scenario, steps="96"):
        return next(r for r in self.aud["numerical_convergence"] if r["scenario"] == scenario)["h2_uptake_by_steps"][steps]

    def comp(self, c):
        return self.ev[self.ev.compartment == c]

    def conv_max_rel_pct(self):
        return max(r["relative_error_96_vs_384"] for r in self.aud["numerical_convergence"]) * 100

    def contrast(self, rx):
        return next(c for c in self.aud["matched_compartment_contrasts"] if re.fullmatch(rx, c["endpoint"] + "|" + c["contrast"]))

    def aud_contrast(self, endpoint, contrast):
        return next(c for c in self.aud["matched_compartment_contrasts"]
                    if c["endpoint"] == endpoint and c["contrast"] == contrast)

    def fit_values(self, key):
        return [f["tuned"][key] for f in self.ext["sonoran"]["fits"]]

    def sonoran_model(self, name):
        return next(m for m in self.ext["sonoran"]["field_metrics"] if m["model"] == name)

    def namib(self, site):
        return self.nam_c[self.nam_c.site == site].iloc[0]

    def sonoran_range(self, lo, hi):
        d = self.son_h[(self.son_h.dry_gwc.round(3) == lo) & (self.son_h.wet_gwc.round(3) == hi)]
        return d.difference.min(), d.difference.max()

    def grid(self):
        m = re.search(r"product\((\[[^\]]*\]),\s*(\[[^\]]*\]),\s*(\[[^\]]*\]),\s*(\[[^\]]*\])\)", self.src["reproducibility/v05/external_transfer.py"])
        return [ast.literal_eval(g) for g in m.groups()]

    def four_pulse_days(self):
        assert "int(i * 14 / n) for i in range(n)" in self.src["reproducibility/v06/analyse.py"]
        return sorted({int(i * 14 / 4) for i in range(4)})

    def matched_contrast(self, a, b, col):
        """Deep-minus-root style contrast: replicate means within site x campaign x compartment, matched pairs."""
        ps = self.ps.copy()
        g = ps.groupby(["site", "trip", "compartment"])[col].mean().unstack("compartment")
        d = (g[a] - g[b]).dropna()
        return d.mean(), len(d), d.index.get_level_values("site").nunique()


# ----------------------------------------------------------------------------------------------- comparison
def dec_places(lit):
    core = lit.replace(",", "").rstrip("%")
    return len(core.split(".")[1]) if "." in core else 0


MACRO_VALUES = {}


def resolve(lit):
    """A placeholder may hold a value macro (\\Name{}); it is replaced by its definition in paper/tables/*.tex."""
    if lit.startswith("\\"):
        name = lit[1:].split("{")[0]
        assert name in MACRO_VALUES, f"macro \\{name} is not defined in paper/tables/*.tex"
        return MACRO_VALUES[name].replace("\\%", "%")
    return lit


def parse(lit):
    core = lit.replace(",", "").rstrip("%")
    return Decimal(core.replace("E", "e"))


def rounded(value, places):
    q = Decimal(1).scaleb(-places)
    return Decimal(repr(round(float(value), 9))).quantize(q, rounding=ROUND_HALF_UP)


def literal_matches(lit, value):
    if value is None or (isinstance(value, float) and not np.isfinite(value)):
        return False
    if "E" in lit:
        return float(parse(lit)) == float(value)
    return rounded(value, dec_places(lit)) == parse(lit)


class OneOf:
    """The literal must equal (after rounding) one of a set computed from the outputs."""
    def __init__(self, fn):
        self.fn = fn


def template_regex(t):
    parts = re.split(r"\[\[(\w+)\]\]", t)
    rx = ""
    for i, p in enumerate(parts):
        if i % 2 == 0:
            rx += re.escape(p).replace(r"\ ", r"\s+")
        else:
            rx += f"(?P<{p}>{NUMP}|\\\\[A-Za-z]+\\{{\\}})"
    return re.compile(rx)


REGISTRY = []   # (template, {group: expr}, min_matches)
ALLOW = []      # (template, reason, check, min_matches)


def chk(template, n=1, **exprs):
    REGISTRY.append((template, exprs, n))


def allow(template, reason, check=None, n=1):
    ALLOW.append((template, reason, check, n))


DESIGN = "design input, see process/contract"


# ----------------------------------------------------------------------------------------------- registry
def build_registry():
    # ---- abstract
    chk("from [[a]] primary soil metagenomes collected across five campaigns and [[b]] reference genomes",
        a=lambda d: d.pc["samples"], b=lambda d: d.res["reference_MAGs"])
    chk("resolves ambiguous hydrogenase--carbon-fixation overlap in [[a]] samples and supports [[b]] candidate genomes",
        a=lambda d: d.pc["ambiguous_marginals"], b=lambda d: d.ont["entailed_class_counts"]["H2AutotrophicCandidateMAG"])
    chk("reduces held-campaign error from [[a]] to [[b]] nmol", a=lambda d: d.sonoran_model("frozen_reference")["rmse"],
        b=lambda d: d.sonoran_model("calibrated_retention")["rmse"])
    chk("compared with [[a]] for a training-mean baseline", a=lambda d: d.sonoran_model("training_mean")["rmse"])
    chk("a fixed [[a]] mm water addition", a=lambda d: d.v6["treatment_versions"]["total_water_mm"])

    # ---- results: primary observations
    chk("The twin represents [[a]] soil metagenomes", a=lambda d: d.pc["samples"])
    chk("cover [[a]] site--campaign visits at [[b]] locations:",
        a=lambda d: d.pc["site_campaigns"], b=lambda d: d.pc["sites"])
    chk("extends from [[a]] to [[b]]$^\\circ$E and from [[c]] to [[d]]$^\\circ$N",
        a=lambda d: d.ev.lon.min(), b=lambda d: d.ev.lon.max(), c=lambda d: d.ev.lat.min(), d=lambda d: d.ev.lat.max())
    chk("contribute [[a]], [[b]] and [[c]] metagenomes, respectively",
        a=lambda d: d.pc["compartments"]["Surface"], b=lambda d: d.pc["compartments"]["Deep"],
        c=lambda d: d.pc["compartments"]["PlantRoot"])
    chk("Of the [[a]] specimens, [[b]] have usable mapped compositions",
        a=lambda d: d.pc["samples"], b=lambda d: d.pc["usable_mapping_samples"])

    # ---- results: deduction
    chk("including [[a]] supported genomes with transaldolase",
        a=lambda d: d.aud["pathway_audit"]["transaldolase_route_without_sbpase"])
    chk("positive width in [[a]] of [[b]] mapped samples, with median width [[c]] percentage points",
        a=lambda d: d.pc["ambiguous_marginals"], b=lambda d: d.pc["usable_mapping_samples"],
        c=lambda d: d.pc["median_overlap_interval_percentage_points"])
    chk("A single Rubisco screen includes [[a]] MAGs, and Rubisco plus phosphoribulokinase includes [[b]]; complete module support narrows the admitted set to [[c]]",
        a=lambda d: d.ont["entailed_class_counts"]["RuBisCOMAG"], b=lambda d: d.ont["entailed_class_counts"]["RuBisCOPRKMAG"],
        c=lambda d: d.ont["entailed_class_counts"]["CBBCompleteMAG"])
    chk("a carbohydrate-carbon ceiling of [[a]] nmol",
        a=lambda d: d.trace["ontology_supported"]["fluxes"]["carbohydrate_c"] * NM)
    chk("raises the computed ceiling to [[a]], a [[b]]-fold change",
        a=lambda d: d.trace["unsupported_CO_admission"]["fluxes"]["carbohydrate_c"] * NM,
        b=lambda d: d.trace["unsupported_CO_admission"]["fluxes"]["carbohydrate_c"] /
        d.trace["ontology_supported"]["fluxes"]["carbohydrate_c"])
    chk("genome-resolved intersection, for [[a]] samples", a=lambda d: d.pc["mapping_overlap_samples"])

    # ---- results: compartments
    chk("is [[a]] in subsurface sand, [[b]] in surface sand and [[c]] in root-adjacent soil",
        a=lambda d: d.comp("Deep").cbb_all_steps_detected_share_mapped.median() * 100,
        b=lambda d: d.comp("Surface").cbb_all_steps_detected_share_mapped.median() * 100,
        c=lambda d: d.comp("PlantRoot").cbb_all_steps_detected_share_mapped.median() * 100)
    chk("Subsurface sand exceeded roots by [[a]] percentage points across [[b]] matched groups, with a site-bootstrap 95% interval of [[c]]--[[d]] points",
        a=lambda d: d.aud_contrast("cbb_all_steps_detected_share_mapped", "Deep minus PlantRoot")["mean_difference_percentage_points"],
        b=lambda d: d.aud_contrast("cbb_all_steps_detected_share_mapped", "Deep minus PlantRoot")["matched_site_campaigns"],
        c=lambda d: d.aud_contrast("cbb_all_steps_detected_share_mapped", "Deep minus PlantRoot")["site_bootstrap_95_percentile"][0],
        d=lambda d: d.aud_contrast("cbb_all_steps_detected_share_mapped", "Deep minus PlantRoot")["site_bootstrap_95_percentile"][1])
    chk("resample sites with all paired visits ([[a]] draws)", a=lambda d: d.aud["bootstrap_draws"])
    chk("visits at [[a]] sites, geographic and seasonal ridge regression reduced RMSE relative to persistence by [[b]] on the",
        a=lambda d: next(c for c in d.aud["forecast_paired_comparisons"] if c["baseline"] == "persistence")["sites"],
        b=lambda d: -next(c for c in d.aud["forecast_paired_comparisons"] if c["baseline"] == "persistence")["rmse_difference_ridge_minus_baseline"])
    chk("with a site-bootstrap interval of [[a]]--[[b]]. Its difference",
        a=lambda d: -next(c for c in d.aud["forecast_paired_comparisons"] if c["baseline"] == "persistence")["site_bootstrap_95_percentile"][1],
        b=lambda d: -next(c for c in d.aud["forecast_paired_comparisons"] if c["baseline"] == "persistence")["site_bootstrap_95_percentile"][0])
    chk("and nominal [[a]] residual intervals covered", a=lambda d: d.cl("F3", r"(\d+)-percent"))

    # ---- results: cross-desert tests
    chk("Mean wet-to-dry ratios were [[a]] (within-site bootstrap 95% interval [[b]]--[[c]]) at C14-2 and [[e]] ([[f]]--[[g]]) at C14-4",
        a=lambda d: d.namib("C14-2").fold, b=lambda d: d.namib("C14-2").fold_lcl, c=lambda d: d.namib("C14-2").fold_ucl,
        e=lambda d: d.namib("C14-4").fold, f=lambda d: d.namib("C14-4").fold_lcl, g=lambda d: d.namib("C14-4").fold_ucl)
    chk("at bulk densities of [[a]] and [[b]] g cm$^{-3}$", a=OneOf(lambda d: set(d.nam_p.assumed_bulk_density_g_cm3)),
        b=OneOf(lambda d: set(d.nam_p.assumed_bulk_density_g_cm3)))
    chk("At [[a]] g cm$^{-3}$, hydrated C14-2 soil nearly fills", a=OneOf(lambda d: set(d.nam_p.assumed_bulk_density_g_cm3)))
    chk("increasing gravimetric water content from [[a]] to [[b]] g g$^{-1}$ increased mean hydrogen uptake by [[c]]--[[d]] nmol",
        a=lambda d: 0.02, b=lambda d: 0.20, c=lambda d: d.sonoran_range(0.02, 0.2)[0], d=lambda d: d.sonoran_range(0.02, 0.2)[1])
    chk("Uptake also increased from [[a]] to [[b]] g g$^{-1}$ in every plot, by [[c]]--[[d]] nmol",
        a=lambda d: 0.10, b=lambda d: 0.20, c=lambda d: d.sonoran_range(0.1, 0.2)[0], d=lambda d: d.sonoran_range(0.1, 0.2)[1])
    chk("an assumed bulk density of [[a]] g cm$^{-3}$, the twin instead", a=lambda d: d.ext["sonoran"]["lab_default_assumed_density_g_cm3"])
    chk("against [[a]] Sonoran field measurements, aggregated into [[b]] plot--campaign means",
        a=lambda d: d.ext["sonoran"]["field_collars"], b=lambda d: d.ext["sonoran"]["field_plot_campaign_means"])
    chk("the frozen reference model had RMSE [[a]] nmol", a=lambda d: d.sonoran_model("frozen_reference")["rmse"])
    chk("gave RMSE [[a]]. Allowing", a=lambda d: d.sonoran_model("calibrated_kinetics")["rmse"])
    chk("reduced error to [[a]], compared with [[b]] for the training mean and [[c]] for a linear",
        a=lambda d: d.sonoran_model("calibrated_retention")["rmse"], b=lambda d: d.sonoran_model("training_mean")["rmse"],
        c=lambda d: d.sonoran_model("linear_moisture_temperature")["rmse"])
    chk("no-temperature-amplification member, $Q_{10}=[[a]]$", a=lambda d: _unique(d.fit_values("q10")))
    chk("computed on [[a]] plot--campaign means", a=lambda d: d.ext["sonoran"]["field_plot_campaign_means"])

    # ---- results: interventions
    chk("after a [[a]] mm pulse, $\\mathrm{do}(\\psi=-[[b]] \\mathrm{MPa})$ can",
        a=lambda d: d.v6["treatment_versions"]["total_water_mm"], b=lambda d: abs(d.vr["targets"][0]["psi_mpa"]))
    chk("at $-[[a]]$ MPa, where both versions are valid", a=lambda d: abs(d.vr["targets"][1]["psi_mpa"]))
    chk("Daily additions totalling [[a]] mm yield [[b]] $\\mu$mol H$_2$ m$^{-2}$ at the reference parameters",
        a=lambda d: d.v6["treatment_versions"]["total_water_mm"], b=lambda d: d.nc("daily_pulses") * UM)
    chk("increasing integration resolution from [[a]] to [[b]] steps per day changes reference-schedule uptake by at most [[c]]",
        a=OneOf(lambda d: set(d.aud["numerical_convergence"][1]["h2_uptake_by_steps"])),
        b=OneOf(lambda d: set(d.aud["numerical_convergence"][1]["h2_uptake_by_steps"])),
        c=lambda d: d.conv_max_rel_pct())
    chk("``Add [[a]] mm of water'' is therefore", a=lambda d: d.v6["treatment_versions"]["total_water_mm"])
    chk("each add [[a]] mm of water and for the untreated control, across [[b]] parameter members",
        a=lambda d: d.v6["treatment_versions"]["total_water_mm"], b=lambda d: d.v6["treatment_versions"]["members"])

    # ---- methods: evidence inputs
    chk("the [[a]] core-site reference coordinates", a=lambda d: d.pc["core_sites"])
    chk("molecular analysis starts from [[a]] metagenomes", a=lambda d: d.pc["samples"])
    chk("with a [[a]] identity floor", a=lambda d: d.cl("M4", r"(\d+)-percent identity"))
    chk("implements the [[a]] conjunctive requirements of KEGG", a=lambda d: d.cl("M2", r"the (\d+) KEGG M00165"))
    chk("requires completeness at least [[a]] and contamination below [[b]], without",
        a=lambda d: d.cl("M5", r"at least (\d+) percent"), b=lambda d: d.cl("M5", r"below (\d+) percent"))
    # ---- methods: process model
    chk("The day is divided into $n=[[a]]$ substeps", a=lambda d: d.defaults["integration_steps_per_day"])
    chk("We compare [[a]], [[b]] and [[c]] substeps for the four reference schedules",
        a=OneOf(lambda d: set(d.aud["numerical_convergence"][1]["h2_uptake_by_steps"])),
        b=OneOf(lambda d: set(d.aud["numerical_convergence"][1]["h2_uptake_by_steps"])),
        c=OneOf(lambda d: set(d.aud["numerical_convergence"][1]["h2_uptake_by_steps"])))
    chk("dilute domain, default [[a]] mol L$^{-1}$", a=lambda d: d.defaults["osmotic_domain_mol_l"])
    chk("Scenario comparison crosses $b\\in\\{[[a]],[[b]]\\}$, $\\psi_d\\in\\{-[[c]],-[[d]]\\}$ MPa, $a\\in\\{[[e]],[[f]]\\}$ and kinetic multipliers $\\{[[g]],[[h]],[[i]]\\}$",
        a=OneOf(lambda d: set(d.tvt.retention_b)), b=OneOf(lambda d: set(d.tvt.retention_b)),
        c=OneOf(lambda d: set(-d.tvt.dry_threshold_mpa)), d=OneOf(lambda d: set(-d.tvt.dry_threshold_mpa)),
        e=OneOf(lambda d: set(d.tvt.activation_shape)), f=OneOf(lambda d: set(d.tvt.activation_shape)),
        g=OneOf(lambda d: set(d.tvt.k_h2_s / 1e-3)), h=OneOf(lambda d: set(d.tvt.k_h2_s / 1e-3)),
        i=OneOf(lambda d: set(d.tvt.k_h2_s / 1e-3)))
    chk("The four schedules use a [[a]]-day horizon, initial water [[b]] mm, temperature $[[c]] ^\\circ$C and evaporation [[e]] mm day$^{-1}$",
        a=lambda d: _src_num(d, "analyse.py", r"for day in range\((\d+)\)"),
        b=lambda d: _src_num(d, "analyse.py", r"State\(([\d.]+)\), s\[\"forcing\"\]"),
        c=lambda d: _src_num(d, "analyse.py", r"\"temp_c\": ([\d.]+)"),
        e=lambda d: _src_num(d, "analyse.py", r"\"evaporation_mm\": ([\d.]+)"))
    chk("Three schedules add [[a]] mm as one pulse, four pulses or daily pulses", a=lambda d: d.v6["treatment_versions"]["total_water_mm"])
    chk("Pulses fall on days [[a]], [[b]], [[c]] and [[e]] for the four-pulse schedule and on day [[f]] for the single pulse",
        a=OneOf(lambda d: set(d.four_pulse_days())), b=OneOf(lambda d: set(d.four_pulse_days())),
        c=OneOf(lambda d: set(d.four_pulse_days())), e=OneOf(lambda d: set(d.four_pulse_days())),
        f=lambda d: 0 if _src_num(d, "analyse.py", r"\"water_mm\": 12\. if d == (\d)") == 0 else None)
    chk("prescribed disturbed area, set to [[a]] m$^2$ for each wetting schedule",
        a=lambda d: _src_num(d, "analyse.py", r"0\. if n == 0 else ([\d.]+)"))
    chk("Ridge regression with penalty [[a]] predicts", a=lambda d: d.cl("F2", r"alpha=(\d+)"))
    chk("Nominal [[a]] intervals use absolute residuals from up to five site-held-out training folds",
        a=lambda d: d.cl("F3", r"(\d+)-percent"))
    chk("using [[a]] percentile-bootstrap draws and seed [[b]]", a=lambda d: d.aud["bootstrap_draws"], b=lambda d: d.aud["seed"])
    # ---- methods: external datasets
    chk("at three moisture levels ([[a]] bottles)", a=lambda d: d.ext["sonoran"]["lab_bottles"])
    chk("sensitivity values [[a]], [[b]] and [[c]] g cm$^{-3}$", a=OneOf(lambda d: set(d.nam_p.assumed_bulk_density_g_cm3)),
        b=OneOf(lambda d: set(d.nam_p.assumed_bulk_density_g_cm3)), c=OneOf(lambda d: set(d.nam_p.assumed_bulk_density_g_cm3)))
    chk("are evaluated against [[a]] plot--campaign means", a=lambda d: d.ext["sonoran"]["field_plot_campaign_means"])
    chk("additionally searches porosity $\\{[[a]],[[b]],[[c]]\\}$, retention exponent $\\{[[e]],[[f]],[[g]]\\}$, entry potentials of $[[h]]$, $[[i]]$ and $[[j]]$ MPa and $Q_{10}\\in\\{[[k]],[[l]]\\}$",
        a=lambda d: d.grid()[0][0], b=lambda d: d.grid()[0][1], c=lambda d: d.grid()[0][2],
        e=lambda d: d.grid()[1][0], f=lambda d: d.grid()[1][1], g=lambda d: d.grid()[1][2],
        h=lambda d: d.grid()[2][0], i=lambda d: d.grid()[2][1], j=lambda d: d.grid()[2][2],
        k=lambda d: d.grid()[3][0], l=lambda d: d.grid()[3][1])
    chk("The kinetic search is bounded between $[[a]]$ and $[[b]]$ s$^{-1}$",
        a=lambda d: 10.0 ** _src_tuple(d, "external_transfer.py", r"bounds=\((-?\d+),(-?\d+)\)")[0],
        b=lambda d: 10.0 ** _src_tuple(d, "external_transfer.py", r"bounds=\((-?\d+),(-?\d+)\)")[1])
    chk("provide [[a]] wet and ten dry atmospheric oxidation rates", a=lambda d: int(((d.nam_o.state == "wet")).sum()))
    chk("to each of [[a]] dilution-corrected bottle time series", a=lambda d: d.res["external_bottles"])
    chk("at $\\psi_d=-[[a]]$ MPa", a=lambda d: -d.defaults["dry_threshold_mpa"])

    # ---- supplement
    chk("Gas uptake integrates midpoint rates over [[a]] within-day steps by default", a=lambda d: d.defaults["integration_steps_per_day"])
    chk("compares [[a]], [[b]] and [[c]] steps for the four reference schedules",
        a=OneOf(lambda d: set(d.aud["numerical_convergence"][1]["h2_uptake_by_steps"])),
        b=OneOf(lambda d: set(d.aud["numerical_convergence"][1]["h2_uptake_by_steps"])),
        c=OneOf(lambda d: set(d.aud["numerical_convergence"][1]["h2_uptake_by_steps"])))
    chk("across the [[a]]-member grid at [[b]] and [[c]] steps", a=lambda d: d.v6["treatment_versions"]["members"],
        b=OneOf(lambda d: set(d.aud["numerical_convergence"][1]["h2_uptake_by_steps"])),
        c=OneOf(lambda d: set(d.aud["numerical_convergence"][1]["h2_uptake_by_steps"])))
    chk("module support as a conjunction of [[a]] requirements", a=lambda d: d.cl("M2", r"the (\d+) KEGG M00165"))
    chk("for all [[a]] hydrogenase-plus-CBB candidates and all [[b]] CO-carrier MAGs",
        a=lambda d: d.jus["counts"]["H2AutotrophicCandidateMAG"], b=lambda d: d.jus["counts"]["COCarrierMAG"])
    chk("maintenance divided by [[a]], and $[[b]]$. Constraint", a=lambda d: d.cl("C7", r"(\d+) and \d+ kJ/mol") * 1000,
        b=lambda d: 1e-20 if "1e-20" in d.src["src/eq_biodt/v03/resources.py"] else None)
    chk("checked at $[[a]]$ times the flux scale and resource violations at $[[b]]$ times that scale",
        a=lambda d: 1e-10 if "1e-10*scale" in d.src["src/eq_biodt/v03/resources.py"] else None,
        b=lambda d: 1e-7 if "1e-7*scale" in d.src["src/eq_biodt/v03/resources.py"] else None)
    chk("The released cohort contains [[a]] metagenomes at [[b]] sites, including [[c]] of the [[g]] core transect sites and two auxiliary locations (sites [[e]] and [[f]]). Core site [[h]] belongs",
        a=lambda d: d.pc["samples"], b=lambda d: d.pc["sites"], c=lambda d: d.pc["core_sites"],
        g=lambda d: len(_ledger_core_sites()), h=lambda d: _missing_core_sites(d)[0],
        e=OneOf(lambda d: set(d.pc["auxiliary_sites"])), f=OneOf(lambda d: set(d.pc["auxiliary_sites"])))
    chk("[[a]] of the [[b]] core sites on the repeated Rub' al-Khali transect", a=lambda d: d.pc["core_sites"],
        b=lambda d: len(_ledger_core_sites()))
    chk("covers all [[a]] core sites; no metagenome is available from site [[b]]", a=lambda d: len(_ledger_core_sites()),
        b=lambda d: _missing_core_sites(d)[0])
    chk("the Rub' al-Khali entry summarizes [[a]] samples", a=lambda d: d.pc["usable_mapping_samples"])
    chk("holds for every sample and in every regime; the release checks nine such inequalities over all [[a]] mapped samples",
        a=lambda d: d.psum["samples"])
    chk("with nominal [[a]] interval coverage for ridge regression", a=lambda d: d.cl("F3", r"(\d+)-percent"))
    chk("The Sonoran field grid contains [[a]] combinations", a=lambda d: len(list(product(*d.grid()))))
    chk("have porosity [[a]], retention exponent [[b]] and $Q_{10}=[[c]]$; entry potential is $[[e]]$ MPa when June is held out and $[[f]]$ MPa in the other folds",
        a=lambda d: _unique(d.fit_values("phi")), b=lambda d: _unique(d.fit_values("b")), c=lambda d: _unique(d.fit_values("q10")),
        e=lambda d: _unique([f["tuned"]["entry"] for f in d.ext["sonoran"]["fits"] if f["heldout_campaign"].startswith("2018-06")]),
        f=lambda d: _unique([f["tuned"]["entry"] for f in d.ext["sonoran"]["fits"] if not f["heldout_campaign"].startswith("2018-06")]))
    chk("At [[a]] and [[b]] g cm$^{-3}$ the reference model predicts higher hydrated than dry uptake",
        a=OneOf(lambda d: set(d.nam_p.assumed_bulk_density_g_cm3)), b=OneOf(lambda d: set(d.nam_p.assumed_bulk_density_g_cm3)))
    chk("At [[a]] g cm$^{-3}$, the hydrated C14-2 soil (gravimetric water content [[b]] g g$^{-1}$) reaches a volumetric water content of [[c]] against porosity [[e]]; the predicted dry flux ([[f]]) then exceeds the hydrated flux ([[g]], both",
        a=OneOf(lambda d: set(d.nam_p.assumed_bulk_density_g_cm3)),
        b=lambda d: _namib_row(d).wet_gwc, c=lambda d: _namib_row(d).wet_gwc * _namib_row(d).assumed_bulk_density_g_cm3,
        e=lambda d: d.defaults["porosity"], f=lambda d: _namib_row(d).dry_reference_flux, g=lambda d: _namib_row(d).wet_reference_flux)
    chk("recomputes all [[a]] Sonoran field predictions", a=lambda d: len(d.son_p))


def _unique(vals):
    vals = {round(float(v), 9) for v in vals}
    assert len(vals) == 1, f"expected a single value, got {sorted(vals)}"
    return vals.pop()


def _namib_row(d):
    r = d.nam_p[(d.nam_p.site == "C14-2") & (d.nam_p.assumed_bulk_density_g_cm3 == 1.8)]
    return r.iloc[0]


def _src_num(d, fname, rx):
    path = {"analyse.py": "reproducibility/v06/analyse.py", "external_transfer.py": "reproducibility/v05/external_transfer.py"}[fname]
    m = re.search(rx, d.src[path])
    assert m, f"{path}: {rx!r} not found"
    return float(m.group(1))


def _src_tuple(d, fname, rx):
    path = {"analyse.py": "reproducibility/v06/analyse.py", "external_transfer.py": "reproducibility/v05/external_transfer.py"}[fname]
    m = re.search(rx, d.src[path])
    assert m, f"{path}: {rx!r} not found"
    return tuple(int(g) for g in m.groups())


# ----------------------------------------------------------------------------------------------- allow-list
def build_allow():
    # Literals that are definitions, arithmetic constants or fixed inputs of the cited experiments or the sampling programme.
    allow("propagates zeros: $q_s(C_g)=[[a]]$ entails $\\gamma_g=[[b]]$, $k_g=[[c]]$ and $J_g=[[d]]$, and $\\psi\\leq\\psi_d$ entails $A=[[e]]$",
          "definition: zero propagation rules of the deductive closure (identification.closure)")
    allow("licenses the deletion of the intervention (Rule [[a]])", "name of do-calculus rule 3, not a quantity")
    allow("site-bootstrap [[a]] interval", "interval level; key site_bootstrap_95_percentile in the audit outputs",
          lambda d: all("site_bootstrap_95_percentile" in c for c in d.aud["matched_compartment_contrasts"]))
    allow("within-site bootstrap [[a]] interval", "interval level; percentile bootstrap, see reviews/v05/namib_hydration_contrasts.tsv")
    allow("Endpoint & Contrast & Pairs & Mean difference & [[a]] interval", "table header: interval level; percentile bootstrap, see audit.json",
          lambda d: True)
    allow("Ridge-minus-baseline RMSE differences with site-bootstrap [[a]] intervals", "interval level; percentile bootstrap, see audit.json",
          lambda d: True)
    allow("In Namib Desert microcosms, [[a]] $\\mu$L of added water", "design input of the cited Namib experiment (Tribbia et al.), see reviews/v05/external-evaluation-protocol.md")
    allow("over [[a]] days, using a shared initial state", "design input, see process/contract: 14-day horizon fixed in reproducibility/v06/analyse.py",
          lambda d: _src_num(d, "analyse.py", r"for day in range\((\d+)\)") == 14)
    allow("bulk subsurface soil at [[a]]--[[b]] cm", "sampling protocol of the primary programme; not an analysis output")
    allow("percentages to $[[a]]\\pm[[b]]%$", "closure tolerance of the evidence builder, see src/eq_biodt/v03/evidence.py")
    allow("and [[a]] g soil mass", "design input of the cited Namib experiment (Tribbia et al.), see reviews/v05/external-evaluation-protocol.md")
    allow("Hydration contrasts match the source batch, and uncertainty resamples the four source batches within each plot ([[a]] draws)",
          "design input, see process/contract: bootstrap draws in reproducibility/v05/external_transfer.py",
          lambda d: "(5000," in d.src["reproducibility/v05/external_transfer.py"].replace(" ", ""))
    allow("translates uniquely to $\\mathrm{do}(W=[[a]]L\\theta^*)$", "definition: water content to storage, theta = W/(1000 L)")
    allow("The reference scenarios use $\\gamma_g=[[a]]$", "definition: reference community factor of the model (Methods eq. for gamma)")
    allow("Third, if $\\gamma_g=[[a]]$ for a community", "definition: closed deductive gate, see Methods")
    allow("sets $v_{CO}=[[a]]$ without CO-carrier entailment", "definition of the admission gate, see Methods and src/eq_biodt/v03/resources.py")
    allow("If $\\gamma_g=[[a]]$, then $k_g=[[b]]$ and $J_g=[[c]]$ in every regime", "definition: consequences of the closed deductive gate, see Methods")
    allow("need not sum to [[a]]", "arithmetic constant (percentages)")


# ----------------------------------------------------------------------------------------------- macros
def macro_expectations():
    """Independent recomputation of every value macro that the text uses (name -> callable giving the value)."""
    E = {}
    E["ReferenceMAGs"] = lambda d: d.res["reference_MAGs"]
    E["CBBFull"] = lambda d: d.ont["entailed_class_counts"]["CBBCompleteMAG"]
    E["HtwoCBB"] = lambda d: d.ont["entailed_class_counts"]["H2AutotrophicCandidateMAG"]
    E["HtwoCBBStrict"] = lambda d: d.ont["entailed_class_counts"]["HighQualityH2AutotrophicCandidateMAG"]
    E["DualHtwoCO"] = lambda d: d.ont["entailed_class_counts"]["DualH2COAutotrophicCandidateMAG"]
    E["AmbiguousSamples"] = lambda d: d.pc["ambiguous_marginals"]
    E["TestVisits"] = lambda d: d.res["test_visits"]
    E["ConformalCoverage"] = lambda d: d.res["conformal_90_coverage_ridge"] * 100
    E["MappedSamples"] = lambda d: d.pc["usable_mapping_samples"]
    E["ClosedCOGates"] = lambda d: int((d.ps.q_co == 0).sum())
    E["ClosedCOZeroEffects"] = lambda d: int((d.ps.q_co == 0).sum()) * 3 * 2   # samples x schedules x couplings
    E["VersionMembers"] = lambda d: len(d.tvt)
    E["VersionSingleFirst"] = lambda d: int((d.tvt.uptake_single_pulse > d.tvt.uptake_four_pulses).sum())
    E["VersionFourFirst"] = lambda d: int((d.tvt.uptake_four_pulses > d.tvt.uptake_single_pulse).sum())
    E["VersionEffectSingleLow"] = lambda d: (d.tvt.uptake_single_pulse - d.tvt.control).min() * UM
    E["VersionEffectSingleHigh"] = lambda d: (d.tvt.uptake_single_pulse - d.tvt.control).max() * UM
    E["RefHtwoSingle"] = lambda d: d.nc("single_pulse") * UM
    E["RefHtwoFour"] = lambda d: d.nc("four_pulses") * UM
    for comp, name in (("Surface", "Surface"), ("Deep", "Deep"), ("PlantRoot", "Root")):
        E[f"HtwoPropSingleMedian{name}"] = (lambda c: lambda d: d.ps[d.ps.compartment == c].h2eff_proportional_single_pulse.median() * UM)(comp)
    E["GraphNodes"] = lambda d: len(d.adm)
    E["AdmissibleNodes"] = lambda d: int((d.adm.verdict == "admissible").sum())
    E["AdmissibleInputs"] = lambda d: int(((d.adm.verdict == "admissible") & (d.adm.zone == "exogenous_manipulable")).sum())
    E["AdmissibleStates"] = lambda d: int(((d.adm.verdict == "admissible") & (d.adm.zone == "causal")).sum())
    # matched contrast, proportional coupling, single pulse: Deep minus root
    E["ContrastDeepRootEst"] = lambda d: d.matched_contrast("Deep", "PlantRoot", "h2eff_proportional_single_pulse")[0] * UM
    E["ContrastDeepRootVisits"] = lambda d: d.matched_contrast("Deep", "PlantRoot", "h2eff_proportional_single_pulse")[1]
    ci = lambda d: next(c for c in d.v6["paired"]["contrasts"] if c["coupling"] == "proportional" and c["contrast"] == "Deep minus PlantRoot")
    E["ContrastDeepRootLow"] = lambda d: ci(d)["site_bootstrap_95_percentile"][0] * UM
    E["ContrastDeepRootHigh"] = lambda d: ci(d)["site_bootstrap_95_percentile"][1] * UM
    # version relevance demonstration
    def vt(d, psi):
        return next(t for t in d.vr["targets"] if t["psi_mpa"] == psi)
    for tag, psi in (("Strong", -1.0), ("Mild", -0.4)):
        for ver, vt_ in (("matric", "Matric"), ("osmotic", "Osmotic")):
            E[f"Psi{tag}{vt_}A"] = (lambda p, v: lambda d: vt(d, p)["versions"][v]["A"])(psi, ver)
            E[f"Psi{tag}{vt_}JHtwo"] = (lambda p, v: lambda d: vt(d, p)["versions"][v]["J_h2_mol_m2_s"] * NM)(psi, ver)
            E[f"Psi{tag}{vt_}Water"] = (lambda p, v: lambda d: vt(d, p)["versions"][v]["water_mm"])(psi, ver)
            E[f"Psi{tag}{vt_}Salt"] = (lambda p, v: lambda d: vt(d, p)["versions"][v]["salt_mol_m2"] / vt(d, p)["versions"][v]["water_mm"])(psi, ver)
    E["PsiMildUptakeMatric"] = lambda d: vt(d, -0.4)["versions"]["matric"]["h2_uptake_3d_mol_m2"] * UM
    E["PsiMildUptakeOsmotic"] = lambda d: vt(d, -0.4)["versions"]["osmotic"]["h2_uptake_3d_mol_m2"] * UM
    return {k: v for k, v in E.items() if v is not None}


def table_checks(d):
    """Generated row tables against their sources (the text cites these through \\ContrastRows etc.)."""
    out, bad = 0, []

    def rows(name):
        return [[c.strip() for c in ln.rstrip("\\ \n").split("&")] for ln in (PAPER / "tables" / name).read_text().splitlines() if "&" in ln]

    def cmp(label, lit, value):
        nonlocal out
        out += 1
        if not literal_matches(lit, value):
            bad.append(f"{label}: table has {lit}, source gives {value}")

    for r in rows("v04_contrasts.tex"):
        ep = {"CBB": "cbb_all_steps_detected_share_mapped"}.get(r[0], None) or "h2_cbb_all_steps_share_mapped"
        c = d.aud_contrast(ep, r[1].replace("Root", "PlantRoot"))
        cmp(f"contrast {r[0]} {r[1]} pairs", r[2], c["matched_site_campaigns"])
        cmp(f"contrast {r[0]} {r[1]} mean", r[3], c["mean_difference_percentage_points"])
        lo, hi = re.findall(NUMP, r[4])
        cmp(f"contrast {r[0]} {r[1]} lo", lo, c["site_bootstrap_95_percentile"][0])
        cmp(f"contrast {r[0]} {r[1]} hi", hi, c["site_bootstrap_95_percentile"][1])
    names = {"Frozen reference": "frozen_reference", "Kinetic coefficient calibrated": "calibrated_kinetics",
             "Diffusion-only coefficient calibrated": "calibrated_diffusion_only", "Training mean": "training_mean",
             "Linear moisture and temperature": "linear_moisture_temperature", "Soil-response parameters calibrated": "calibrated_retention"}
    for r in rows("v05_external.tex"):
        m = d.sonoran_model(names[r[0]])
        cmp(f"external {r[0]} RMSE", r[1], m["rmse"])
        cmp(f"external {r[0]} MAE", r[2], m["mae"])
    for r in rows("v05_campaigns.tex"):
        c = d.camp[d.camp.iloc[:, 0] == int(r[0])].iloc[0]
        cmp(f"campaign {r[0]} samples", r[1], c.iloc[1])
        cmp(f"campaign {r[0]} sites", r[2], c.iloc[2])
    base = {"persistence": "persistence", "compartment mean": "compartment_mean", "global mean": "global_mean"}
    for r in rows("v04_paired.tex"):
        c = next(c for c in d.aud["forecast_paired_comparisons"] if c["baseline"] == base[r[0]])
        cmp(f"paired {r[0]} diff", r[1], c["rmse_difference_ridge_minus_baseline"])
        lo, hi = re.findall(NUMP, r[2])
        cmp(f"paired {r[0]} lo", lo, c["site_bootstrap_95_percentile"][0])
        cmp(f"paired {r[0]} hi", hi, c["site_bootstrap_95_percentile"][1])
    adm = rows("v06_admissibility.tex")
    out += 1
    if len(adm) != len(d.adm):
        bad.append(f"admissibility table has {len(adm)} rows, admissibility.tsv has {len(d.adm)}")
    return out, bad


# ----------------------------------------------------------------------------------------------- driver
def macro_defs():
    defs = {}
    for p in sorted((PAPER / "tables").glob("*.tex")):
        for m in re.finditer(r"\\newcommand\{\\([A-Za-z]+)\}\{([^\n]*)\}\s*$", p.read_text(), re.M):
            defs[m.group(1)] = m.group(2)
    return defs


BUILTIN_MACROS = {"FloatBarrier", "Bigl", "Bigr", "Phi"}


def main():
    global PAPER
    ap = argparse.ArgumentParser()
    ap.add_argument("-v", "--verbose", action="store_true")
    ap.add_argument("--paper", default=str(PAPER), help="directory with main.tex, supplementary.tex and tables/ (default: paper/)")
    ap.add_argument("--list-literals", action="store_true", help="print every numeric literal with its context and exit")
    args = ap.parse_args()
    PAPER = Path(args.paper)

    texts, masks = {}, {}
    for f in FILES:
        texts[f] = normalise((PAPER / f).read_text())
        masks[f] = mask(texts[f])
    if args.list_literals:
        for f in FILES:
            for m in NUM.finditer(masks[f]):
                print(f, m.group(), "|", texts[f][max(0, m.start() - 40):m.end() + 25].replace("\n", " "))
        return 0

    MACRO_VALUES.update(macro_defs())
    d = D()
    build_registry()
    build_allow()
    errors, verified, allowed = [], [], []
    covered = {f: [] for f in FILES}
    REGISTRY[:] = [r for r in REGISTRY if r[2] != 0]

    for template, exprs, n in REGISTRY:
        rx = template_regex(template)
        hits = [(f, m) for f in FILES for m in rx.finditer(texts[f])]
        if len(hits) < n:
            errors.append(f"STALE registry entry (no match in the text): {template!r}")
            continue
        names = set(rx.groupindex)
        if names != set(exprs):
            errors.append(f"registry entry {template!r} has groups {sorted(names)} but expressions {sorted(exprs)}")
            continue
        for f, m in hits:
            for g, expr in exprs.items():
                raw_lit = m.group(g)
                lit = resolve(raw_lit)
                covered[f].append(m.span(g))
                try:
                    exp = expr.fn(d) if isinstance(expr, OneOf) else expr(d)
                except Exception as e:  # noqa: BLE001
                    errors.append(f"{f}: cannot compute expected value for {lit!r} in {template!r}: {e!r}")
                    continue
                if isinstance(expr, OneOf):
                    ok = any(literal_matches(lit, v) for v in exp)
                    expd = f"one of {sorted(exp)}"
                else:
                    ok = literal_matches(lit, exp)
                    expd = f"{exp!r} rounded to {dec_places(lit)} decimals"
                (verified if ok else errors).append(
                    f"{f}: {raw_lit}  [{template[:70]}...]" if ok else
                    f"MISMATCH {f}: text says {raw_lit} ({lit}), source gives {expd}  [context: {template[:90]}]")

    for template, reason, check, n in ALLOW:
        rx = template_regex(template)
        hits = [(f, m) for f in FILES for m in rx.finditer(texts[f])]
        if n == 0:
            continue
        if len(hits) < n:
            errors.append(f"STALE allow-list entry (no match in the text): {template!r}")
            continue
        if check is not None and not check(d):
            errors.append(f"allow-list cross-check failed for {template!r} ({reason})")
        for f, m in hits:
            for g in rx.groupindex:
                covered[f].append(m.span(g))
                allowed.append(f"{f}: {m.group(g)}  [{reason}]")

    # literals covered by neither list
    for f in FILES:
        for m in NUM.finditer(masks[f]):
            if not any(a <= m.start() and m.end() <= b for a, b in covered[f]):
                errors.append(f"UNREGISTERED literal {m.group()!r} in {f}: ...{texts[f][max(0, m.start() - 50):m.end() + 30]!r}".replace("\\n", " "))

    # macros: defined, and value macros tied to sources
    defs = macro_defs()
    preamble_defs = set()
    for f in FILES:
        raw = (PAPER / f).read_text()
        preamble_defs |= set(re.findall(r"\\(?:re)?newcommand\{\\([A-Za-z]+)\}", raw))
    used = set()
    for f in FILES:
        used |= macros_used((PAPER / f).read_text())
    custom = {m for m in used if m[0].isupper() and m not in BUILTIN_MACROS}
    for m in sorted(custom):
        if m not in defs and m not in preamble_defs:
            errors.append(f"UNDEFINED macro \\{m} (not defined in paper/tables/*.tex or the preamble)")
    for f in FILES:
        for inc in re.findall(r"\\input\{([^}]+)\}", (PAPER / f).read_text()):
            if not (PAPER / (inc if inc.endswith(".tex") else inc + ".tex")).exists():
                errors.append(f"{f}: \\input{{{inc}}} does not exist")
    expect = macro_expectations()
    macro_ok = 0
    for m in sorted(custom & set(defs)):
        if m not in expect:
            errors.append(f"macro \\{m} is used in the text but has no independent source computation in MACROS")
            continue
        try:
            val = expect[m](d)
        except Exception as e:  # noqa: BLE001
            errors.append(f"macro \\{m}: cannot compute expected value: {e!r}")
            continue
        lit = defs[m].replace("\\%", "%")
        if literal_matches(lit, val):
            macro_ok += 1
        else:
            errors.append(f"MISMATCH macro \\{m}: tables/*.tex defines {lit}, source gives {val!r}")

    tcount, tbad = table_checks(d)
    errors += [f"MISMATCH generated table {b}" for b in tbad]

    if args.verbose:
        print("\n".join("ok " + v for v in verified))
        print("\n".join("allowed " + a for a in allowed))
    print(f"prose/caption literals verified against outputs: {len(verified)}")
    print(f"literals allow-listed (design inputs, definitions, constants): {len(allowed)}")
    print(f"value macros used in the text recomputed from outputs: {macro_ok} of {len(custom & set(defs))}")
    print(f"generated table cells cross-checked: {tcount}")
    if errors:
        print(f"\nFAILED: {len(errors)} problem(s)")
        for e in errors:
            print(" -", e)
        return 1
    print("manuscript numbers verified")
    return 0


if __name__ == "__main__":
    sys.exit(main())
