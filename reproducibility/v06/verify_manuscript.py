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
            "reproducibility/v06/analyse.py", "reproducibility/v05/external_transfer.py", "src/eq_biodt/v03/resources.py",
            "reproducibility/v07/digitise_fig1.py", "reproducibility/v07/composition.py")}
        self.v7_sum = J("reviews/v07/composition/summary.json")
        self.v7 = {r["id"]: J(f"reviews/v07/composition/{r['id']}.json") for r in self.v7_sum["experiments"]}
        self.v7_reg = {c["id"]: c for c in J("data/external/composition_constants.json")["constants"]}
        self.v7_fig1 = J("data/external/de_zelicourt_2018/fig1_digitised.json")

    def v7_lvl(self, exp, **kw):
        return next(l for l in self.v7[exp]["levels"] if all(l.get(k, l.get("action", {}).get(k)) == v for k, v in kw.items()))

    def v7_head(self, exp):
        return self.v7[exp]["headline"][0]

    def v7_const(self, cid, field="value"):
        return self.v7_reg[cid][field]

    def v7_status_count(self, status):
        return sum(1 for c in self.v7_reg.values() if c["status"] == status)

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

    # ---- v0.7 composition: main.tex Results & Discussion
    chk("across [[a]] chained experiments on the [[b]]-day single-pulse baseline ([[c]] $\\mu$mol H$_2$ m$^{-2}$ and [[e]] $\\mu$mol CO m$^{-2}$",
        a=lambda d: len(d.v7_sum["experiments"]),
        b=lambda d: d.v7["drip_vs_flood"]["baseline"]["period_days"],
        c=lambda d: d.v7["drip_vs_flood"]["levels"][0]["twin_baseline"]["h2_uptake_mol_m2"] * UM,
        e=lambda d: d.v7["drip_vs_flood"]["levels"][0]["twin_baseline"]["co_uptake_mol_m2"] * UM)
    chk("Adding irrigation water ([[a]] mm day$^{-1}$ seasonal demand",
        a=lambda d: d.v7_const("crop_etc_mm_per_day"))
    chk("Daily drip irrigation changes [[a]]-day H$_2$ uptake by $[[b]]$ ($[[c]]$ to $[[e]]$), compared with $[[f]]$ for [[g]]-day pulses and $[[h]]$ for [[i]]-day flood pulses",
        a=lambda d: d.v7["drip_vs_flood"]["baseline"]["period_days"],
        b=lambda d: d.v7_lvl("drip_vs_flood", irrigation_interval_days=1)["difference_vs_baseline"]["h2_uptake_mol_m2"]["percent"]["value"],
        c=lambda d: d.v7_lvl("drip_vs_flood", irrigation_interval_days=1)["difference_vs_baseline"]["h2_uptake_mol_m2"]["percent"]["low"],
        e=lambda d: d.v7_lvl("drip_vs_flood", irrigation_interval_days=1)["difference_vs_baseline"]["h2_uptake_mol_m2"]["percent"]["high"],
        f=lambda d: d.v7_lvl("drip_vs_flood", irrigation_interval_days=3)["difference_vs_baseline"]["h2_uptake_mol_m2"]["percent"]["value"],
        g=lambda d: d.v7_lvl("drip_vs_flood", irrigation_interval_days=3)["action"]["irrigation_interval_days"],
        h=lambda d: d.v7_lvl("drip_vs_flood", irrigation_interval_days=10)["difference_vs_baseline"]["h2_uptake_mol_m2"]["percent"]["value"],
        i=lambda d: d.v7_lvl("drip_vs_flood", irrigation_interval_days=10)["action"]["irrigation_interval_days"])
    chk("Under daily brackish irrigation from [[a]] to [[b]] dS m$^{-1}$, twin H$_2$ uptake declines by [[c]] to [[e]] as accumulated salt lowers osmotic potential (with high-salinity envelope runs at [[f]] dS m$^{-1}$ exceeding the [[g]] mol L$^{-1}$ dilute domain), while downstream Maas--Hoffman relative yield at [[h]] dS m$^{-1}$ is [[i]] ([[j]] to [[k]]) for alfalfa, [[l]] for barley, [[m]] ([[o]] to [[p]]) for date palm and [[q]] ([[r]] to [[s]]) for tomato",
        a=lambda d: d.v7_lvl("brackish_irrigation", ec_w_dS_m=0.5)["action"]["ec_w_dS_m"],
        b=lambda d: d.v7_lvl("brackish_irrigation", ec_w_dS_m=8)["action"]["ec_w_dS_m"],
        c=lambda d: -d.v7_lvl("brackish_irrigation", ec_w_dS_m=0.5)["difference_vs_baseline"]["h2_uptake_mol_m2"]["percent"]["value"],
        e=lambda d: -d.v7_lvl("brackish_irrigation", ec_w_dS_m=8)["difference_vs_baseline"]["h2_uptake_mol_m2"]["percent"]["value"],
        f=lambda d: d.v7_lvl("brackish_irrigation", ec_w_dS_m=8)["action"]["ec_w_dS_m"],
        g=lambda d: d.defaults["osmotic_domain_mol_l"],
        h=lambda d: d.v7_lvl("brackish_irrigation", ec_w_dS_m=8)["action"]["ec_w_dS_m"],
        i=lambda d: d.v7_lvl("brackish_irrigation", ec_w_dS_m=8)["downstream"]["maas_hoffman_alfalfa"]["relative_yield_alfalfa"]["value"],
        j=lambda d: d.v7_lvl("brackish_irrigation", ec_w_dS_m=8)["downstream"]["maas_hoffman_alfalfa"]["relative_yield_alfalfa"]["low"],
        k=lambda d: d.v7_lvl("brackish_irrigation", ec_w_dS_m=8)["downstream"]["maas_hoffman_alfalfa"]["relative_yield_alfalfa"]["high"],
        l=lambda d: d.v7_lvl("brackish_irrigation", ec_w_dS_m=8)["downstream"]["maas_hoffman_barley"]["relative_yield_barley"]["value"],
        m=lambda d: d.v7_lvl("brackish_irrigation", ec_w_dS_m=8)["downstream"]["maas_hoffman_date_palm"]["relative_yield_date_palm"]["value"],
        o=lambda d: d.v7_lvl("brackish_irrigation", ec_w_dS_m=8)["downstream"]["maas_hoffman_date_palm"]["relative_yield_date_palm"]["low"],
        p=lambda d: d.v7_lvl("brackish_irrigation", ec_w_dS_m=8)["downstream"]["maas_hoffman_date_palm"]["relative_yield_date_palm"]["high"],
        q=lambda d: d.v7_lvl("brackish_irrigation", ec_w_dS_m=8)["downstream"]["maas_hoffman_tomato"]["relative_yield_tomato"]["value"],
        r=lambda d: d.v7_lvl("brackish_irrigation", ec_w_dS_m=8)["downstream"]["maas_hoffman_tomato"]["relative_yield_tomato"]["low"],
        s=lambda d: d.v7_lvl("brackish_irrigation", ec_w_dS_m=8)["downstream"]["maas_hoffman_tomato"]["relative_yield_tomato"]["high"])
    chk("full organic mulch changes [[a]]-day H$_2$ uptake by [[b]], solar-panel cover with [[c]]-day washing by [[e]], cloud seeding (+[[f]] rain on rain days) by [[g]], +[[h]] $^\\circ$C warming at unchanged rainfall by [[i]], shade-net cooling ($-[[j]] ^\\circ$C) by $[[k]]$, [[l]] w/w biochar by $[[m]]$ and [[o]] off-road vehicle passes by $[[p]]$",
        a=lambda d: d.v7["mulch_shade"]["baseline"]["period_days"],
        b=lambda d: d.v7_head("mulch_shade")["value"],
        c=lambda d: d.v7_const("panel_wash_interval_days"),
        e=lambda d: d.v7_head("solar_park")["value"],
        f=lambda d: d.v7_const("cloud_seeding_rain_increase_frac") * 100,
        g=lambda d: d.v7_head("cloud_seeding")["value"],
        h=lambda d: d.v7_lvl("warming", warming_c=4, rain_change_frac=0.0)["action"]["warming_c"],
        i=lambda d: d.v7_head("warming")["value"],
        j=lambda d: d.v7_const("shade_net_soil_temp_reduction_c"),
        k=lambda d: d.v7_lvl("mulch_shade", mulch_cover_fraction=0, shade_net_cover_fraction=1)["difference_vs_baseline"]["h2_uptake_mol_m2"]["percent"]["value"],
        l=lambda d: d.v7_lvl("biochar", biochar_dose_pct_w_w=5)["action"]["biochar_dose_pct_w_w"],
        m=lambda d: d.v7_head("biochar")["value"],
        o=lambda d: d.v7_lvl("compaction", offroad_passes=50)["action"]["offroad_passes"],
        p=lambda d: d.v7_head("compaction")["value"])
    chk("Perturbing atmospheric H$_2$ above the [[a]] ppm background scales the [[b]]-day soil H$_2$ sink linearly across [[c]], [[e]] and [[f]] times background (+[[g]], +[[h]] and +[[i]] uptake, reaching [[j]] $\\mu$mol m$^{-2}$ at [[k]] ppm), because first-order reaction--diffusion is proportional to surface concentration and peak O$_2$ demand reaches only [[l]] of diffusive O$_2$ supply, leaving oxygen unlimiting and CO uptake unchanged at [[m]] $\\mu$mol m$^{-2}$",
        a=lambda d: d.v7_const("h2_background_ppm"),
        b=lambda d: d.v7["h2_leakage"]["baseline"]["period_days"],
        c=lambda d: d.v7_lvl("h2_leakage", h2_multiplier=1.5)["action"]["h2_multiplier"],
        e=lambda d: d.v7_lvl("h2_leakage", h2_multiplier=2)["action"]["h2_multiplier"],
        f=lambda d: d.v7_lvl("h2_leakage", h2_multiplier=10)["action"]["h2_multiplier"],
        g=lambda d: d.v7_lvl("h2_leakage", h2_multiplier=1.5)["difference_vs_baseline"]["h2_uptake_mol_m2"]["percent"]["value"],
        h=lambda d: d.v7_lvl("h2_leakage", h2_multiplier=2)["difference_vs_baseline"]["h2_uptake_mol_m2"]["percent"]["value"],
        i=lambda d: d.v7_lvl("h2_leakage", h2_multiplier=10)["difference_vs_baseline"]["h2_uptake_mol_m2"]["percent"]["value"],
        j=lambda d: d.v7_lvl("h2_leakage", h2_multiplier=10)["twin"]["h2_uptake_mol_m2"]["value"] * UM,
        k=lambda d: d.v7_lvl("h2_leakage", h2_multiplier=10)["upstream_ppm"]["value"],
        l=lambda d: d.v7_lvl("h2_leakage", h2_multiplier=10)["oxygen_used_over_supply_max"] * 100,
        m=lambda d: d.v7_lvl("h2_leakage", h2_multiplier=10)["twin"]["co_uptake_mol_m2"]["value"] * UM)
    chk("knallgas CBB stoichiometry ([[a]]/[[b]] mol CO$_2$ per mol H$_2$ ) under the deductive gate gives an annualised carbon-fixation upper bound in high-autotroph sample ([[c]] autotroph share) of [[e]] g C m$^{-2}$ yr$^{-1}$ (envelope [[f]] to [[g]] g C m$^{-2}$ yr$^{-1}$; [[h]] $\\mu$mol CO$_2$ m$^{-2}$ over [[i]] days), equivalent to [[j]] of a [[k]] g C m$^{-2}$ yr$^{-1}$ dryland biocrust reference flux",
        a=lambda d: 1,
        b=lambda d: 1.0 / d.v7_const("knallgas_co2_per_h2"),
        c=lambda d: d.v7_lvl("carbon_fixation_bound", sample="F14Dr1", scenario="baseline (12 mm pulse)")["h2_autotroph_share"] * 100,
        e=lambda d: d.v7_head("carbon_fixation_bound")["value"] / 1e6,
        f=lambda d: d.v7_head("carbon_fixation_bound")["low"] / 1e6,
        g=lambda d: d.v7_head("carbon_fixation_bound")["high"] / 1e6,
        h=lambda d: d.v7_lvl("carbon_fixation_bound", sample="F14Dr1", scenario="baseline (12 mm pulse)")["co2_fixed_upper_bound"]["per_period_umol_m2"]["value"],
        i=lambda d: d.v7["carbon_fixation_bound"]["baseline"]["period_days"],
        j=lambda d: d.v7_lvl("carbon_fixation_bound", sample="F14Dr1", scenario="baseline (12 mm pulse)")["bound_over_context"] * 100,
        k=lambda d: d.v7_const("biocrust_net_c_fixation_g_c_m2_yr"))
    chk("where the entailed share is [[a]], the downstream deductive gate returns [[b]] fixed carbon with provenance while twin H$_2$ oxidation continues ([[c]] $\\mu$mol m$^{-2}$ under the gate). Coupled CO oxidation converts CO to CO$_2$ in equimolar stoichiometry with [[e]] net change in atmospheric carbon atoms",
        a=lambda d: d.v7_lvl("carbon_fixation_bound", sample="14PRr2", scenario="baseline (12 mm pulse)")["h2_autotroph_share"] * 100,
        b=lambda d: d.v7_lvl("carbon_fixation_bound", sample="14PRr2", scenario="baseline (12 mm pulse)")["c_fixed_upper_bound"]["per_year_ug_c_m2"]["value"],
        c=lambda d: d.v7_lvl("carbon_fixation_bound", sample="14PRr2", scenario="baseline (12 mm pulse)")["twin"]["h2_uptake_mol_m2"]["value"] * UM,
        e=lambda d: d.v7_lvl("carbon_fixation_bound", sample="F14Dr1", scenario="baseline (12 mm pulse)")["co_oxidation"]["net_atmospheric_c_change_mol_m2"])
    chk("which has [[a]] mapped CO-carrier share and a closed CO gate ([[b]] $\\mu$mol m$^{-2}$ uptake), introducing CO-carrier MAG at [[c]], [[e]] or [[f]] share opens the gate and increases [[g]]-day CO uptake to [[h]] $\\mu$mol m$^{-2}$ under gate coupling, or to [[i]], [[j]] and [[k]] $\\mu$mol m$^{-2}$ under proportional coupling",
        a=lambda d: float(d.ps[d.ps["sample"] == "51PRr3"].q_co.iloc[0]) * 100,
        b=lambda d: d.v7["inoculation"]["gate"]["co_uptake_before_umol_m2"],
        c=lambda d: d.v7_lvl("inoculation", sample="51PRr3", coupling="gate", share=0.001)["share"] * 100,
        e=lambda d: d.v7_lvl("inoculation", sample="51PRr3", coupling="gate", share=0.01)["share"] * 100,
        f=lambda d: d.v7_lvl("inoculation", sample="51PRr3", coupling="gate", share=0.05)["share"] * 100,
        g=lambda d: d.v7["inoculation"]["baseline"]["period_days"],
        h=lambda d: d.v7["inoculation"]["gate"]["co_uptake_after_umol_m2"],
        i=lambda d: d.v7_lvl("inoculation", sample="51PRr3", coupling="proportional", share=0.001, mag_id="11Dr3__single_SemiBin_57")["twin"]["co_uptake_mol_m2"]["value"] * UM,
        j=lambda d: d.v7_lvl("inoculation", sample="51PRr3", coupling="proportional", share=0.01, mag_id="11Dr3__single_SemiBin_57")["twin"]["co_uptake_mol_m2"]["value"] * UM,
        k=lambda d: d.v7_lvl("inoculation", sample="51PRr3", coupling="proportional", share=0.05, mag_id="11Dr3__single_SemiBin_57")["twin"]["co_uptake_mol_m2"]["value"] * UM)
    chk("keeps CO uptake at [[a]] $\\mu$mol m$^{-2}$, and introducing the CO carrier into sample , whose gate is already open ([[b]] initial share), gives [[c]] change under gate coupling while raising proportional CO uptake from [[e]] to [[f]], [[g]] and [[h]] $\\mu$mol m$^{-2}$",
        a=lambda d: d.v7["inoculation"]["gate"]["negative_control_co_change_gate_umol_m2"],
        b=lambda d: float(d.ps[d.ps["sample"] == "V53PRr2"].q_co.iloc[0]) * 100,
        c=lambda d: d.v7["inoculation"]["gate"]["open_sample_co_change_gate_umol_m2"],
        e=lambda d: d.v7_lvl("inoculation", sample="V53PRr2", coupling="proportional", share=0.001)["difference_vs_baseline"]["co_uptake_mol_m2"]["baseline"] * UM,
        f=lambda d: d.v7_lvl("inoculation", sample="V53PRr2", coupling="proportional", share=0.001)["twin"]["co_uptake_mol_m2"]["value"] * UM,
        g=lambda d: d.v7_lvl("inoculation", sample="V53PRr2", coupling="proportional", share=0.01)["twin"]["co_uptake_mol_m2"]["value"] * UM,
        h=lambda d: d.v7_lvl("inoculation", sample="V53PRr2", coupling="proportional", share=0.05)["twin"]["co_uptake_mol_m2"]["value"] * UM)
    chk("Hada Al-Sham alfalfa field trial irrigated at [[a]] and [[b]] dS m$^{-1}$",
        n=2,
        a=lambda d: d.v7_const("sa187_ecw_low_dS_m"),
        b=lambda d: d.v7_const("sa187_ecw_high_dS_m"))
    chk("Digitised season-[[a]] mock high-to-low salinity yield ratios are [[b]] (digitisation envelope [[c]] to [[e]]) for fresh biomass and [[f]] ([[g]] to [[h]]) for dry biomass",
        a=lambda d: 1,
        b=lambda d: d.v7["inoculated_alfalfa_salinity"]["observed"]["fresh"]["value"],
        c=lambda d: d.v7["inoculated_alfalfa_salinity"]["observed"]["fresh"]["low"],
        e=lambda d: d.v7["inoculated_alfalfa_salinity"]["observed"]["fresh"]["high"],
        f=lambda d: d.v7["inoculated_alfalfa_salinity"]["observed"]["dry"]["value"],
        g=lambda d: d.v7["inoculated_alfalfa_salinity"]["observed"]["dry"]["low"],
        h=lambda d: d.v7["inoculated_alfalfa_salinity"]["observed"]["dry"]["high"])
    chk("standard fixed leaching factor ($f=[[a]]$, $\\mathrm{EC_e}=[[b]] \\mathrm{EC_w}$ ) with the Maas--Hoffman alfalfa response (threshold [[c]] dS m$^{-1}$, slope [[e]] dS$^{-1}$ m ) predicts a mock ratio of [[f]] (envelope [[g]] to [[h]]), whereas the twin's bare-sand salt balance predicts [[i]] ([[j]] to [[k]])",
        a=lambda d: d.v7_const("fao29_ece_per_ecw"),
        b=lambda d: d.v7_const("fao29_ece_per_ecw"),
        c=lambda d: d.v7_const("mh_alfalfa_threshold_dS_m"),
        e=lambda d: d.v7_const("mh_alfalfa_slope_pct_per_dS_m"),
        f=lambda d: d.v7["inoculated_alfalfa_salinity"]["validation_mock_ratio"]["A_fao_fixed_leaching"]["predicted_ratio"]["value"],
        g=lambda d: d.v7["inoculated_alfalfa_salinity"]["validation_mock_ratio"]["A_fao_fixed_leaching"]["predicted_ratio"]["low"],
        h=lambda d: d.v7["inoculated_alfalfa_salinity"]["validation_mock_ratio"]["A_fao_fixed_leaching"]["predicted_ratio"]["high"],
        i=lambda d: d.v7["inoculated_alfalfa_salinity"]["validation_mock_ratio"]["B_twin_physics"]["predicted_ratio"]["value"],
        j=lambda d: d.v7["inoculated_alfalfa_salinity"]["validation_mock_ratio"]["B_twin_physics"]["predicted_ratio"]["low"],
        k=lambda d: d.v7["inoculated_alfalfa_salinity"]["validation_mock_ratio"]["B_twin_physics"]["predicted_ratio"]["high"])
    chk("observed mock yield ratio gives $f=[[a]]$ (envelope [[b]] to [[c]]) for fresh biomass and $f=[[e]]$ ([[f]] to [[g]]) for dry biomass, implying a low-salinity root-zone $\\mathrm{EC_e}$ of [[h]] to [[i]] dS m$^{-1}$ ([[j]] dS m$^{-1}$ dry, [[k]] dS m$^{-1}$ fresh) compared with the measured site soil conductivity of [[l]] dS m$^{-1}$",
        a=lambda d: d.v7["inoculated_alfalfa_salinity"]["fitted_leaching"]["fresh"]["f_central"],
        b=lambda d: d.v7["inoculated_alfalfa_salinity"]["fitted_leaching"]["fresh"]["f_range_from_digitisation_error"][0],
        c=lambda d: d.v7["inoculated_alfalfa_salinity"]["fitted_leaching"]["fresh"]["f_range_from_digitisation_error"][1],
        e=lambda d: d.v7["inoculated_alfalfa_salinity"]["fitted_leaching"]["dry"]["f_central"],
        f=lambda d: d.v7["inoculated_alfalfa_salinity"]["fitted_leaching"]["dry"]["f_range_from_digitisation_error"][0],
        g=lambda d: d.v7["inoculated_alfalfa_salinity"]["fitted_leaching"]["dry"]["f_range_from_digitisation_error"][1],
        h=lambda d: d.v7["inoculated_alfalfa_salinity"]["fitted_leaching"]["dry"]["ece_at_3_12_dS_m"],
        i=lambda d: d.v7["inoculated_alfalfa_salinity"]["fitted_leaching"]["fresh"]["ece_at_3_12_dS_m"],
        j=lambda d: d.v7["inoculated_alfalfa_salinity"]["fitted_leaching"]["dry"]["ece_at_3_12_dS_m"],
        k=lambda d: d.v7["inoculated_alfalfa_salinity"]["fitted_leaching"]["fresh"]["ece_at_3_12_dS_m"],
        l=lambda d: d.v7_const("sa187_soil_ec_dS_m"))
    chk("With the fitted leaching factor and the [[a]] ([[b]] to [[c]]) SA187 yield multiplier, predicted inoculated relative yield at [[e]] dS m$^{-1}$ is [[f]] ([[g]] to [[h]], fresh biomass); from a single trial at two salinities, the SA187 benefit is identifiable only as a constant multiplier across the tested [[i]] to [[j]] dS m$^{-1}$ range",
        a=lambda d: d.v7_const("sa187_benefit_b") * 100,
        b=lambda d: d.v7_const("sa187_benefit_b", "low") * 100,
        c=lambda d: d.v7_const("sa187_benefit_b", "high") * 100,
        e=lambda d: d.v7_const("sa187_ecw_high_dS_m"),
        f=lambda d: d.v7_head("inoculated_alfalfa_salinity")["value"],
        g=lambda d: d.v7_head("inoculated_alfalfa_salinity")["low"],
        h=lambda d: d.v7_head("inoculated_alfalfa_salinity")["high"],
        i=lambda d: d.v7_const("sa187_ecw_low_dS_m"),
        j=lambda d: d.v7_const("sa187_ecw_high_dS_m"))
    chk("evaluated against a [[a]]-day single-pulse reference baseline, and [[b]] of the [[c]] external registry constants represent modeller assumptions ([[e]]) or values whose source could not be checked ([[f]])",
        a=lambda d: d.v7["drip_vs_flood"]["baseline"]["period_days"],
        b=lambda d: d.v7_status_count("modeller_assumption") + d.v7_status_count("source_not_checked"),
        c=lambda d: len(d.v7_reg),
        e=lambda d: d.v7_status_count("modeller_assumption"),
        f=lambda d: d.v7_status_count("source_not_checked"))

    # ---- v0.7 composition: supplementary.tex Section S8
    chk("evaluates [[a]] management what-if chains against the [[b]]-day single-pulse rain baseline ([[c]] mm rain on day [[e]], initial water storage [[f]] mm, soil temperature $[[g]] ^\\circ$C, potential evaporation [[h]] mm day$^{-1}$; reference uptake [[i]] $\\mu$mol H$_2$ m$^{-2}$ and [[j]] $\\mu$mol CO m$^{-2}$). Table and Fig. summarize all [[k]] experiments",
        a=lambda d: len(d.v7_sum["experiments"]),
        b=lambda d: d.v7["drip_vs_flood"]["baseline"]["period_days"],
        c=lambda d: d.v6["treatment_versions"]["total_water_mm"],
        e=lambda d: 0,
        f=lambda d: d.v7["drip_vs_flood"]["baseline"]["initial_water_mm"],
        g=lambda d: _src_num(d, "analyse.py", r"\"temp_c\": ([\d.]+)"),
        h=lambda d: _src_num(d, "analyse.py", r"\"evaporation_mm\": ([\d.]+)"),
        i=lambda d: d.v7["drip_vs_flood"]["levels"][0]["twin_baseline"]["h2_uptake_mol_m2"] * UM,
        j=lambda d: d.v7["drip_vs_flood"]["levels"][0]["twin_baseline"]["co_uptake_mol_m2"] * UM,
        k=lambda d: len(d.v7_sum["experiments"]))
    chk("Summary of the [[a]] management what-if compositions on the [[b]]-day single-pulse reference baseline",
        a=lambda d: len(d.v7_sum["experiments"]),
        b=lambda d: d.v7["drip_vs_flood"]["baseline"]["period_days"])
    chk("Small-multiples summary of the [[a]] management what-if compositions on the [[b]]-day single-pulse reference baseline",
        a=lambda d: len(d.v7_sum["experiments"]),
        b=lambda d: d.v7["drip_vs_flood"]["baseline"]["period_days"])
    chk("Upstream , and ([[a]] mm day$^{-1}$) add daily water",
        a=lambda d: d.v7_const("crop_etc_mm_per_day"))
    chk("with high-salinity constant-envelope runs at $\\mathrm{EC_w}=[[a]]$ dS m$^{-1}$ exceeding the [[b]] mol L$^{-1}$ dilute osmotic domain",
        a=lambda d: d.v7_lvl("brackish_irrigation", ec_w_dS_m=8)["action"]["ec_w_dS_m"],
        b=lambda d: d.defaults["osmotic_domain_mol_l"])
    chk("Upstream ([[a]] mm day$^{-1}$) and deliver equal seasonal water demand on [[b]]-, [[c]]- or [[e]]-day intervals on top of the [[f]] mm rain pulse, and the twin integrates [[g]]-day H$_2$ and CO uptake",
        a=lambda d: d.v7_const("crop_etc_mm_per_day"),
        b=lambda d: d.v7_lvl("drip_vs_flood", irrigation_interval_days=1)["action"]["irrigation_interval_days"],
        c=lambda d: d.v7_lvl("drip_vs_flood", irrigation_interval_days=3)["action"]["irrigation_interval_days"],
        e=lambda d: d.v7_lvl("drip_vs_flood", irrigation_interval_days=10)["action"]["irrigation_interval_days"],
        f=lambda d: d.v6["treatment_versions"]["total_water_mm"],
        g=lambda d: d.v7["drip_vs_flood"]["baseline"]["period_days"])
    chk("full mulch with shade increases [[a]]-day H$_2$ uptake by [[b]] ($[[c]]$ to [[e]])",
        a=lambda d: d.v7["mulch_shade"]["baseline"]["period_days"],
        b=lambda d: d.v7_lvl("mulch_shade", mulch_cover_fraction=1, shade_net_cover_fraction=1)["difference_vs_baseline"]["h2_uptake_mol_m2"]["percent"]["value"],
        c=lambda d: d.v7_lvl("mulch_shade", mulch_cover_fraction=1, shade_net_cover_fraction=1)["difference_vs_baseline"]["h2_uptake_mol_m2"]["percent"]["low"],
        e=lambda d: d.v7_lvl("mulch_shade", mulch_cover_fraction=1, shade_net_cover_fraction=1)["difference_vs_baseline"]["h2_uptake_mol_m2"]["percent"]["high"])
    chk("adds [[a]] to soil porosity per percentage point w/w across $\\{0,1,2,5\\}%$ biochar (raising $\\phi$ from [[b]] to [[c]] at [[e]] w/w)",
        a=lambda d: d.v7_const("biochar_porosity_increase_per_pct"),
        b=lambda d: d.defaults["porosity"],
        c=lambda d: d.v7_lvl("biochar", biochar_dose_pct_w_w=5)["upstream"]["biochar_porosity"]["porosity"]["value"],
        e=lambda d: d.v7_lvl("biochar", biochar_dose_pct_w_w=5)["action"]["biochar_dose_pct_w_w"])
    chk("Twin [[a]]-day H$_2$ and CO uptake feed downstream (gated by , stoichiometry [[b]]/[[c]] mol CO$_2$ per mol H$_2$",
        a=lambda d: d.v7["carbon_fixation_bound"]["baseline"]["period_days"],
        b=lambda d: 1,
        c=lambda d: 1.0 / d.v7_const("knallgas_co2_per_h2"))
    chk("across $\\{1,1.5,2,10\\}$ times the [[a]] ppm background",
        a=lambda d: d.v7_const("h2_background_ppm"))
    chk("evaporation reduction and [[a]]-day wash pulses across panel cover fractions",
        a=lambda d: d.v7_const("panel_wash_interval_days"))
    chk("increases bulk density from [[a]] g cm$^{-3}$ toward an asymptotic +[[b]] g cm$^{-3}$ increment across $\\{0,1,5,10,50\\}$ vehicle passes",
        a=lambda d: (1.0 - d.defaults["porosity"]) * d.v7_const("particle_density_g_cm3"),
        b=lambda d: d.v7_const("offroad_bd_max_increase_g_cm3"))
    chk("converts bulk density to porosity at quartz particle density [[a]] g cm$^{-3}$",
        a=lambda d: d.v7_const("particle_density_g_cm3"))
    chk("scaling potential evaporation by [[a]] $^\\circ$C$^{-1}$",
        a=lambda d: d.v7_const("evap_scaling_per_c") * 100)
    chk("In each panel, the script detects the horizontal y-axis tick marks from [[a]] to [[b]] t ha$^{-1}$ (step [[c]] t ha$^{-1}$) at pixel column $x=[[e]]$ and fits a linear pixel-to-yield calibration by ordinary least squares, with a maximum calibration residual of [[f]] t ha$^{-1}$ in the low-salinity panel and [[g]] t ha$^{-1}$ in the high-salinity panel",
        a=lambda d: 0,
        b=lambda d: 14,
        c=lambda d: 2,
        e=lambda d: int(re.search(r"x=(\d+)", d.v7_fig1["method"]).group(1)),
        f=lambda d: d.v7_fig1["values"]["low_3.12/calibration_max_residual_t_ha"],
        g=lambda d: d.v7_fig1["values"]["high_7.81/calibration_max_residual_t_ha"])
    chk("against all [[a]] SA187 percentage increases printed on the figure (+[[b]] and +[[c]] fresh and +[[e]] and +[[f]] dry at [[g]] dS m$^{-1}$; +[[h]] and +[[i]] fresh and +[[j]] and +[[k]] dry at [[l]] dS m$^{-1}$) reproduces all [[m]] printed percentages within [[o]] percentage points (maximum relative difference [[p]])",
        a=lambda d: len(d.v7_fig1["check_against_printed_percentages"]),
        b=lambda d: d.v7_fig1["check_against_printed_percentages"]["low_3.12/fresh_s1"]["printed"] * 100,
        c=lambda d: d.v7_fig1["check_against_printed_percentages"]["low_3.12/fresh_s2"]["printed"] * 100,
        e=lambda d: d.v7_fig1["check_against_printed_percentages"]["low_3.12/dry_s1"]["printed"] * 100,
        f=lambda d: d.v7_fig1["check_against_printed_percentages"]["low_3.12/dry_s2"]["printed"] * 100,
        g=lambda d: d.v7_const("sa187_ecw_low_dS_m"),
        h=lambda d: d.v7_fig1["check_against_printed_percentages"]["high_7.81/fresh_s1"]["printed"] * 100,
        i=lambda d: d.v7_fig1["check_against_printed_percentages"]["high_7.81/fresh_s2"]["printed"] * 100,
        j=lambda d: d.v7_fig1["check_against_printed_percentages"]["high_7.81/dry_s1"]["printed"] * 100,
        k=lambda d: d.v7_fig1["check_against_printed_percentages"]["high_7.81/dry_s2"]["printed"] * 100,
        l=lambda d: d.v7_const("sa187_ecw_high_dS_m"),
        m=lambda d: len(d.v7_fig1["check_against_printed_percentages"]),
        o=lambda d: d.v7["inoculated_alfalfa_salinity"]["digitisation_error"]["relative"] * 100,
        p=lambda d: d.v7["inoculated_alfalfa_salinity"]["digitisation_error"]["relative"])
    chk("Digitised season-[[a]] mock yields at [[b]] and [[c]] dS m$^{-1}$ are [[e]] and [[f]] t ha$^{-1}$ for fresh weight (high-to-low ratio [[g]], digitisation envelope [[h]] to [[i]]) and [[j]] and [[k]] t ha$^{-1}$ for dry biomass (ratio [[l]], [[m]] to [[o]])",
        a=lambda d: 1,
        b=lambda d: d.v7_const("sa187_ecw_low_dS_m"),
        c=lambda d: d.v7_const("sa187_ecw_high_dS_m"),
        e=lambda d: d.v7["inoculated_alfalfa_salinity"]["observed"]["fresh"]["yields_t_ha"]["low"],
        f=lambda d: d.v7["inoculated_alfalfa_salinity"]["observed"]["fresh"]["yields_t_ha"]["high"],
        g=lambda d: d.v7["inoculated_alfalfa_salinity"]["observed"]["fresh"]["value"],
        h=lambda d: d.v7["inoculated_alfalfa_salinity"]["observed"]["fresh"]["low"],
        i=lambda d: d.v7["inoculated_alfalfa_salinity"]["observed"]["fresh"]["high"],
        j=lambda d: d.v7["inoculated_alfalfa_salinity"]["observed"]["dry"]["yields_t_ha"]["low"],
        k=lambda d: d.v7["inoculated_alfalfa_salinity"]["observed"]["dry"]["yields_t_ha"]["high"],
        l=lambda d: d.v7["inoculated_alfalfa_salinity"]["observed"]["dry"]["value"],
        m=lambda d: d.v7["inoculated_alfalfa_salinity"]["observed"]["dry"]["low"],
        o=lambda d: d.v7["inoculated_alfalfa_salinity"]["observed"]["dry"]["high"])
    chk("Three root-zone salinity formulations were chained to the Maas--Hoffman alfalfa response (threshold [[a]] dS m$^{-1}$, slope [[b]] dS$^{-1}$ m) and the SA187 yield multiplier ($b=[[c]]$, range [[e]] to [[f]])",
        a=lambda d: d.v7_const("mh_alfalfa_threshold_dS_m"),
        b=lambda d: d.v7_const("mh_alfalfa_slope_pct_per_dS_m"),
        c=lambda d: d.v7_const("sa187_benefit_b"),
        e=lambda d: d.v7_const("sa187_benefit_b", "low"),
        f=lambda d: d.v7_const("sa187_benefit_b", "high"))
    chk("Model A ( ) uses the standard fixed leaching factor $\\mathrm{EC_e}=[[a]] \\mathrm{EC_w}$ (range [[b]] to [[c]] across [[e]] to [[f]] leaching fractions). At [[g]] and [[h]] dS m$^{-1}$, implied $\\mathrm{EC_e}$ is [[i]] and [[j]] dS m$^{-1}$, predicting a mock high-to-low yield ratio of [[k]] ([[l]] to [[m]]; gap $[[o]]$ fresh and $[[p]]$ dry) and an inoculated relative yield at [[q]] dS m$^{-1}$ of [[r]] ([[s]] to [[t]])",
        a=lambda d: d.v7_const("fao29_ece_per_ecw"),
        b=lambda d: d.v7_const("fao29_ece_per_ecw", "low"),
        c=lambda d: d.v7_const("fao29_ece_per_ecw", "high"),
        e=lambda d: 30,
        f=lambda d: 10,
        g=lambda d: d.v7_const("sa187_ecw_low_dS_m"),
        h=lambda d: d.v7_const("sa187_ecw_high_dS_m"),
        i=lambda d: d.v7["inoculated_alfalfa_salinity"]["fitted_leaching"]["comparison_with_measured_soil_ec"]["A_1.5x3.12"],
        j=lambda d: d.v7_const("fao29_ece_per_ecw") * d.v7_const("sa187_ecw_high_dS_m"),
        k=lambda d: d.v7["inoculated_alfalfa_salinity"]["validation_mock_ratio"]["A_fao_fixed_leaching"]["predicted_ratio"]["value"],
        l=lambda d: d.v7["inoculated_alfalfa_salinity"]["validation_mock_ratio"]["A_fao_fixed_leaching"]["predicted_ratio"]["low"],
        m=lambda d: d.v7["inoculated_alfalfa_salinity"]["validation_mock_ratio"]["A_fao_fixed_leaching"]["predicted_ratio"]["high"],
        o=lambda d: d.v7["inoculated_alfalfa_salinity"]["validation_mock_ratio"]["A_fao_fixed_leaching"]["gap_to_observed_fresh"],
        p=lambda d: d.v7["inoculated_alfalfa_salinity"]["validation_mock_ratio"]["A_fao_fixed_leaching"]["gap_to_observed_dry"],
        q=lambda d: d.v7_const("sa187_ecw_high_dS_m"),
        r=lambda d: d.v7["inoculated_alfalfa_salinity"]["headline_other_models"]["A_fao_fixed_leaching"]["value"],
        s=lambda d: d.v7["inoculated_alfalfa_salinity"]["headline_other_models"]["A_fao_fixed_leaching"]["low"],
        t=lambda d: d.v7["inoculated_alfalfa_salinity"]["headline_other_models"]["A_fao_fixed_leaching"]["high"])
    chk("Model B ( ) runs the twin's bare-sand water and salt balance under daily irrigation at [[a]] mm day$^{-1}$ and converts mean soil-water conductivity via $\\mathrm{EC_e}=\\mathrm{EC_{sw}}/2$. With [[b]] mm day$^{-1}$ bare-soil evaporation and no crop transpiration, drainage keeps implied $\\mathrm{EC_e}$ at [[c]] dS m$^{-1}$ for $\\mathrm{EC_w}=[[e]]$ dS m$^{-1}$, predicting a mock yield ratio of [[f]] ([[g]] to [[h]]; gap +[[i]] fresh and +[[j]] dry) and an inoculated relative yield at [[k]] dS m$^{-1}$ of [[l]] ([[m]] to [[o]]), while twin gas uptake abstains at $\\mathrm{EC_w}=[[p]]$ dS m$^{-1}$ outside the [[q]] mol L$^{-1}$ dilute NaCl domain",
        a=lambda d: d.v7_const("alfalfa_etc_mm_per_day"),
        b=lambda d: _src_num(d, "analyse.py", r"\"evaporation_mm\": ([\d.]+)"),
        c=lambda d: d.v7["inoculated_alfalfa_salinity"]["fitted_leaching"]["comparison_with_measured_soil_ec"]["B_twin_ece_at_3.12"],
        e=lambda d: d.v7_const("sa187_ecw_low_dS_m"),
        f=lambda d: d.v7["inoculated_alfalfa_salinity"]["validation_mock_ratio"]["B_twin_physics"]["predicted_ratio"]["value"],
        g=lambda d: d.v7["inoculated_alfalfa_salinity"]["validation_mock_ratio"]["B_twin_physics"]["predicted_ratio"]["low"],
        h=lambda d: d.v7["inoculated_alfalfa_salinity"]["validation_mock_ratio"]["B_twin_physics"]["predicted_ratio"]["high"],
        i=lambda d: d.v7["inoculated_alfalfa_salinity"]["validation_mock_ratio"]["B_twin_physics"]["gap_to_observed_fresh"],
        j=lambda d: d.v7["inoculated_alfalfa_salinity"]["validation_mock_ratio"]["B_twin_physics"]["gap_to_observed_dry"],
        k=lambda d: d.v7_const("sa187_ecw_high_dS_m"),
        l=lambda d: d.v7["inoculated_alfalfa_salinity"]["headline_other_models"]["B_twin_physics"]["value"],
        m=lambda d: d.v7["inoculated_alfalfa_salinity"]["headline_other_models"]["B_twin_physics"]["low"],
        o=lambda d: d.v7["inoculated_alfalfa_salinity"]["headline_other_models"]["B_twin_physics"]["high"],
        p=lambda d: d.v7_lvl("inoculated_alfalfa_salinity", ec_w_dS_m=12)["action"]["ec_w_dS_m"],
        q=lambda d: d.defaults["osmotic_domain_mol_l"])
    chk("Model C ( and ) fits $f$ in $\\mathrm{EC_e}=f \\mathrm{EC_w}$ so that the Maas--Hoffman mock yield ratio between [[a]] and [[b]] dS m$^{-1}$ matches the digitised first-season mock ratio, yielding $f=[[c]]$ ([[e]] to [[f]]) for fresh weight and $f=[[g]]$ ([[h]] to [[i]]) for dry biomass. At $\\mathrm{EC_w}=[[j]]$ dS m$^{-1}$, these fitted factors imply root-zone $\\mathrm{EC_e}=[[k]]$ dS m$^{-1}$ ([[l]] to [[m]] dS m$^{-1}$, fresh) and $[[o]]$ dS m$^{-1}$ ([[p]] to [[q]] dS m$^{-1}$, dry), matching the measured site soil conductivity of [[r]] dS m$^{-1}$, and predict inoculated relative yields at [[s]] dS m$^{-1}$ of [[t]] ([[u]] to [[v]]) for fresh weight and [[w]] ([[x]] to [[y]]) for dry biomass",
        a=lambda d: d.v7_const("sa187_ecw_high_dS_m"),
        b=lambda d: d.v7_const("sa187_ecw_low_dS_m"),
        c=lambda d: d.v7["inoculated_alfalfa_salinity"]["fitted_leaching"]["fresh"]["f_central"],
        e=lambda d: d.v7["inoculated_alfalfa_salinity"]["fitted_leaching"]["fresh"]["f_range_from_digitisation_error"][0],
        f=lambda d: d.v7["inoculated_alfalfa_salinity"]["fitted_leaching"]["fresh"]["f_range_from_digitisation_error"][1],
        g=lambda d: d.v7["inoculated_alfalfa_salinity"]["fitted_leaching"]["dry"]["f_central"],
        h=lambda d: d.v7["inoculated_alfalfa_salinity"]["fitted_leaching"]["dry"]["f_range_from_digitisation_error"][0],
        i=lambda d: d.v7["inoculated_alfalfa_salinity"]["fitted_leaching"]["dry"]["f_range_from_digitisation_error"][1],
        j=lambda d: d.v7_const("sa187_ecw_low_dS_m"),
        k=lambda d: d.v7["inoculated_alfalfa_salinity"]["fitted_leaching"]["fresh"]["ece_at_3_12_dS_m"],
        l=lambda d: d.v7["inoculated_alfalfa_salinity"]["fitted_leaching"]["fresh"]["ece_at_3_12_range"][0],
        m=lambda d: d.v7["inoculated_alfalfa_salinity"]["fitted_leaching"]["fresh"]["ece_at_3_12_range"][1],
        o=lambda d: d.v7["inoculated_alfalfa_salinity"]["fitted_leaching"]["dry"]["ece_at_3_12_dS_m"],
        p=lambda d: d.v7["inoculated_alfalfa_salinity"]["fitted_leaching"]["dry"]["ece_at_3_12_range"][0],
        q=lambda d: d.v7["inoculated_alfalfa_salinity"]["fitted_leaching"]["dry"]["ece_at_3_12_range"][1],
        r=lambda d: d.v7_const("sa187_soil_ec_dS_m"),
        s=lambda d: d.v7_const("sa187_ecw_high_dS_m"),
        t=lambda d: d.v7_head("inoculated_alfalfa_salinity")["value"],
        u=lambda d: d.v7_head("inoculated_alfalfa_salinity")["low"],
        v=lambda d: d.v7_head("inoculated_alfalfa_salinity")["high"],
        w=lambda d: d.v7["inoculated_alfalfa_salinity"]["headline_other_models"]["C_fitted_dry"]["value"],
        x=lambda d: d.v7["inoculated_alfalfa_salinity"]["headline_other_models"]["C_fitted_dry"]["low"],
        y=lambda d: d.v7["inoculated_alfalfa_salinity"]["headline_other_models"]["C_fitted_dry"]["high"])
    chk("All [[a]] external constants used in the composition layer reside in (SHA256 digest recorded in every output JSON). Each entry records a central value, low and high bounds, unit, source citation and a verification status from a manual audit against the cited sources: [[b]] constants have status (the two leaching and conductivity conversion factors, the eight Maas--Hoffman crop thresholds and slopes, the Julian year length, the tropospheric H$_2$ background, the potential-evaporation warming scaling, and the ten Hada Al-Sham SA187 field-trial constants and digitised ratios); [[c]] constants have status (the two crop water demands, the full-cover organic-mulch evaporation reduction, the biocrust carbon-fixation context value, the solar-panel soil cooling, evaporation reduction and washing interval, the cloud-seeding rainfall enhancement, the off-road compaction $e$-folding pass scale, and the soil-to-air warming pass-through); and [[e]] constants have status",
        a=lambda d: len(d.v7_reg),
        b=lambda d: d.v7_status_count("verified_in_source"),
        c=lambda d: d.v7_status_count("modeller_assumption"),
        e=lambda d: d.v7_status_count("source_not_checked"))
    chk("for ([[a]] ppm) was corrected to ; (ii) was re-centred to [[b]] dS$^{-1}$ m with bounds [[c]] to [[e]] dS$^{-1}$ m (the lower bound from the FAO evapotranspiration paper; the upper bound is the value usually attributed to Maas and Hoffman, not checked here); (iii) the source citations for the Maas--Hoffman slopes of alfalfa ([[f]] dS$^{-1}$ m), date palm ([[g]] dS$^{-1}$ m) and barley ([[h]] dS$^{-1}$ m) were updated to the FAO evapotranspiration paper ( ); (iv) the low and high bounds of (central [[i]]) were set to [[j]] and [[k]] from the concentration factors for [[l]] and [[m]] leaching fractions; (v) the upper bound of (central [[o]] $^\\circ$C$^{-1}$) was adjusted to [[p]] $^\\circ$C$^{-1}$ to match the Clausius--Clapeyron ceiling in Scheff and Frierson (2014); and (vi) the citation for ([[q]]) was corrected to of the FAO evapotranspiration paper and its status set to as a linear extrapolation of the [[r]] evaporation reduction per [[s]] organic-mulch cover rule",
        a=lambda d: d.v7_const("h2_background_ppm"),
        b=lambda d: d.v7_const("mh_tomato_slope_pct_per_dS_m"),
        c=lambda d: d.v7_const("mh_tomato_slope_pct_per_dS_m", "low"),
        e=lambda d: d.v7_const("mh_tomato_slope_pct_per_dS_m", "high"),
        f=lambda d: d.v7_const("mh_alfalfa_slope_pct_per_dS_m"),
        g=lambda d: d.v7_const("mh_date_palm_slope_pct_per_dS_m"),
        h=lambda d: d.v7_const("mh_barley_slope_pct_per_dS_m"),
        i=lambda d: d.v7_const("fao29_ece_per_ecw"),
        j=lambda d: d.v7_const("fao29_ece_per_ecw", "low"),
        k=lambda d: d.v7_const("fao29_ece_per_ecw", "high"),
        l=lambda d: 30,
        m=lambda d: 10,
        o=lambda d: d.v7_const("evap_scaling_per_c"),
        p=lambda d: d.v7_const("evap_scaling_per_c", "high"),
        q=lambda d: d.v7_const("mulch_evap_reduction_frac"),
        r=lambda d: 5,
        s=lambda d: 10)


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
    v7_rows = [r for r in rows("v07_composition.tex") if r[0].isdigit()]
    out += 1
    if len(v7_rows) != len(d.v7_sum["experiments"]):
        bad.append(f"v07_composition.tex has {len(v7_rows)} rows, summary.json has {len(d.v7_sum['experiments'])}")
    for r, s in zip(v7_rows, d.v7_sum["experiments"]):
        cmp(f"v07 row {s['id']} n", r[0], s["n"])
        h = s["headline"]
        for tag, lit, val in zip(("value", "low", "high"), r[3].replace("[", ",").replace("]", "").split(","), (h["value"], h["low"], h["high"])):
            m = re.search(r"\\times10\^\{(-?\d+)\}", lit)
            if m:
                val = val / (10 ** int(m.group(1)))
                lit = lit[:m.start()]
            cmp(f"v07 row {s['id']} {tag}", lit.replace("$", "").strip(), val)
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
