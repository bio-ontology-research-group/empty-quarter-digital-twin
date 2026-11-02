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

7. Statements that depend on the outputs without being numeric literals (which interventions have a direction shared by at
   least 95% of draws, "in every draw" claims, the Sonoran comparison against the training mean) are checked by
   `derived_checks`; the generated row tables of the main text (Table 2) and of the supplement (constants, envelopes) are
   cross-checked against their sources by `table_checks`.

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
            "reproducibility/v07/digitise_fig1.py", "reproducibility/v07/composition.py",
            "src/eq_biodt/v03/uq.py", "reproducibility/v08/gather.py", "reproducibility/v08/paper_figures.py",
            "reproducibility/v08/paper_numbers.py")}
        self.v7_sum = J("reviews/v07/composition/summary.json")
        self.v7 = {r["id"]: J(f"reviews/v07/composition/{r['id']}.json") for r in self.v7_sum["experiments"]}
        self.v7_reg = {c["id"]: c for c in J("data/external/composition_constants.json")["constants"]}
        self.v7_fig1 = J("data/external/de_zelicourt_2018/fig1_digitised.json")
        self.v8_params_list = J("data/external/twin_parameter_ranges.json")["parameters"]
        self.v8_params = {p["id"]: p for p in self.v8_params_list}
        self.v8_run = J("reviews/v08/uq/run_record.json")
        self.v8_idx = J("reviews/v08/questions/index.json")
        self.v8_q = {q["slug"]: J(f"reviews/v08/questions/{q['slug']}.json") for q in self.v8_idx["questions"]}
        self.v8_uq = {q["slug"]: J(f"reviews/v08/uq/{q['slug']}.json") for q in self.v8_idx["questions"]}

    def v7_lvl(self, exp, **kw):
        return next(l for l in self.v7[exp]["levels"] if all(l.get(k, l.get("action", {}).get(k)) == v for k, v in kw.items()))

    def v7_head(self, exp):
        return self.v7[exp]["headline"][0]

    def v7_const(self, cid, field="value"):
        return self.v7_reg[cid][field]

    def v7_status_count(self, status):
        return sum(1 for c in self.v7_reg.values() if c["status"] == status)

    def v8_param(self, pid, field="value"):
        return self.v8_params[pid][field]

    def v8_q_lvl(self, slug, series_id, level_idx):
        ser = next(s for s in self.v8_q[slug]["series"] if s["id"] == series_id)
        return ser["levels"][level_idx]

    def v8_q_out(self, slug, series_id, level_idx, out_key):
        return self.v8_q_lvl(slug, series_id, level_idx)["outcomes"][out_key]

    def v8_uq_lvl(self, slug, series_id, level_idx):
        ser = next(s for s in self.v8_uq[slug]["series"] if s["id"] == series_id)
        return next(l for l in ser["levels"] if l["level"] == level_idx)

    def v8_uq_out(self, slug, series_id, level_idx, out_key):
        return self.v8_uq_lvl(slug, series_id, level_idx)["outputs"][out_key]

    def v8_total_levels(self):
        return sum(len(s["levels"]) for u in self.v8_uq.values() for s in u["series"])

    def v8_all_outputs(self):
        return [o for u in self.v8_uq.values() for s in u["series"] for l in s["levels"] for o in l["outputs"].values()]

    def v8_total_outputs(self):
        return len(self.v8_all_outputs())

    def v8_flagged_shifts_pct(self):
        return [o["convergence"]["half_n_shift_frac"] * 100 for o in self.v8_all_outputs() if o.get("convergence", {}).get("flagged")]

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
        return MACRO_VALUES[name].replace("\\%", "%").replace("{,}", ",")
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

    # ---- results: primary observations
    chk("The twin represents [[a]] soil metagenomes", a=lambda d: d.pc["samples"])
    chk("contribute [[a]], [[b]] and [[c]] metagenomes, respectively",
        a=lambda d: d.pc["compartments"]["Surface"], b=lambda d: d.pc["compartments"]["Deep"],
        c=lambda d: d.pc["compartments"]["PlantRoot"])

    # ---- results: deduction
    chk("a carbohydrate-carbon ceiling of [[a]] nmol",
        a=lambda d: d.trace["ontology_supported"]["fluxes"]["carbohydrate_c"] * NM)
    chk("genome-resolved intersection, for [[a]] samples", a=lambda d: d.pc["mapping_overlap_samples"])

    # ---- results: compartments
    chk("is [[a]] in subsurface sand, [[b]] in surface sand and [[c]] in root-adjacent soil",
        a=lambda d: d.comp("Deep").cbb_all_steps_detected_share_mapped.median() * 100,
        b=lambda d: d.comp("Surface").cbb_all_steps_detected_share_mapped.median() * 100,
        c=lambda d: d.comp("PlantRoot").cbb_all_steps_detected_share_mapped.median() * 100)
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
    chk("no-temperature-amplification member, $Q_{10}=[[a]]$", a=lambda d: _unique(d.fit_values("q10")))
    chk("computed on [[a]] plot--campaign means", a=lambda d: d.ext["sonoran"]["field_plot_campaign_means"])

    # ---- results: interventions
    chk("Daily additions totalling [[a]] mm yield [[b]] $\\mu$mol H$_2$ m$^{-2}$ at the reference parameters",
        a=lambda d: d.v6["treatment_versions"]["total_water_mm"], b=lambda d: d.nc("daily_pulses") * UM)
    chk("increasing integration resolution from [[a]] to [[b]] steps per day changes reference-schedule uptake by at most [[c]]",
        a=OneOf(lambda d: set(d.aud["numerical_convergence"][1]["h2_uptake_by_steps"])),
        b=OneOf(lambda d: set(d.aud["numerical_convergence"][1]["h2_uptake_by_steps"])),
        c=lambda d: d.conv_max_rel_pct())
    chk("each add [[a]] mm of water and for the untreated control, across [[b]] parameter members",
        a=lambda d: d.v6["treatment_versions"]["total_water_mm"], b=lambda d: d.v6["treatment_versions"]["members"])

    # ---- methods: evidence inputs
    chk("the [[a]] core-site reference coordinates", a=lambda d: d.pc["core_sites"])
    chk("molecular analysis starts from [[a]] metagenomes", a=lambda d: d.pc["samples"])
    chk("with a [[a]] identity floor", a=lambda d: d.cl("M4", r"(\d+)-percent identity"))
    chk("implements the [[a]] conjunctive requirements of KEGG", a=lambda d: d.cl("M2", r"the (\d+) KEGG M00165"))
    # ---- methods: process model
    chk("The day is divided into $n=[[a]]$ substeps", a=lambda d: d.defaults["integration_steps_per_day"])
    chk("dilute domain, default [[a]] mol L$^{-1}$", a=lambda d: d.defaults["osmotic_domain_mol_l"])
    chk("Three schedules add [[a]] mm as one pulse, four pulses or daily pulses", a=lambda d: d.v6["treatment_versions"]["total_water_mm"])
    chk("Pulses fall on days [[a]], [[b]], [[c]] and [[e]] for the four-pulse schedule and on day [[f]] for the single pulse",
        a=OneOf(lambda d: set(d.four_pulse_days())), b=OneOf(lambda d: set(d.four_pulse_days())),
        c=OneOf(lambda d: set(d.four_pulse_days())), e=OneOf(lambda d: set(d.four_pulse_days())),
        f=lambda d: 0 if _src_num(d, "analyse.py", r"\"water_mm\": 12\. if d == (\d)") == 0 else None)
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
    chk("Adding irrigation water ([[a]] mm day$^{-1}$ seasonal demand",
        a=lambda d: d.v7_const("crop_etc_mm_per_day"))

    # ---- v0.7 composition: supplementary.tex Section S8
    chk("Summary of the [[a]] management what-if compositions on the [[b]]-day single-pulse reference baseline",
        a=lambda d: len(d.v7_sum["experiments"]),
        b=lambda d: d.v7["drip_vs_flood"]["baseline"]["period_days"])
    chk("Small-multiples summary of the [[a]] management what-if compositions on the [[b]]-day single-pulse reference baseline",
        a=lambda d: len(d.v7_sum["experiments"]),
        b=lambda d: d.v7["drip_vs_flood"]["baseline"]["period_days"])
    chk("evaporation reduction and [[a]]-day wash pulses across panel cover fractions",
        a=lambda d: d.v7_const("panel_wash_interval_days"))
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

    # ---- v0.8 uncertainty quantification: main.tex Results, Methods, Discussion

    # ---- v0.8 uncertainty quantification: supplementary.tex Section S6
    chk("Headline joint parameter-uncertainty summary across the [[a]] questions",
        a=lambda d: len(d.v8_idx["questions"]))
    chk("In a [[a]]-task cross-machine verification check ([[b]] scalar outputs), [[c]] raw floating-point values differed between workstations at a maximum relative difference of $[[e]]\\times [[f]]$, and all published summaries are rounded to [[g]] significant digits.",
        a=lambda d: d.v8_run["pilot"]["tasks"],
        b=lambda d: int(re.search(r"of (\d+) raw values", d.v8_run["pilot"]["comparison_with_ws"]).group(1)),
        c=lambda d: int(re.search(r"(\d+) of \d+ raw values", d.v8_run["pilot"]["comparison_with_ws"]).group(1)),
        e=lambda d: float(re.search(r"difference ([\d.e-]+)", d.v8_run["pilot"]["comparison_with_ws"]).group(1)) / 1e-11,
        f=lambda d: 1e-11,
        g=lambda d: int(re.search(r"SIG_DIGITS(?::\s*int)?\s*=\s*(\d+)", d.src["reproducibility/v08/gather.py"]).group(1)))
    chk("Global sensitivity of the [[a]] headline question outputs under joint parameter uncertainty.",
        a=lambda d: len(d.v8_idx["questions"]))


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


# ----------------------------------------------------------------------------------------------- registry, restructured manuscript
def _tw(d):
    return d.v6["treatment_versions"]["total_water_mm"]


def _days(d):
    return d.v7["drip_vs_flood"]["baseline"]["period_days"]


def _bd(d):
    return set(d.nam_p.assumed_bulk_density_g_cm3)


def _fresh_alf(d):
    return d.v7["inoculated_alfalfa_salinity"]


def _lvl_x(d, slug, series, level):
    return d.v8_q_lvl(slug, series, level)["x"]


def build_registry_restructured():
    """Numbers of the restructured manuscript (main text and supplement), tied to generated outputs."""
    nq = lambda d: len(d.v8_idx["questions"])
    # ---- abstract and introduction
    chk(r"from [[a]] primary soil metagenomes that separates", a=lambda d: d.pc["samples"])
    chk(r"uncertainty through [[a]] management questions", a=nq)
    # ---- results 1: primary observations and deduction
    chk(r"They cover [[a]] site--campaign visits at [[b]] of the [[c]] core sites",
        a=lambda d: d.pc["site_campaigns"], b=lambda d: d.pc["core_sites"], c=lambda d: len(_ledger_core_sites()))
    chk(r"a screen for Rubisco alone admits [[a]] MAGs", a=lambda d: d.ont["entailed_class_counts"]["RuBisCOMAG"])
    chk(r"(median width [[a]] percentage points)", a=lambda d: d.pc["median_overlap_interval_percentage_points"])
    chk(r"raises the ceiling [[a]]-fold, to [[b]]",
        a=lambda d: d.trace["unsupported_CO_admission"]["fluxes"]["carbohydrate_c"] / d.trace["ontology_supported"]["fluxes"]["carbohydrate_c"],
        b=lambda d: d.trace["unsupported_CO_admission"]["fluxes"]["carbohydrate_c"] * NM)
    chk(r"in mapped CBB carrier share by [[a]] percentage points across [[b]] matched site--campaign groups (site-bootstrap 95% interval [[c]]--[[d]];",
        a=lambda d: d.aud_contrast("cbb_all_steps_detected_share_mapped", "Deep minus PlantRoot")["mean_difference_percentage_points"],
        b=lambda d: d.aud_contrast("cbb_all_steps_detected_share_mapped", "Deep minus PlantRoot")["matched_site_campaigns"],
        c=lambda d: d.aud_contrast("cbb_all_steps_detected_share_mapped", "Deep minus PlantRoot")["site_bootstrap_95_percentile"][0],
        d=lambda d: d.aud_contrast("cbb_all_steps_detected_share_mapped", "Deep minus PlantRoot")["site_bootstrap_95_percentile"][1])
    chk(r"Site reference coordinates between [[a]] and [[b]]$^\circ$E and [[c]] and [[d]]$^\circ$N",
        a=lambda d: d.ev.lon.min(), b=lambda d: d.ev.lon.max(), c=lambda d: d.ev.lat.min(), d=lambda d: d.ev.lat.max())
    # ---- results 2: interventions
    chk(r"After a [[a]] mm pulse, both versions of $\mathrm{do}(\psi=-[[b]] \mathrm{MPa})$ give the same activation",
        a=_tw, b=lambda d: abs(d.vr["targets"][1]["psi_mpa"]))
    chk(r"such as adding [[a]] mm of water is ambiguous", a=_tw)
    # ---- results 3: external desert measurements
    chk(r"mean wet-to-dry ratios of [[a]] at C14-2 and [[b]] at C14-4",
        a=lambda d: d.namib("C14-2").fold, b=lambda d: d.namib("C14-4").fold)
    chk(r"uptake increased from [[a]] to [[b]] g g$^{-1}$ gravimetric water content in every plot, by [[c]]--[[d]] nmol",
        a=lambda d: 0.10, b=lambda d: 0.20, c=lambda d: d.sonoran_range(0.1, 0.2)[0], d=lambda d: d.sonoran_range(0.1, 0.2)[1])
    chk(r"The frozen reference model had RMSE [[a]] nmol", a=lambda d: d.sonoran_model("frozen_reference")["rmse"])
    chk(r"kinetic coefficient gave [[a]]. Calibrating", a=lambda d: d.sonoran_model("calibrated_kinetics")["rmse"])
    chk(r"reduced the error to [[a]], compared with [[b]] for the training mean and [[c]] for a linear",
        a=lambda d: d.sonoran_model("calibrated_retention")["rmse"], b=lambda d: d.sonoran_model("training_mean")["rmse"],
        c=lambda d: d.sonoran_model("linear_moisture_temperature")["rmse"])
    # ---- results 4: composition and the Hada Al-Sham trial
    chk(r"The [[a]] constants of the external models are recorded with their sources; we verified [[b]] against the cited source, [[c]] are modeller assumptions and [[e]] could not be checked",
        a=lambda d: len(d.v7_reg), b=lambda d: d.v7_status_count("verified_in_source"),
        c=lambda d: d.v7_status_count("modeller_assumption"), e=lambda d: d.v7_status_count("source_not_checked"))
    chk(r"alfalfa was irrigated at [[a]] and [[b]] dS m$^{-1}$ with and without",
        a=lambda d: d.v7_const("sa187_ecw_low_dS_m"), b=lambda d: d.v7_const("sa187_ecw_high_dS_m"))
    chk(r"the high-salinity treatment yielded [[a]] of the low-salinity yield for fresh biomass and [[b]] for dry biomass",
        a=lambda d: _fresh_alf(d)["observed"]["fresh"]["value"], b=lambda d: _fresh_alf(d)["observed"]["dry"]["value"])
    chk(r"$\mathrm{EC_e}=[[a]] \mathrm{EC_w}$ , combined with the Maas--Hoffman alfalfa response", a=lambda d: d.v7_const("fao29_ece_per_ecw"))
    chk(r"predicts a ratio of [[a]] (Fig", a=lambda d: _fresh_alf(d)["validation_mock_ratio"]["A_fao_fixed_leaching"]["predicted_ratio"]["value"])
    chk(r"salt balance predicts [[a]], because drainage",
        a=lambda d: _fresh_alf(d)["validation_mock_ratio"]["B_twin_physics"]["predicted_ratio"]["value"])
    chk(r"gives $f=[[a]]$ for fresh and $f=[[b]]$ for dry biomass",
        a=lambda d: _fresh_alf(d)["fitted_leaching"]["fresh"]["f_central"], b=lambda d: _fresh_alf(d)["fitted_leaching"]["dry"]["f_central"])
    chk(r"root-zone $\mathrm{EC_e}$ of [[a]] to [[b]] dS m$^{-1}$ under the low-salinity treatment",
        a=lambda d: _fresh_alf(d)["fitted_leaching"]["dry"]["ece_at_3_12_dS_m"], b=lambda d: _fresh_alf(d)["fitted_leaching"]["fresh"]["ece_at_3_12_dS_m"])
    chk(r"measured site soil conductivity of [[a]] dS m$^{-1}$, a quantity independent of the fit", a=lambda d: d.v7_const("sa187_soil_ec_dS_m"))
    chk(r"The standard rule implies [[a]] dS m$^{-1}$", a=lambda d: _fresh_alf(d)["fitted_leaching"]["comparison_with_measured_soil_ec"]["A_1.5x3.12"])
    chk(r"relative yield at [[a]] dS m$^{-1}$ has a median of", a=lambda d: d.v7_const("sa187_ecw_high_dS_m"))
    chk(r"across the tested range of [[a]] to [[b]] dS m$^{-1}$",
        a=lambda d: d.v7_const("sa187_ecw_low_dS_m"), b=lambda d: d.v7_const("sa187_ecw_high_dS_m"))
    chk(r"in the Hada Al-Sham trial ([[a]] versus [[b]] dS m$^{-1}$)",
        a=lambda d: d.v7_const("sa187_ecw_high_dS_m"), b=lambda d: d.v7_const("sa187_ecw_low_dS_m"))
    # ---- results 5-7: joint uncertainty
    chk(r"We posed [[a]] questions about desert land management", a=nq)
    chk(r"common [[a]]-day baseline of a single [[b]] mm rain event. We sampled the [[c]] twin process parameters",
        a=_days, b=_tw, c=lambda d: len(d.v8_params))
    chk(r"Partitioning [[a]] mm of rain into four [[b]] mm events", a=_tw, b=lambda d: _tw(d) / 4)
    chk(r"Cloud seeding that adds [[a]] to each rain event", a=lambda d: d.v7_const("cloud_seeding_rain_increase_frac") * 100)
    chk(r"warming the soil by [[a]] $^\circ$C at fixed moisture", a=lambda d: _lvl_x(d, "hotter", "warmer-soil-only", 3))
    chk(r"Climate warming of [[a]] $^\circ$C with unchanged rainfall", a=lambda d: _lvl_x(d, "hotter", "warmer-climate", 3))
    chk(r"Biochar at [[a]] w/w reduced uptake", a=lambda d: _lvl_x(d, "biochar", "biochar-dose", 3))
    chk(r"and [[a]] off-road vehicle passes in", a=lambda d: _lvl_x(d, "compaction", "offroad-passes", 3))
    chk(r"seasonal demand of [[a]] mm day$^{-1}$ reduced uptake", a=lambda d: d.v7_const("crop_etc_mm_per_day"))
    chk(r"Irrigating daily with brackish water at [[a]] dS m$^{-1}$ changed uptake", a=lambda d: _lvl_x(d, "salty-water", "daily-irrigation", 3))
    chk(r"for the same irrigation at [[a]] dS m$^{-1}$; the added water", a=lambda d: _lvl_x(d, "salty-water", "daily-irrigation", 0))
    chk(r"Brackish irrigation at [[a]] dS m$^{-1}$ acts on microbial uptake", a=lambda d: _lvl_x(d, "salty-water", "daily-irrigation", 3))
    chk(r"relative to a single [[a]] mm rain event for each intervention", a=_tw)
    chk(r"Percentage change in [[a]]-day hydrogen uptake", a=_days)
    chk(r"Answers to the [[a]] management and community questions", a=nq)
    chk(r"relative to a single [[a]] mm rain event over [[b]] days", a=_tw, b=_days)
    chk(r"In [[a]] of the [[b]] headline outputs, the inputs that we classified as modeller assumptions",
        a=lambda d: sum(1 for q in d.v8_q.values() if q["uncertainty"]["group_ST"]["modeller_assumption"][0] >= 0.85), b=nq)
    chk(r"In the samples with mapped carriers, [[a]]-day CO uptake has a median", a=_days)
    chk(r"a [[a]] inoculum has no effect", a=lambda d: _lvl_x(d, "adding-microbes", "V53PRr2-proportional", 1))
    chk(r"in a sample with a [[a]] initial carrier share", a=lambda d: float(d.ps[d.ps["sample"] == "V53PRr2"].q_co.iloc[0]) * 100)
    chk(r"knallgas CBB stoichiometry of [[a]]/[[b]] mol CO$_2$ per mol H$_2$",
        a=lambda d: 1, b=lambda d: 1.0 / d.v7_const("knallgas_co2_per_h2"))
    chk(r"highest autotroph share ([[a]]), the annualized bound",
        a=lambda d: d.v7_lvl("carbon_fixation_bound", sample="F14Dr1", scenario="baseline (12 mm pulse)")["h2_autotroph_share"] * 100)
    chk(r"of a [[a]] g C m$^{-2}$ yr$^{-1}$ dryland biocrust reference flux", a=lambda d: d.v7_const("biocrust_net_c_fixation_g_c_m2_yr"))
    chk(r"above the [[a]] ppm background", a=lambda d: d.v7_const("h2_background_ppm"))
    chk(r"Of the [[a]] external constants, [[b]] are modeller assumptions or values whose source could not be checked",
        a=lambda d: len(d.v7_reg), b=lambda d: d.v7_status_count("modeller_assumption") + d.v7_status_count("source_not_checked"))
    # ---- methods
    chk(r"; [[a]] supported genomes use the transaldolase route", a=lambda d: d.aud["pathway_audit"]["transaldolase_route_without_sbpase"])
    chk(r"requires completeness at least [[a]] and contamination below [[b]], a subset of the MIMAG",
        a=lambda d: d.cl("M5", r"at least (\d+) percent"), b=lambda d: d.cl("M5", r"below (\d+) percent"))
    chk(r"Increasing the resolution from [[a]] to [[b]] substeps per day changes reference-schedule uptake by at most [[c]]",
        a=OneOf(lambda d: set(d.aud["numerical_convergence"][1]["h2_uptake_by_steps"])),
        b=OneOf(lambda d: set(d.aud["numerical_convergence"][1]["h2_uptake_by_steps"])), c=lambda d: d.conv_max_rel_pct())
    chk(r"The schedule comparison crosses $b\in\{[[a]],[[b]]\}$, $\psi_d\in\{-[[c]],-[[d]]\}$ MPa, $a\in\{[[e]],[[f]]\}$ and kinetic multipliers $\{[[g]],[[h]],[[i]]\}$",
        a=OneOf(lambda d: set(d.tvt.retention_b)), b=OneOf(lambda d: set(d.tvt.retention_b)),
        c=OneOf(lambda d: set(-d.tvt.dry_threshold_mpa)), d=OneOf(lambda d: set(-d.tvt.dry_threshold_mpa)),
        e=OneOf(lambda d: set(d.tvt.activation_shape)), f=OneOf(lambda d: set(d.tvt.activation_shape)),
        g=OneOf(lambda d: set(d.tvt.k_h2_s / 1e-3)), h=OneOf(lambda d: set(d.tvt.k_h2_s / 1e-3)), i=OneOf(lambda d: set(d.tvt.k_h2_s / 1e-3)))
    chk(r"The schedules use a [[a]]-day horizon, initial water [[b]] mm, temperature $[[c]] ^\circ$C and evaporation [[e]] mm day$^{-1}$",
        a=lambda d: _src_num(d, "analyse.py", r"for day in range\((\d+)\)"),
        b=lambda d: _src_num(d, "analyse.py", r"State\(([\d.]+)\), s\[\"forcing\"\]"),
        c=lambda d: _src_num(d, "analyse.py", r"\"temp_c\": ([\d.]+)"),
        e=lambda d: _src_num(d, "analyse.py", r"\"evaporation_mm\": ([\d.]+)"))
    chk(r"Three schedules add [[a]] mm as one pulse, four pulses or daily pulses", a=_tw)
    chk(r"reproduces the printed percentages to within [[a]] percentage points",
        a=lambda d: d.v7["inoculated_alfalfa_salinity"]["digitisation_error"]["relative"] * 100)
    chk(r"Each of the [[a]] questions is a set of composed chains evaluated against a common baseline: a single [[b]] mm rain event on day [[c]], initial water storage [[e]] mm, soil temperature [[f]] $^\circ$C and evaporation demand [[g]] mm day$^{-1}$ over [[h]] days",
        a=nq, b=_tw, c=lambda d: 0, e=lambda d: d.v7["drip_vs_flood"]["baseline"]["initial_water_mm"],
        f=lambda d: _src_num(d, "analyse.py", r"\"temp_c\": ([\d.]+)"), g=lambda d: _src_num(d, "analyse.py", r"\"evaporation_mm\": ([\d.]+)"), h=_days)
    chk(r"irrigation intervals of [[a]], [[b]] and [[c]] days",
        a=OneOf(lambda d: {l["x"] for l in d.v8_q["drip-or-flood"]["series"][0]["levels"]}),
        b=OneOf(lambda d: {l["x"] for l in d.v8_q["drip-or-flood"]["series"][0]["levels"]}),
        c=OneOf(lambda d: {l["x"] for l in d.v8_q["drip-or-flood"]["series"][0]["levels"]}))
    chk(r"and the [[a]] twin process parameters, all of which are modeller assumptions", a=lambda d: len(d.v8_params))
    chk(r"follow log-uniform distributions on $[[[a]],[[b]]]$ s$^{-1}$ around $[[c]]$ s$^{-1}$; $\psi_e$ and $\psi_w$ follow log-uniform distributions on magnitude over $[-[[e]],-[[f]]]$ and $[-[[g]],-[[h]]]$ MPa; $Q_{10}$ is uniform on $[[[i]],[[j]]]$; and $b\in[[[k]],[[l]]]$, $\psi_d\in[-[[m]],-[[w]]]$ MPa, $a\in[[[o]],[[p]]]$, $D_{\mathrm{H_2}}$, $D_{\mathrm{CO}}$, $d\in[[[q]],[[r]]]$ day$^{-1}$, $f\in[[[s]],[[t]]]$ and $\phi\in[[[u]],[[v]]]$ are triangular",
        a=lambda d: d.v8_param("k_h2_s", "low"), b=lambda d: d.v8_param("k_h2_s", "high"), c=lambda d: d.v8_param("k_h2_s", "value"),
        e=lambda d: -d.v8_param("psi_entry_mpa", "low"), f=lambda d: -d.v8_param("psi_entry_mpa", "high"),
        g=lambda d: -d.v8_param("wet_activation_mpa", "low"), h=lambda d: -d.v8_param("wet_activation_mpa", "high"),
        i=lambda d: d.v8_param("q10", "low"), j=lambda d: d.v8_param("q10", "high"),
        k=lambda d: d.v8_param("retention_b", "low"), l=lambda d: d.v8_param("retention_b", "high"),
        m=lambda d: -d.v8_param("dry_threshold_mpa", "low"), w=lambda d: -d.v8_param("dry_threshold_mpa", "high"),
        o=lambda d: d.v8_param("activation_shape", "low"), p=lambda d: d.v8_param("activation_shape", "high"),
        q=lambda d: d.v8_param("drainage_fraction_day", "low"), r=lambda d: d.v8_param("drainage_fraction_day", "high"),
        s=lambda d: d.v8_param("field_capacity_fraction", "low"), t=lambda d: d.v8_param("field_capacity_fraction", "high"),
        u=lambda d: d.v8_param("porosity", "low"), v=lambda d: d.v8_param("porosity", "high"))
    chk(r"Saltelli design with base sample size $N=[[a]]$ generates $N(d+[[b]])$ runs per level",
        a=lambda d: d.v8_run["full_run"]["n"], b=lambda d: 2 + len(d.v8_q["hotter"]["uncertainty"]["group_ST"]))
    chk(r"with [[a]] percentile-bootstrap draws (95% intervals; seed derived from the question, series and level)",
        a=lambda d: int(re.search(r"boot:\s*int\s*=\s*(\d+)", d.src["src/eq_biodt/v03/uq.py"]).group(1)))
    chk(r"The [[a]] non-discrete levels whose interval endpoints shifted by more than [[b]] of the full-$N$ interval width were rerun at $N=[[c]]$ ([[e]] array tasks), leaving [[f]] of [[g]] outputs flagged (endpoint shifts of [[h]] to [[i]])",
        a=lambda d: len(d.v8_run["convergence_rerun"]["levels"]),
        b=lambda d: float(re.search(r"shift\s*>\s*([\d.]+)", d.src["src/eq_biodt/v03/uq.py"]).group(1)) * 100,
        c=lambda d: d.v8_run["convergence_rerun"]["n"], e=lambda d: d.v8_run["convergence_rerun"]["tasks"],
        f=lambda d: d.v8_run["convergence_rerun"]["still_flagged_after_rerun"], g=lambda d: d.v8_total_outputs(),
        h=lambda d: min(d.v8_flagged_shifts_pct()), i=lambda d: max(d.v8_flagged_shifts_pct()))
    chk(r"When fewer than [[a]] of draws violate a validity or osmotic domain, quantiles and indices are computed on admitted draws and the refusal fraction is reported; when at least [[b]] are refused",
        a=lambda d: d.v8_q["hotter"]["methods"]["admissibility"]["refusal_threshold_fraction"] * 100,
        b=lambda d: d.v8_q["hotter"]["methods"]["admissibility"]["refusal_threshold_fraction"] * 100)
    chk(r"The full design ([[a]] tasks, [[b]] evaluations)",
        a=lambda d: d.v8_run["full_run"]["tasks"], b=lambda d: d.v8_run["full_run"]["evaluations"])
    chk(r"In a [[a]]-task cross-machine check ([[b]] scalar outputs), raw floating-point values agreed to within $[[c]]\times [[e]]$ relative difference, and published summaries are rounded to [[f]] significant digits",
        a=lambda d: d.v8_run["pilot"]["tasks"],
        b=lambda d: int(re.search(r"of (\d+) raw values", d.v8_run["pilot"]["comparison_with_ws"]).group(1)),
        c=lambda d: float(re.search(r"difference ([\d.e-]+)", d.v8_run["pilot"]["comparison_with_ws"]).group(1)) / 1e-11,
        e=lambda d: 1e-11,
        f=lambda d: int(re.search(r"SIG_DIGITS(?::\s*int)?\s*=\s*(\d+)", d.src["reproducibility/v08/gather.py"]).group(1)))
    chk(r"and the [[a]] questions can be explored interactively", a=nq)


def build_allow_restructured():
    thr = lambda d: "frac_neg >= 0.95" in d.src["reproducibility/v08/paper_figures.py"] and "frac_pos >= 0.95" in d.src["reproducibility/v08/paper_figures.py"] \
        or "fn >= 0.95" in d.src["reproducibility/v08/paper_figures.py"]
    allow("at least [[a]] of parameter draws", "direction-reporting threshold: a direction is reported when at least 95% of draws share it (fn >= 0.95 in reproducibility/v08/paper_figures.py)", thr, n=3)
    allow("at least [[a]] of draws", "direction-reporting threshold: a direction is reported when at least 95% of draws share it (fn >= 0.95 in reproducibility/v08/paper_figures.py)", thr, n=3)
    allow("the intervention by Rule [[a]].", "name of do-calculus rule 3, not a quantity")
    allow("We digitised Fig. [[a]] of de Z", "figure number in the cited article, see data/external/de_zelicourt_2018/fig1_digitised.json",
          lambda d: "Fig 1" in d.v7_fig1["source"])
    allow("group total-order index of at least [[a]]", "threshold of the count UQHeadlinesAssumGeEightyFive, defined in reproducibility/v08/paper_numbers.py",
          lambda d: "0.85" in d.src["reproducibility/v08/paper_numbers.py"])
    allow("exact zero with interval $[[[a]],[[b]]]$ and deduced provenance", "definition: closed deductive gate returns exact zero [0, 0] (see src/eq_biodt/v03/uq.py)")




# ----------------------------------------------------------------------------------------------- registry, supplement
def _axis(d, exp, key):
    return {l["action"][key] for l in d.v7[exp]["levels"] if key in l["action"]}


def _pc(d, exp, **kw):
    return d.v7_lvl(exp, **kw)["difference_vs_baseline"]["h2_uptake_mol_m2"]["percent"]


def build_registry_supplement():
    """Numbers of the supplementary information, tied to generated outputs."""
    A = lambda d: d.v7["inoculated_alfalfa_salinity"]
    nq = lambda d: len(d.v8_idx["questions"])
    # ---- genome evidence and compartments
    chk("A screen for Rubisco alone admits [[a]] MAGs, and Rubisco together with phosphoribulokinase admits [[b]].",
        a=lambda d: d.ont["entailed_class_counts"]["RuBisCOMAG"], b=lambda d: d.ont["entailed_class_counts"]["RuBisCOPRKMAG"])
    chk("Subsurface sand exceeds root-adjacent soil by [[a]] percentage points across [[b]] matched groups, with a site-bootstrap 95% interval of [[c]]--[[d]] points",
        a=lambda d: d.aud_contrast("cbb_all_steps_detected_share_mapped", "Deep minus PlantRoot")["mean_difference_percentage_points"],
        b=lambda d: d.aud_contrast("cbb_all_steps_detected_share_mapped", "Deep minus PlantRoot")["matched_site_campaigns"],
        c=lambda d: d.aud_contrast("cbb_all_steps_detected_share_mapped", "Deep minus PlantRoot")["site_bootstrap_95_percentile"][0],
        d=lambda d: d.aud_contrast("cbb_all_steps_detected_share_mapped", "Deep minus PlantRoot")["site_bootstrap_95_percentile"][1])
    # ---- invariance, versions and schedules
    chk(r"After a [[a]] mm pulse, $\mathrm{do}(\psi=-[[b]] \mathrm{MPa})$ is realized by", a=_tw, b=lambda d: abs(d.vr["targets"][0]["psi_mpa"]))
    chk("At $-[[a]]$ MPa, where both versions are valid", a=lambda d: abs(d.vr["targets"][1]["psi_mpa"]))
    chk("Adding [[a]] mm of water therefore names several interventions", a=_tw)
    chk("which is set to [[a]] m$^2$ for each wetting schedule", a=lambda d: _src_num(d, "analyse.py", r"0\. if n == 0 else ([\d.]+)"))
    # ---- composition: overview and registry
    chk(r"The composition runner evaluates [[a]] management what-if chains against the [[b]]-day single-pulse rain baseline ([[c]] mm rain on day [[e]], initial water storage [[f]] mm, soil temperature $[[g]] ^\circ$C, potential evaporation [[h]] mm day$^{-1}$; reference uptake [[i]] $\mu$mol H$_2$ m$^{-2}$ and [[j]] $\mu$mol CO m$^{-2}$)",
        a=lambda d: len(d.v7_sum["experiments"]), b=_days, c=_tw, e=lambda d: 0,
        f=lambda d: d.v7["drip_vs_flood"]["baseline"]["initial_water_mm"],
        g=lambda d: _src_num(d, "analyse.py", r"\"temp_c\": ([\d.]+)"), h=lambda d: _src_num(d, "analyse.py", r"\"evaporation_mm\": ([\d.]+)"),
        i=lambda d: d.v7["drip_vs_flood"]["levels"][0]["twin_baseline"]["h2_uptake_mol_m2"] * UM,
        j=lambda d: d.v7["drip_vs_flood"]["levels"][0]["twin_baseline"]["co_uptake_mol_m2"] * UM)
    chk("summarize all [[a]] experiments with their one-at-a-time envelopes", a=lambda d: len(d.v7_sum["experiments"]))
    chk("Summary of the [[a]] management what-if compositions on the [[b]]-day single-pulse reference baseline",
        a=lambda d: len(d.v7_sum["experiments"]), b=_days)
    chk("Small-multiples summary of the [[a]] management what-if compositions on the [[b]]-day single-pulse reference baseline",
        a=lambda d: len(d.v7_sum["experiments"]), b=_days)
    chk("All [[a]] external constants used in the composition layer are recorded in a registry",
        a=lambda d: len(d.v7_reg))
    chk("[[a]] constants are verified in source, [[b]] are modeller assumptions and [[c]] have a source that was not checked",
        a=lambda d: d.v7_status_count("verified_in_source"), b=lambda d: d.v7_status_count("modeller_assumption"),
        c=lambda d: d.v7_status_count("source_not_checked"))
    chk("for the background hydrogen mixing ratio ([[a]] ppm) was corrected to", a=lambda d: d.v7_const("h2_background_ppm"))
    chk("The tomato salinity slope was re-centred to [[a]] dS$^{-1}$ m with bounds [[b]] to [[c]] dS$^{-1}$ m",
        a=lambda d: d.v7_const("mh_tomato_slope_pct_per_dS_m"), b=lambda d: d.v7_const("mh_tomato_slope_pct_per_dS_m", "low"),
        c=lambda d: d.v7_const("mh_tomato_slope_pct_per_dS_m", "high"))
    chk("slopes of alfalfa ([[a]] dS$^{-1}$ m), date palm ([[b]] dS$^{-1}$ m) and barley ([[c]] dS$^{-1}$ m)",
        a=lambda d: d.v7_const("mh_alfalfa_slope_pct_per_dS_m"), b=lambda d: d.v7_const("mh_date_palm_slope_pct_per_dS_m"),
        c=lambda d: d.v7_const("mh_barley_slope_pct_per_dS_m"))
    chk("The low and high bounds of the leaching concentration factor (central [[a]]) were set to [[b]] and [[c]] from the concentration factors for [[e]] and [[f]] leaching fractions",
        a=lambda d: d.v7_const("fao29_ece_per_ecw"), b=lambda d: d.v7_const("fao29_ece_per_ecw", "low"),
        c=lambda d: d.v7_const("fao29_ece_per_ecw", "high"), e=lambda d: 30, f=lambda d: 10)
    chk("The upper bound of the evaporation scaling per degree (central [[a]] $^\\circ$C$^{-1}$) was adjusted to [[b]] $^\\circ$C$^{-1}$ to match",
        a=lambda d: d.v7_const("evap_scaling_per_c"), b=lambda d: d.v7_const("evap_scaling_per_c", "high"))
    chk("The citation for the mulch evaporation reduction ([[a]]) was corrected to", a=lambda d: d.v7_const("mulch_evap_reduction_frac"))
    # ---- composition: chain specifications
    chk("add a generic crop water demand ([[a]] mm day$^{-1}$)", a=lambda d: d.v7_const("crop_etc_mm_per_day"))
    chk("across irrigation conductivities of [[a]], [[b]], [[c]] and [[e]] dS m$^{-1}$. Downstream",
        a=OneOf(lambda d: _axis(d, "brackish_irrigation", "ec_w_dS_m")), b=OneOf(lambda d: _axis(d, "brackish_irrigation", "ec_w_dS_m")),
        c=OneOf(lambda d: _axis(d, "brackish_irrigation", "ec_w_dS_m")), e=OneOf(lambda d: _axis(d, "brackish_irrigation", "ec_w_dS_m")))
    chk("The high-salinity constant-envelope runs at [[a]] dS m$^{-1}$ exceed the [[b]] mol L$^{-1}$ dilute osmotic domain",
        a=lambda d: d.v7_lvl("brackish_irrigation", ec_w_dS_m=8)["action"]["ec_w_dS_m"], b=lambda d: d.defaults["osmotic_domain_mol_l"])
    chk("a generic crop water demand ([[a]] mm day$^{-1}$) and an irrigation schedule deliver equal seasonal water demand on [[b]]-, [[c]]- or [[e]]-day intervals on top of the [[f]] mm rain pulse, and the twin integrates [[g]]-day hydrogen and CO uptake",
        a=lambda d: d.v7_const("crop_etc_mm_per_day"),
        b=OneOf(lambda d: _axis(d, "drip_vs_flood", "irrigation_interval_days")), c=OneOf(lambda d: _axis(d, "drip_vs_flood", "irrigation_interval_days")),
        e=OneOf(lambda d: _axis(d, "drip_vs_flood", "irrigation_interval_days")), f=_tw, g=_days)
    chk("across mulch fractions of [[a]], [[b]] and [[c]] and shade levels of [[e]], [[f]] and [[g]]; full mulch with shade increases [[h]]-day hydrogen uptake by [[i]] ($[[j]]$ to [[k]])",
        a=OneOf(lambda d: _axis(d, "mulch_shade", "mulch_cover_fraction")), b=OneOf(lambda d: _axis(d, "mulch_shade", "mulch_cover_fraction")),
        c=OneOf(lambda d: _axis(d, "mulch_shade", "mulch_cover_fraction")),
        e=OneOf(lambda d: _axis(d, "mulch_shade", "shade_net_cover_fraction")), f=OneOf(lambda d: _axis(d, "mulch_shade", "shade_net_cover_fraction")),
        g=OneOf(lambda d: _axis(d, "mulch_shade", "shade_net_cover_fraction")), h=lambda d: d.v7["mulch_shade"]["baseline"]["period_days"],
        i=lambda d: _pc(d, "mulch_shade", mulch_cover_fraction=1, shade_net_cover_fraction=1)["value"],
        j=lambda d: _pc(d, "mulch_shade", mulch_cover_fraction=1, shade_net_cover_fraction=1)["low"],
        k=lambda d: _pc(d, "mulch_shade", mulch_cover_fraction=1, shade_net_cover_fraction=1)["high"])
    chk("an upstream model adds [[a]] to soil porosity per percentage point w/w across [[b]], [[c]], [[e]] and [[f]] biochar, which raises $\\phi$ from [[g]] to [[h]] at [[i]] w/w",
        a=lambda d: d.v7_const("biochar_porosity_increase_per_pct"),
        b=OneOf(lambda d: _axis(d, "biochar", "biochar_dose_pct_w_w")), c=OneOf(lambda d: _axis(d, "biochar", "biochar_dose_pct_w_w")),
        e=OneOf(lambda d: _axis(d, "biochar", "biochar_dose_pct_w_w")), f=OneOf(lambda d: _axis(d, "biochar", "biochar_dose_pct_w_w")),
        g=lambda d: d.defaults["porosity"], h=lambda d: d.v7_lvl("biochar", biochar_dose_pct_w_w=5)["upstream"]["biochar_porosity"]["porosity"]["value"],
        i=lambda d: d.v7_lvl("biochar", biochar_dose_pct_w_w=5)["action"]["biochar_dose_pct_w_w"])
    chk("a non-carrier MAG at shares of [[a]], [[b]] and [[c]] into a sample",
        a=OneOf(lambda d: {s * 100 for s in _axis(d, "inoculation", "inoculum_share")}), b=OneOf(lambda d: {s * 100 for s in _axis(d, "inoculation", "inoculum_share")}),
        c=OneOf(lambda d: {s * 100 for s in _axis(d, "inoculation", "inoculum_share")}))
    chk("In the carbon fixation experiment, twin [[a]]-day hydrogen and CO uptake feed", a=lambda d: d.v7["carbon_fixation_bound"]["baseline"]["period_days"])
    chk("uses the stoichiometry of [[a]]/[[b]] mol CO$_2$ per mol H$_2$ of Bongers", a=lambda d: 1, b=lambda d: 1.0 / d.v7_const("knallgas_co2_per_h2"))
    chk("across [[a]], [[b]], [[c]] and [[e]] times the [[f]] ppm background",
        a=OneOf(lambda d: _axis(d, "h2_leakage", "h2_multiplier")), b=OneOf(lambda d: _axis(d, "h2_leakage", "h2_multiplier")),
        c=OneOf(lambda d: _axis(d, "h2_leakage", "h2_multiplier")), e=OneOf(lambda d: _axis(d, "h2_leakage", "h2_multiplier")),
        f=lambda d: d.v7_const("h2_background_ppm"))
    chk("and [[a]]-day wash pulses across panel cover fractions of [[b]], [[c]] and [[e]] before the twin runs",
        a=lambda d: d.v7_const("panel_wash_interval_days"),
        b=OneOf(lambda d: _axis(d, "solar_park", "panel_cover_fraction")), c=OneOf(lambda d: _axis(d, "solar_park", "panel_cover_fraction")),
        e=OneOf(lambda d: _axis(d, "solar_park", "panel_cover_fraction")))
    chk("increases bulk density from [[a]] g cm$^{-3}$ toward an asymptotic increment of [[b]] g cm$^{-3}$ across [[c]], [[e]], [[f]] and [[g]] vehicle passes",
        a=lambda d: (1.0 - d.defaults["porosity"]) * d.v7_const("particle_density_g_cm3"), b=lambda d: d.v7_const("offroad_bd_max_increase_g_cm3"),
        c=OneOf(lambda d: _axis(d, "compaction", "offroad_passes")), e=OneOf(lambda d: _axis(d, "compaction", "offroad_passes")),
        f=OneOf(lambda d: _axis(d, "compaction", "offroad_passes")), g=OneOf(lambda d: _axis(d, "compaction", "offroad_passes")))
    chk("at a quartz particle density of [[a]] g cm$^{-3}$", a=lambda d: d.v7_const("particle_density_g_cm3"))
    chk("crosses air warming offsets of [[a]], [[b]], [[c]] and [[e]] $^\\circ$C with rainfall scaling of $[[f]]$, [[g]] and $+[[h]]$, and scales potential evaporation by [[i]] $^\\circ$C$^{-1}$",
        a=OneOf(lambda d: _axis(d, "warming", "warming_c")), b=OneOf(lambda d: _axis(d, "warming", "warming_c")),
        c=OneOf(lambda d: _axis(d, "warming", "warming_c")), e=OneOf(lambda d: _axis(d, "warming", "warming_c")),
        f=OneOf(lambda d: {r * 100 for r in _axis(d, "warming", "rain_change_frac")}), g=OneOf(lambda d: {r * 100 for r in _axis(d, "warming", "rain_change_frac")}),
        h=OneOf(lambda d: {r * 100 for r in _axis(d, "warming", "rain_change_frac")}), i=lambda d: d.v7_const("evap_scaling_per_c") * 100)
    chk("across irrigation conductivities of [[a]], [[b]], [[c]], [[e]], [[f]], [[g]] and [[h]] dS m$^{-1}$.",
        a=OneOf(lambda d: _axis(d, "inoculated_alfalfa_salinity", "ec_w_dS_m")), b=OneOf(lambda d: _axis(d, "inoculated_alfalfa_salinity", "ec_w_dS_m")),
        c=OneOf(lambda d: _axis(d, "inoculated_alfalfa_salinity", "ec_w_dS_m")), e=OneOf(lambda d: _axis(d, "inoculated_alfalfa_salinity", "ec_w_dS_m")),
        f=OneOf(lambda d: _axis(d, "inoculated_alfalfa_salinity", "ec_w_dS_m")), g=OneOf(lambda d: _axis(d, "inoculated_alfalfa_salinity", "ec_w_dS_m")),
        h=OneOf(lambda d: _axis(d, "inoculated_alfalfa_salinity", "ec_w_dS_m")))
    # ---- composition: one-at-a-time results
    chk("Adding irrigation water ([[a]] mm day$^{-1}$ seasonal demand)", a=lambda d: d.v7_const("crop_etc_mm_per_day"))
    chk("daily drip irrigation changes [[a]]-day hydrogen uptake by $[[b]]$ (one-at-a-time constant envelope $[[c]]$ to $[[e]]$), compared with $[[f]]$ for [[g]]-day pulses and $[[h]]$ for [[i]]-day flood pulses",
        a=_days, b=lambda d: _pc(d, "drip_vs_flood", irrigation_interval_days=1)["value"],
        c=lambda d: _pc(d, "drip_vs_flood", irrigation_interval_days=1)["low"], e=lambda d: _pc(d, "drip_vs_flood", irrigation_interval_days=1)["high"],
        f=lambda d: _pc(d, "drip_vs_flood", irrigation_interval_days=3)["value"], g=lambda d: d.v7_lvl("drip_vs_flood", irrigation_interval_days=3)["action"]["irrigation_interval_days"],
        h=lambda d: _pc(d, "drip_vs_flood", irrigation_interval_days=10)["value"], i=lambda d: d.v7_lvl("drip_vs_flood", irrigation_interval_days=10)["action"]["irrigation_interval_days"])
    chk("Under daily brackish irrigation from [[a]] to [[b]] dS m$^{-1}$, twin hydrogen uptake declines by [[c]] to [[e]] as accumulated salt lowers osmotic potential, and the high-salinity envelope runs at [[f]] dS m$^{-1}$ exceed the [[g]] mol L$^{-1}$ dilute domain. Downstream Maas--Hoffman relative yield at [[h]] dS m$^{-1}$ is [[i]] ([[j]] to [[k]]) for alfalfa, [[l]] for barley, [[m]] ([[o]] to [[p]]) for date palm and [[q]] ([[r]] to [[s]]) for tomato",
        a=lambda d: d.v7_lvl("brackish_irrigation", ec_w_dS_m=0.5)["action"]["ec_w_dS_m"], b=lambda d: d.v7_lvl("brackish_irrigation", ec_w_dS_m=8)["action"]["ec_w_dS_m"],
        c=lambda d: -_pc(d, "brackish_irrigation", ec_w_dS_m=0.5)["value"], e=lambda d: -_pc(d, "brackish_irrigation", ec_w_dS_m=8)["value"],
        f=lambda d: d.v7_lvl("brackish_irrigation", ec_w_dS_m=8)["action"]["ec_w_dS_m"], g=lambda d: d.defaults["osmotic_domain_mol_l"],
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
    chk("full organic mulch changes [[a]]-day hydrogen uptake by [[b]], solar-panel cover with [[c]]-day washing by [[e]], cloud seeding (+[[f]] rain on rain days) by [[g]], +[[h]] $^\\circ$C warming at unchanged rainfall by [[i]], shade-net cooling ($-[[j]] ^\\circ$C) by $[[k]]$, [[l]] w/w biochar by $[[m]]$ and [[o]] off-road vehicle passes by $[[p]]$",
        a=lambda d: d.v7["mulch_shade"]["baseline"]["period_days"], b=lambda d: d.v7_head("mulch_shade")["value"],
        c=lambda d: d.v7_const("panel_wash_interval_days"), e=lambda d: d.v7_head("solar_park")["value"],
        f=lambda d: d.v7_const("cloud_seeding_rain_increase_frac") * 100, g=lambda d: d.v7_head("cloud_seeding")["value"],
        h=lambda d: d.v7_lvl("warming", warming_c=4, rain_change_frac=0.0)["action"]["warming_c"], i=lambda d: d.v7_head("warming")["value"],
        j=lambda d: d.v7_const("shade_net_soil_temp_reduction_c"),
        k=lambda d: _pc(d, "mulch_shade", mulch_cover_fraction=0, shade_net_cover_fraction=1)["value"],
        l=lambda d: d.v7_lvl("biochar", biochar_dose_pct_w_w=5)["action"]["biochar_dose_pct_w_w"], m=lambda d: d.v7_head("biochar")["value"],
        o=lambda d: d.v7_lvl("compaction", offroad_passes=50)["action"]["offroad_passes"], p=lambda d: d.v7_head("compaction")["value"])
    chk("Perturbing atmospheric hydrogen above the [[a]] ppm background scales the [[b]]-day soil hydrogen sink linearly across [[c]], [[e]] and [[f]] times background (+[[g]], +[[h]] and +[[i]] uptake, reaching [[j]] $\\mu$mol m$^{-2}$ at [[k]] ppm), because first-order reaction--diffusion is proportional to surface concentration. Peak oxygen demand reaches only [[l]] of diffusive oxygen supply, which leaves oxygen unlimiting and CO uptake unchanged at [[m]] $\\mu$mol m$^{-2}$",
        a=lambda d: d.v7_const("h2_background_ppm"), b=lambda d: d.v7["h2_leakage"]["baseline"]["period_days"],
        c=lambda d: d.v7_lvl("h2_leakage", h2_multiplier=1.5)["action"]["h2_multiplier"], e=lambda d: d.v7_lvl("h2_leakage", h2_multiplier=2)["action"]["h2_multiplier"],
        f=lambda d: d.v7_lvl("h2_leakage", h2_multiplier=10)["action"]["h2_multiplier"],
        g=lambda d: _pc(d, "h2_leakage", h2_multiplier=1.5)["value"], h=lambda d: _pc(d, "h2_leakage", h2_multiplier=2)["value"],
        i=lambda d: _pc(d, "h2_leakage", h2_multiplier=10)["value"],
        j=lambda d: d.v7_lvl("h2_leakage", h2_multiplier=10)["twin"]["h2_uptake_mol_m2"]["value"] * UM,
        k=lambda d: d.v7_lvl("h2_leakage", h2_multiplier=10)["upstream_ppm"]["value"],
        l=lambda d: d.v7_lvl("h2_leakage", h2_multiplier=10)["oxygen_used_over_supply_max"] * 100,
        m=lambda d: d.v7_lvl("h2_leakage", h2_multiplier=10)["twin"]["co_uptake_mol_m2"]["value"] * UM)
    chk("knallgas CBB stoichiometry ([[a]]/[[b]] mol CO$_2$ per mol H$_2$) under the deductive gate gives an annualised carbon-fixation upper bound in the sample with the highest autotroph share ([[c]] autotroph share) of [[e]] g C m$^{-2}$ yr$^{-1}$ (envelope [[f]] to [[g]] g C m$^{-2}$ yr$^{-1}$; [[h]] $\\mu$mol CO$_2$ m$^{-2}$ over [[i]] days), equivalent to [[j]] of a [[k]] g C m$^{-2}$ yr$^{-1}$ dryland biocrust reference flux",
        a=lambda d: 1, b=lambda d: 1.0 / d.v7_const("knallgas_co2_per_h2"),
        c=lambda d: d.v7_lvl("carbon_fixation_bound", sample="F14Dr1", scenario="baseline (12 mm pulse)")["h2_autotroph_share"] * 100,
        e=lambda d: d.v7_head("carbon_fixation_bound")["value"] / 1e6, f=lambda d: d.v7_head("carbon_fixation_bound")["low"] / 1e6,
        g=lambda d: d.v7_head("carbon_fixation_bound")["high"] / 1e6,
        h=lambda d: d.v7_lvl("carbon_fixation_bound", sample="F14Dr1", scenario="baseline (12 mm pulse)")["co2_fixed_upper_bound"]["per_period_umol_m2"]["value"],
        i=lambda d: d.v7["carbon_fixation_bound"]["baseline"]["period_days"],
        j=lambda d: d.v7_lvl("carbon_fixation_bound", sample="F14Dr1", scenario="baseline (12 mm pulse)")["bound_over_context"] * 100,
        k=lambda d: d.v7_const("biocrust_net_c_fixation_g_c_m2_yr"))
    chk("In a sample where the entailed hydrogen-autotroph share is [[a]], the downstream deductive gate returns [[b]] fixed carbon with deduced provenance while twin hydrogen oxidation continues ([[c]] $\\mu$mol m$^{-2}$ under the gate). Coupled CO oxidation converts CO to CO$_2$ in equimolar stoichiometry with [[e]] net change in atmospheric carbon atoms",
        a=lambda d: d.v7_lvl("carbon_fixation_bound", sample="14PRr2", scenario="baseline (12 mm pulse)")["h2_autotroph_share"] * 100,
        b=lambda d: d.v7_lvl("carbon_fixation_bound", sample="14PRr2", scenario="baseline (12 mm pulse)")["c_fixed_upper_bound"]["per_year_ug_c_m2"]["value"],
        c=lambda d: d.v7_lvl("carbon_fixation_bound", sample="14PRr2", scenario="baseline (12 mm pulse)")["twin"]["h2_uptake_mol_m2"]["value"] * UM,
        e=lambda d: d.v7_lvl("carbon_fixation_bound", sample="F14Dr1", scenario="baseline (12 mm pulse)")["co_oxidation"]["net_atmospheric_c_change_mol_m2"])
    chk("In a sample with [[a]] mapped CO-carrier share and a closed CO gate ([[b]] $\\mu$mol m$^{-2}$ uptake), introducing a CO-carrier MAG at [[c]], [[e]] or [[f]] share opens the gate and increases [[g]]-day CO uptake to [[h]] $\\mu$mol m$^{-2}$ under gate coupling, or to [[i]], [[j]] and [[k]] $\\mu$mol m$^{-2}$ under proportional coupling",
        a=lambda d: float(d.ps[d.ps["sample"] == "51PRr3"].q_co.iloc[0]) * 100, b=lambda d: d.v7["inoculation"]["gate"]["co_uptake_before_umol_m2"],
        c=lambda d: d.v7_lvl("inoculation", sample="51PRr3", coupling="gate", share=0.001)["share"] * 100,
        e=lambda d: d.v7_lvl("inoculation", sample="51PRr3", coupling="gate", share=0.01)["share"] * 100,
        f=lambda d: d.v7_lvl("inoculation", sample="51PRr3", coupling="gate", share=0.05)["share"] * 100,
        g=lambda d: d.v7["inoculation"]["baseline"]["period_days"], h=lambda d: d.v7["inoculation"]["gate"]["co_uptake_after_umol_m2"],
        i=lambda d: d.v7_lvl("inoculation", sample="51PRr3", coupling="proportional", share=0.001, mag_id="11Dr3__single_SemiBin_57")["twin"]["co_uptake_mol_m2"]["value"] * UM,
        j=lambda d: d.v7_lvl("inoculation", sample="51PRr3", coupling="proportional", share=0.01, mag_id="11Dr3__single_SemiBin_57")["twin"]["co_uptake_mol_m2"]["value"] * UM,
        k=lambda d: d.v7_lvl("inoculation", sample="51PRr3", coupling="proportional", share=0.05, mag_id="11Dr3__single_SemiBin_57")["twin"]["co_uptake_mol_m2"]["value"] * UM)
    chk("keeps CO uptake at [[a]] $\\mu$mol m$^{-2}$. Introducing the CO carrier into a sample whose gate is already open ([[b]] initial share) gives [[c]] change under gate coupling and raises proportional CO uptake from [[e]] to [[f]], [[g]] and [[h]] $\\mu$mol m$^{-2}$ across the same shares",
        a=lambda d: d.v7["inoculation"]["gate"]["negative_control_co_change_gate_umol_m2"],
        b=lambda d: float(d.ps[d.ps["sample"] == "V53PRr2"].q_co.iloc[0]) * 100, c=lambda d: d.v7["inoculation"]["gate"]["open_sample_co_change_gate_umol_m2"],
        e=lambda d: d.v7_lvl("inoculation", sample="V53PRr2", coupling="proportional", share=0.001)["difference_vs_baseline"]["co_uptake_mol_m2"]["baseline"] * UM,
        f=lambda d: d.v7_lvl("inoculation", sample="V53PRr2", coupling="proportional", share=0.001)["twin"]["co_uptake_mol_m2"]["value"] * UM,
        g=lambda d: d.v7_lvl("inoculation", sample="V53PRr2", coupling="proportional", share=0.01)["twin"]["co_uptake_mol_m2"]["value"] * UM,
        h=lambda d: d.v7_lvl("inoculation", sample="V53PRr2", coupling="proportional", share=0.05)["twin"]["co_uptake_mol_m2"]["value"] * UM)
    # ---- composition: Hada Al-Sham
    chk("Hada Al-Sham alfalfa field trial irrigated at [[a]] and [[b]] dS m$^{-1}$",
        a=lambda d: d.v7_const("sa187_ecw_low_dS_m"), b=lambda d: d.v7_const("sa187_ecw_high_dS_m"))
    chk("With the fitted fresh-weight leaching factor and the [[a]] ([[b]] to [[c]]) SA187 yield multiplier, the predicted inoculated relative yield at [[e]] dS m$^{-1}$ is [[f]] (one-at-a-time envelope [[g]] to [[h]]). A single trial at two salinities identifies the SA187 benefit only as a constant multiplier across the tested [[i]] to [[j]] dS m$^{-1}$ range",
        a=lambda d: d.v7_const("sa187_benefit_b") * 100, b=lambda d: d.v7_const("sa187_benefit_b", "low") * 100,
        c=lambda d: d.v7_const("sa187_benefit_b", "high") * 100, e=lambda d: d.v7_const("sa187_ecw_high_dS_m"),
        f=lambda d: d.v7_head("inoculated_alfalfa_salinity")["value"], g=lambda d: d.v7_head("inoculated_alfalfa_salinity")["low"],
        h=lambda d: d.v7_head("inoculated_alfalfa_salinity")["high"], i=lambda d: d.v7_const("sa187_ecw_low_dS_m"),
        j=lambda d: d.v7_const("sa187_ecw_high_dS_m"))
    # ---- uncertainty quantification
    chk("across all [[a]] management and community questions (", a=nq)
    chk("evaluate [[a]] management and community questions across [[b]] scenario levels and [[c]] output series", a=nq,
        b=lambda d: d.v8_total_levels(), c=lambda d: d.v8_total_outputs())
    chk("lists all [[a]] uncertain twin process parameters", a=lambda d: len(d.v8_params))
    chk("[[a]] constants in total, each assigned a triangular distribution", a=lambda d: len(d.v7_reg))
    chk("All [[a]] parameters carry the verification status", a=lambda d: len(d.v8_params))
    chk("for the headline level and output of each of the [[a]] questions. For the doubling of atmospheric hydrogen ($[[b]]\\times$ background H$_2$), the relative H$_2$ uptake change is $+[[c]]$ across every parameter draw (5--95% interval $[[e]]$ to $[[f]]$) because the first-order reaction--diffusion equation is linear in surface concentration. The Sobol' indices in Table and Fig. for this question are therefore reported on the paired [[g]]-day H$_2$ uptake contrast (central $[[h]] \\mu$mol m$^{-2}$, median $[[i]]$, 5--95% interval $[[j]]$ to $[[k]] \\mu$mol m$^{-2}$, where the H$_2$ rate constant has $S_T=[[l]]$ and the retention exponent has $S_T=[[m]]$)",
        a=nq, b=lambda d: d.v8_q_lvl("air-hydrogen", "h2-multiplier", 2)["x"],
        c=lambda d: d.v8_q_out("air-hydrogen", "h2-multiplier", 2, "h2_uptake_change_pct")["value"],
        e=lambda d: d.v8_q_out("air-hydrogen", "h2-multiplier", 2, "h2_uptake_change_pct")["q05"],
        f=lambda d: d.v8_q_out("air-hydrogen", "h2-multiplier", 2, "h2_uptake_change_pct")["q95"],
        g=lambda d: d.v7["h2_leakage"]["baseline"]["period_days"],
        h=lambda d: d.v8_q_out("air-hydrogen", "h2-multiplier", 2, "h2_uptake_contrast_umol_m2")["reference"],
        i=lambda d: d.v8_q_out("air-hydrogen", "h2-multiplier", 2, "h2_uptake_contrast_umol_m2")["median"],
        j=lambda d: d.v8_q_out("air-hydrogen", "h2-multiplier", 2, "h2_uptake_contrast_umol_m2")["q05"],
        k=lambda d: d.v8_q_out("air-hydrogen", "h2-multiplier", 2, "h2_uptake_contrast_umol_m2")["q95"],
        l=lambda d: d.v8_q["air-hydrogen"]["uncertainty"]["ST"]["k_h2_s"][0],
        m=lambda d: d.v8_q["air-hydrogen"]["uncertainty"]["ST"]["retention_b"][0])
    chk("The full Saltelli design ($N=[[a]]$ base samples, [[b]] bootstrap draws for 95% Sobol' confidence intervals, [[c]] Slurm array tasks and [[e]] model evaluations) was checked for sample-size convergence by comparing the 5--95% interval endpoints between the first $N/2=512$ base samples and all $N=[[f]]$ base samples. Seven non-discrete levels exceeded the [[g]] interval-width shift threshold and were rerun at $N=[[h]]$ ([[i]] array tasks). After the $N=[[j]]$ rerun, [[k]] of the [[l]] evaluated outputs remain flagged with half-$N$ endpoint shifts of [[m]] to [[o]] of the interval width",
        a=lambda d: d.v8_run["full_run"]["n"], b=lambda d: int(re.search(r"boot:\s*int\s*=\s*(\d+)", d.src["src/eq_biodt/v03/uq.py"]).group(1)),
        c=lambda d: d.v8_run["full_run"]["tasks"], e=lambda d: d.v8_run["full_run"]["evaluations"], f=lambda d: d.v8_run["full_run"]["n"],
        g=lambda d: float(re.search(r"shift\s*>\s*([\d.]+)", d.src["src/eq_biodt/v03/uq.py"]).group(1)) * 100,
        h=lambda d: d.v8_run["convergence_rerun"]["n"], i=lambda d: d.v8_run["convergence_rerun"]["tasks"], j=lambda d: d.v8_run["convergence_rerun"]["n"],
        k=lambda d: d.v8_run["convergence_rerun"]["still_flagged_after_rerun"], l=lambda d: d.v8_total_outputs(),
        m=lambda d: min(d.v8_flagged_shifts_pct()), o=lambda d: max(d.v8_flagged_shifts_pct()))
    chk("at base sample size $N=[[a]]$) shows that absolute [[b]]-day hydrogen uptake on the single-pulse baseline spans an order of magnitude",
        a=lambda d: d.v8_run["full_run"]["n"], b=_days)
    chk("Four [[a]] mm showers", a=lambda d: _tw(d) / 4)
    chk("daily [[a]] mm drizzle", a=lambda d: _tw(d) / 14)
    chk("soil-only heating by [[a]] $^\\circ$C at fixed moisture", a=lambda d: _lvl_x(d, "hotter", "warmer-soil-only", 3))
    chk("negative and positive percentage changes for [[a]] $^\\circ$C climate warming at unchanged rainfall", a=lambda d: _lvl_x(d, "hotter", "warmer-climate", 3))
    chk("), [[a]] w/w biochar (", a=lambda d: _lvl_x(d, "biochar", "biochar-dose", 3))
    chk("), [[a]] vehicle passes (", a=lambda d: _lvl_x(d, "compaction", "offroad-passes", 3))
    chk("), [[a]]-day flood irrigation (", a=lambda d: _lvl_x(d, "drip-or-flood", "irrigation-interval", 2))
    chk("daily brackish irrigation at [[a]] dS m$^{-1}$ (", a=lambda d: _lvl_x(d, "salty-water", "daily-irrigation", 3))
    chk("refused outside the [[a]] mol L$^{-1}$ dilute osmotic domain", a=lambda d: d.defaults["osmotic_domain_mol_l"])
    chk("Open-gate [[a]]-day CO uptake has median", a=lambda d: d.v7["inoculation"]["baseline"]["period_days"])
    chk("while closed-gate samples remain identically [[a]] $\\mu$mol m$^{-2}$ (interval [[b]] to [[c]])",
        a=lambda d: d.v8_q_out("missing-co-genes", "co-gate-samples", 0, "co_uptake_umol_m2")["value"],
        b=lambda d: d.v8_q_out("missing-co-genes", "co-gate-samples", 0, "co_uptake_umol_m2")["q05"],
        c=lambda d: d.v8_q_out("missing-co-genes", "co-gate-samples", 0, "co_uptake_umol_m2")["q95"])
    chk("alfalfa relative yield at [[a]] dS m$^{-1}$ under the fitted fresh-weight leaching model", a=lambda d: d.v7_const("sa187_ecw_high_dS_m"))
    chk("contributed by constants verified in source and [[a]] by modeller assumptions",
        a=lambda d: d.v8_q["alfalfa-bacterium"]["uncertainty"]["group_ST"]["modeller_assumption"][0] * 100)
    # ---- Hada Al-Sham digitisation and salinity models (copied from the earlier composition text, unchanged numbers)
    chk("In each panel, the digitiser detects the horizontal y-axis tick marks from [[a]] to [[b]] t ha$^{-1}$ (step [[c]] t ha$^{-1}$) at pixel column $x=[[e]]$ and fits a linear pixel-to-yield calibration by ordinary least squares, with a maximum calibration residual of [[f]] t ha$^{-1}$ in the low-salinity panel and [[g]] t ha$^{-1}$ in the high-salinity panel",
            a=lambda d: 0,
            b=lambda d: 14,
            c=lambda d: 2,
            e=lambda d: int(re.search(r"x=(\d+)", d.v7_fig1["method"]).group(1)),
            f=lambda d: d.v7_fig1["values"]["low_3.12/calibration_max_residual_t_ha"],
            g=lambda d: d.v7_fig1["values"]["high_7.81/calibration_max_residual_t_ha"])
    chk("Model A uses the standard fixed leaching factor $\\mathrm{EC_e}=[[a]] \\mathrm{EC_w}$ (range [[b]] to [[c]] across [[e]] to [[f]] leaching fractions). At [[g]] and [[h]] dS m$^{-1}$, implied $\\mathrm{EC_e}$ is [[i]] and [[j]] dS m$^{-1}$, predicting a mock high-to-low yield ratio of [[k]] ([[l]] to [[m]]; gap $[[o]]$ fresh and $[[p]]$ dry) and an inoculated relative yield at [[q]] dS m$^{-1}$ of [[r]] ([[s]] to [[t]])",
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
    chk("Model B runs the twin's bare-sand water and salt balance under daily irrigation at [[a]] mm day$^{-1}$ and converts mean soil-water conductivity via $\\mathrm{EC_e}=\\mathrm{EC_{sw}}/2$. With [[b]] mm day$^{-1}$ bare-soil evaporation and no crop transpiration, drainage keeps implied $\\mathrm{EC_e}$ at [[c]] dS m$^{-1}$ for $\\mathrm{EC_w}=[[e]]$ dS m$^{-1}$, predicting a mock yield ratio of [[f]] ([[g]] to [[h]]; gap +[[i]] fresh and +[[j]] dry) and an inoculated relative yield at [[k]] dS m$^{-1}$ of [[l]] ([[m]] to [[o]]), while twin gas uptake abstains at $\\mathrm{EC_w}=[[p]]$ dS m$^{-1}$ outside the [[q]] mol L$^{-1}$ dilute NaCl domain",
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
    chk("Model C fits $f$ in $\\mathrm{EC_e}=f \\mathrm{EC_w}$ so that the Maas--Hoffman mock yield ratio between [[a]] and [[b]] dS m$^{-1}$ matches the digitised first-season mock ratio, yielding $f=[[c]]$ ([[e]] to [[f]]) for fresh weight and $f=[[g]]$ ([[h]] to [[i]]) for dry biomass. At $\\mathrm{EC_w}=[[j]]$ dS m$^{-1}$, these fitted factors imply root-zone $\\mathrm{EC_e}=[[k]]$ dS m$^{-1}$ ([[l]] to [[m]] dS m$^{-1}$, fresh) and $[[o]]$ dS m$^{-1}$ ([[p]] to [[q]] dS m$^{-1}$, dry), which matches the measured site soil conductivity of [[r]] dS m$^{-1}$, and predict inoculated relative yields at [[s]] dS m$^{-1}$ of [[t]] ([[u]] to [[v]]) for fresh weight and [[w]] ([[x]] to [[y]]) for dry biomass",
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


def build_allow_supplement():
    allow("(Table [[a]])", "table number in the cited FAO evapotranspiration paper (source page of the alfalfa slope in the constants registry)",
          lambda d: "Table 23" in d.v7_reg["mh_alfalfa_slope_pct_per_dS_m"]["source"]["authors"], n=1)
    allow("corrected to Chapter [[a]] of the FAO", "chapter number in the cited FAO evapotranspiration paper (source page of the mulch constant)",
          lambda d: "Chapter 9" in d.v7_reg["mulch_evap_reduction_frac"]["source"]["page_or_table"])
    allow("linear extrapolation of the [[a]] evaporation reduction per [[b]] organic-mulch cover rule",
          "the FAO rule of 5% evaporation reduction per 10% organic-mulch cover, quoted in the registry source page",
          lambda d: "by about 5% for each 10%" in d.v7_reg["mulch_evap_reduction_frac"]["source"]["page_or_table"])
    allow("we digitised Fig. [[a]] of the published article", "figure number in the cited article, see data/external/de_zelicourt_2018/fig1_digitised.json",
          lambda d: "Fig 1" in d.v7_fig1["source"])
    allow("return exact [[a]] with interval $[[[b]],[[c]]]$ and deduced provenance", "definition: closed deductive gate returns exact zero [0, 0] (see src/eq_biodt/v03/uq.py)")


# ----------------------------------------------------------------------------------------------- allow-list
def build_allow():
    # Literals that are definitions, arithmetic constants or fixed inputs of the cited experiments or the sampling programme.
    allow("propagates zeros: $q_s(C_g)=[[a]]$ entails $\\gamma_g=[[b]]$, $k_g=[[c]]$ and $J_g=[[d]]$, and $\\psi\\leq\\psi_d$ entails $A=[[e]]$",
          "definition: zero propagation rules of the deductive closure (identification.closure)")
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
    allow("so their sum is free to differ from [[a]] in", "arithmetic constant (percentages)")
    allow("[[a]]--[[b]] interval", "quantile interval level (5th to 95th percentile, q05 and q95 in reviews/v08/uq/*.json)",
          lambda d: all("q05" in o and "q95" in o for o in d.v8_all_outputs()))
    allow("[[a]]--[[b]] parameter-uncertainty interval", "quantile interval level (5th to 95th percentile, q05 and q95 in reviews/v08/uq/*.json)",
          lambda d: all("q05" in o and "q95" in o for o in d.v8_all_outputs()))
    allow("[[a]] intervals; seed derived", "bootstrap confidence interval level (95%, see src/eq_biodt/v03/uq.py)")
    allow("[[a]] Sobol' confidence intervals", "bootstrap confidence interval level (95%, see src/eq_biodt/v03/uq.py)")
    allow("with [[a]] bootstrap intervals", "bootstrap confidence interval level (95%, see src/eq_biodt/v03/uq.py)", n=2)


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
    E = v8_macro_expectations(E)
    return {k: v for k, v in E.items() if v is not None}


# ----------------------------------------------------------------------------------------------- v0.8 macros (reviews/v08/questions)
def sig(value, digits):
    """Round to a number of significant digits (the precision of a generated macro)."""
    return 0.0 if value == 0 else float(f"{float(value):.{digits - 1}e}")


def pct_sig(value):
    """Percentage-change macros carry 3 significant digits from 10 upwards and 2 below."""
    return sig(value, 3 if abs(sig(value, 2)) >= 10 else 2)


# macro prefix -> (question, series, level index, output): the sampled levels that the text reports
V8_PCT_LEVELS = {
    "QRainFour": ("rain-pattern", "rain-schedule", 2, "h2_uptake_change_pct"),
    "QRainDrizzle": ("rain-pattern", "rain-schedule", 3, "h2_uptake_change_pct"),
    "QDripDaily": ("drip-or-flood", "irrigation-interval", 0, "h2_uptake_change_pct"),
    "QDripThreeDay": ("drip-or-flood", "irrigation-interval", 1, "h2_uptake_change_pct"),
    "QDripTenDay": ("drip-or-flood", "irrigation-interval", 2, "h2_uptake_change_pct"),
    "QCloudFive": ("cloud-seeding", "seeding-intensity", 1, "h2_uptake_change_pct"),
    "QCloudTen": ("cloud-seeding", "seeding-intensity", 2, "h2_uptake_change_pct"),
    "QCloudTwenty": ("cloud-seeding", "seeding-intensity", 3, "h2_uptake_change_pct"),
    "QSaltyFreshDaily": ("salty-water", "daily-irrigation", 0, "h2_uptake_change_pct"),
    "QSaltyTwo": ("salty-water", "daily-irrigation", 1, "h2_uptake_change_pct"),
    "QSaltyFour": ("salty-water", "daily-irrigation", 2, "h2_uptake_change_pct"),
    "QSaltyEight": ("salty-water", "daily-irrigation", 3, "h2_uptake_change_pct"),
    "QHotterFour": ("hotter", "warmer-climate", 3, "h2_uptake_change_pct"),
    "QHotterSoilFour": ("hotter", "warmer-soil-only", 3, "h2_uptake_change_pct"),
    "QHotterSoilTen": ("hotter", "warmer-soil-only", 4, "h2_uptake_change_pct"),
    "QMulchHalf": ("mulch-shade", "mulch-only", 1, "h2_uptake_change_pct"),
    "QMulchFull": ("mulch-shade", "mulch-only", 2, "h2_uptake_change_pct"),
    "QShadeHalf": ("mulch-shade", "shade-only", 1, "h2_uptake_change_pct"),
    "QShadeFull": ("mulch-shade", "shade-only", 2, "h2_uptake_change_pct"),
    "QMulchShadeFull": ("mulch-shade", "mulch-and-shade", 2, "h2_uptake_change_pct"),
    "QSolarHalf": ("solar-park", "panel-cover", 1, "h2_uptake_change_pct"),
    "QSolarFull": ("solar-park", "panel-cover", 2, "h2_uptake_change_pct"),
    "QCompactionOne": ("compaction", "offroad-passes", 1, "h2_uptake_change_pct"),
    "QCompactionTen": ("compaction", "offroad-passes", 2, "h2_uptake_change_pct"),
    "QCompactionFifty": ("compaction", "offroad-passes", 3, "h2_uptake_change_pct"),
    "QCompactionDirect": ("compaction", "direct-porosity-v06", 0, "h2_uptake_change_pct"),
    "QBiocharOne": ("biochar", "biochar-dose", 1, "h2_uptake_change_pct"),
    "QBiocharTwo": ("biochar", "biochar-dose", 2, "h2_uptake_change_pct"),
    "QBiocharFive": ("biochar", "biochar-dose", 3, "h2_uptake_change_pct"),
    "QAirHtwoOneHalf": ("air-hydrogen", "h2-multiplier", 1, "h2_uptake_change_pct"),
    "QAirHtwoDouble": ("air-hydrogen", "h2-multiplier", 2, "h2_uptake_change_pct"),
    "QAirHtwoTen": ("air-hydrogen", "h2-multiplier", 3, "h2_uptake_change_pct"),
    "QAddingOnePct": ("adding-microbes", "V53PRr2-proportional", 1, "co_uptake_change_pct"),
    "QAddingFivePct": ("adding-microbes", "V53PRr2-proportional", 2, "co_uptake_change_pct"),
}
# macro prefix -> question whose headline inputs rank the S_T macros
V8_HEADLINES = {"QRainFour": "rain-pattern", "QDripTenDay": "drip-or-flood", "QCloudTen": "cloud-seeding",
                "QSaltyEight": "salty-water", "QAlfalfaSa": "alfalfa-bacterium", "QHotterFour": "hotter",
                "QMulchFull": "mulch-shade", "QSolarFull": "solar-park", "QCompactionFifty": "compaction",
                "QBiocharFive": "biochar", "QMissingCOOpen": "missing-co-genes", "QAddingOnePct": "adding-microbes",
                "QAirHtwoDouble": "air-hydrogen", "QCarbonBound": "carbon-fixation"}


def v8_macro_expectations(E):
    """Independent recomputation of the macros of paper/tables/v08_numbers.tex from reviews/v08/questions/*.json."""
    for prefix, (slug, series, lvl, key) in V8_PCT_LEVELS.items():
        o = lambda d, s=slug, r=series, l=lvl, k=key: d.v8_q_out(s, r, l, k)
        E[prefix + "Median"] = lambda d, o=o: pct_sig(o(d)["median"])
        E[prefix + "Low"] = lambda d, o=o: pct_sig(o(d)["q05"])
        E[prefix + "High"] = lambda d, o=o: pct_sig(o(d)["q95"])
        E[prefix + "Ref"] = lambda d, o=o: pct_sig(o(d)["reference"])
        E[prefix + "Pneg"] = lambda d, o=o: o(d)["frac_neg"] * 100
        E[prefix + "Ppos"] = lambda d, o=o: o(d)["frac_pos"] * 100
    E["UQQuestions"] = lambda d: len(d.v8_idx["questions"])
    E["UQLevels"] = lambda d: d.v8_total_levels()
    E["UQOutputs"] = lambda d: d.v8_total_outputs()
    E["UQEvaluations"] = lambda d: d.v8_run["full_run"]["evaluations"]
    E["UQHeadlinesAssumGeEightyFive"] = lambda d: sum(1 for q in d.v8_q.values() if q["uncertainty"]["group_ST"]["modeller_assumption"][0] >= 0.85)
    base = lambda d: d.v8_q_out("rain-pattern", "rain-schedule", 1, "h2_uptake_umol_m2")
    E["QBaseHtwoMedian"] = lambda d: sig(base(d)["median"], 3)
    E["QBaseHtwoLow"] = lambda d: sig(base(d)["q05"], 3)
    E["QBaseHtwoHigh"] = lambda d: sig(base(d)["q95"], 3)
    E["QBaseHtwoRef"] = lambda d: sig(base(d)["reference"], 4)
    E["QBaseHtwoSTkHtwo"] = lambda d: d.v8_q["rain-pattern"]["uncertainty"]["baseline_h2_ST"]["k_h2_s"][0]
    E["QBaseHtwoSTb"] = lambda d: d.v8_q["rain-pattern"]["uncertainty"]["baseline_h2_ST"]["retention_b"][0]
    E["QSaltyEightRefusalPct"] = lambda d: sig(d.v8_q_out("salty-water", "daily-irrigation", 3, "h2_uptake_change_pct")["refusal_fraction"] * 100, 2)
    cf = lambda d: d.v8_q_out("carbon-fixation", "baseline-12mm", 1, "c_fixed_upper_bound_g_m2_yr")
    E["QCarbonBoundMedian"] = lambda d: sig(cf(d)["median"], 2)
    E["QCarbonBoundLow"] = lambda d: sig(cf(d)["q05"], 2)
    E["QCarbonBoundHigh"] = lambda d: sig(cf(d)["q95"], 2)
    E["QCarbonBoundRef"] = lambda d: sig(cf(d)["reference"], 2)
    E["QCarbonBoundPctTen"] = lambda d: sig(100.0 * cf(d)["median"] / 10.0, 2)
    co = lambda d: d.v8_q_out("missing-co-genes", "co-gate-samples", 6, "co_uptake_umol_m2")
    E["QMissingCOOpenMedian"] = lambda d: sig(co(d)["median"], 3)
    E["QMissingCOOpenLow"] = lambda d: sig(co(d)["q05"], 3)
    E["QMissingCOOpenHigh"] = lambda d: sig(co(d)["q95"], 3)
    E["QMissingCOClosedSamples"] = lambda d: sum(
        1 for l in d.v8_q["missing-co-genes"]["series"][0]["levels"]
        if float(d.ps[d.ps["sample"] == l["x"]].q_co.iloc[0]) == 0.0)
    alf = lambda d: d.v8_q_out("alfalfa-bacterium", "C_fitted_fresh", 4, "relative_yield_alfalfa_sa187")
    E["QAlfalfaSaMedian"] = lambda d: sig(alf(d)["median"], 3)
    E["QAlfalfaSaLow"] = lambda d: sig(alf(d)["q05"], 3)
    E["QAlfalfaSaHigh"] = lambda d: sig(alf(d)["q95"], 3)
    E["QAlfalfaVerifiedGroupST"] = lambda d: d.v8_q["alfalfa-bacterium"]["uncertainty"]["group_ST"]["verified_in_source"][0]
    E["QAlfalfaVerifiedGroupSTPct"] = lambda d: pct_sig(100.0 * d.v8_q["alfalfa-bacterium"]["uncertainty"]["group_ST"]["verified_in_source"][0])
    for prefix, slug in V8_HEADLINES.items():
        def ranked(d, s=slug):
            return sorted(d.v8_q[s]["uncertainty"]["ST"].items(), key=lambda kv: (-kv[1][0], kv[0]))
        E[prefix + "STOne"] = lambda d, r=ranked: r(d)[0][1][0]
        E[prefix + "STTwo"] = lambda d, r=ranked: r(d)[1][1][0]
        E[prefix + "GroupSTAssum"] = lambda d, s=slug: d.v8_q[s]["uncertainty"]["group_ST"]["modeller_assumption"][0]
    return E


def derived_checks(d):
    """Statements of the text that are not numeric literals but depend on the outputs (directions, counts, orderings)."""
    bad, n = [], 0

    def out(slug, series, lvl, key="h2_uptake_change_pct"):
        return d.v8_q_out(slug, series, lvl, key)

    def claim(ok, text):
        nonlocal n
        n += 1
        if not ok:
            bad.append(f"DERIVED claim not supported by the outputs: {text}")

    determined = {  # intervention -> sign of the effect that the text states for it
        "four 3 mm showers": (out("rain-pattern", "rain-schedule", 2), -1),
        "daily drizzle": (out("rain-pattern", "rain-schedule", 3), -1),
        "full organic mulch": (out("mulch-shade", "mulch-only", 2), +1),
        "full solar-panel cover": (out("solar-park", "panel-cover", 2), +1),
        "cloud seeding of 10%": (out("cloud-seeding", "seeding-intensity", 2), +1),
        "full shade net": (out("mulch-shade", "shade-only", 2), -1),
        "soil warming of 4 C": (out("hotter", "warmer-soil-only", 3), +1)}
    undetermined = {
        "climate warming of 4 C": out("hotter", "warmer-climate", 3), "biochar 5%": out("biochar", "biochar-dose", 3),
        "50 vehicle passes": out("compaction", "offroad-passes", 3), "daily irrigation": out("drip-or-flood", "irrigation-interval", 0),
        "irrigation every ten days": out("drip-or-flood", "irrigation-interval", 2),
        "brackish daily irrigation at 8 dS/m": out("salty-water", "daily-irrigation", 3)}
    claim(len(determined) == 7, "seven interventions are listed as determined")
    for name, (o, sign) in determined.items():
        share = o["frac_pos"] if sign > 0 else o["frac_neg"]
        claim(share >= 0.95, f"{name} changes uptake in the stated direction in at least 95% of draws (share {share:.4f})")
    for name, o in undetermined.items():
        claim(max(o["frac_pos"], o["frac_neg"]) < 0.95, f"{name} has no direction shared by 95% of draws (shares {o['frac_neg']:.3f}, {o['frac_pos']:.3f})")
    for name in ("cloud seeding of 10%", "full shade net", "soil warming of 4 C"):
        o, sign = determined[name]
        claim((o["frac_pos"] if sign > 0 else o["frac_neg"]) == 1.0, f"{name}: the text states 'in every draw'")
    claim(out("air-hydrogen", "h2-multiplier", 2)["q05"] == out("air-hydrogen", "h2-multiplier", 2)["q95"] == 100.0,
          "doubling atmospheric hydrogen doubles uptake in every draw")
    claim(bool((d.ps.q_h2 > 0).all()), "every sampled community contains mapped hydrogenase carriers")
    claim(int((d.ps.q_co == 0).sum()) > 0 and all(o == 0 for o in d.ps[d.ps.q_co == 0].co_gate_single_pulse),
          "communities without CO carriers have zero CO response")
    # Sonoran: calibrated model against the training mean in each held campaign
    p = d.son_p
    wins, losses = [], []
    for c, g in p.groupby("campaign"):
        r = {m: float(np.sqrt(((x.observed - x.predicted) ** 2).mean())) for m, x in g.groupby("model")}
        (wins if r["calibrated_retention"] < r["training_mean"] else losses).append(c)
    claim(len(wins) == 3 and len(losses) == 1 and losses[0].startswith("2018-08"),
          "the calibrated model improves on the training mean in three of four held campaigns, with a larger error in August")
    # Namib: ten sites, two with quantified dry endpoint
    claim(d.nam_o[d.nam_o.state == "dry"].site.nunique() == 2 and len(d.nam_o[d.nam_o.state == "wet"].site.unique()) == 10,
          "dry rates are quantified for two of ten sites, wet rates for ten")
    # exception of the inoculated alfalfa yield among the headline outputs
    low = [s for s, q in d.v8_q.items() if q["uncertainty"]["group_ST"]["modeller_assumption"][0] < 0.85]
    claim(low == ["alfalfa-bacterium"], "the alfalfa yield is the only headline output with a modeller-assumption group index below 0.85")
    claim(d.v8_q["alfalfa-bacterium"]["uncertainty"]["group_ST"]["verified_in_source"][0] > 0.85, "alfalfa variance comes from constants verified in the trial publication")
    claim(sum(1 for l in d.v8_q["missing-co-genes"]["series"][0]["levels"] if l["outcomes"]["co_uptake_umol_m2"].get("mode") == "deduced") == 5,
          "five samples of the CO question have a closed gate")
    return n, bad


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
    env_rows = [r for r in rows("v08_envelopes.tex") if r[0].isdigit()]
    out += 1
    if len(env_rows) != len(d.v7_sum["experiments"]):
        bad.append(f"v08_envelopes.tex has {len(env_rows)} rows, summary.json has {len(d.v7_sum['experiments'])}")
    for r, s in zip(env_rows, d.v7_sum["experiments"]):
        cmp(f"envelope row {s['id']} n", r[0], s["n"])
        h = s["headline"]
        scale = 1e6 if s["id"] == "carbon_fixation_bound" else 1.0   # table reports g C m-2 yr-1, summary.json ug C m-2 yr-1
        for tag, lit, val in zip(("value", "low", "high"), r[3].replace("[", ",").replace("]", "").split(","), (h["value"], h["low"], h["high"])):
            cmp(f"envelope row {s['id']} {tag}", lit.replace("$", "").strip(), val / scale)
    const_rows = [r for r in rows("v08_constants_rows.tex") if r[0].isdigit()]
    out += 1
    if len(const_rows) != len(d.v7_reg):
        bad.append(f"v08_constants_rows.tex has {len(const_rows)} rows, composition_constants.json has {len(d.v7_reg)}")
    status_text = {"verified_in_source": "verified in source", "modeller_assumption": "modeller assumption", "source_not_checked": "source not checked"}
    for idx, (r, c) in enumerate(zip(const_rows, d.v7_reg.values()), start=1):
        cmp(f"constants row {c['id']} n", r[0], idx)
        for tag, lit, val in zip(("value", "low", "high"), r[2].replace("[", ",").replace("]", "").split(","), (c["value"], c["low"], c["high"])):
            cmp(f"constants row {c['id']} {tag}", lit.replace("$", "").strip(), val)
        out += 1
        if r[5] != status_text[c["status"]]:
            bad.append(f"constants row {c['id']}: status {r[5]!r} differs from the registry ({c['status']})")
    # Table 2 of the main text: the 14 answers
    answer_slugs = ["rain-pattern", "drip-or-flood", "cloud-seeding", "salty-water", "mulch-shade", "solar-park", "hotter",
                    "compaction", "biochar", "air-hydrogen", "missing-co-genes", "adding-microbes", "carbon-fixation", "alfalfa-bacterium"]
    ans = [r for r in rows("v08_answers.tex") if r[0] != "Question"]
    out += 1
    if len(ans) != len(answer_slugs):
        bad.append(f"v08_answers.tex has {len(ans)} rows, expected {len(answer_slugs)}")

    def clean(t):
        return t.replace("$-$", "-").replace("$", "").replace("\\%", "%").strip()

    def triple(cell):
        m = re.match(r"\s*([^\[]+)\[([^,\]]+),\s*([^\]]+)\]", clean(cell))
        return [x.strip().lstrip("+") for x in m.groups()]

    for r, slug in zip(ans, answer_slugs):
        u = d.v8_q[slug]["uncertainty"]
        if slug == "missing-co-genes":
            closed = d.v8_q_out(slug, "co-gate-samples", 0, "co_uptake_umol_m2")
            opn = d.v8_q_out(slug, "co-gate-samples", 6, "co_uptake_umol_m2")
            parts = r[3].split(";")
            for tag, lit, val in zip(("closed median", "closed q05", "closed q95"), triple(parts[0]), (closed["value"], closed["q05"], closed["q95"])):
                cmp(f"answers {slug} {tag}", lit, val)
            for tag, lit, val in zip(("open median", "open q05", "open q95"), triple(parts[1]), (opn["median"], opn["q05"], opn["q95"])):
                cmp(f"answers {slug} {tag}", lit, val)
            n_closed = sum(1 for l in d.v8_q[slug]["series"][0]["levels"] if float(d.ps[d.ps["sample"] == l["x"]].q_co.iloc[0]) == 0.0)
            cmp(f"answers {slug} closed samples", re.search(r"\((\d+) closed\)", r[3]).group(1), n_closed)
            cmp(f"answers {slug} closed samples (direction)", re.search(r"\((\d+) samples\)", r[4]).group(1), n_closed)
            out += 1
            if "all draws" not in r[4] or opn["frac_pos"] != 1.0 or opn["frac_neg"] != 0.0:
                bad.append(f"answers {slug}: direction text {r[4]!r} not supported by the open-gate output")
        else:
            for tag, lit, val in zip(("median", "q05", "q95"), triple(r[3]), (u["median"], u["q05"], u["q95"])):
                cmp(f"answers {slug} {tag}", lit, val)
            o = d.v8_q_out(slug, u["headline_series"], u["headline_level"], u["headline_output"])
            dtext = clean(r[4])
            m = re.search(r"(\d+(?:\.\d+)?)% (?:of draws|(increase|decrease))|(decrease|increase), (\d+(?:\.\d+)?)% of draws", dtext)
            out += 1
            if "all draws" in dtext and "constant" not in dtext:
                sign = "frac_pos" if ("increase" in dtext or "positive" in dtext) else "frac_neg"
                if o[sign] != 1.0:
                    bad.append(f"answers {slug}: {dtext!r} but {sign}={o[sign]}")
            elif "constant" in dtext:
                if not (o["q05"] == o["q95"] == o["median"] == 100.0):
                    bad.append(f"answers {slug}: {dtext!r} but the interval is not exactly +100%")
            else:
                mm = re.search(r"(\d+(?:\.\d+)?)% increase|(\d+(?:\.\d+)?)% decrease|(increase|decrease), (\d+(?:\.\d+)?)% of draws", dtext)
                pct = next(g for g in (mm.group(1), mm.group(2), mm.group(4)) if g)
                inc = ("increase" in dtext)
                cmp(f"answers {slug} direction share", pct, o["frac_pos" if inc else "frac_neg"] * 100)
        ranked = sorted(u["ST"].items(), key=lambda kv: (-kv[1][0], kv[0]))[:2]
        for (pname, st), lit in zip(ranked, re.findall(r"\((\d\.\d+)\)", r[5])):
            cmp(f"answers {slug} ST {pname}", lit, st[0])
    v8p_rows = [r for r in rows("v08_parameters.tex") if r[0].isdigit()]
    out += 1
    if len(v8p_rows) != len(d.v8_params_list):
        bad.append(f"v08_parameters.tex has {len(v8p_rows)} rows, twin_parameter_ranges.json has {len(d.v8_params_list)}")
    for idx, (r, p) in enumerate(zip(v8p_rows, d.v8_params_list), start=1):
        cmp(f"v08_parameters row {p['id']} n", r[0], idx)
        for tag, lit, val in zip(("value", "low", "high"), r[3].replace("[", ",").replace("]", "").split(","), (p["value"], p["low"], p["high"])):
            m = re.search(r"\\times10\^\{(-?\d+)\}", lit)
            if m:
                val = val / (10 ** int(m.group(1)))
                lit = lit[:m.start()]
            cmp(f"v08_parameters row {p['id']} {tag}", lit.replace("$", "").strip(), val)
    v8u_rows = [r for r in rows("v08_uq.tex") if r[0].isdigit()]
    out += 1
    if len(v8u_rows) != len(d.v8_idx["questions"]):
        bad.append(f"v08_uq.tex has {len(v8u_rows)} rows, index.json has {len(d.v8_idx['questions'])}")
    for idx, (r, qinfo) in enumerate(zip(v8u_rows, d.v8_idx["questions"]), start=1):
        slug = qinfo["slug"]
        doc = d.v8_q[slug]
        u = doc["uncertainty"]
        out_obj = d.v8_q_out(slug, u["headline_series"], u["headline_level"], u["headline_output"])
        cmp(f"v08_uq row {slug} n", r[0], idx)
        cmp(f"v08_uq row {slug} central", r[3].replace("$", "").strip(), out_obj["reference"])
        for tag, lit, val in zip(("median", "q05", "q95"), r[4].replace("[", ",").replace("]", "").split(","), (u["median"], u["q05"], u["q95"])):
            cmp(f"v08_uq row {slug} {tag}", lit.replace("$", "").strip(), val)
        ranked = sorted(u["ST"].items(), key=lambda kv: (-kv[1][0], kv[0]))[:2]
        top_lits = re.findall(r"\(([-+]?\d[^)]*)\)", r[5])
        for (pname, st), lit in zip(ranked, top_lits):
            cmp(f"v08_uq row {slug} ST {pname}", lit.replace("$", "").strip(), st[0])
        cmp(f"v08_uq row {slug} group_ST assum", r[6].replace("$", "").strip(), u["group_ST"]["modeller_assumption"][0])
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
    build_registry_restructured()
    build_allow_restructured()
    build_registry_supplement()
    build_allow_supplement()
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
        lit = defs[m].replace("\\%", "%").replace("{,}", ",")
        if literal_matches(lit, val):
            macro_ok += 1
        else:
            errors.append(f"MISMATCH macro \\{m}: tables/*.tex defines {lit}, source gives {val!r}")

    tcount, tbad = table_checks(d)
    errors += [f"MISMATCH generated table {b}" for b in tbad]
    dcount, dbad = derived_checks(d)
    errors += dbad

    if args.verbose:
        print("\n".join("ok " + v for v in verified))
        print("\n".join("allowed " + a for a in allowed))
    print(f"prose/caption literals verified against outputs: {len(verified)}")
    print(f"literals allow-listed (design inputs, definitions, constants): {len(allowed)}")
    print(f"value macros used in the text recomputed from outputs: {macro_ok} of {len(custom & set(defs))}")
    print(f"generated table cells cross-checked: {tcount}")
    print(f"derived claims (directions, counts, orderings) checked: {dcount}")
    if errors:
        print(f"\nFAILED: {len(errors)} problem(s)")
        for e in errors:
            print(" -", e)
        return 1
    print("manuscript numbers verified")
    return 0


if __name__ == "__main__":
    sys.exit(main())
