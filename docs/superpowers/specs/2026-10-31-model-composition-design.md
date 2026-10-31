# EQ-BioDT v0.7 — composing the twin with external models ("chained what-ifs")

Status: approved by R. Hoehndorf 2026-10-31. Branch `v0.7-composition` (worktree `~/tmp/eq-twin-v07`).

## Intent

The twin models one thing: soil water/salt/heat → microbial activity → H2 and CO uptake, with community
factors γ deduced from the ontology. Users ask management questions it cannot answer alone (crops,
carbon, solar parks, cloud seeding). We answer them by **composing** the twin with other published models
or causal theories:

* **upstream** models translate a management action into twin inputs (forcing, prescribed soil parameters,
  state W/N, atmospheric gas concentration, community);
* **downstream** models translate twin outputs (H2/CO uptake, activity windows, final W/N/salinity) into an
  endpoint people care about (crop relative yield, CO2 fixed, …).

The composition must keep the v0.6 guarantees: every number has provenance, interventions stay admissible,
deductive nulls carry through chains, and nothing is reported outside the domain where its model is valid.

## Goals

1. A small composition layer: `src/eq_biodt/v03/composition.py` (contract + chain runner) and
   `src/eq_biodt/v03/external_models.py` (the external models). Pure Python, no new heavy dependencies,
   deterministic.
2. A registry of every literature constant: `data/external/composition_constants.json`. Each entry:
   `id, value, low, high, unit, meaning, source {authors, year, title, doi_or_url, page_or_table},
   verified: false`. **No literature constant may appear anywhere else in the code.** Use central value plus
   low/high for scenario ranges. If you are not sure of a value, still enter your best value with
   `"note": "uncertain"`; the supervisor verifies every entry against the source afterwards.
3. Eleven experiments (below), run by `reproducibility/v07/composition.py`, writing
   `reviews/v07/composition/<experiment>.json` plus one summary `reviews/v07/composition/summary.json`,
   a LaTeX table `paper/tables/v07_composition.tex` and one figure `paper/figures/v07_composition.pdf`
   (matplotlib, small multiples, one panel per experiment, central value with low–high band).
4. `reproducibility/v07/reproduce.py`: runs everything, records sha256 of inputs (constants registry,
   ontology digest, process.py, composition.py, external_models.py) and of all outputs into
   `reviews/v07/reproducibility.json`; a second run must produce byte-identical JSON (no timestamps inside
   result files; put timings only in reproducibility.json).
5. Tests `tests/test_composition.py` (pytest): contract checks, refusal cases, deductive gate carried
   through, determinism, and one numerical sanity test per experiment. The existing 107 tests must still pass.
6. Theory note: new section in `docs/causal-ontological-theory.md` "7. Composition with external models"
   (≤ 1.5 pages) and a CHANGELOG entry `0.7.0-dev`.

## Non-goals / anti-goals

* Do **not** claim plant growth, biomass, N cycling or net carbon removal from the twin. External-model
  endpoints are labelled with that model as their source, never as twin predictions.
* Do not change v0.6 numbers. All existing reproducibility outputs and the manuscript verifier must be
  unchanged (run `uv run python reproducibility/v06/verify_manuscript.py` and the test suite at the end).
* Do not touch the website or the paper text (only the new table/figure files).
* No randomness except where stated; no network access at run time.

## The contract (composition.py)

```python
@dataclass(frozen=True)
class Port:            # a typed variable
    name: str          # e.g. "nacl_mol_l", "h2_uptake_mol_m2", "crop_relative_yield"
    term: str          # ontology CURIE already used by the twin (ChEBI/UO/ENVO/PATO/SIO) or rak: term
    unit: str          # UO CURIE or plain SI string, compared exactly
@dataclass(frozen=True)
class ExternalModel:
    id: str; name: str; kind: Literal["upstream", "downstream"]
    inputs: tuple[Port, ...]; outputs: tuple[Port, ...]
    validity: dict[str, tuple[float, float]]   # per input port, closed interval; outside → refuse
    constants: tuple[str, ...]                  # ids into composition_constants.json
    fn: Callable[[dict, dict], dict]            # (inputs, constants) -> outputs
    gate: str | None = None                     # ontology class whose reasoner-entailed share must be > 0
```

Chain runner `run_chain(action, upstream: list, twin_request, downstream: list, scenario="central")`:

1. Run upstream models; map their outputs onto the twin intervention dict accepted by
   `interventions.intervene()` (keys `forcing`, `parameters`, `state`, `community`, `coupling`). Port
   matching is by `(term, unit)` equality; mismatch → `CompositionError`.
2. Validity: any input outside its model's interval → refusal object `{refused: true, model, port, value,
   interval}` (not an exception in the experiment runner; it is a result to report).
3. Admissibility: the twin's existing `theory.check_admissible` / `intervene` checks apply unchanged.
4. Downstream gates: if a downstream model has `gate=C` and the reasoner-entailed share of `C` in the chosen
   community is 0, its output is `0` with provenance `deduced` and the derivation string, without
   calling `fn` (mirror `interventions.py`'s `deduced` block). Absence of entailment ≠ biological absence;
   say so in the derivation.
5. Acyclicity: the composed graph (upstream → twin → downstream) must be a DAG; a downstream output that
   is also an upstream input is refused (feedback is out of scope for v0.7).
6. Output: every value is wrapped `{value, low, high, unit, provenance: "twin" | "external:<id>" |
   "deduced", constants: [...]}`; low/high come from rerunning the chain with the constants' low/high
   (one-at-a-time for each uncertain constant, then report the envelope), not from random sampling.
7. Include the twin's `trace` (digests) in each chain result.

## Two small twin extensions (needed by experiments)

* **Atmospheric H2/CO as exogenous forcing.** `h2_mol_m3` and `co_mol_m3` are boundary conditions. Add them
  as admissible **forcing** targets (`h2_mol_m3`, `co_mol_m3` per day, defaulting to the Parameters value)
  so do(C_H2) is an intervention on an exogenous node. Add the node to `theory.graph()` /
  `identification.py` as exogenous parent of J_h2 (resp. J_co), check `check_admissible` accepts it, and add
  a test. Do not change defaults; v0.6 outputs must stay identical.
* **Inoculation = constitutive intervention on the community.** v0.6 treats choosing a community as
  evidence. Inoculation adds organisms. Implement `community.inoculate(samples, mag_id, share)`: build a
  hypothetical community = chosen sample(s) plus MAG `mag_id` at relative share `share` (others rescaled),
  take memberships of `mag_id` from the existing reasoner classification (`load_classification`), and
  recompute shares and γ through the existing code path. Report it as a *new kind* of change in the theory
  note: the TBox K and the evidence E are unchanged; the deduction D_K is recomputed for a counterfactual
  unit. Show the deductive gate opening: sample `51PRr3` (CO gate closed in v0.6) inoculated with a
  `COCarrierMAG` MAG gives γ_co > 0.

## Baseline

Use the same baseline forcing, initial state, parameters and default community that
`reproducibility/v06/analyse.py` uses for its primary scenario (read it; reuse its helpers, don't copy
numbers). State the baseline in each result JSON.

## The experiments (eleven + experiment 12 below)

Each experiment = one function in `reproducibility/v07/composition.py`, with a short docstring giving
the question in plain words. Levels are central values; every literature constant has low/high.

| # | id | question | upstream | twin | downstream |
|---|----|----------|----------|------|------------|
| 1 | `brackish_irrigation` | What if we irrigate with brackish groundwater or treated wastewater? | EC_w (0.5, 2, 4, 8 dS/m) → NaCl mol/L (≈ 10 mmol/L per dS/m, registry) | irrigation pulses with that salinity | soil-salinity ECe = f·EC_w (FAO-29 leaching relation, registry) → Maas–Hoffman relative yield for date palm, barley, alfalfa, tomato (thresholds/slopes from FAO-29 / Maas & Hoffman 1977, registry) |
| 2 | `drip_vs_flood` | Same water, delivered as daily drip or occasional flooding? | equal seasonal total (registry: FAO-56 crop water demand order of magnitude, used only to set the total) split as daily, every 3 days, every 10 days | water pulses | — (contrast with the v0.6 finding that daily pulses give the lowest uptake) |
| 3 | `mulch_shade` | What if we mulch or put up shade nets? | mulch → evaporation reduction fraction; shade net → soil temperature reduction (°C) | evaporation, temperature | — |
| 4 | `biochar` | What if we add biochar? | dose (0, 1, 2, 5 % w/w) → porosity change (meta-analysis, registry); retention unchanged unless registry gives a value | do(porosity) | — |
| 5 | `inoculation` | What if we add microbes that oxidise CO or H2? | inoculum MAG + share (0.1, 1, 5 %) | `community.inoculate` → γ | — ; report gate opening for 51PRr3 and the effect on a sample whose gate is already open |
| 6 | `carbon_fixation_bound` | How much CO2 could the H2 oxidisers fix? | — | baseline and wet scenarios | CO2_fixed ≤ J_H2 · y, y = mol CO2 per mol H2 for knallgas autotrophy (registry, e.g. ~1/6); gate: reasoner-entailed `H2AutotrophicCandidateMAG` share > 0 else 0 by deduction. Also report CO oxidation as CO→CO2 (no net change in atmospheric C atoms; say so). Express µmol and µg C m⁻² per period and per year; put next to one literature context value (e.g. desert biocrust net C fixation, registry, labelled context) |
| 7 | `h2_leakage` | What if a hydrogen economy raises atmospheric H2? | mixing ratio scenarios 1×, 1.5×, 2×, 10× background (registry: background ~0.53 ppm; global scenarios from hydrogen-economy literature; 10× = near-source plume, labelled hypothetical) → mol/m3 | do(C_H2) forcing (new) | — ; report sink per m² and check linearity / oxygen limitation |
| 8 | `solar_park` | What if we build a solar park? | under-panel temperature reduction and evaporation reduction (registry), panel-washing water pulses (L/m² per wash, interval; registry) | temperature, evaporation, water pulses | — |
| 9 | `cloud_seeding` | What if cloud seeding adds rain? | +0, +5, +10, +20 % rain on existing rain days (registry: efficacy uncertain, cite WMO statement) | rain pulses | — ; if the baseline has few rain days, the honest result is a small effect; say so |
| 10 | `compaction` | What if off-road traffic compacts the soil? | passes (0, 1, 10, 50) → bulk-density increase (registry, Mojave off-road vehicle studies) → porosity = 1 − ρ_b/ρ_s | do(porosity) via `process.compact` | — |
| 11 | `warming` | What if the climate warms? | +1.5, +2, +4 °C on temperature; evaporation scaled (registry, % per °C); rain ±20 % | all forcing | — |

### Experiment 12 `inoculated_alfalfa_salinity` — the only plant experiment with field data (KSA)

Question: *What if we irrigate alfalfa with brackish water and coat the seed with the desert bacterium
Enterobacter sp. SA187?* Source: de Zélicourt A, … Hirt H (2018) PLoS Genet 14:e1007273,
doi:10.1371/journal.pgen.1007273 — field trial at Hada Al-Sham near Jeddah, KSA; irrigation groundwater
EC 3.12 (low) and 7.81 (high) dS/m; soil EC 1.95 dS/m; two seasons (2015–16, 2016–17); n = 4 / 3 plots.
Registry entries (all with this source, `page_or_table` as given):

* `sa187_effect_fresh_low_s1` = 0.16, `…_s2` = 0.12; `sa187_effect_dry_low_s1` = 0.14, `…_s2` = 0.17
  (Results text, relative increase at 3.12 dS/m; high salinity "similar"). Use central 0.15, low 0.12,
  high 0.17 as the multiplicative benefit `b` (yield × (1+b)).
* Digitised mock yields, **season 1 only** (season 2 confounded by exceptional rain, per the authors),
  t/ha, from Fig 1, read by eye ±0.15, `"note": "digitised from figure, not a published table"`:
  fresh low 5.4, fresh high 4.3; dry low 2.1, dry high 1.75. Observed mock ratio high/low:
  fresh ≈ 0.80, dry ≈ 0.83.
* Maas–Hoffman alfalfa threshold 2.0 dS/m ECe and slope 7.3 % per dS/m (FAO-29, Ayers & Westcot 1985,
  Table 4) — verify.
* FAO-29 soil-water to saturation-extract relation EC_sw ≈ 2 · ECe (soil water at field capacity) — verify.

Chain (three alternative root-zone salinity models, report all three):

* **A (FAO fixed leaching):** ECe = 1.5 · EC_w (FAO-29, 15–20 % leaching fraction).
* **B (twin physics):** run the twin's salt/water bucket with irrigation at EC_w (NaCl mol/L via registry
  conversion), using an irrigation schedule representative of drip alfalfa (registry; state it), and convert
  the twin's mean soil-solution concentration over the season to EC_sw, then ECe = EC_sw / 2. Label clearly
  that the twin's soil and depth are Rub' al Khali sand, not Hada Al-Sham soil.
* **C (fitted):** ECe = f · EC_w with one free factor f fitted so that the Maas–Hoffman ratio
  RY(7.81)/RY(3.12) matches the observed season-1 mock ratio (fit separately for fresh and dry; closed-form
  or bisection, deterministic). Report f and its range from the ±0.15 t/ha reading error (propagate to the
  ratio by taking the extreme combinations).

Then relative yield RY = 1 − slope · max(0, ECe − threshold), floored at 0; inoculated RY = RY · (1+b).
Report for EC_w = 0.5, 2, 3.12, 4, 7.81, 8, 12 dS/m: RY mock and inoculated, under A, B, C; and the
predicted vs observed high/low ratio for A and B (a validation: how far is each from the field result).
Plain-language text must say: the bacterium's +12–17 % comes from one field trial at two salinities; whether
the benefit grows with salinity cannot be told from these data; the twin contributes only the soil salt
balance; plant growth is the crop model's and the field trial's, not the twin's.

Each experiment JSON: `question`, `chain` (list of model ids with kind), `levels`, per level the
wrapped outputs (`h2_uptake`, `co_uptake`, active days, final water/salinity, any downstream endpoints),
refusals, deduced nulls, `plain_language` (2–3 sentences a member of the public can read, no jargon,
stating what is the twin's and what is another model's), and `caveats` (list).

## Done means

* `uv run pytest -q` all green (107 old + new).
* `uv run python reproducibility/v07/reproduce.py` run twice → identical result files (check sha256).
* `uv run python reproducibility/v06/verify_manuscript.py` passes; v0.6 outputs unchanged
  (`git diff --stat reviews/v06 paper/tables/v06_* ` empty after rerunning v0.6 analysis is not required —
  just do not modify them).
* Commit on branch `v0.7-composition` in small commits (registry; twin extensions; composition layer;
  experiments; tests; docs). Do not push.
* Write `reviews/v07/composition/IMPLEMENTATION_NOTES.md`: what you built, every judgement call you made,
  anything you could not do, and for each registry constant how confident you are.

## Left to the builder

Internal structure, plotting style (follow `reproducibility/v06/figures.py`), names of helper functions.
