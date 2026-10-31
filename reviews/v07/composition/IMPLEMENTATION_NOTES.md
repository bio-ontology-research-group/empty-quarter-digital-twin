# v0.7 composition experiments: implementation notes (part 2)

Branch `v0.7-composition`. Part 1 (composition layer, external models, registry, twin extensions) is unchanged except for the SA187 mock-yield entries (below).

## What was built

* `reproducibility/v07/composition.py`: twelve experiment functions (docstring = the question). Writes `reviews/v07/composition/<id>.json` and `summary.json`. Deterministic: `sort_keys`, no timestamps, floats rounded to 6 significant digits. Each JSON has `question, chain, baseline, levels` (per level: wrapped `twin`, `twin_baseline`, `difference_vs_baseline` with percent, `upstream`, `downstream`, refusals, deduced nulls, envelope refusals, trace digests), `refusals, deduced_nulls, plain_language, caveats, headline, figure` (plot data), `constants_used, registry_sha256`.
* `reproducibility/v07/figures.py`: `paper/figures/v07_composition.{pdf,png}` (3x4 small multiples, central line with low-high band, v0.6 style) and `paper/tables/v07_composition.tex` (booktabs, one row per experiment). `main.tex` / `supplementary.tex` untouched.
* `reproducibility/v07/reproduce.py`: runs `digitise_fig1.py`, `composition.py`, `figures.py`; writes `reviews/v07/reproducibility.json` (sha256 of inputs and outputs; timings only there).
* `tests/test_v07_experiments.py`: one numerical sanity test per experiment (run into a temp dir), plus a test that rerunning reproduces the committed JSON exactly, plus the registry-mirror test. Total 162 tests pass (150 before).
* v0.6 compatibility: `reproducibility/v06/v06_graph.py` (see below).

## Baseline (all experiments)

The v0.6 primary scenario: `reproducibility/v06/analyse.py` `reference_schedules()[1]` ("single_pulse": 14 days, 12 mm on day 0, evaporation 1 mm/d, 25 C, no NaCl), `Parameters(k_h2_s=.001, k_co_s=.001)`, initial water 1.8 mm, no community (gamma = 1) unless the experiment names one. The "season" of the irrigation experiments is therefore these 14 days.

## Judgement calls

1. **Experiment 12 reads the pixel digitisation.** `external_models.digitised_trial()`, `mock_yield()`, `digitisation_error()`, `observed_mock_ratio()` read `data/external/de_zelicourt_2018/fig1_digitised.json`. The four eye-read `sa187_mock_*_t_ha` registry entries are removed (not just marked) and replaced by `sa187_mock_ratio_fresh` / `_dry` (0.7958 / 0.8282, low-high = digitisation error). The fit (`fit_leaching_factor`) needs the ratio as a registry constant so that the chain's one-at-a-time envelope can vary it; a test asserts the registry mirror equals the value computed from the digitised file (to 6e-5), and `observed_mock_ratio` and the experiment's observed yields read the file directly.
2. **Error choice.** The eye-read +-0.15 t/ha is dropped: it is 10x the digitiser's measured calibration residual and cannot be justified against the pixel data. I use absolute 0.0148 t/ha (largest y-axis calibration residual of the two panels) plus relative 0.0083 (largest disagreement between the digitised SA187 effects and the percentages printed on the figure), and take extreme combinations for the ratio: fresh 0.7958 [0.7779, 0.8141], dry 0.8282 [0.8018, 0.8554]. This covers reading error only, not the sampling error shown by the figure's error bars (n = 3-4 plots), which cannot be recovered from the published data; so the f range is a lower bound on its uncertainty.
3. **Gate vs proportional in experiment 5** (supervisor request): both modes are reported. Gate = whether the process can run (any share > 0 gives gamma = 1, so no dose dependence); proportional = how much (gamma = q/median(q)). Survival/establishment of the inoculum is stated as not modelled.
4. **Experiment 5 samples and MAGs**: closed gate `51PRr3`; already-open sample = the sample with the largest COCarrierMAG share (`V53PRr2`); CO carrier MAG = first (sorted) reference MAG entailed to COCarrierMAG (`11Dr3__single_SemiBin_57`); negative-control MAG = first MAG with no entailed class (`10Dr2__single_SemiBin_0`).
5. **Experiment 6**: samples `14PRr2` (H2AutotrophicCandidateMAG share 0, deduced zero) and `F14Dr1` (share 0.0264, closest to the median share). "Wet" = a single 36 mm pulse (3x the baseline). Annualised linearly (x 365.25/14), labelled extrapolation. Context value (biocrust 10 g C m-2 yr-1, registry) is labelled context.
6. **Experiment 7** oxygen check: maximum over days of oxygen used / `oxygen_supply_mol_m2_s` from the intervened trajectory.
7. **Experiment 12 chain B** uses the twin's baseline forcing (bare sand, 1 mm/d evaporation, 14 days) with daily irrigation at the registry alfalfa demand (8 mm/d). I tried evaporation = ETc = 8 mm/d: then irrigation equals evaporation, nothing drains, and the soil solution reaches > 1 mol/L (twin abstains / chain refuses). With 1 mm/d the bucket leaches ~85 % of the water, so B strongly under-estimates root-zone salinity relative to a real stand; B's ECe at 3.12 dS/m is 1.44. Stated in the caveats.
8. **Headline definitions** (summary.json): twin experiments use the percent change in 14-day H2 uptake relative to the baseline at a stated reference level (central, envelope low/high); experiments 1 and 12 use the crop endpoint (external model); 5 uses CO uptake; 6 the carbon bound.
9. **v0.6 outputs kept identical** (see below).

## Results (headline per experiment; central [low, high])

| # | experiment | headline | central [low, high] | provenance |
|---|---|---|---|---|
| 1 | brackish_irrigation | alfalfa relative yield, EC_w = 8 dS/m, daily irrigation | 0.887 [0.800, 0.938] | external:maas_hoffman_alfalfa (salinity from the twin) |
| 2 | drip_vs_flood | H2 uptake vs unirrigated baseline, every 10 d | -5.55 % [-5.55, -5.55] (daily: -38.4 % [-63.4, -17.7]; every 3 d: -37.3 %) | twin |
| 3 | mulch_shade | H2 uptake, full mulch | +50.1 % [36.3, 53.9]; full shade net -24.8 % | twin |
| 4 | biochar | H2 uptake, 5 % w/w | -2.26 % [-4.24, -0.64] | twin |
| 5 | inoculation | CO uptake of 51PRr3 + 1 % CO carrier (gate) | 144.4 umol m-2 / 14 d (baseline exactly 0) | twin |
| 6 | carbon_fixation_bound | upper bound C fixed, F14Dr1, baseline | 5.76e4 ug C m-2 yr-1 [3.45e4, 8.64e4] (= 183.7 umol CO2 m-2 per 14 d); 0 deduced for 14PRr2 | external:knallgas_carbon_fixation |
| 7 | h2_leakage | H2 sink at 2x background | +100 % (exactly linear; 10x: +900 %) | twin |
| 8 | solar_park | H2 uptake, full panel cover | +16.6 % [-1.9, 30.5] | twin |
| 9 | cloud_seeding | H2 uptake, +10 % rain on rain days | +12.7 % [0, 24.7] | twin |
| 10 | compaction | H2 uptake, 50 passes | -1.49 % [-7.58, +1.00] | twin |
| 11 | warming | H2 uptake, +4 C, rain unchanged | +7.16 % [-2.89, 15.1] | twin |
| 12 | inoculated_alfalfa_salinity | inoculated relative yield at 7.81 dS/m, fitted salinity (fresh) | 0.915 [0.891, 0.936] | external:sa187_benefit |

Details:

* **Exp 5.** 51PRr3: CO gate closed in v0.6 (CO uptake 0); with any carrier share the gate opens (gamma_co 0 -> 1; CO uptake 144.4 umol m-2 per 14 d, the reference value). `V53PRr2` (gate already open): no change in gate mode (`gate_opened` false). Non-carrier MAG: exactly no change in gate mode (CO and H2). Proportional mode, 51PRr3 + CO carrier: CO uptake 10.4 / 77.4 / 221 umol m-2 at 0.1 / 1 / 5 %.
* **Exp 6.** Per level the JSON gives umol CO2 per period and per year and ug C per period and per year, the gate result for `14PRr2` (provenance `deduced`, derivation string) and `F14Dr1`, the context value labelled context (the bound is 0.58 % of the biocrust figure at baseline), and CO oxidation as CO -> CO2 with net atmospheric C change 0.
* **Exp 7.** uptake / multiplier is the same at all four levels to 2e-16; peak oxygen used / supply = 7.3e-6, so oxygen does not limit. Linearity is a property of the first-order uptake law.
* **Exp 12.** Fitted ECe = f x EC_w: f(fresh) = 0.614 (0.582-0.648), f(dry) = 0.557 (0.510-0.604). f x 3.12 = 1.92 (1.82-2.02) fresh, 1.74 (1.59-1.88) dry, versus the measured soil EC 1.95 dS/m (independent; the registry does not say which EC definition) and model A's 1.5 x 3.12 = 4.68, model B's 1.44. Predicted vs observed mock ratio (7.81 / 3.12): A 0.362 [0.008, 0.627]; B 0.883 [0.795, 0.935]; observed 0.796 fresh / 0.828 dry. A (FAO 1.5 leaching factor) is far from the field result; B (twin bucket) is within 0.09 but for the wrong reason (see call 7); C matches by construction. Inoculated RY (A / B / C-fresh) at 7.81 dS/m: 0.33 / 1.02 / 0.92. Whether the bacterium's benefit depends on salinity cannot be told (one trial, two salinities).

## Surprises and how I checked them

* **Irrigation lowers H2 uptake by 38-44 % (exp 1, 2).** Not a bug: 84-96 mm of water keeps the sand wet and gas diffusion limits uptake (consistent with the v0.6 finding that daily pulses are worst; the baseline dries out after one 12 mm pulse). Final water 24.6 mm daily vs 19.9 mm every 10 d.
* **Exp 2, every 10 days: low = high.** The two big pulses (7 x demand each) saturate the soil, so the registry range of the water demand does not change the result.
* **Compaction is non-monotone and biochar lowers uptake (exp 4, 10).** Checked by a direct porosity scan with `simulate` (no chain): uptake peaks near porosity 0.36 (1112 umol m-2) and falls on both sides (0.30: 1079; 0.46: 1044) because lower porosity wets the same water volume (more activation) while raising diffusion resistance. Effects are small and the compaction envelope spans zero.
* **Mulch +50 % (exp 3).** The baseline soil dries out after the pulse (12 of 14 active days); halving evaporation keeps it active all 14 days. Specific to the single-pulse baseline.
* **Non-carrier inoculum in proportional mode lowers H2 uptake by ~4.9 % at 5 %**: renormalisation of the shares (the inoculum displaces the community), a model artefact, noted in the caveats; under the gate it is exactly zero.
* **Exp 1 low/high of H2 uptake is None at 8 dS/m**: with the registry NaCl factor at its high end the soil solution leaves the ideal-dilute domain on some day and the twin abstains.
* **Pre-existing nondeterminism found**: `identification.d_connecting_path` iterated sets, so the witness open back-door path in `reviews/v06/identification.json` flipped between `D_h2@1` and `D_co@1` with PYTHONHASHSEED. I sorted the neighbour order (both are valid d-connecting paths). The committed file had `D_h2@1`; the regenerated, now deterministic, file has `D_co@1` (one line; conclusions unchanged; not quoted in the paper; `verify_manuscript.py` passes).

## Reproducibility

* `uv run python reproducibility/v07/reproduce.py` run twice: sha256 of all 12 experiment JSONs, summary.json, the .tex table, **and the PDF and PNG** are byte-identical (PDF made deterministic with `CreationDate/ModDate/Producer = None`).
* v0.6: `reproducibility/v06/reproduce.py --skip-v05` re-run. `git status --short reviews/v06 paper` shows no change except: `results.json` (contains `runtime_seconds`, differs on every run) and `manifest.json` (timings and input digests; the digests were already stale after part 1), which I reverted rather than committed; and the one-line `identification.json` witness above (committed). `admissibility.tsv`, all `paper/` files unchanged. This needed `reproducibility/v06/v06_graph.py` (the v0.6 view of the graph without `C_h2`/`C_co` and with the old `C_a` label), used by `v06/analyse.py` and `v06/figures.py`; otherwise figures.py raised `KeyError: 'C_h2'` and the admissibility macros would have changed.
* `reproducibility/v06/verify_manuscript.py` passes; `uv run pytest -q`: 162 passed.

## Open problems

* All registry constants are still `verified: false`; every number inherits that. The supervisor/reviewer corrections should only need a rerun of `reproduce.py` (the mirrored `sa187_mock_ratio_*` entries are the only ones tied to a data file; the test will flag drift).
* The 14-day single-pulse window is short for "seasonal" questions; annualisation (exp 6) is linear scaling of one wetting event.
* Model B of exp 12 is a bare-sand bucket without root uptake; a proper test needs an evapotranspiration-driven salt balance with a leaching fraction, outside the twin's scope.
* The inoculum's survival and the true dose are not modelled; gate mode gives no dose-response by construction.
* Digitisation error excludes the figure's sampling error bars.
* `reviews/v06/manifest.json` is stale relative to the part-1 source changes; regenerate with the full v0.6 pipeline when the release is cut.
* Table column widths in `v07_composition.tex` assume a full-width environment; not compiled into the manuscript (main.tex untouched).
