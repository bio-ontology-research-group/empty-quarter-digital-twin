# v0.8.0 — Full uncertainty analysis and question catalogue — 1 November 2026

- Add the twin process-parameter uncertainty registry `data/external/twin_parameter_ranges.json` (13 parameters,
  all `modeller_assumption`) and `src/eq_biodt/v03/uq_inputs.py` to merge external constants (`triangular` over
  `[low, value, high]`) and non-clamped twin parameters into a deterministic input vector per scenario level.
- Add the unified 14-question catalogue `src/eq_biodt/v03/questions.py` (111 scenario levels across 12 composed
  what-if chains plus rainfall splitting and the deductive CO gate).
- Add Saltelli (2010) / Jansen (1999) quasi-Monte Carlo uncertainty and global sensitivity analysis in
  `src/eq_biodt/v03/uq.py`: scrambled Sobol' design with individual and verification-status group re-sampling
  blocks (`verified_in_source`, `modeller_assumption`, `source_not_checked`), 500-draw bootstrap CIs, paired
  per-row contrasts and percentage changes, refusal indicator sensitivity when refusals reach 5 %, exact
  deductive zeros (`[0, 0]`) for closed ontological gates, and $N/2$ vs $N$ quantile convergence checks.
- Add cluster array task runner and deterministic 6-significant-digit gatherer (`reproducibility/v08/{tasks,uq_task,gather}.py`)
  with full $N=1,024$ and $N=2,048$ rerun summaries in `reviews/v08/uq/` and cluster execution metadata in
  `reviews/v08/uq/run_record.json` (2,285,568 base evaluations + 171 rerun tasks).
- Add paired plain-language and scientific question builders (`reproducibility/v08/{texts,questions}.py`,
  `reviews/v08/questions/*.json`), paper sensitivity figure and tables (`reproducibility/v08/figures.py`,
  `paper/figures/v08_sensitivity.pdf`, `paper/tables/v08_{parameters,uq}.tex`), and `reproducibility/v08/reproduce.py`
  (`reviews/v08/reproducibility.json`).
- Update `scripts/web/export_web_bundle.py` to copy `reviews/v08/questions/*.json` into `data/questions/` of the
  web bundle, and extend `reproducibility/v06/verify_manuscript.py` to verify all v0.8 manuscript numbers and tables.

# v0.7.0 — Chained what-ifs: composition with external models — 31 October 2026

- Add `composition.py` and `external_models.py`: chained what-ifs `action -> upstream -> twin -> downstream`
  with typed ports, validity domains as refusals, deductive gates carried through, acyclicity check and
  one-at-a-time low/high envelopes. Values carry provenance (`twin`, `external:<id>`, `deduced`).
- Add the registry of literature constants `data/external/composition_constants.json`; each entry was checked
  independently and carries a status (23 verified in source, 10 modeller assumptions, 9 source not checked).
  Corrections: H2 background DOI, tomato slope (9.5, range 9.0–9.9), alfalfa slope cited to FAO-56, FAO-29
  leaching range 1.0–2.1, evaporation scaling bound.
- Add twelve experiments (`reproducibility/v07/composition.py`, `reviews/v07/composition/`), figure, table and
  `reproducibility/v07/reproduce.py` with input/output digests; two runs are byte-identical.
- Add a reproducible digitisation of the SA187 alfalfa field trial (de Zélicourt et al. 2018, Fig 1), pinned by
  hash and checked against the printed percentages; fit of the root-zone/irrigation salinity factor.
- v0.6 scripts read a v0.6 view of the graph (`reproducibility/v06/v06_graph.py`) so v0.6 outputs are unchanged;
  the d-connecting witness path is now deterministic.
- Atmospheric H2 and CO concentrations (`h2_mol_m3`, `co_mol_m3`) are exogenous forcing targets, per day,
  defaulting to the Parameters value; new graph nodes `C_h2`, `C_co`. v0.6 results are unchanged when the keys
  are absent.
- Add `community.inoculate` and an `inoculation` intervention key: a constitutive intervention on the community
  (K and E unchanged, D_K recomputed for a counterfactual unit).
- Theory note section 7 "Composition with external models" (scope is now section 8).

# v0.6.0 — Stratified causal–ontological twin — 30 October 2026

- Add the stratified causal–ontological graph (`theory.py`, `GET /v3/graph`,
  `eq-biodt graph`): every node carries its zone, manipulability and equation.
  The theory is stated in `docs/causal-ontological-theory.md`.
- Couple rate coefficients to the selected community through reasoner-entailed
  carrier shares (`community.py`; `reference`, `gate`, `proportional`;
  `GET /v3/community/{sample}`, `eq-biodt community`).
- Check admissibility of interventions, with zone-specific errors for ontology,
  evidence, entailment, constituted and outcome targets and for infeasible regimes.
- Add state interventions (`water_mm`, `theta`, `psi_mpa` with a matric or osmotic
  version, `nacl_mol_m2`), version relevance (`POST /v3/version-relevance`) and a
  treatment-version block in `compare`.
- Replace hard-coded invariance flags by runtime checks of memberships, shares,
  constitutive identities and balances in both arms; report deductive null effects.
- `eq-biodt intervene` now accepts top-level `community`, `coupling` and `state`
  as the API does.
- Ontology: verified GO/ChEBI/ENVO/PATO/UO/SIO annotation mappings
  (`data/v03/ontology_mappings.tsv`), measurement endpoint classes and declared
  observation maps (`observation_maps.json`), runtime verification of ontology and
  reasoner digests, and justifications for all 136 candidates and 198 CO carriers.
- Add do-calculus identification with deductive value closure
  (`identification.py`, `POST /v3/identify`, `GET|POST /v3/closure`,
  `eq-biodt identify`). It derives the deductive CO null as Rule 3, the
  independence of hydrogen uptake from the CO pathway where no CO carrier is
  mapped, and back-door identification of the water effect under a latent season.
- Rebuild the local console around the pipeline. It adds paired regimes,
  treatment versions, a sampling map, reasoner justifications, a ψ-version
  example, an admissibility explorer, transfer tests and a glossary. New
  endpoints are `/v3/mags/{id}`, `/v3/transfer` and `/v3/admissibility`.
- Check every number in the manuscript against the outputs
  (`reproducibility/v06/verify_manuscript.py`, also run in CI).
- Export the public web page bundle (`scripts/web/export_web_bundle.py`). The
  page at https://rubalkhali.science/model has a plain-language "What if…?"
  view, a research view, and a guide with a glossary at `/model/guide`.
- Move the v0.2 exploratory model to `legacy/` and remove it from the package.
- Add `reproducibility/v06` analyses and figures; revise the manuscript.
- Model version reported by the API, contract and outputs now comes from package
  metadata. `/v3` paths, `eq_biodt.v03` and `data/v03` keep their names.

## v0.5.0 — Causal–ontological desert twin and external transfer

- Reframe the manuscript around ontological deduction, causal intervention
  mechanisms, primary Rub’ al-Khali data and independent environmental tests.
- Add verified Namib and Sonoran numerical datasets, campaign-held-out field
  prediction, held-soil temperate tests and descriptive cross-desert marker data.
- Preserve failures of the frozen parameterization and distinguish environmental
  calibration from transfer without calibration. Publish all held-out predictions.
- Trace primary sampling coverage, genome-linked overlap bounds and a worked
  deductive reaction-admission calculation, with a same-information flat control.
- Add three publication figures, expanded supplementary methods, a single offline
  reproduction command and an independent numerical audit.
- Confirm ENA PRJEB125153; raw-read public release is scheduled for late November
  2026. Derived inputs remain available for reproduction.
- Retain v03 Python and API paths. Process equations and default parameters are
  unchanged from v0.4; release metadata advances to v0.5.

# Changes

## v0.4.0 — 28 October 2026

- Implement the full KEGG M00165 alternative regeneration branch. Counts change
  from 83 to 232 CBB candidates, 23 to 136 hydrogenase-plus-CBB candidates, and
  20 to 100 in the declared quality subset. An independent expression parser
  verifies all 5,438 assignments against the pinned module definition.
- Replace unsupported biological false-positive and archaeal Rubisco-function
  claims with annotation-support comparisons. Two archaeal Rubisco-positive
  MAGs also have phosphoribulokinase annotations.
- Use verified exported OWL memberships for resource admission and independently
  gate CO oxidation. Reject stale or modified classification artifacts.
- Integrate gas uptake at subdaily midpoints. Short pulses have positive uptake;
  single-versus-four-pulse rankings reverse across parameter settings. Verify
  the reversal at four times the default temporal resolution.
- Add paired, version-specific intervention queries with immutable evidence and
  definitions, input validation, and same-input flat-calculation checks.
- Report empirical residual intervals without a distribution-free guarantee.
  Add paired site-bootstrap forecast comparisons, compartment coverage and
  matched site--campaign ecological contrasts.
- Bundle all current analysis inputs; pin the reasoner binary and Python lock;
  add CI, package release metadata, installation instructions and source attribution.
- Regenerate the manuscript, tables and figures and state the experimental
  validation and analysis-reproduction boundaries explicitly.

The earlier v0.2 implementations and their reproduction material are kept under `legacy/` (unpackaged, untested).
The local v0.3 manuscript and implementation are preserved in commit `06347de`.
