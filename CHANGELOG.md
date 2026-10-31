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
