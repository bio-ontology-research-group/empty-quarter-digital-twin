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

Earlier v0.2 implementations and their reproduction material remain historical.
The local v0.3 manuscript and implementation are preserved in commit `06347de`.
