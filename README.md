# EQ-BioDT v0.6.0

EQ-BioDT is a causal–ontological digital twin of microbial processes in the
Rub’ al-Khali, the world’s largest continuous sand desert. It connects primary
genome-resolved observations, deductive pathway definitions and causal
intervention calculations.
It provides an interactive console, HTTP API and command-line interface.

The release contains 551 sample records and 5,438 reference MAGs. The complete
KEGG M00165 annotation rule, including its transaldolase alternative, supports
232 CBB candidates; 136 also carry screened hydrogenase markers. These are
annotation-supported candidates, not measurements of activity or carbon removal.


The manuscript is **A causal–ontological digital twin of the world’s largest
continuous sand desert**. Its principal external tests use independent Namib
and Sonoran observations and contrasting temperate soils. The frozen response
parameters do not transfer uniformly: environmental calibration improves
Sonoran held-campaign RMSE to 0.69 nmol H₂ m⁻² s⁻¹, compared with 1.65 for the
reference model and 0.79 for a training mean. Dry Namib activity and the Sonoran
high-moisture response expose failures of the reference assumptions.

Primary shotgun reads are deposited under [PRJEB125153](https://www.ebi.ac.uk/ena/browser/view/PRJEB125153),
with public release scheduled for late November 2026. The derived inputs needed
for the analyses are bundled, independently of that raw-read release schedule.
See [external data and measurement definitions](data/external/README.md),
[the retrospective protocol](reviews/v05/external-evaluation-protocol.md),
and [the independent numerical audit](reviews/v05/independent_audit.json).

## Use it online

The twin runs in the browser at <https://rubalkhali.science/model>, with no
installation. The released `process.py` and `theory.py` run there through
self-hosted Pyodide, and the page checks them against the release manifest.

- **What if…?** (default) shows interventions in plain terms: one heavy rain
  versus drizzle with the same total, salty water, heat, compacted sand, and
  places where the genes for carbon-monoxide uptake were not found.
- **Research view** (`/model?view=research`) shows the sampling footprint, the
  causal–ontological graph, live paired `do()` regimes, treatment versions,
  community coupling, reasoner justifications, identification with deductive
  closure, transfer tests and ontology mappings.
- **How to read the twin** (`/model/guide`) explains every term, with an A–Z
  glossary. The same definitions appear in the info buttons of the research view.

The page data are exported by `scripts/web/export_web_bundle.py`, which records
the source commit; `tests/test_web_bundle.py` checks the export against the
analysis outputs and the manuscript.

## Install and open the tool

Install [uv](https://docs.astral.sh/uv/getting-started/installation/), then:

```bash
git clone https://github.com/bio-ontology-research-group/empty-quarter-digital-twin.git
cd empty-quarter-digital-twin
git checkout v0.6.0
uv sync --locked
uv run eq-biodt serve --port 8000
```

Open <http://127.0.0.1:8000> for scenarios and evidence inspection, or
<http://127.0.0.1:8000/docs> for the interactive API. The default server binds
to localhost. Python 3.13 is the tested reproduction environment; package
metadata permits Python 3.10 or later.

The [GitHub release](https://github.com/bio-ontology-research-group/empty-quarter-digital-twin/releases/tag/v0.6.0)
also provides an installable wheel and manuscript PDFs. Install the downloaded
wheel with `uv tool install path/to/eq_biodt-0.6.0-py3-none-any.whl`, then run
`eq-biodt serve`. Normal tool use needs neither a separate ontology reasoner
nor a sibling repository.

## Run calculations

```bash
uv run eq-biodt contract
uv run eq-biodt simulate --input data/v03/example_scenario.json
uv run eq-biodt compare --input data/v03/example_comparison.json
uv run eq-biodt intervene --input data/v03/example_intervention.json
uv run eq-biodt resource-ceiling --input data/v03/example_resource.json
uv run eq-biodt backtest
```

`simulate` integrates a supplied water/NaCl history and conditional gas uptake.
`compare` evaluates four explicit objectives across a 24-member sensitivity
grid. `intervene` evaluates paired regimes with shared initial state and fixed
definitions; see [intervention semantics](docs/causal-design.md).
`resource-ceiling` uses a MAG's verified reasoner entailments to admit each
reaction. It requires independent substrate and energy bounds.

The console allows selection of every reference MAG and downloads complete
scenario inputs, assumptions and outputs. The API exposes the same calculations:
`/v3/simulate`, `/v3/compare`, `/v3/intervene`, `/v3/resource-ceiling`,
`/v3/evidence`, `/v3/mags`, `/v3/ontology` and `/v3/contract`.
The API additionally exposes `/v3/graph`, `/v3/community/{sample}`,
`/v3/version-relevance`, `/v3/experiment` and `/v3/carbon-account`.
The `/v3` API paths, the `eq_biodt.v03` Python package and the `data/v03`
directory keep their names across v0.4 to v0.6 for compatibility; they do not
denote the model version.

## Causal–ontological design

The full statement is in [docs/causal-ontological-theory.md](docs/causal-ontological-theory.md);
the endpoint reference is [docs/causal-design.md](docs/causal-design.md).

- **Three kinds of change.** An observation update changes evidence, a model
  revision changes the ontology or a mechanism, and an intervention changes only
  exogenous inputs or designated physical states. Choosing a community (a sample
  or site visit) selects evidence; it conditions the twin and is not an intervention.
- **Two zones.** The ontological zone holds the OWL classes, evidence, reasoner
  entailments and constitutive definitions (for example `theta = W/(1000 L)`). The
  causal zone holds the mechanisms for water, salt, activation and uptake. An
  intervention replaces causal-zone assignments only.
- **Admissible targets.** Forcing inputs, prescribed soil parameters
  (`porosity`, `psi_entry_mpa`, `retention_b`) and the states `water_mm` and
  `nacl_mol_m2`. `theta` is translated to `water_mm`; `psi_mpa` needs a named
  version (`matric` or `osmotic`). Ontology, evidence, entailments, activation,
  rate coefficients and uptake cannot be set, and infeasible regimes (for example
  stored water above pore capacity) are rejected with a zone-specific error.
- **Community coupling.** `reference` (default) leaves rate coefficients unchanged;
  `gate` sets a gas's coefficient to zero when the selected community has no
  reasoner-entailed carrier (`H2CarrierMAG`, `COCarrierMAG`); `proportional` scales
  it by the carrier's read-weighted share relative to the median sample. The
  proportional scaling is an assumption, not a validated law.
- **Version relevance.** Whether a constituted target such as `psi` has a
  version-dependent effect on an outcome is decided from the graph: it does if some
  physical handle reaches the outcome by a path that avoids the target.
- **Deductive null effects.** If a gate closes, uptake of that gas is zero in every
  admissible regime. The response reports this as `zero_by_deduction` with its
  derivation, checks that memberships and shares are identical in both arms, and
  does not read absent entailment as biological absence.

Worked examples (inputs in [docs/examples](docs/examples)):

```bash
# stratified graph: nodes with zone, manipulability and equation; edges with kind
uv run eq-biodt graph | head
curl -s http://127.0.0.1:8000/v3/graph | head -c 400

# reasoner-entailed carrier shares and coupling factors for one sample
uv run eq-biodt community --sample 51PRr3
curl -s http://127.0.0.1:8000/v3/community/51PRr3

# is the effect of do(psi) on H2 uptake version-dependent?
curl -s -X POST http://127.0.0.1:8000/v3/version-relevance \
  -H 'content-type: application/json' -d @docs/examples/version_relevance.json

# paired intervention on a selected community with gate coupling
uv run eq-biodt intervene --input docs/examples/intervention_community_gate.json
curl -s -X POST http://127.0.0.1:8000/v3/intervene \
  -H 'content-type: application/json' -d @docs/examples/intervention_community_gate.json
```

State interventions and a ψ target use the same request, for example
`"state": {"0": {"psi_mpa": -2.0, "version": "matric"}}` (also `water_mm`, `theta`,
`nacl_mol_m2`). `POST /v3/compare` accepts `community`, `coupling` and reports the
treatment-version block.

For new observations, issue a forecast with a future target date:

```bash
uv run eq-biodt freeze --target-date YYYY-MM-DD --output forecast.json
uv run eq-biodt score --forecast forecast.json --observations new_observations.tsv
```

Replace `YYYY-MM-DD` with a date after issuance. A forecast cannot overwrite
an existing file. Validation observations require `site`, `compartment`,
`collection_date`, and `total_h2_marker_ge`; matching observations must occur
after issuance and within the forecast's target window. These predictions
concern marker genome equivalents, not gas flux. The nominal 90% intervals
are empirical training-residual intervals, with no guaranteed coverage under
campaign shift.

## Reproduce the paper

The analysis starts from the released MAG annotations, mapping summaries,
marker counts and curated metadata. All these inputs, plus the KEGG definition,
published bottle observations and external desert supplements, are in
[data/sources](data/sources/README.md) and [data/external](data/external/README.md).
It does not rerun raw sequencing, assembly or primary annotation.

On Linux x86_64:

```bash
uv sync --locked
uv run python scripts/install_reasoner.py
uv run python reproducibility/v06/reproduce.py
uv run pytest -q
```

The installer downloads Kobayashi-MaRust v1.5.3 and verifies its pinned SHA-256.
After dependencies and this binary are installed, the analysis runs offline.
It rebuilds evidence from the source snapshot, classifies the ontology,
exports verified memberships, runs the independent pathway audit, matched
compartment comparisons, campaign backtests, bottle analysis, numerical
convergence checks, scenario ensemble, primary-provenance trace, cross-desert
transfer evaluations and independent numerical audit, and regenerates paper
tables and figures. The v0.6 command runs the v0.5 reproduction first (`--skip-v05` reuses its outputs)
and then the v0.6 admissibility, treatment-version and paired-regime analyses. Outputs are in
[reviews/v06](reviews/v06/manifest.json),
[reviews/v05](reviews/v05/reproduction_manifest.json) and
[reviews/v03](reviews/v03/report.md).

For other platforms, use the Linux GitHub Actions workflow or build the same
[km release](https://github.com/bio-ontology-research-group/kobayashi-marust/releases/tag/v1.5.3)
and set `KM_BIN` to its executable. The report records the binary digest.

With `latexmk` and TeX Live installed:

```bash
cd paper
latexmk -pdf -interaction=nonstopmode -halt-on-error main.tex
latexmk -pdf -interaction=nonstopmode -halt-on-error supplementary.tex
```

For a separate output directory, `eq-biodt rebuild-evidence --output DIR`
recreates the evidence tables from bundled inputs without overwriting the
installed package. `eq-biodt reason --ofn FILE` regenerates and classifies the
ontology in a user-selected location. For alternative datasets,
`eq-biodt build-evidence --help` lists all explicit required inputs.

## Scientific scope

- Logical classification establishes consequences of supplied annotations
  and axioms. A separately implemented KEGG-expression parser gives identical
  pathway calls with the same information.
- Physical kinetics and retention parameters are supplied assumptions until
  locally measured. The software is a digital-twin prototype without an
  automated field-state assimilation loop.
- Hydrogen uptake is integrated within each day. Wetting rankings can reverse
  across the parameter grid; there is no universally preferred schedule.
- A published bottle dataset tests a moisture envelope, with its threshold
  informed by that same study. It does not validate the complete soil simulator.
- High-salinity gas outputs and unsupported genomic resource ceilings return
  unavailable. Water and salt ledgers remain explicit.
- Net carbon removal requires independently measured stocks, sources,
  counterfactual change and lifecycle emissions.

See [the scientific audit](reviews/v03/report.md),
[changes in v0.6](CHANGELOG.md) and [the manuscript](paper/main.tex).
The earlier v0.2 exploratory model is retained for provenance under
[legacy/](legacy/README.md). It is neither packaged nor tested, and its heuristic
intervention endpoints return HTTP 410 from the current server.

Code is [MIT licensed](LICENSE). External data and standards retain their
original terms; source attribution is in the snapshot README. Cite this release
using [CITATION.cff](CITATION.cff).
