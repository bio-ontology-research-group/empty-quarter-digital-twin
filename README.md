# EQ-BioDT v0.5.0

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

## Install and open the tool

Install [uv](https://docs.astral.sh/uv/getting-started/installation/), then:

```bash
git clone https://github.com/bio-ontology-research-group/empty-quarter-digital-twin.git
cd empty-quarter-digital-twin
git checkout v0.5.0
uv sync --locked
uv run eq-biodt serve --port 8000
```

Open <http://127.0.0.1:8000> for scenarios and evidence inspection, or
<http://127.0.0.1:8000/docs> for the interactive API. The default server binds
to localhost. Python 3.13 is the tested reproduction environment; package
metadata permits Python 3.10 or later.

The [GitHub release](https://github.com/bio-ontology-research-group/empty-quarter-digital-twin/releases/tag/v0.5.0)
also provides an installable wheel and manuscript PDFs. Install the downloaded
wheel with `uv tool install path/to/eq_biodt-0.5.0-py3-none-any.whl`, then run
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
The `/v3` API and `eq_biodt.v03` Python paths are retained for compatibility.

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
uv run python reproducibility/v05/reproduce.py
uv run pytest -q
```

The installer downloads Kobayashi-MaRust v1.5.3 and verifies its pinned SHA-256.
After dependencies and this binary are installed, the analysis runs offline.
It rebuilds evidence from the source snapshot, classifies the ontology,
exports verified memberships, runs the independent pathway audit, matched
compartment comparisons, campaign backtests, bottle analysis, numerical
convergence checks, scenario ensemble, primary-provenance trace, cross-desert
transfer evaluations and independent numerical audit, and regenerates paper
tables and figures. Outputs are in
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
[changes in v0.5](CHANGELOG.md) and [the manuscript](paper/main.tex).
The old v0.2 modules and bundled historical analyses are retained for provenance;
their heuristic intervention endpoints return HTTP 410 from the current server.

Code is [MIT licensed](LICENSE). External data and standards retain their
original terms; source attribution is in the snapshot README. Cite this release
using [CITATION.cff](CITATION.cff).
