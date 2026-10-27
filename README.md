# Empty Quarter Digital Twin v0.3

EQ-BioDT connects metagenomic evidence to explicit ontological rules, conditional soil-process scenarios and testable forecasts. Its scientific contract distinguishes observations, constitutive relations, causal assumptions and decision objectives.

The release uses **551 soil metagenomes and 5,438 corrected reference MAG annotations**. An eleven-step CBB screen finds **83** MAGs with every step detected; **23** also carry a screened hydrogenase, and **20** pass completeness ≥90%, contamination <5%. These are annotation-based capacity measurements. Local activity, absolute biomass and gas kinetics need independent measurements.

[Seven-concern disposition](reviews/v03/report.md) · [Complete executable contract](data/v03/model_contract.json) · [Results](reviews/v03/results.json) · [Manuscript](https://www.overleaf.com/project/6addceec3f96e54f5b1ddd02)

## Run

~~~bash
uv sync
uv run eq-biodt serve --port 8000
~~~

Open http://127.0.0.1:8000/ for the process console or /docs for the API. The console runs the same Python model as the command line and exports inputs, assumptions and daily mass balances.

~~~bash
uv run eq-biodt contract
uv run eq-biodt simulate --input data/v03/example_scenario.json
uv run eq-biodt compare --input data/v03/example_comparison.json
uv run eq-biodt experiment --input data/v03/example_scenario.json
uv run eq-biodt backtest
~~~

Kinetic coefficients are required inputs. The examples supply 0.001 s⁻¹ for H₂ and CO; these are **assumptions**, not fitted Empty Quarter rates.

## Supported uses

- Compare repeated water pulses with a common water and salt history.
- Examine prescribed porosity, temperature, evaporation and atmospheric trace-gas changes.
- Screen pathways and hydrogenase carriers from full MAG annotations, with quality and missingness exposed.
- Plan a randomized moisture × porosity experiment with a withheld condition.
- Compare Pareto sets and uptake rankings across 24 explicit parameter alternatives.
- Freeze and score genomic-capacity forecasts against persistence and mean baselines.

The gas model couples dry-end activation, retention, diffusion, concentration, oxygen and energy. Water and salt are conserved. A separate three-reaction ceiling balances atoms and substrate/energy supplies. Incomplete pathway evidence causes abstention.

High-salinity states retain their water/salt ledger but abstain from gas prediction beyond the dilute NaCl domain. Species-specific growth, cloud-seeding efficacy, nitrogen cycling, sulfate reduction, mineral equilibria and community turnover require additional mechanisms and measurements. Carbon removal is unavailable from simulated gas uptake; a separate account requires measured stocks, carbon sources, a counterfactual and lifecycle emissions.

## Forecast registration

No next sampling date has been scheduled, so **no live prospective forecast is issued**. Once a date is known:

~~~bash
uv run eq-biodt freeze --target-date YYYY-MM-DD --output forecast.json
uv run eq-biodt score --forecast forecast.json --observations new_observations.tsv
~~~

The target must be after issuance. Existing files cannot be overwritten. Scoring verifies integrity and matches observations after issuance within ±14 days of the target. New observations require site, compartment, collection_date and total_h2_marker_ge columns. There are no calibrated predictive intervals. A supplied issue-date override in the Python test interface labels the artifact as a fixture.

## Scientific evaluation

The campaign-blocked backtest has 409 visit means and 278 held-out visits. Geography/season ridge improves RMSE over persistence in all four test campaigns; compartment means perform better in two. The corrected reference catalogue is fixed across this retrospective analysis.

Reanalysis of 77 published bottles supports dry-end activation in a directional moisture-envelope comparison. Its threshold comes from the same study. Local rate calibration and independent field intervention validation remain outstanding.

The full inventory is generated from [contract.py](src/eq_biodt/v03/contract.py). It covers causal assumptions, ontology rules, observation mappings, all default parameters, pathway rules, prediction and planning choices. The manuscript's supplement is generated from that same inventory.

## Reproduce

~~~bash
uv run python reproducibility/v03/reproduce.py --paper paper
uv run pytest --basetemp="$HOME/tmp/eq-biodt-pytest"
~~~

The paper argument is optional; the manuscript is a separate Overleaf repository. The script regenerates backtests, external bottle fits, example scenarios, the experiment schedule, figures, tables and assumptions. External data are retrieved from a pinned commit and checked by SHA-256.

Compact evidence is bundled in the wheel and source checkout. Its [manifest](data/v03/manifest.json) identifies annotation, quality, mapping, marker and metadata sources. Rebuild with:

~~~bash
uv run eq-biodt build-evidence \
  --full-ko /path/to/rep_ko_v2.tsv \
  --quality /path/to/refined_mag_manifest.tsv \
  --metadata-panel data/processed/eq_wgs_validation_panel.parquet \
  --output data/v03
~~~

The full KO source is at ibex:/ibex/scratch/projects/c2014/rob/metagenomics/04_functional/repko_clean/results/rep_ko_v2.tsv. Its [reconstruction notes](reviews/v03/repko-source-README.md) explain retained-contig corrections. No raw reads or annotation reruns are needed for the bundled analysis.

The website export contains computed Python results without browser scenario multipliers:

~~~bash
uv run python reproducibility/v03/export_website.py /path/to/website/frontend/public/twin-data/eq_biodt_v03.json
~~~

The website source update does not deploy the live service. Deployments can set VITE_TWIN_API_URL to link a live process console.

## Version boundary

The default CLI and API serve v0.3. Historical /v1/* heuristic routes return HTTP 410 with migration guidance. Old engines, state cubes and reproduction fixtures remain available as historical v0.2 artifacts, with labelled tests. Current scientific tests are in tests/test_v03.py.

The manuscript is **“Ontological-causal reasoning for metagenomic evidence and testable soil-process scenarios in a digital twin of the Rub' al-Khali”**, by Robert Hoehndorf and Magnus Rueping. This release is a research model with conditional process parameters; absolute field predictions require calibration.
