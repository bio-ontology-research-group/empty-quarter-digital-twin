# Digital Twin scientific and implementation audit — 26 October 2026

**Verdict: a reproducible observational predictor, ontology-constraint demonstrator and conditional scenario tool. The available evidence does not validate its management interventions, portfolio recommendations or regional carbon removal. Its inputs are partly measured, partly imputed; many response functions are prescribed.** The revised paper and website make these boundaries explicit. Model equations and numerical outputs were preserved during this audit.

## What was established

| Question | Evidence and conclusion |
|---|---|
| Does the checked implementation reproduce its selected results? | All **9 existing pytest tests** pass, including **205 numerical fixtures**. Fixtures compare computed values with constants; they do not inspect manuscript prose or independently validate biology. |
| Are the registered source files identifiable? | All **23 registered files exist** and their recomputed hashes match the bundled manifest. This establishes file identity, not correctness of measurements, exhaustive provenance or causal support. |
| Does ontology reconciliation enforce the specified hierarchy? | The restricted compiler loads **17 variables, 13 subsumptions and 5 partitions**. Existing tests establish feasibility for their evaluated states. This is an additive constraint system, with no general OWL entailment service. |
| Is predictive performance established? | The supplied evaluation reports moderate held-out skill: group 1l spatial Spearman rho 0.543; Trip 5 R² 0.299. The environmental comparator reaches 0.290 and the same within-site Pearson correlation, 0.409. Upstream pooling/imputation limits a strictly prospective interpretation. |
| Are intervention effects causally identified? | No intervention experiment or identification analysis establishes these numerical effects. The engines evaluate stipulated structural functions. Causal interpretation assumes invariant responses, defined doses, transportability and the absence of omitted interactions. |
| Is carbon removal demonstrated? | No. Annual uptake calibration, vegetation production, mineral accumulation and retention assumptions drive the totals. **Vegetation contributes 98.3046% of P2's corridor total** (8,278.868 of 8,421.650 kt CO₂-equivalent/year). This is neither measured net removal nor a permanence assessment. |
| Are the maps the Python mechanistic model? | Several distinct products coexist: fitted statistical rasters, interpolation with browser offsets, exported Python scenarios and a procedural 3D illustration. Slider responses are not uniformly ontology-constrained causal predictions. |

Evidence: [test log](pytest.log), [numerical fixture audit](../../reproducibility/outputs/paper_numbers_audit.json), [source comparisons](source-hash-comparison.json), [independent contract checks](model-contracts.json).

## Findings requiring model or validation work

1. **Upstream information enters the validation covariates.** Working-state pH, DNA and other summaries can be imputed or pooled across campaigns before fold-specific preprocessing. Sequence reference catalogues also use all campaigns. Training-only median replacement of remaining missing values cannot undo this. Freeze the reference catalogue, environmental summaries and imputation within training partitions and rerun site/campaign holdouts before claiming prospective forecasting. The revised paper discloses this limitation; the reported metrics have not been recomputed under that stricter design. See [validator](../../src/eq_biodt/validation/twin_validator.py) and [cube builder](../../src/eq_biodt/data/build_twin_state_cube.py).

2. **Ontology consistency has a narrower meaning than scientific validity.** Subsumption inequalities require nonnegative additive measures with compatible denominators; partitions require genuinely disjoint categories. GE gene categories and organism-carrier sets must be distinguished. OCS counts states satisfying all compiled algebraic constraints. It does not test predictive accuracy, causal direction, ATP feasibility or empirical truth. Scenario preprocessing already constructs most equalities before projection. The mathematical projection guarantee applies to the supplied vector and exact metric projection, whereas the fallback performs a different repair. The full-state unit-weight objective also mixes fractions, GE and power on their numerical scales. See [constraint engine](../../src/eq_biodt/models/causal_ontological_twin.py).

3. **Input and solver contracts fail additional checks.** Projection and constitutive evaluation accept missing coordinates by filling them; an indefinite negative identity weighting matrix is accepted. A mocked numerical HiGHS failure receives the biological label `MAINTENANCE_DEFICIT_DORMANCY`, although the native status is retained. Require complete inputs or explicit imputation metadata, validate the projection metric, and separate infeasibility from numerical failure. The LP is a reduced ten-flux resource-allocation model; a chemically balanced genome-scale reconstruction would require additional reactions and balances. See [LP](../../src/eq_biodt/models/metabolic_community_engine.py).

4. **The virome model uses heuristic response equations.** The 0.254 Mantel correlation is used as a slope without a fitted response model. Viral/plasmid richness floors turn a zero baseline into one/two under a null intervention. Supplied climate, pH and evaporite inputs do not affect the calculation. ARG scores and the 0.25 flag are unvalidated screening indices. A class name containing JSDM does not establish a fitted joint species distribution model. See [virome engine](../../src/eq_biodt/models/spatial_jsdm_virome.py).

5. **Portfolio evaluation does not couple the interventions.** Six engines start independently from the same baseline. Their scores combine endpoints and time horizons; a joint resource-conserving state transition is absent. Monte Carlo draws add prescribed Gaussian noise to scores, rather than propagating sampling and parameter uncertainty through mechanisms. P1 includes current-use activities; it is not a zero-intervention control. Rankings depend on the six chosen candidates and utility weights. See [optimizer](../../src/eq_biodt/interventions/pareto_optimizer.py).

6. **Regional carbon is a separate, selectively coupled calculation.** Changing P2 rainfall from 6.5 to 100 mm produces identical regional cell outputs. Other omitted inputs include wild-ungulate density, refugium integrity and hydrocarbon remediation. The regional calculator uses different mixing and maintenance calculations from the six site engines. Literature-context ranges partly overlap calibration and compare unlike estimands: mapped reads versus cells, stock versus production, and assay versus annual rates. Low/high totals combine parameter endpoints; they are sensitivity envelopes. Establishment, water delivery, operating emissions, respiration and permanence are missing from net-removal accounting. See [carbon engine](../../src/eq_biodt/models/carbon_budget_engine.py).

7. **Geographical support is not wholly observed or empirically validated.** Of 128 supported cells, **35 inherit imputed surface guild states** from their assigned primary site. The 45 km mask uses approximate distance and interpolated environmental bounds; the atlas assignment is a fixed rule. Constant cell area and mean background beyond the corridor add representativeness assumptions. This mask does not measure prediction error or establish support for intervention doses.

8. **Provenance needs stronger execution binding.** The default simulation source list hashes processed products even when raw-source rebuilding is possible; cached digests need invalidation if files change. A fallback manifest can report frozen entries without checking the corresponding files. The independent hash comparison in this audit avoids that ambiguity for the present snapshot. A content-derived identifier does not by itself establish persistent resolution or complete standards compliance.

The [contract script](../../scripts/audit_model_contracts.py) records seven failed claim/implementation checks and four passes. The rainfall and imputed-grid findings are scope mismatches, not necessarily programming errors. Existing fixes remain effective: unconsumed DOC no longer produces fallback ATP, zero CBB co-carriage preserves zero microbial fixation, and pathway totals add correctly. Scientific response coefficients were not adjusted merely to make tests pass.

## Paper revisions

The [main manuscript](https://www.overleaf.com/project/6addceec3f96e54f5b1ddd02) now explicitly defines **ontological-causal reasoning**:

- Ontological operations translate asserted class subsumptions and exhaustive partitions into additive constraints and aggregates. They express constitutive relationships.
- Causal operations change designated inputs in hypothesized structural response functions. Direction and effect size come from those functions and their assumptions.
- Statistical fits estimate associations and predictive performance under the stated validation design.

The title, abstract, Methods and Discussion were revised accordingly. Carbon Methods now describe the implemented uptake references, conversion units, co-carriage allocation, biomass scaling, NPP, mineral gates, retention assumptions and geographical transfer. Existing numerical tables were retained. Calibration comparisons are distinguished from independent validation. The erroneous claim that the rainfall lag is established by a universal 1–4-day resuscitation window was removed.

[Supplementary Note 4](https://www.overleaf.com/project/6addceec3f96e54f5b1ddd02) documents measurement and preprocessing assumptions, each intervention engine, ontology scope and numerical weighting, reduced metabolism, virome response floors, portfolio uncertainty, carbon transfer and browser interpretation. Table S7 distinguishes fitted, transferred and assumed coefficients. An [executable inventory](../../scripts/inventory_model_assumptions.py) records **1,500 source locations across 28 Python files**, with hashes, in [JSON](../../reproducibility/outputs/model_coefficient_inventory.json). It includes structural/index constants and deliberately does not infer empirical support from their presence. Data values, ontology axioms and browser code remain separate audit surfaces. This is a substantially fuller disclosure, not a claim that every biological assumption can be proved or that every literature link was exhaustively adjudicated.

The DOI audit checked 43 DOI-bearing entries from the original 47-entry bibliography. Two did not resolve and were corrected against primary records:

- Imminger et al. (2024): title, authors, article number and DOI corrected to [*Survival and rapid resuscitation permit limited productivity in desert microbial communities*](https://www.nature.com/articles/s41467-024-46920-6). Rapid activation can occur within minutes; that observation does not estimate the model's 2.5-day diversity kernel.
- Oldoni et al. (2020): wrong F1000Research venue and DOI replaced with the matching [bioRxiv preprint](https://doi.org/10.1101/2020.03.23.983601).
- Added [Sand et al. (2023)](https://www.nature.com/articles/s43247-023-00857-8) for hydrogen GWP 11.6 ± 2.8 (one standard deviation). The model's [6,16] interval is explicitly a scenario range. [Warwick et al.](https://acp.copernicus.org/articles/23/13451/2023/) report a different estimate and should not be the sole attribution for 11.6.

All three corrected/new DOI records resolve. DOI identity checks do not establish support for each sentence. See [original audit](doi-audit.log) and [corrected records](doi-corrections.log).

## Website changes

The public website source is in the neighbouring `empty-quarter/website/frontend` repository. Its existing unrelated changes were preserved. Changes are local and have not been deployed.

- Added a persistent evidence guide showing observed versus imputed states, ontology consistency and conditional mechanisms, with selected-map explanations and keyboard-accessible detail.
- Corrected the interpolated group 1l map to display the log₂(GE + 0.01) quantity stated by its legend; cells outside the atlas envelope no longer inherit a selected portfolio's colour.
- Replaced stale hard-coded carbon card results with values from the exported scenario bundle. Clarified assumed retention, current-use P1, regional background and the separate browser equations.
- Corrected the carbon illustration's unsupported permanence claims, labelled retention percentages as assumptions, and made its small-screen text horizontally scrollable. Replaced the API console's unavailable Carto basemap with OpenStreetMap tiles.
- Corrected metric explanations including CLR interpretation, bootstrap count, spatial blocking, OCS, objective definitions, units and Pareto dominance.
- Improved the local API console's responsive layout, provenance disclosure and error handling. Missing results display `Unavailable`; genuine OCS zero remains zero. Failed requests clear stale outputs, and opening the page does not silently run an intervention.

See [desktop screenshot](website-desktop.png), [carbon view](website-carbon.png), [mobile screenshot](website-mobile.png), [API console](console-desktop.png), and [browser checks](browser-checks.json).

## Scope and next acceptance criteria

This audit inspects model code, processed data, registered-file identity, published-source identity, reported validation and the rendered interfaces. It does not independently repeat sequencing QC, genome assembly/annotation or laboratory assays, and it supplies no new field experiments. The predictive, causal and deployment limitations above remain substantive.

Before making empirical management or carbon-removal claims: correct solver/input contracts; rerun validation with training-only preprocessing; estimate response functions from reserved intervention experiments; propagate parameter and structural uncertainty through coupled scenarios; reconcile regional/site calculations; and measure gas fluxes, production, mineral sources, lifecycle costs and retention. The current system can support transparent hypothesis exploration while that evidence is collected.

Validation of the revisions: production frontend build passed; main and supplementary LaTeX builds passed; ten browser checks passed, including a real API intervention, keyboard operation, mobile overflow and HTTP error handling. Rendered PDF and browser screenshots were inspected. The [reviewable diff](reviewed-changes.diff) includes modifications to pre-existing files; new audit scripts, the evidence-guide component and assumptions note are separate additions.
