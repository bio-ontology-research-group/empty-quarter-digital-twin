# Digital Twin v0.3: disposition of the seven concerns

27 October 2026. This report separates implemented corrections from remaining empirical requirements. The preceding audit is preserved in reviews/validation-2026-10-26/. Current code is src/eq_biodt/v03/; the default API/CLI and website use it.

| Concern | Implemented correction | Remaining requirement |
|---|---|---|
| 1. Incorrect moisture response | Water-potential activation, declared retention and Millington–Quirk diffusion produce dry shutdown, an interior optimum and saturated-pore limitation. Reanalysed published concentration time series. | Local retention, activation shape and gas flux. The external comparison is not independent threshold validation. |
| 2. Disconnected physical quantities | Atmospheric concentration and diffusivity enter analytical molar flux. Shared oxygen limits oxidation. Energy derives from flux; maintenance demand is independent or unavailable. | Local kinetics, active biomass, maintenance and gas production; nonlinear kinetics outside the trace regime. |
| 3. Ontology mixed with priors | Endpoint/unit/denominator typing, disjoint named gene categories, carrier set unions, exact complements and propagated missingness. Closure never shifts causal leaves. Historical projector rejects missing coordinates, invalid metrics and solver failures. | Annotation correctness remains empirical. No ontology-induced predictive gain is claimed. |
| 4. Conflated observations | All 5,438 corrected full MAG KO sets, exact quality joins and read mappings. Eleven CBB steps, 60%-identity marker GE, quality sensitivity, mapped-MAG entropy and unmapped fractions. No 16S functional inference. | Expression, flux, biomass and viability remain distinct. Missing steps do not establish biological absence. |
| 5. Coarse metabolism and chemistry | Three element/electron-balanced reactions with donor, oxygen, carbon and energy bounds. Independent maintenance. Solver infeasibility differs from numerical failure. Conservative NaCl state; nonideal chemistry abstains. | Genome-scale metabolism, nitrogen species, sulfate reduction, ion activities, precipitation and multicomponent brines need additional modules. |
| 6. Cosmetic interventions and missing coupling | Physical input levers; repeated water and salt memory; porosity affects transport and storage. Unsupported species names and efficacy claims are rejected. | Species growth, biomass memory, microbial turnover, dose-to-porosity and cloud-seeding efficacy require evidence. |
| 7. Carbon and objectives | Net removal unavailable from gas uptake. Separate measured stock/source/counterfactual/lifecycle account. No fixed NPP or retention. Pareto objectives preserve intact soil, water, salt and supplied disturbance. | Complete carbon boundaries, permanence and ecological objectives beyond the represented process. |

## Evidence

The full CBB screen detects 83 MAGs, versus 1,078 with RuBisCO, 416 with RuBisCO plus PRK and 315 with the additional small-subunit marker. Twenty-three full-screen MAGs also carry a screened hydrogenase; twenty meet completeness ≥90%, contamination <5%. These use corrected post-refinement annotations and describe genomic capacity.

The mapping source has undefined weights for F5Dr1 and the non-soil control SPositive. F5Dr1 remains missing in soil composition. All 551 soil samples are retained, with 550 marker totals and 550 mapping compositions.

Published Reji data contain 746 concentrations from 77 bottles. A fixed-threshold activation × diffusion envelope reduces held-replicate RMSE relative to diffusion-only in all three soils. Five supersaturated bottles remain in evaluation with zero air-pore transport. Negative fitted decay rates and measured air-dried trajectories are retained. Empty Quarter kinetics remain uncalibrated.

## Next-version workflows

**Prediction.** Four future-campaign splits score 278 visits; 409 means are available. Ridge uses geography, compartment and bounded seasonal coordinates, with training-only fitting. It beats persistence RMSE in four campaigns; compartment means win in two. Reference discovery is outside the temporal holdout. No live forecast is issued without a sampling date.

**Experiment.** A randomized moisture × porosity design has 18 cells and 54 initial allocations. Every -3 MPa condition is withheld. The protocol specifies flux, water potential, retention, abundance, activity, blanks and killed controls. Pilot variance is needed for power and final replication.

**Comparison.** One, four and fourteen pulses each supply 12 mm over 14 days, alongside an intact control. Twenty-four kinetic, retention and activation settings yield sensitivity ranges and per-member Pareto sets. This grid has no uptake rank reversals. Intact soil and the single-pulse example remain non-dominated in every member. This conditional result is not a field recommendation.

**Other uses.** Implemented levers support pulse timing, porosity disturbances, temperature/evaporation scenarios, dilute salt accumulation, atmospheric trace-gas sensitivity, and conservation-versus-water-addition comparisons. Vegetation, solar shade with biocrust photosynthesis, erosion, multi-ion brines, inoculation and dust suppression require the additional mechanisms named in the paper.

## Manuscript and interfaces

The paper describes actual v0.3 equations, data processing, semantics, forecast splits, external-test scope, resource constraints and accounting. The supplement includes every executable contract clause and default, full KO alternatives and numerical conventions. Superseded portfolio, NPP and basin-removal claims are removed from the current manuscript.

Default API endpoints use /v3; heuristic /v1 routes return 410. Historical engines and fixtures remain labelled. The console runs Python scenarios, downloads inputs and balances, exposes missingness and clears results after failures. The website displays exported v0.3 trajectories and evidence without old browser offsets.

## Verification

- Nine historical tests passed in the integrated gate; all 19 final current scientific tests pass.
- Reproduction script: current results, figures, tables and full generated assumptions.
- Browser checks: live execution/download, missing observation, high-salinity abstention, stale-result clearing, keyboard contract, mobile layouts and HTTP failures.
- Website: production Vite build. Packaged wheel installs outside the checkout with bundled evidence and all four forecast splits available.
- Manuscript: pdfLaTeX/BibTeX main manuscript and supplement, with rendered figures/pages inspected.

Source identities and result files accompany this report. Tests establish implemented contracts. Local calibration and genuinely prospective field validation remain necessary.
