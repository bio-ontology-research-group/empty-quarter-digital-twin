1. Candid editorial verdict on the six upgrades

I would still recommend rejection at Nature Sustainability or Nature Ecology & Evolution in the manuscript’s present form. The revision repairs several implementation errors and describes its assumptions more honestly, but it does not yet establish a major ecological discovery or an empirically supported management decision.

The six upgrades are not equally valuable. Upgrades 4 and 6 are essential foundations; selected parts of Upgrade 1 could produce substantive science; Upgrade 3 could establish a useful computational contribution. Upgrade 2 needs targeted model curation before integration. Upgrade 5 should come last, after demonstrating that the intervention responses being mapped are defensible.

The audit also changes several premises of the proposed plan:

- The 2,480 figure denotes core specimen records, not 2,480 sequenced 16S profiles.
- The five molecular sampling campaigns begin in March 2023 in the inspected records, not 2022.
- The environmental panel has campaign-level date substitutions, including an August 2024 date inconsistent with the curated field dates.
- The current code confuses numeric auxiliary sample-site codes with ontology identifier suffixes.
- A substantial existing GSMM pilot explicitly concludes that the available reconstructions do not support the proposed hydrogen-metabolism inference.
- I did not find the claimed metaproteomic expression dataset.
- The inspected ontology modules do not contain the five named disjoint-union axioms described in the Supplement.

These are scientific prerequisites, not presentation details.

Nature Ecology & Evolution is the better initial target for an analysis using the existing evidence. Its scope includes microbial ecology, community ecology and biogeochemistry. A compelling paper would establish a general ecological result about environmental filtering, community replacement or temporal persistence, with the computational framework supporting that result. [NEE aims and scope](https://www.nature.com/natecolevol/aims).

Nature Sustainability becomes preferable only if the study establishes a consequential management trade-off using defensible outcomes, realistic intervention definitions and explicit resource constraints. Adding a Saudi Green Initiative framing to a model of gene abundance is insufficient. Its scope emphasizes sustainability problems, policy dimensions and possible solutions. [Nature Sustainability aims and scope](https://www.nature.com/natsustain/aims).

My preferred scientific question is:

“Do the environmental relationships observed across desert habitats predict changes within habitats through time, and which biological quantities and management conclusions remain defensible when measurement, causal and mechanistic uncertainty are propagated?”

That is falsifiable. It could reveal that spatial gradients are poor substitutes for temporal responses, which would itself matter for dryland restoration. It does not presume the answer.

Upgrade 1: empirical causal identification

The direction is right; the proposed terminology is premature. Replacing prescribed coefficients with fitted coefficients produces an empirical model. It produces causal identification only when the design and assumptions identify a specified intervention effect. DML addresses nuisance-estimation bias under an identification framework; it does not supply that framework. [Chernozhukov et al.](https://academic.oup.com/ectj/article/21/1/C1/5056401).

(a) Rainfall

Within-site analysis is worthwhile, but five unevenly sampled campaigns are not a densely observed longitudinal experiment. There are only seven WGS sites represented in all five trips. Replicate specimens do not create independent rainfall exposures.

The existing lifeonair analysis already contains within-site rainfall models. Its 90-day ERA5 and POWER anomalies agree weakly, with within-site Spearman correlation approximately 0.193. Some apparent ERA5 responses weaken or change direction with POWER. A new estimator must confront this exposure uncertainty rather than select the product producing the preferred effect.

A credible analysis should:

- Reconstruct exposure windows relative to actual site-specific collection dates.
- Separate short pulse responses from seasonal and annual associations.
- Quantify the independent storms and weather-grid cells supporting each contrast.
- Use site effects and parsimonious campaign adjustment, then inspect how much exposure variation remains.
- Treat SMAP moisture as a possible mediator of rainfall, rather than automatically adjusting for it when estimating rainfall’s total effect.
- Test alternative weather products, sampling-route effects and prespecified timing falsifications.

P7, P30, P90 and P365 overlap substantially. Fitting all four as unrestricted “distributed lags” invites unstable coefficients. Use a small, prespecified set of non-overlapping windows or a constrained lag basis, conditional on adequate exposure support.

The 1–4-day resuscitation claim cannot be calibrated from 30-day SMAP and 90-day precipitation averages. Daily weather exists, but daily microbial observations do not.

(b) Revegetation

The 106 root-versus-bulk visits at 47 sites are valuable paired observations. They estimate a habitat contrast. They do not directly estimate the effect of planting shrubs on previously bare soil.

Root presence, plant identity, establishment history, substrate supply, shade and microsite selection are bundled together. NDVI/EVI are remotely sensed vegetation indices, not calibrated shrub-cover treatments. A 60-day pixel average cannot be substituted directly for the simulator’s “30% Haloxylon cover.”

Start with a paired hierarchical model and explicit estimands. Use DML as a comparison if the effective sample size and overlap support it. Cross-fit by site. For continuous exposure, define a dose-response or supported exposure-shift estimand; a binary propensity-score workflow is insufficient.

Cinelli–Hazlett analysis would make omitted-confounding assumptions visible, but its sensitivity region is not an unconditional causal bound. Use an implementation appropriate to the fitted estimator and preserve clustered uncertainty. [Cinelli and Hazlett](https://escholarship.org/uc/item/0pg6471b).

Most importantly, DNA extract concentration is not calibrated microbial biomass. Plant DNA, extraction efficiency and extracellular DNA undermine the current conversion from relative depletion to a 1.38-fold increase in absolute scavenger density.

(c) Salinity and ionic composition

A geochemical response surface is feasible as an observational analysis. Calling it the causal effect of Ca-SO4 versus Na-Cl is not justified.

XRF measures bulk elemental composition. It does not directly measure dissolved ion activities, sulfate speciation, electrical conductivity or porewater water activity. Ratios such as Ca/S and Na/Cl also discard total concentration. A soil containing little salt and a hypersaline soil can share a ratio.

The aquatic samples provide actual dissolved-ion measurements, but from three water samples at two sites. The associated hundreds of genomes do not supply hundreds of independently replicated salinity treatments.

A common-baseline factorial simulation would repair the existing confounded scenario comparison. It would isolate factors within the model, not experimentally identify their effects in nature.

(d) Thermal refugium

The paired depth comparison is useful, but depth changes moisture, UV exposure, substrate availability and DNA preservation as well as temperature. Removing a horizon is also a different intervention from observing an existing surface community.

The existing temperature-support audit is particularly important: campaign explains approximately 95.8% of temperature variation across the 68 paired visits. Residual temperature variation after campaign adjustment is much smaller. A temperature interaction could therefore be largely a campaign comparison.

Test the depth contrast and its modification by measured exposure, but do not claim that this identifies an 11°C cooling effect or a causal benefit of preserving a specified horizon fraction.

(e) Grazing and oil spills

This component currently lacks the proposed empirical foundation.

I found named-site environmental records and notes about additional samples, but no linked Camel1, Camel2 or named oil-spill molecular specimens in the canonical specimen/WGS tables. The apparent WGS samples from “site 64” are campground samples, not Camel1.

Vanadium and nickel are not unique indicators of petroleum exposure. Neither provides stocking density, compaction history, manure loading or hydrocarbon concentration. ARG annotations measure genetic potential, not clinical risk or transfer rates.

Remove quantitative grazing, seismic and oil-remediation recommendations from the main paper unless archived exposure and assay records can be linked reliably. These modules can remain explicitly illustrative.

Upgrade 2: real GSMMs

Do not implement this as a wholesale replacement of the small LP with thousands of uncurated models.

The existing [metabolic pilot report]( /home/leechuck/Public/software/empty-quarter-lifeonair/lifeonair/results/metabolic_pilot/REPORT.txt) examined 72 MAGs selected from 5,438 reconstructions. Among 48 hydrogenase-marker carriers, only one had a native internal hydrogen reaction whose gene–protein–reaction rule included the exactly matched marker. No group-1l or group-1h model gained ATP from hydrogen in the tested primary reconstructions. Alternative reconstructions disagreed on hydrogen phenotypes.

This is evidence that model availability and biological readiness differ substantially.

The defensible upgrade is a small, curated panel of representative lineages with:

- Verified hydrogenase identity, electron-transfer partners and reaction coupling.
- Explicit oxygen, carbon, nitrogen and sulfur balances.
- Separated native and gap-filled reactions.
- Closed-exchange, ATP-cycle, mass/charge and positive-control checks.
- Maintenance requirements that remain requirements when infeasible.
- Flux-variability and reconstruction sensitivity analyses.

Standardized model testing is necessary, but cannot establish physiological validity by itself. [MEMOTE](https://www.nature.com/articles/s41587-020-0446-y).

Atmospheric gas concentration must be connected to uptake through a physical transport and kinetic model. An arbitrary exchange bound is not that connection. Genome abundance must also be linked to biomass before model fluxes can become ecosystem fluxes.

If these conditions cannot be met using existing evidence, retain a smaller, balanced feasibility model and restrict its conclusions. A large GSMM ensemble would otherwise amplify unsupported assumptions.

Upgrade 3: ontology compiler and weighted reconciliation

This is a worthwhile software contribution, but three claims need separating:

- Compiling valid quantitative constraints from an ontology.
- Improving predictions through reconciliation.
- Identifying causal effects.

Success at the first does not establish the other two.

The compiler should support a declared, measurement-grounded fragment of OWL, not claim to translate arbitrary OWL 2 DL into a polyhedron. Subclass and disjoint-union semantics concern sets. Quantitative inequalities additionally require compatible measures, denominators, contexts and units. Disjointness alone does not establish an exhaustive partition. [OWL 2 Direct Semantics](https://www.w3.org/TR/owl2-direct-semantics/).

For weighted reconciliation, use

\[
\hat z=\arg\min_{z\in C}(z-\hat y)^\top \widehat W^{-1}(z-\hat y),
\]

with shrinkage covariance estimated from inner, training-only out-of-fold residuals. Compare identity, diagonal and shrinkage covariance estimates. The corresponding projection guarantee concerns the weighted loss and assumes the true target lies in the constraint set; it does not promise improvement in every taxon, ordinary MSE or downstream ecological decision.

MinT and nonnegative reconciliation are established methods. The novel contribution would be correct semantic compilation, provenance and reuse, not the weighted projection itself. [MinT](https://1www.robjhyndman.com/publications/mint/), [nonnegative reconciliation](https://1www.robjhyndman.com/publications/nnmint/index.html).

An indispensable comparator is a flat implementation given exactly the same constraints and inputs. Numerical equivalence would demonstrate correct compilation. It would not demonstrate a special predictive advantage of ontological representation.

Upgrade 4: empirical state cube

This is the highest-priority upgrade.

However, “100% empirical” should mean that observed entries have traceable empirical provenance. It cannot mean filling every site–trip–compartment–assay combination with measurements that do not exist.

Maintain separate observation and prediction layers. Missing assays remain missing. Record sample identity, collection time, assay time, compartment, measurement units, catalogue version, detection limits and provenance.

Keep 16S read proportions, MAG mapped-read abundance, inferred carrier abundance and gene copies per genome equivalent distinct. They have different observation processes. Retain unmapped and unclassified fractions rather than forcing the recovered catalogue to represent the entire community.

Upgrade 5: spatial atlas and optimization

A 248-cell map is a display resolution, not 248 independent observations. The current grid is interpolated climate over a corridor. It is not a basin-wide empirical sample of the Rub’ al-Khali.

The immediate deliverable should be a corridor prediction atlas with explicit support and uncertainty masks. Cells lacking suitable geochemistry, habitat representation or treatment overlap should carry an abstention or extrapolation label.

Continuous optimization is premature while the engines independently modify copies of the same baseline. The current objectives average diversity from incompatible scenarios and sum fluxes that can count baseline production repeatedly.

Before optimization:

- Define shared state variables and a common spatial and temporal accounting basis.
- Couple water, salt, porosity, vegetation and gas transport.
- Include an actual no-action option.
- Use meaningful resource constraints and outcomes.
- Separate natural saline-community conservation from anthropogenic salinization.
- Propagate uncertainty through fitting, reconciliation, simulation and optimization.

Bootstrap intervals are confidence intervals; Bayesian posterior intervals are credible intervals. Monte Carlo sampling of assumed coefficients is neither automatically. Report the interpretation explicitly.

For management, probabilities of dominance, regret and rank stability are more informative than a four-decimal utility score.

Upgrade 6: forward-in-time validation

Essential, but rename it “temporally ordered forecasting validation.” Time ordering does not establish causality.

The inspected WGS metadata contain 91 Trip-5 samples at 52 sites, representing 89 distinct site–compartment combinations. Eighty-five of those combinations have earlier WGS observations before further assay/QC exclusions. This is a useful test set.

Separate:

- Forecasting before Trip 5, using information available at a fixed forecast origin.
- Predicting microbial state after observing environmental conditions immediately before collection.
- Retrospective prediction using contemporaneous laboratory assays.

These are different tasks. Realized Trip-5 rainfall cannot be used in a forecast issued before that rainfall unless treated as a conditional scenario.

All preprocessing, climate summaries, covariance estimation, feature selection and model selection must respect the cutoff. Test whether MAG/reference catalogues themselves incorporate Trip-5 discoveries. Report a frozen-reference forecast separately from a retrospective analysis using the complete catalogue.

Trip 5 has already been inspected in this project and companion analyses. It can support a transparently specified retrospective temporal holdout, but cannot honestly be described as a previously untouched confirmatory test.

2. Audit of what actually exists on disk

I inspected the manuscript sources, previous report, relevant implementation, local tables, manifests, ontology files, model summaries and selected SBML locations. I independently counted key tables and reconstructed selected joins. I did not rerun the complete validation suite or verify the remote IBEX payloads.

Several directory names in the proposal do not exist at the specified `~/Documents/papers/` locations: `empty-quarter-map`, `empty-quarter-viruses`, `empty-quarter-saline`, `empty-quarter-models`, `empty-quarter-models-aquatic`, `empty-quarter-models-thermophiles` and `empty-quarter-rak-pr`. Relevant material exists under different repositories, described below. Code comments naming the missing directories are not provenance.

A. Sampling universe, dates and identity

The central files are:

- [sample_ledger.tsv](/home/leechuck/Public/software/empty-quarter-ecology-data/evidence/release/sample_ledger.tsv)
- [release_evidence.json](/home/leechuck/Public/software/empty-quarter-ecology-data/evidence/release/release_evidence.json)
- [site_aliases.tsv](/home/leechuck/Public/software/empty-quarter-ecology-data/metadata/samples/site_aliases.tsv)
- [environmental_measurements_curated.tsv](/home/leechuck/Public/software/empty-quarter-ecology-data/metadata/environmental/environmental_measurements_curated.tsv)

The ledger has 2,550 records: 2,480 core-site records, 36 Trip-1-only auxiliary records and 34 laboratory-control records. The 2,480 include plant-material records.

The release evidence distinguishes:

- 1,271 feature-table profiles.
- 1,242 unique normalized field identifiers in that feature table.
- 1,237 ecology-analysis profiles.
- 1,227 primary-site ecology profiles.

It also records repeated normalized identifiers whose laboratory explanation remains undocumented. These profiles must not automatically become independent biological replicates.

The ledger has no collection-date field. Dates must be joined from curated field records. The primary-site campaigns are:

- Trip 1: 17–20 March 2023.
- Trip 2: 14 July 2023, eight sites.
- Trip 3: 17–21 February 2024.
- Trip 4: 29–31 August 2024.
- Trip 5: 12–15 October 2025.

The existing 2022–2025 export dates should therefore be corrected for these observations.

B. Rainfall, vegetation and temperature

The exact [env_site_trip.tsv](/home/leechuck/Public/software/empty-quarter-lifeonair/lifeonair/results/inputs/env_site_trip.tsv) contains 244 rows: 60, 8, 60, 60 and 56 for Trips 1–5.

Its exposure columns are:

`ndvi_W60d`, `evi_W60d`, `smap_sm_W30d`, `P90_era5`, `P365_era5`, `P365_bestmatch`, `P90_power`, `P365_power`.

There is no P7 or P30 rainfall column. `P365_bestmatch` is missing in 56 rows.

More seriously, `trip_date` contains one date per campaign. Trip 4 is assigned 20 August 2024, whereas the curated records place sampling on 29–31 August. Rebuild every antecedent exposure affected by these substitutions.

Useful daily source files already exist:

- [era5_daily_1995_2025.tsv](/home/leechuck/Public/software/empty-quarter-lifeonair/lifeonair/results/inputs/era5_daily_1995_2025.tsv): 675,720 rows, 60 sites, through 31 October 2025.
- [power_daily_2021_2025.tsv](/home/leechuck/Public/software/empty-quarter-lifeonair/lifeonair/results/inputs/power_daily_2021_2025.tsv): 105,900 rows, 60 sites.
- [daily_weather_canonical.tsv](/home/leechuck/Public/software/empty-quarter-ecology-data/metadata/climate/daily_weather_canonical.tsv).
- [rain_event_product_exposures.tsv](/home/leechuck/Public/software/empty-quarter-ecology-data/metadata/climate/rain_event_product_exposures.tsv): 300 rows covering five products at 60 sites for the Trip-1, preceding-days-1–4 comparison.

The last file provides a useful exposure-product sensitivity test, not a five-trip daily microbial panel.

Temperature support exists in `lifeonair/results/inputs/per_trip_site_temperature.tsv` and the diagnostic [A3_temperature_residual_support.tsv](/home/leechuck/Public/software/empty-quarter-lifeonair/lifeonair/results/h123_v2/A3_temperature_residual_support.tsv).

C. Paired observations and shotgun abundances

[WGS sample_metadata.tsv](/home/leechuck/Public/software/empty-quarter-workbench/data/metadata/sample_metadata.tsv) contains 551 rows: 117 Surface, 166 Deep and 268 PlantRoot.

The paired contrast files are already available:

- [h1_visit_log2fc_root_vs_bulk.tsv](/home/leechuck/Public/software/empty-quarter-lifeonair/lifeonair/results/h123_v2/h1_visit_log2fc_root_vs_bulk.tsv): 106 visits, 47 sites.
- [h1_visit_log2fc_surface_vs_bulk.tsv](/home/leechuck/Public/software/empty-quarter-lifeonair/lifeonair/results/h123_v2/h1_visit_log2fc_surface_vs_bulk.tsv): 68 visits, 36 sites.

These are paired visit summaries, not 106 or 68 independent sites.

[guild_genome_equivalents_v2_r3noflagged.tsv](/home/leechuck/Public/software/empty-quarter-lifeonair/lifeonair/results/reads/guild_genome_equivalents_v2_r3noflagged.tsv) has 1,719 rows covering 573 libraries at three identity thresholds. Select the declared soil subset and harmonize the threshold used for estimation, validation and simulation.

A full local abundance source exists at [coverm_rep_abundance_long.tsv.gz](/home/leechuck/Documents/papers/empty-quarter-mags/lifeonair/results/inputs/coverm_rep_abundance_long.tsv.gz). It is headerless and contains 487,639 rows, 573 libraries, 5,438 MAG identifiers and an unmapped category. Unmapped fractions range from approximately 7.3% to 100%. Exclude controls appropriately and preserve that substantial unrepresented fraction.

The revised manuscript still incorrectly says that all compartments at Sites 59 and 60 lack metagenomes. The source tables contain:

- Site 59: `59PRr1`.
- Site 60: `60PRr1`, `S60PRr1`, `V60PRr2`.

These populate root-adjacent states. Their surface and subsurface states remain unobserved in this WGS panel.

D. Auxiliary sites: a consequential identity error

The [alias table](/home/leechuck/Public/software/empty-quarter-ecology-data/metadata/samples/site_aliases.tsv) establishes that Trip-1 numeric sites 61–64 map to ontology sites 67–70:

- 61 → water well, location 2.
- 62 → road, location 1.
- 63 → road, location 2.
- 64 → campground.

The builder instead derives auxiliary numeric keys from IRI suffixes, and the [validation panel joins on those numbers](/home/leechuck/Public/software/empty-quarter-digital-twin/src/eq_biodt/validation/twin_validator.py:177). Consequently, `62PRr2` is assigned the hot-spring site and the three `64*` samples are assigned Camel1.

This directly undermines the proposed grazing calibration and alters spatial covariates for those four WGS samples.

In the canonical tables inspected:

- Camel1 has a named environmental observation, but zero linked specimen/WGS records.
- Camel2 has no canonical named-site entry or linked molecular records.
- The named oil-spill site has environmental records from Trips 3 and 5, but zero linked canonical specimen/WGS records.
- Trip-5 notes mention additional oil and camel samples, including A1–A4. Notes establish collection evidence, not completed assays or usable quantitative endpoints.

Recover any additional sample manifests before claiming a small calibration cohort exists.

E. 16S, XRF and pH

The actual 16S sources include:

- [feature-table-trips1-5.tsv](/home/leechuck/Public/software/empty-quarter-ecology-data/metadata/taxonomy/feature-table-trips1-5.tsv): 351,472 feature rows and 1,271 profile columns.
- [taxonomy-trips1-5.tsv](/home/leechuck/Public/software/empty-quarter-ecology-data/metadata/taxonomy/taxonomy-trips1-5.tsv).
- `feature_table_filtered_raw_w_reps.tsv`: 63,859 feature rows, 975 sample columns.
- `feature_table_filtered_normalized.tsv`: 63,859 feature rows, 531 columns.

These are distinct processing products. Select the release and filtering policy explicitly; do not combine their denominators.

The current twin consumes [xrf_site_compartment_panel.tsv](/home/leechuck/Public/software/empty-quarter-lifeonair/lifeonair/results/inputs/xrf_site_compartment_panel.tsv), a 180-row panel with geochemistry and Shannon diversity. It is not a specimen-by-trip environmental panel. The lab-filtered source has 178 rows and elemental/oxide columns; campaign-specific measurement files also exist.

[ph_accepted_measurements.tsv](/home/leechuck/Public/software/empty-quarter-ecology-data/metadata/ph/ph_accepted_measurements.tsv) contains the stated 712 admitted measurements, with specimen IDs, trip and assay dates. Some assays were performed in 2026 on archived material. They may support retrospective inference, but cannot silently become information available to a forecast issued before October 2025.

F. GSMMs and flux outputs

Three useful bodies of material exist:

- Soil: [lifeonair metabolic_pilot](/home/leechuck/Public/software/empty-quarter-lifeonair/lifeonair/results/metabolic_pilot). It contains the 5,438-model inventory, 72-MAG panel, 10,584 primary scenario results, 216 model-variant audits and a 25-genome alternative-reconstruction comparison. XML payloads are documented on IBEX; the local compact evidence is substantive and explicitly negative about readiness.
- Aquatic: [empty-quarter-tiannyu](/home/leechuck/Documents/papers/empty-quarter-tiannyu). Its provenance records 211 genomes and 210 model files on IBEX. The [flux supplement](/home/leechuck/Documents/papers/empty-quarter-tiannyu/sections/13_supplement_flux.tex) still requests archived media, bounds, diagnostic reactions, gene support and loopless outputs. It reports model-feasibility classifications, not measured environmental flux.
- Thermal: [metabolic-pilot-20261012](/home/leechuck/Documents/papers/empty-quarter-thermal/results/metabolic-pilot-20261012). Actual local SBML files exist, alongside `model_summary.tsv`, `growth_scenarios.tsv`, `maintenance_scenarios.tsv`, `mass_balance_flags.tsv`, reaction audits and curated-model comparisons. These are isolate-model experiments, not a calibrated desert-community temperature engine.

There is also a [16-community cross-feeding pilot](/home/leechuck/Public/software/empty-quarter-proteome/results/smetana/CROSSFEED.md). It reports no support for greater donor potential among atmospheric-guild MAGs under its tested conditions. Preserve this negative evidence.

G. Aquatic chemistry, mobilome, resistome and claimed proteomics

The [aquatic supplementary workbook](</home/leechuck/Documents/papers/empty-quarter-tiannyu/source/tiannyu-2026-10-18/Empty quarter project supplementary table-latest-9.xlsx>) contains actual chemistry:

- Table S1 lists eight environmental sample identities.
- Table S2 contains three surface-water samples, with conductivities 216.2, 231.4 and 10.65 mS/cm, and measured Ca, Na, Cl and sulfate.
- Table S3 contains solid-sample chemistry.

These are better anchors for brine chemistry than XRF ratios alone, but their regional replication is limited. Existing 70 versus 150 g/L enrichment material also appears in the aquatic analysis; audit its replication and treatment history before using it as an independent contrast.

The twin reads [sample_mobilome_load.tsv](/home/leechuck/Documents/papers/empty-quarter-mobilome/results/sample_mobilome_load.tsv), with 552 rows. The current companion manuscript reports 19,300 corrected vOTUs and 322,755 pOTUs. Neither the proposed 2,340/2,213 counts nor the twin’s older catalogue totals should be carried forward without a version crosswalk.

A newer [temporal/environmental analysis](/home/leechuck/Documents/papers/empty-quarter-mobilome/plans/2026-10-14-temporal-environment-analysis/README.txt) already provides frozen environmental inputs, matched-cohort flags and temporal analyses. Reuse its provenance work.

[amr_long.tsv](/home/leechuck/Public/software/empty-quarter-proteome/results/amr_long.tsv) contains 12,331 annotation rows. These are MAG-associated gene annotations, not calibrated ARG exposure or risk measurements.

I found predicted-protein catalogues, structural models and annotation products, but no identifiable sample-linked mass-spectrometry matrix supporting “53.4% expressed dark proteome.” The [current proteome manuscript](/home/leechuck/Documents/papers/empty-quarter-proteome/main.tex:113) explicitly describes coding potential and says proteomics would establish translation. The twin’s expression fraction must therefore be removed or supported by a separately located assay dataset.

H. Ontology and current implementation

[rubalkhali.owl](/home/leechuck/Public/software/empty-quarter-ecology-data/ontology/rubalkhali.owl) is a real ontology with imported measurement, sample, site and taxonomy modules. However, the inspected local modules contain neither the named hydrogenase/DNA partitions nor explicit disjoint-union axioms described as A1–A5.

The engine constructs its 17 variables and constraints in Python. Supplementary Note 1 should not state that those exact axioms reside in the released ontology without supplying the missing versioned module.

I also checked the current consistency routine directly: an empty state dictionary receives OCS = 1 because missing coordinates default to zero. Explicit NaNs are rejected, but absent measurements are still treated as zeros. This must be fixed before extending the compiler.

The [CKG research guidance](/home/leechuck/Documents/papers/causal-knowledge-graphs/CKG_USEFULNESS.md) already states the correct standard: distinguish information, representation, inference and prediction benefits; compare against equally informed flat implementations; never treat simulated outputs as independent biological validation.

3. Additional high-impact analyses

A. Test spatial-to-temporal transport explicitly

This is the most promising addition.

Estimate environmental relationships between sites and within sites separately. Test whether the between-site relationship predicts the direction and magnitude of within-site change. Compare models containing site history, environmental anomalies and both.

A reliable failure of transport would challenge the use of spatial desert gradients to predict responses to changing rainfall or vegetation. That is potentially more consequential than a slightly better overall correlation.

B. Separate ecological replacement from changes within lineages

Decompose changes in community functional capacity into lineage turnover, carrier prevalence and gene-copy contribution. Ask whether the same decomposition explains temporal change and root/depth contrasts.

Some of this already exists in lifeonair, including `A2_*` turnover analyses and `A6_copy_number_vs_carrier_representation.tsv`. The novel contribution must be prediction or transport across contexts, not republication of the same decomposition.

Also separate hydrogen oxidation potential from carbon fixation potential. The existing CBB table records CBB co-carriage in 101 of 558 group-1l-only representatives. Group-1l abundance cannot simply stand for chemosynthetic primary production.

C. Quantify measurement sensitivity before causal sensitivity

Construct an observation model incorporating:

- Read depth and mapping threshold.
- Unmapped MAG abundance.
- Catalogue and taxonomy versions.
- Extraction and sequencing batches.
- Repeated profile identifiers.
- Detection limits and contamination flags.
- Relative versus absolute abundance uncertainty.

A sophisticated causal estimator cannot recover an effect if differential measurement error creates the apparent contrast.

D. Use existing PMA observations narrowly

[PMA_ASV_table.tsv](/home/leechuck/Public/software/empty-quarter-ecology-data/metadata/relic-dna/PMA_ASV_table.tsv) contains nine treated/untreated pairs plus controls. These can assess whether particular diversity conclusions are sensitive to PMA treatment in those aliquots.

They cannot identify a basin-wide viable fraction, a two-pool DNA turnover rate or a surface–subsurface viability difference. The companion paper already makes these limitations clear.

E. Test the ontology’s actual contribution

A stronger evaluation would demonstrate concrete errors prevented by typed identities and measurement semantics:

- Numeric site code versus ontology IRI suffix.
- Collection date versus assay date.
- Gene copies versus organisms.
- Bulk elemental concentration versus dissolved-ion activity.
- Carbon-fixation potential versus realized flux.

Measure compilation correctness, unsupported mappings rejected, provenance completeness and maintenance effort. Compare against a competent flat implementation. This is more convincing than presenting zero arithmetic violations as ecological validation.

F. Analyze decision robustness and the value of existing information

Before continuous optimization, ask which decisions remain unchanged across defensible causal assumptions, weather products, reconstruction variants and outcome weights.

Where no decision is stable, identify which already available assay or metadata recovery would reduce uncertainty most. Abstention can be scientifically useful if it is tied to a concrete limitation.

G. Seek external validation from existing datasets

Use archived independent desert datasets where compatible measurements exist. Harmonize outcomes and sampling depths before testing transfer. The local Atacama comparator files provide a starting point, but check whether they have already informed model choices.

External replication would contribute more to broad ecological relevance than increasing the number of interpolated cells.

4. Execution plan and manuscript restructuring

Step 1: freeze the scientific question and publication boundaries

Create a claim-overlap matrix against the amplicon/ecology, lifeonair, mobilome, aquatic, proteome and thermal papers.

Several proposed “new” analyses already exist: rainfall-pulse associations, within-site precipitation models, paired root effects, temperature interactions, lineage turnover and mobilome temporal structure.

Decide whether EQ-BioDT will establish a distinct synthesis result or whether the strongest route is a consolidated flagship paper. Do not build a new manuscript whose main evidence is a collection of companion-paper headlines.

Step 2: rebuild identities, chronology and denominators

Produce one versioned sample–visit–assay crosswalk. Resolve auxiliary aliases, repeated profile IDs, real collection dates and assay dates. Recompute rainfall and vegetation windows.

Required outputs are an assay-coverage table, site–trip sampling matrix, catalogue-version manifest and explicit exclusion ledger.

Stop causal analysis if treatment timing or sample identity cannot be resolved.

Step 3: build the observation cube before the forecasting cube

Replace phylum formulas with observed taxonomic profiles. Integrate WGS guild estimates and CoverM abundances while preserving measurement distinctions and missingness.

Use site × visit × compartment as the principal environmental unit, with specimens and assay replicates nested underneath. Do not pool campaigns before temporal validation.

Retire inferred occurrences and formula-derived evenness from observational exports.

Step 4: specify estimands and assess identification

Choose one primary contrast and a small number of secondary contrasts. For each, record:

- Intervention or exposure.
- Population and support.
- Outcome and time horizon.
- Confounders, mediators and selection mechanisms.
- Identification assumptions.
- Measurement limitations.
- Criteria for refusing causal interpretation.

Start with paired root/depth contrasts and rainfall anomalies. Defer grazing, seismic and inoculation effects unless stronger archived evidence emerges.

Step 5: establish forecasting performance

Use rolling-origin development within Trips 1–4, then evaluate the prespecified retrospective Trip-5 holdout.

Include persistence, site/compartment history, environmental ridge, a flexible environmental model and the proposed method on identical inputs.

Score within-site changes explicitly. Report error, calibration, change-direction accuracy and performance by compartment. Assess uncertainty using the actual dependence structure, not specimen counts alone.

If the model cannot beat a useful persistence baseline for temporal change, do not present a dynamic intervention atlas as validated.

Step 6: implement the restricted semantic compiler

Create the missing ontology module and measurement mappings. Generate traceable constraints and explicit rejection reasons for unsupported mappings.

Compare bottom-up aggregation, ordinary projection, diagonal weighting, shrinkage MinT, nonnegative reconciliation and a matched flat implementation. Estimate residual covariance inside training only.

Separate numerical correctness from predictive benefit.

Step 7: curate the smallest mechanistic component needed

Choose representative GSMMs only after reaction-level audits. Preserve the existing negative pilot results. Require a shared gas-supply calculation and balanced conservation equations before reporting ecosystem energy or carbon flux.

If the mechanism remains unidentifiable, report feasible ranges or conditional scenarios. Do not manufacture point estimates by relaxing maintenance demands or prescribing expression fractions.

Step 8: build a supported corridor atlas

Map empirical predictions first, then supported conditional responses. Include uncertainty, data support, extrapolation and abstention.

Only subsequently optimize a small number of well-defined interventions under common accounting, realistic constraints and joint uncertainty. Continuous optimization is a computational choice, not evidence of ecological validity.

Step 9: restructure the manuscript around findings

For NEE, I would use five principal displays:

- Figure 1: sampling design, actual assay overlap and independent exposure support.
- Figure 2: paired ecological contrasts and separation of turnover from within-lineage effects.
- Figure 3: between-site versus within-site relationships, including exposure-product sensitivity.
- Figure 4: temporal forecasting and the incremental contribution of reconciliation.
- Figure 5: spatial domain of supported inference, with a narrowly justified ecological implication.

Move software architecture, theorem proofs, standards exports, exhaustive parameter tables and unsupported intervention demonstrations to the Supplement or a separate methods paper.

For Nature Sustainability, the final two figures would instead need to establish a concrete decision: feasible actions, meaningful outcomes, water or land constraints, uncertainty, opportunity costs and robust preferred choices. Without those results, retain the NEE framing.

The title should name the finding once it is known. A defensible working title is “Spatial and temporal controls of microbial energy-acquisition potential across a hyperarid desert.” Stronger titles about restoration trade-offs should wait for stronger evidence.

Step 10: apply an editorial stop rule

Proceed toward either target only if the reanalysis establishes a result that survives measurement and exposure sensitivity, adds something materially distinct from the companion papers, and supports a broad ecological or sustainability conclusion.

A repaired pipeline, larger model, automated ontology compiler and attractive atlas would make EQ-BioDT a better system. The Nature-level contribution must be the scientific conclusion those components allow the authors to establish.