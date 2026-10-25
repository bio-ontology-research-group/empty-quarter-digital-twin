**EQ-BioDT is not ready for submission to Nature Sustainability in its present form. I recommend keeping Life on Air separate.** The strongest defensible EQ-BioDT contribution is a reproducible forecasting and constraint-reconciliation framework with explicitly conditional management scenarios. The current claims of validated intervention effects, regional carbon sequestration, and geographically supported prescriptions exceed what the implementation and evidence establish.

I reran the complete validation workflow without writing outputs. Its numerical results reproduced the stored report; two citation strings differed. Four selected tests also passed. I additionally traced cohort membership, recalculated CBB shares, tested failure cases, and checked the supplied figure PDFs against their captions. Thus, most problems below concern **what reproducible calculations mean**, rather than unexplained numerical drift. I did not regenerate upstream sequencing analyses.

**A. Nature Sustainability readiness: the principal scientific and editorial problems**

The additions improve the manuscript: temporal holdout evaluation, WLS/MinT comparisons, explicit abstention, and dimensional accounting address legitimate weaknesses. However, several additions currently create stronger claims than their supporting analyses warrant.

Nature Sustainability seeks research connecting environmental understanding with consequential sustainability problems and solutions. A model-only contribution can fit; new experiments are not categorically necessary. Here, the missing step is a defensible demonstration that the framework improves a sustainability decision under realistic uncertainty and constraints. The present portfolio rankings largely reflect prescribed response functions and preferences. [Journal scope](https://www.nature.com/natsustain/aims)

1. **The carbon headline is not an estimate of net atmospheric carbon removal.**

   The manuscript calls the baseline 1.49 Mt CO₂-eq yr⁻¹ and restoration 9.06 Mt CO₂-eq yr⁻¹ “binding” or “sequestration.” The accounting adds quantities with different interpretations:

   - Vegetation NPP is not net ecosystem carbon accumulation: heterotrophic respiration, turnover, mortality and losses remain outside the calculation.
   - Microbial carbon fixation is not necessarily persistent storage.
   - Carbonate accumulation is not automatically net atmospheric CO₂ removal. The manuscript’s own precipitation reaction releases CO₂; the net balance depends on the origin of calcium, alkalinity and dissolved inorganic carbon.
   - Existing soil H₂ uptake is not automatically an additional, creditable climate benefit. A marginal atmospheric counterfactual is needed before interpreting uptake × GWP as intervention mitigation.
   - Water production, pumping, establishment and maintenance emissions are absent.

   The distinction matters numerically. Vegetation contributes about **94.4% of P2’s basin CO₂-eq total**; direct microbial fixation contributes less than 0.1%. Consequently, the flagship mitigation result is driven predominantly by transferred vegetation assumptions, rather than the measured microbial data.

   **Reconcile:** report separate gross biological fixation, possible stock accumulation, carbonate accumulation and H₂ climate-equivalent service. Present intervention-minus-baseline differences. Do not combine them into “net sequestration” without an explicit carbon balance.

   Evidence: [carbon results and table](/home/leechuck/Public/software/empty-quarter-digital-twin/paper/main.tex:275), [carbon implementation](/home/leechuck/Public/software/empty-quarter-digital-twin/src/eq_biodt/models/carbon_budget_engine.py:249).

2. **The carbon Methods describe a materially different model from the executable code.**

   This is a submission blocker, even though the reported totals reproduce.

   | Component | Manuscript | Current implementation |
   |---|---|---|
   | Gas uptake | Power × depth × time / ΔG | Power scaled to prescribed annual uptake references |
   | H₂ uptake reference | Thermodynamic conversion | 12.5 mg m⁻² yr⁻¹ at 28 µW m⁻³ |
   | CO uptake reference | Thermodynamic conversion | 50 mg m⁻² yr⁻¹ at 1.95 µW m⁻³ |
   | Vegetation coefficient | 340 [220,480] | 295 [190,420] t CO₂ km⁻² yr⁻¹ |
   | Rainfall response | Clipped linear function | Clipped power function |
   | Inoculation multipliers | 0.28 and 0.18 | 0.25 and 0.15 |
   | H₂ GWP₁₀₀ | 11.6 | 11.0 |
   | CBB allocation | Stated conditional share and FBA flux ratio | State-dependent clipped share and prescribed multipliers |
   | Regional calculation | Each supported cell evaluated at its climate | Mean of 60 surface-state budgets multiplied by area |

   The power-to-moles equation also omits the **10⁻⁶ conversion from µW to W**. More fundamentally, correctly applying the stated equation gives:

   - 28 µW m⁻³ over 0.10 m: **0.878 mg H₂ m⁻² yr⁻¹**, versus the code’s 12.5.
   - 1.95 µW m⁻³ over 0.10 m: **0.718 mg CO m⁻² yr⁻¹**, versus the code’s 50.

   These differ by factors of approximately **14.2 and 69.6**. The quantities cannot simultaneously retain their stated physical interpretations without explaining an additional conversion or efficiency.

   **Reconcile:** choose one dimensional model, implement it once, generate the Methods and parameter table from that implementation, and rerun the budgets.

   Evidence: [Methods](/home/leechuck/Public/software/empty-quarter-digital-twin/paper/main.tex:398), [parameters](/home/leechuck/Public/software/empty-quarter-digital-twin/src/eq_biodt/models/carbon_budget_engine.py:48), [uptake calculation](/home/leechuck/Public/software/empty-quarter-digital-twin/src/eq_biodt/models/carbon_budget_engine.py:224).

3. **The purported CBB community share is actually an rbcL-only share.**

   The loader computes an `rbcL AND prk` flag, then ignores it and selects representatives using **rbcL alone**. The main text describes rbcL + prk, while the supplement introduces another definition involving RuBisCO forms or prk.

   Recalculation using the existing annotation and abundance tables, with the same 155-state aggregation, gives:

   | Screen | Catalogue representatives | Mean mapped-read share |
   |---|---:|---:|
   | rbcL alone—the implemented screen | 1,078 | **22.9348%** |
   | rbcL + prk | 416 | **13.2087%** |
   | rbcL + rbcS + prk—Life on Air’s candidate screen | 315 | **11.2340%** |

   These are alternative marker-screen results, **not independently established functional CBB prevalence**. Form-specific annotation, incomplete genomes and pathway completeness still require treatment. The three-marker screen should not automatically become a universal definition of every CBB implementation.

   The engine additionally clips conditional CBB carriage to at least 0.15, so a genuine zero can become a positive fixation allocation.

   **Reconcile:** distinguish rbcL carriage, candidate CBB carriage, and supported complete pathways; preserve zero and missing values; rerun the carbon allocation and benchmark comparisons.

   Evidence: [loader](/home/leechuck/Public/software/empty-quarter-digital-twin/src/eq_biodt/data/manifest_registry.py:696), [supplementary definitions](/home/leechuck/Public/software/empty-quarter-digital-twin/paper/supplementary.tex:375).

4. **“Matched seven cross-desert benchmarks” is not presently defensible as external validation.**

   The validator checks whether seven model outputs lie inside hard-coded intervals. Several intervals overlap the calibration priors, and several comparisons use incompatible quantities.

   | Benchmark | Audit finding | Required correction |
   |---|---|---|
   | Global CBB share | Leung reports approximately **21% primary producers across pathways** and **17% rbcL-bearing cells**. Neither is the implemented mapped-read CBB share. | Match marker definition and denominator; establish the provenance of the claimed 17–28% interval. |
   | Dark CO₂ fixation | The cited 0.01–1 range concerns radiolabel incorporation under several conditions, including processes other than H₂-supported autotrophy. The model’s annual rate is concentrated into an assumed 15-day, 2-cm window. | Compare assay-matched quantities; treat duty cycle and depth as uncertain assumptions. |
   | Namib H₂ oxidation | Rates are **per 16S gene copy**, whereas EQ-BioDT divides by estimated guild cells. The interval mixes wet-soil and dry-soil observations. | Keep wet/dry treatments separate and harmonize normalization. |
   | Annual areal H₂ uptake | The 5–25 range is also a model input envelope. | Separate calibration evidence from independent evaluation; document annualization. |
   | Negev co-carriage | The EQ denominator is an exclusive hydrogenase group; the exact matching population and derivation of the literature interval are insufficiently established. | Extract source-level counts using matching gene and carrier definitions. |
   | Mojave carbonate | The carbonate-to-CO₂ mass conversion can be valid, but does not establish net atmospheric removal; the same range calibrates the model. | Label mineral accumulation and distinguish calibration from validation. |
   | Haloxylon corridor total | Plantation stock accumulation, NPP and a four-pathway corridor total are different estimands, transferred across species, climates and stand ages. | Compare like quantities and retain transfer uncertainty. |

   Two source checks are especially decisive:

   - Tribbia reports dry-site rates of 5.31×10⁻⁸ and 2.47×10⁻⁸ nmol h⁻¹ **16S-copy⁻¹**, wet rates up to 2.20×10⁻⁶, and a wet hyperarid-zone mean of 2.82×10⁻⁸. These do not constitute the manuscript’s single cell-specific interval. Wetting also accelerated oxidation, which must be confronted by the moisture model. [Tribbia et al.](https://journals.asm.org/doi/full/10.1128/aem.00265-26)
   - Leung’s 21% and 17% refer to distinct functional populations; the tracer measurements are not exclusively CBB fixation powered by atmospheric hydrogen. [Author-posted primary manuscript](https://www.researchgate.net/publication/398120611_Aridity_drives_global_convergence_of_desert_microbiomes_and_biogeochemical_activities)

   Ma’s plantation study concerns *Haloxylon ammodendron* and age-dependent carbon-stock accumulation in a substantially different environment, not a directly measured full-cover NPP coefficient for *H. persicum* in the Empty Quarter. [Ma et al.](https://www.sciencedirect.com/science/article/pii/S2351989421001578)

   **Reconcile:** rename this section “external plausibility comparisons” until matched, independent predictions are evaluated. Do not make seven passing range checks a submission-readiness criterion.

5. **The Trip 5 analysis is retrospective temporal holdout, with information leakage vulnerabilities.**

   Training targets precede Trip 5, which is useful. However:

   - Validation features merge pH, DNA and XRF summaries from the **all-campaign compartment cube**. Historical predictors can therefore incorporate later measurements.
   - Some imputation occurs during cube construction, before fold-specific imputation.
   - Predictors include environmental information at collection time; there is no explicit forecast issuance time or lead horizon.
   - Training history features use other training campaigns, including campaigns later than the training observation being predicted.
   - There is no demonstrated pre-Trip-5 frozen forecast.

   These points do not make the results meaningless, but they preclude “strictly prospective forecasting.”

   The incremental evidence is also modest: longitudinal ridge RMSE **1.836**, environmental ridge **1.849**, random forest **1.803**, gradient boosting **1.794**. The longitudinal model does not beat the stronger baselines on squared-error accuracy. Change-prediction R² is approximately **0.055**, with an interval spanning zero. Correlation after subtracting the same historical baseline from observations and predictions should not alone establish temporal skill.

   **Reconcile:** rebuild every feature as of a declared cutoff, use rolling-origin evaluation, and report paired differences in forecast loss against persistence and environmental models. Call the present analysis a retrospective holdout.

   Evidence: [feature merge](/home/leechuck/Public/software/empty-quarter-digital-twin/src/eq_biodt/validation/twin_validator.py:249), [Trip 5 implementation](/home/leechuck/Public/software/empty-quarter-digital-twin/src/eq_biodt/validation/twin_validator.py:612), [reported results](/home/leechuck/Public/software/empty-quarter-digital-twin/paper/main.tex:185).

6. **Reconciliation works, but its present demonstration establishes limited methodological novelty.**

   The supplement appropriately acknowledges that the projection theorem is standard and that observed aggregate targets are constructed from their leaves. Consequently, per-sample non-increase in total squared error is expected.

   Further limitations:

   - Different features and penalties deliberately make the base predictors incoherent.
   - Nonnegative diagonal WLS performs slightly better than Euclidean projection.
   - The covariance estimates use in-sample residuals; inner out-of-fold residuals would provide a more credible comparison.
   - The full module’s subsumption inequalities are redundant under its exhaustive partitions and nonnegativity.
   - The benchmark does not demonstrate an advantage from more expressive ontology constraints.
   - The theorem concerns the chosen vector norm, not individual targets, biological correctness or downstream decisions.

   **Reconcile:** frame the contribution as traceable compilation and reliable enforcement of measurement constraints. Add nested tuning, paired uncertainty on method differences, and an ablation showing what the ontology contributes beyond a manually specified summing matrix.

   Evidence: [reconciliation Methods](/home/leechuck/Public/software/empty-quarter-digital-twin/paper/main.tex:420), [theorem qualification](/home/leechuck/Public/software/empty-quarter-digital-twin/paper/supplementary.tex:139).

7. **Paired field contrasts do not identify intervention effects.**

   Root-adjacent versus bulk soil and surface versus subsurface soil are observational compartment contrasts. They do not identify the effects of planting, inoculation, soil excavation or thermal manipulation.

   Matching a pooled contrast used to choose coefficients is calibration, not independent validation. Plant-associated soil also differs in substrate inputs, host selection, biomass, moisture and history. A surface/depth contrast cannot isolate temperature.

   Particularly vulnerable claims are:

   - Extract DNA concentration → microbial biomass → absolute gene density.
   - Natural roots → newly planted landscape mosaics.
   - Candidate plasmid carriage → successful transfer, expression and increased C/N flux.
   - Higher Shannon diversity → improved ecological condition.
   - Unknown protein function → active beneficial “dark proteome.”

   The current 1.38-fold absolute density increase is therefore a **model-derived proxy**, not an observed increase in hydrogenase copies per soil mass.

8. **The support mask does not validate the decision atlas or basin extrapolation.**

   The implementation uses distance plus MAP/elevation ranges, whereas the manuscript describes a climate hull. Climate values at grid points are interpolated from sampled sites, making range membership a weak test of external environmental support. Distances use degrees ×111 rather than geodesic distance.

   More importantly:

   - Geographic support for an observational predictor is not support for an intervention dose.
   - The atlas flags Sites 59/60 as requiring another habitat module, yet the manuscript assigns sabkha prescriptions.
   - Carbon upscaling averages 60 surface states, including imputed states, rather than evaluating grid cells.
   - The unsampled 557,840 km² receives transferred rates. Holding that region unchanged between scenarios does **not** make its baseline observed.
   - The supplied atlas function returns support metadata; I found no corresponding calculation substantiating the claimed **118 P3 / 8 P2 / 2 P4** assignments.

   **Reconcile:** separate observational prediction support, intervention support and regional transfer assumptions. Use a real geographic boundary and cell areas; report unsupported-region budgets separately, or omit them.

   Evidence: [mask construction](/home/leechuck/Public/software/empty-quarter-digital-twin/src/eq_biodt/data/build_twin_state_cube.py:476), [atlas](/home/leechuck/Public/software/empty-quarter-digital-twin/src/eq_biodt/interventions/pareto_optimizer.py:509), [upscaling](/home/leechuck/Public/software/empty-quarter-digital-twin/src/eq_biodt/models/carbon_budget_engine.py:387).

9. **Portfolio ranking is a preference-conditioned scenario comparison, not a validated optimum.**

   The six engines start from the same baseline independently. Their outputs are combined without jointly updating shared water, vegetation, soil, nutrient or microbial states.

   The objective functions combine Shannon diversity, unknown-protein fractions, power, dimensionless fluxes and risk scores. The “water/energy cost” is a weighted index, not water volume, energy consumption or money. Relative scaling depends on the chosen candidate set.

   The uncertainty analysis adds assumed Gaussian noise to objective values; it does not propagate fitted parameter, climate, measurement or structural uncertainty through the engines. Overlapping intervals do not establish equivalence.

   The ranking narrative also omits **P5’s 13.5% probability of ranking first** at Site 30. Its omission materially simplifies the apparent three-way comparison.

   **Reconcile:** report conditional rankings and regret across explicit preference and model-assumption ranges. Include all six portfolios and a genuine zero-intervention control. Remove “equivalent choices” unless an equivalence margin and appropriate analysis support it.

10. **The sustainability decision itself needs sharper definition.**

   Increasing shrub cover or Shannon diversity is not automatically restoration of a naturally sparse desert ecosystem. Halobacteriota abundance is not intrinsically ecological degradation in a natural sabkha.

   The spatial scale makes feasibility central: a 6-mm pulse across 92,160 km² represents approximately **553 million m³ of precipitation**. This is not necessarily supplied irrigation water, but it shows why a dimensionless water-cost score is inadequate. Cloud-seeding efficacy is also distinct from the ecological response to a prescribed rain pulse.

   A stronger question would be: **Which management choices remain preferable across plausible microbial responses, water constraints and conservation priorities—and where should the model abstain?**

11. **The figures and captions are inconsistent across all three main figures.**

   I checked the supplied PDFs and generating script:

   - **Figure 1d:** caption describes MAG-versus-16S phylum shares; the implementation shows 16S composition and unmapped shotgun-read fraction.
   - **Figure 2:** panel assignments differ from the caption. The supplied figure places precipitation slopes in b, reconciliation in c, and paired contrasts in d. Panel c also says **n=554**, whereas the benchmark uses 550.
   - **Figure 3b:** the caption describes inoculated versus uninoculated curves; the figure compares area-weighted GE, biomass-weighted GE and absolute-density multiplier under one inoculation setting.
   - **Figure 3c/d:** the supplied PDF shows brine responses and the Pareto comparison. It does **not** contain the decision atlas described as panel d.

   These are direct, visible submission defects, not interpretation differences.

   Evidence: [figure generator](/home/leechuck/Public/software/empty-quarter-digital-twin/paper/generate_figures.py:170), [Figure 2 PDF](/home/leechuck/Public/software/empty-quarter-digital-twin/paper/figures/fig2_validation_causal_polytope.pdf), [Figure 3 PDF](/home/leechuck/Public/software/empty-quarter-digital-twin/paper/figures/fig3_interventions_pareto.pdf).

12. **The manuscript requires substantial editorial compression.**

   The title is approximately **27 words/215 characters**; the abstract is approximately **413 words**. There are seven main display items. Nature Sustainability’s Article guidance specifies a short title, a 150-word abstract, approximately 3,500 main-text words and up to six display items. My rough source count is approximately 9,580 words including main-text tables, captions and equations, so a formal count is needed, but substantial reduction is clearly required. Repository/accession/acknowledgement placeholders remain. [Article requirements](https://www.nature.com/natsustain/content)

   The abstract should lead with one substantive contribution and its strongest independent evidence. The present accumulation of prediction, ontology, intervention and carbon claims obscures their different evidential status.

**B. Cross-paper numerical consistency and overlap**

Several apparent discrepancies are legitimate nested denominators. Others are confirmed implementation errors.

| Item | Audit result | Reconciliation |
|---|---|---|
| Sites | 70 registered sites = 60 primary +10 auxiliary; Life on Air’s retained WGS cohort covers 59 primary sites. | State registered, sampled and assay-retained site counts separately. |
| WGS libraries | EQ metadata: **551**. One sample, **F5Dr1**, has zero detected SCGs. Four auxiliary libraries are **62PRr2, 64Dr3, 64PRr2, 64Sr2**. Removing these five gives **546**, matching Life on Air. | Explain 551 metadata /550 finite guild /546 primary analytical libraries. |
| WGS compartments | EQ: **117/166/268** surface/deep/root; Life on Air: **116/164/266**. | Exactly reconciles through the exclusions above. |
| Guild validation | EQ’s 550 includes four auxiliary libraries; Life on Air’s 546 does not. | Use the same primary cohort for the principal comparison; retain auxiliary transfer as a separate test. |
| 16S profiles | Descriptor/ecology: **1,271 −24 controls −10 low-depth =1,237**, including **1,227 primary-site profiles**. EQ uses **1,212**. | Confirmed join/QC bug, detailed below. |
| 16S join defect | Against the ecology canonical cohort, EQ omits **16 retained preparation-suffix profiles**, but includes **e681_T1Dr1**, with **934 reads**. Net difference: −15. | Use canonical preparation-to-specimen mappings and profile-level QC. Do not repair by arbitrary suffix stripping. |
| Paired root/bulk | **106 visits/47 sites**, median log₂ contrast approximately **−2.27**, agrees. | Retain paired-visit and site denominators explicitly. |
| Depletion frequency | Life on Air: **89% of visits**; EQ: **87% of sites**. | Potentially compatible different aggregation levels; regenerate both from one paired table. |
| DNA contrast | log₂ **3.37**, approximately **10.3-fold**; Life on Air’s DNA evidence has **47 paired visits**, not all 106 guild pairs. | State the DNA-specific denominator and avoid interpreting extract concentration as measured biomass. |
| Surface/deep | **68 visits/36 sites**, raw contrast approximately **−1.05**, agrees. | Do not confuse raw contrasts with Life on Air’s marker-control-adjusted effects. |
| Turnover | Approximately **97–97.4%**, consistent in broad conclusion. | EQ’s GE decomposition and Life on Air’s carrier-fraction/copy-number decomposition are different estimands; report scales and aggregation. |
| Long-term precipitation | EQ emphasizes blended-product ρ≈**−0.514**. Life on Air reports approximately **−0.40 ERA5**, **−0.43 POWER**, **−0.51 blended**, with important spatial/compositional sensitivity. | Import the current robustness qualifications, not only the strongest correlation. |
| Within-site precipitation | EQ’s raw slopes and Life on Air’s IQR-scaled slopes need not match numerically. Surface/root ERA5 effects differ from POWER; deep effects are imprecise. | Give predictor transformation, scaling, product and window beside each coefficient. |
| Product agreement | 90-day ρ≈**0.193**, sign agreement **55.7%**; annual ρ≈**0.612**, sign agreement **79.5%**, n=244. | Consistent with companion rounding; retain as exposure uncertainty. |
| Synthetic benchmarks | Counts **18/18/17/19**, mean Pearson **0.9649**, mean Spearman **0.9580** and eight-pair R² **0.9888** agree. | Explain 20 simulated communities versus guild-specific evaluable subsets and eight evaluable contrasts. |
| Synthetic bias | log₂ biases **0.6653/0.8042/0.5959/0.8410** agree. | Propagate calibration uncertainty; full-reference rank recovery does not prove unbiased absolute capacity. Include the companion’s held-out-reference sensitivity. |
| Hydrogenase carriers | **597 group 1l** and **410 group 1h** are inclusive counts. **558** and **371** are exclusive groups; **39** carry both. | The denominators are reconcilable, but EQ labels exclusive groups as all carriers. |
| CBB co-carriage | Exclusive groups: **101/558=18.10%**, **30/371=8.09%**; dual group: **11/39**. Inclusive counts therefore become **112/597=18.76%** and **41/410=10.00%** under that screen. | Label exclusive/inclusive populations and use the matching annotation definition. |
| pH | **712 accepted measurements**, **178 EQ states**, **709 ecology-linked measurements**, and **560 composition/563 alpha-diversity groups** concern different joins and aggregations. | Publish a denominator ledger; do not force equality. |
| XRF | **725 laboratory records** and **71 field sessions** are distinct. Eleven retained laboratory analytes are not the number of measurements or all reported channels. | Keep assay and analyte counts separate. Preserve unknown units. |
| MAG catalogues | Current catalogue **5,438**; ecology’s older validation uses **991 genomes**, a **990-genome KO intersection**, and **975 annotated genomes**. | Legitimate separate releases, but require explicit catalogue IDs and hashes. |
| Regional area | Descriptor says approximately **500,000 km²**; EQ budgets use **650,000 km²**. | Define and archive the geographic polygon; area is an analytical denominator, not a cosmetic wording difference. |
| Climate summaries | Ecology’s monthly, limited-period summaries differ from long-term annual climate products in EQ. | Specify product, dates, aggregation and spatial extraction; do not present them as interchangeable MAP/MAT. |

The 16S omission is traceable to the exact-name merge in [the registry](/home/leechuck/Public/software/empty-quarter-digital-twin/src/eq_biodt/data/manifest_registry.py:658). The companion already supplies the needed preparation mappings in its [profile join audit](/home/leechuck/Public/software/empty-quarter-ecology-data/evidence/ecology-canonical/depth_extraction/profile_extraction_join_audit.tsv).

There are also **source-version risks beyond counts**:

- The descriptor explicitly records a sequence-identity issue involving **1,010 refined `_sub` representatives** versus substituted parent bins, and a rebuilding requirement for annotations and recruitment. EQ still consumes existing Life on Air annotation/coverage files. Identical representative names and a shared total of 5,438 do not demonstrate identical sequences. Resolve this dependency through FASTA, annotation and recruitment hashes before freezing any carrier or abundance claim. [Descriptor release notes in source](/home/leechuck/Documents/papers/empty-quarter-data-current/02_methods.tex:216)
- EQ contains unused macros for **636 specimens, 35,342 ASVs, 56 phyla and 372 MAGs**. These are not currently rendered numerical claims, but should be removed to prevent stale values re-entering the paper. Ecology’s current raw table reports 351,472 ASVs.
- The 72-MAG/10,584-run pilot description should say **one of 48 carrier models had the relevant native GPR wiring**, not one of 48 “keywords.”
- Mobilome catalogue constants differ between the manuscript’s macros and `SpatialJSDMViromeModel` metadata. Catalogue-wide totals, corrected catalogues and sample richness require explicit version labels even though they are outside the three principal companions.

**The IRI audit mostly supports identifier reuse, but reveals measurement-semantic problems.**

The registry correctly preserves the project namespace and auxiliary-site alias distinctions—for example, numeric Trip-1 site aliases 61–64 are not interchangeable with identically numbered canonical site-IRI suffixes. The supplied mapping audit supports the phylum identifiers used by the registry.

However:

- Supplementary ontology semantics define the phylum partition using **16S reads**, while scenario variables and exports commonly use **MAG-mapped read shares**.
- `RAK_0000078` is an abundance **quality**, not a unit.
- `genome_equivalents_per_cell` contradicts “copies per genome equivalent.”
- XRF exports assign **percent** despite the descriptor explicitly withholding undocumented units.
- Annual precipitation’s text unit is mm/year but its unit IRI denotes millimetres alone.
- A zero sequence-based estimate should not automatically become biological absence.
- Functional gene classes should not be represented as ordinary taxonomic ranks without a clearly documented extension.
- B-Cubed occurrence counts currently reuse sample counts; these need a defensible detection/occurrence interpretation.
- The twin constraint module and KG v3.0.1 need separately identified versions and provenance.

Evidence: [export implementation](/home/leechuck/Public/software/empty-quarter-digital-twin/src/eq_biodt/standards/dwc_emof_exporter.py:62), [descriptor measurement semantics](/home/leechuck/Documents/papers/empty-quarter-data-current/03_knowledge_representation.tex:48).

**The main overlap risk is duplication of results and claimed novelty, rather than merely sharing samples.**

| Manuscript | Results it should primarily own | Appropriate use in EQ-BioDT |
|---|---|---|
| Life on Air | Marker typing, relative group-1l depletion, controls, lineage turnover, co-carriage, spatial/temporal precipitation sensitivity, assay benchmarks | Cited biological evidence and calibration inputs; independently held-out model tests |
| Ecology | Community assembly, compartment diversity, distance decay, precipitation associations, PICRUSt reliability | Cited observational context and carefully defined response data |
| Data/KG | Sampling resource, identifier reconciliation, assay provenance, ontology/KG construction and release validation | Dependency description and the genuinely new constraint/decision layer |
| EQ-BioDT | Leakage-controlled forecasts, reconciliation, support assessment, uncertainty-aware decision analysis | Main novel contribution |

EQ currently gives several Life on Air results enough prominence to look like a second discovery paper: the root depletion, turnover explanation, CBB distinction, aridity relationship and synthetic benchmark. Calling these “validation tiers” does not make them independent new findings.

Likewise, the ecology paper’s exploratory rainfall–**richness** relationship does not establish EQ’s prescribed rainfall–**Shannon** response. Ecology also reports lower root-associated alpha diversity, whereas EQ assigns positive revegetation diversity increments. A landscape-mixture interpretation could reconcile these, but it must be calculated rather than assumed.

Related manuscripts should be disclosed and supplied to the editor, with an explicit contribution map. Citation alone does not resolve substantial duplication of principal results. [Nature Portfolio policy](https://www.nature.com/commsenv/editorial-policies/plagiarism)

**C. Keep Life on Air separate, or merge?**

**My recommendation is to keep them separate, while substantially narrowing EQ-BioDT.** Life on Air’s current manuscript and submission notes already present a coherent, cautious ISME Journal story. Its evidence does not become stronger by attaching unvalidated management and carbon claims.

| Consideration | Option 1: Separate | Option 2: Merge |
|---|---|---|
| Scientific argument | Clear distinction between biological discovery and prediction/decision methods | Potentially stronger observation-to-decision narrative |
| Main advantage | Protects Life on Air’s empirical contribution; allows EQ’s methods to mature independently | Removes repeated presentation of the biological anchor results |
| Main disadvantage | EQ must demonstrate substantial novelty beyond repackaging companion findings | Scope becomes very broad; weaker modelling claims can undermine the empirical paper |
| Referee burden | More focused microbial and modelling reviews | Requires confidence across microbial ecology, forecasting, ontology, biogeochemistry and restoration |
| Publication risk | Manageable overlap through clear ownership and disclosure | Merging consumes Life on Air’s standalone publication while not guaranteeing flagship suitability |
| Present recommendation | **Preferred** | Consider only after narrowing to one integrated claim and repairing the modelling foundation |

A strong **separate EQ-BioDT structure** would be:

1. What management-relevant microbial quantities can be predicted?
2. Does temporally and spatially honest validation demonstrate useful skill?
3. What does reconciliation improve beyond established alternatives?
4. Where do predictions and interventions lack support?
5. Which decisions remain preferable under uncertainty and competing objectives?

The companion biological findings would appear briefly as cited inputs. The six interventions would become a broader supplementary scenario library; the main paper would develop one or two defensible decisions.

A viable **merged structure** would instead be:

1. Group-specific microbial responses to plant association and soil position.
2. Evidence distinguishing lineage replacement from copy-number change.
3. Limits of temporal and spatial environmental transport.
4. A narrowly evaluated predictive model of those responses.
5. Conditional implications for vegetation expansion or soil disturbance.

That merged paper should probably **drop most of the six-intervention catalogue and the current basin-sequestration headline**. Otherwise it is a concatenation of manuscripts rather than a stronger scientific argument.

For the journals you named:

- **Nature Sustainability:** strongest fit if the central advance is a robust, consequential management decision involving water, biodiversity and climate trade-offs.
- **Nature Microbiology:** strongest fit for a general microbial mechanism; the current evidence remains largely gene representation and observational association.
- **Nature Ecology & Evolution:** potentially suitable for a general ecological assembly or functional trade-off result; ontology/software breadth alone will not provide that advance.
- **Nature:** neither packaging currently establishes the breadth and strength of advance expected; merging does not repair identification or validation weaknesses.

These are editorial judgments, not predictions of acceptance. Nature’s stated emphasis is exceptional scientific importance and interdisciplinary interest; Nature Microbiology likewise emphasizes the strongest research across microbiology. [Nature criteria](https://www.nature.com/nature/journal-information), [Nature Microbiology](https://www.nature.com/nmicrobiol/journal-information)

**D. Additional code and method defects requiring correction**

The relevant implementations are named `causal_ontological_twin.py`, `moisture_pulse_diffusion.py` and `metabolic_community_engine.py`; the two alternative filenames in the request are not present.

Beyond the major issues above, I found the following concrete problems.

1. **Missing and NaN states can receive OCS=1.**  
   I reproduced:
   - `evaluate_constitutive_layer({})` → default biological values, OCS=1.
   - `evaluate_constitutive_layer({"hydB_1l_ge": NaN})` → group-1l value 0, OCS=1.
   - `polytope.project({})` → completed zero vector, OCS=1.

   Direct validation rejects incomplete states, but public preprocessing/projection routes can conceal the problem. Validate required inputs before defaulting or clipping; keep completeness and consistency as separate statuses. [Code](/home/leechuck/Public/software/empty-quarter-digital-twin/src/eq_biodt/models/causal_ontological_twin.py:216)

2. **Projection fallback is not a certified metric projection.**  
   On solver failure it clips and repairs partitions without establishing optimality or exposing failure. The theorem does not apply to that fallback. Return solver status and residuals; fail explicitly when projection cannot be certified. Validate weighting-matrix dimensions and positive definiteness.

3. **The full projection metric mixes physical scales.**  
   Fractions, GE and µW m⁻³ cannot be treated as scientifically equivalent Euclidean coordinates without a declared scaling or uncertainty model. Changing units can change the solution. Scenario preprocessing usually makes sums coherent before projection, so perfect scenario OCS provides little independent evidence.

4. **The metabolic dormancy fallback can create maintenance energy.**  
   With H₂=CO=external N=sulfate=0, DOC bound=0.1 and EC=216.2, the fallback returns zero DOC and other donor fluxes but an ATP maintenance drain of **0.22**. The theoretical ATP maximum includes DOC that the returned vector does not consume. Solve a separate maintenance-only LP and verify its complete returned flux balance. [Code](/home/leechuck/Public/software/empty-quarter-digital-twin/src/eq_biodt/models/metabolic_community_engine.py:203)

5. **The reduced LP is not a genome-scale, element-balanced metabolic model.**  
   It has ten flux variables and a few inequalities, without explicit carbon/electron/element mass balances. Nitrogen pools are combined, and sulfur energy accounting needs donor bookkeeping. “Thermodynamic invariants” currently test selected imposed relationships, not full thermodynamic validity. Describe it as a reduced resource-allocation model unless those balances are added. [Constraints](/home/leechuck/Public/software/empty-quarter-digital-twin/src/eq_biodt/models/metabolic_community_engine.py:171)

6. **Solver failure is conflated with biological dormancy.**  
   Infeasibility, numerical failure and other unsuccessful solver outcomes should not all become `maintenance_deficit_arrest`. Preserve native solver status and distinguish model failure from predicted biology.

7. **Carbon and metabolic maintenance are disconnected.**  
   Regional carbon uses `EC <150` as its feasibility rule rather than the metabolic solver’s ATP requirement and supply. Several intervention inputs never reach the carbon calculation. Plasmid boosts are applied through overlapping multipliers. Integrate a shared state transition or explicitly label the calculations independent. [Code](/home/leechuck/Public/software/empty-quarter-digital-twin/src/eq_biodt/models/carbon_budget_engine.py:430)

8. **Moisture-response uncertainty excludes major biological mechanisms.**  
   The diffusion response does not adequately establish the balance between reduced gas transport and moisture-induced activation. The Shannon response is prescribed, not fitted to the ecology paper’s reported endpoint. Daily sampling also misses the immediate post-pulse extreme. Test alternative activation–diffusion models and timestep convergence. [Model](/home/leechuck/Public/software/empty-quarter-digital-twin/src/eq_biodt/models/moisture_pulse_diffusion.py)

9. **Revegetation mixes incompatible abundance concepts.**  
   Area-weighted GE, biomass-weighted GE and absolute density are useful distinctions, but DNA extraction concentration does not establish the biomass conversion. Carbon then applies a biomass multiplier raised to 0.35, rather than the manuscript’s absolute-density scaling. Preserve these as explicit alternative assumptions. [Revegetation](/home/leechuck/Public/software/empty-quarter-digital-twin/src/eq_biodt/interventions/revegetation.py:90)

10. **Several ecological mechanisms are asserted more specifically than their evidence supports.**  
    Examples include phylum-level attribution of nitrification/3-HP/4-HB activity, calcium-sulfate “buffering” of osmotic stress, and gene/plasmid inoculation benefits. Map pathway assumptions to actual gene evidence and distinguish pool-water chemistry from soil exposure.

11. **The virome/JSDM component is a heuristic response model.**  
    A Mantel correlation of 0.254 becomes an exponential response coefficient; these are not interchangeable statistical quantities. The class contains fixed response functions rather than a fitted joint species distribution model. Its ARG threshold is not a validated biosafety classifier. [Code](/home/leechuck/Public/software/empty-quarter-digital-twin/src/eq_biodt/models/spatial_jsdm_virome.py:62)

12. **“Active dark proteome” is a renamed coding fraction.**  
    `dark_proteome_active_fraction` is returned equal to `dark_proteome_coding_fraction`, without expression or activity data. Remove “active” and do not treat the quantity as a validated ecosystem-benefit indicator.

13. **Sibling-effect uncertainty ignores shared sampling.**  
    The three hydrogenase contrasts come from overlapping paired visits. Their errors are correlated; independent random-effects weighting does not capture that covariance. Estimate aggregate contrasts directly and bootstrap sites jointly. The supplement acknowledges this limitation but the calculation still needs replacement.

14. **Validation should be strengthened without multiplying superficial tests.**  
    Existing passing tests include checks that outcomes fall within the same ranges used to set parameters, or satisfy monotonicity imposed by formulas. Add tests for the actual failures above: missingness, valid units, donor/energy balance, solver failure, zero-intervention identity and fold/cutoff isolation.

15. **Provenance is incomplete despite the FDO terminology.**  
    Run crates contain empty input-digest dictionaries and an ontology filename where a digest is expected. The registry uses Robert-specific absolute paths. Generated HTTP identifiers do not establish resolvable persistent identifiers. The JSON report also contains NaN tokens, which are not strict JSON. Package resolved inputs, hashes, model versions and a portable rebuild command; serialize unavailable values as null. [Provenance implementation](/home/leechuck/Public/software/empty-quarter-digital-twin/src/eq_biodt/standards/fdo_provenance.py:40)

16. **Citation and parameter provenance needs a systematic audit.**  
    The stored report has stale Schlesinger and plantation citations. The manuscript attributes H₂ GWP 11.6 to Warwick, whereas Warwick reports 12±6 and discusses the 11.6 estimate from other work. Correct the source, uncertainty interpretation and code value together. [Warwick et al.](https://acp.copernicus.org/articles/23/13451/2023/)

**E. Prioritized improvements using existing data**

I would implement the following in order.

1. **Freeze a common evidence release across all four papers.**  
   Resolve representative-sequence identity first. Produce machine-readable cohort manifests linking specimen, preparation, library, profile, visit, site and assay. Record exclusions and input hashes. Generate manuscript counts from these manifests.

   **Completion criterion:** every reported denominator has an executable derivation and an explicit relationship to companion denominators.

2. **Repair the confirmed implementation defects.**  
   Fix the 16S preparation join and read-depth filter; harmonize CBB screens; preserve true zero/missing carriage; repair metabolic fallback; reject incomplete/nonfinite projection inputs; expose solver status.

   **Completion criterion:** focused regression tests demonstrate each failure is corrected.

3. **Rebuild validation with fold-local and time-local inputs.**  
   Recompute summaries and imputation within training folds and forecasting cutoffs. Use rolling-origin evaluation across available campaigns and the primary 546-library cohort, with auxiliary-site transfer reported separately.

   **Completion criterion:** report uncertainty on improvement over persistence and environmental baselines, including within-site change—not just level correlations.

4. **Separate calibration, validation and plausibility evidence.**  
   Reserve independent sites/campaigns for evaluating calibrated field contrasts. Label borrowed biological results as companion findings. Replace “seven benchmarks passed” with source-matched comparisons and explicit transfer limitations.

   **Completion criterion:** each validation target is independent of the data or range used to set its model coefficients.

5. **Rebuild carbon accounting as a transparent conditional balance.**  
   Resolve the power/uptake mismatch, eliminate undocumented physical conversions, and distinguish production from storage. Report separate pathways and additional effects. Treat carbonate net removal and unsupported-region fluxes as unresolved or bounded assumptions.

   **Completion criterion:** every term has compatible dimensions, a stated counterfactual and a traceable source; no sequestration claim rests solely on NPP or carbonate accumulation.

6. **Exploit existing public evidence for meaningful external tests.**  
   Use the desert studies’ source data where available to evaluate matched gas-uptake or fixation predictions across moisture states and soils. Harmonize normalization first. Evaluate errors and calibration, not merely membership in broad intervals.

   **Completion criterion:** at least one externally evaluated process module has predictive performance that was not guaranteed by its parameterization.

7. **Replace objective noise with propagated uncertainty.**  
   Resample sites and paired visits; vary precipitation products, marker assignment/calibration, CBB definitions, biomass assumptions and process parameters; rerun the engines. Separate data uncertainty, parameter uncertainty, structural alternatives and preference uncertainty.

   **Completion criterion:** portfolio conclusions survive documented combinations of plausible assumptions, or are explicitly reported as unresolved.

8. **Make the management problem physically and ecologically meaningful.**  
   Concentrate on one or two decisions. Add water quantities, establishment/maintenance assumptions, intervention extent and habitat-specific conservation objectives. Include a genuine no-action reference and show trade-offs without requiring one universal utility weighting.

   **Completion criterion:** the analysis identifies either a robust action, a robust reason to avoid an action, or information needed before choosing.

9. **Validate abstention empirically and implement the actual atlas.**  
   Use geodesic distances and real environmental layers. Assess prediction error against distance/environmental novelty on held-out sites. Add intervention-dose support. Generate prescriptions only where both prediction and intervention evidence permit them.

   **Completion criterion:** the map is reproducible from per-cell evaluations, and abstention demonstrably changes the accuracy–coverage trade-off.

10. **Rewrite and regenerate the submission package last.**  
    Generate tables, counts, figures and captions from the frozen results. Remove unsupported “prospective,” “causal,” “sequestration,” “biosafety pass,” and “equivalent choices” wording where the revised analysis does not justify it. Disclose the three companions and clearly assign ownership of shared findings.

    A defensible working title is **“Forecasting desert microbial responses under ecological and management uncertainty.”**

I would make the Nature Sustainability submission decision after steps 1–8. If those analyses establish a consequential, robust management conclusion, the journal remains a reasonable target. If they establish useful forecasting and software infrastructure but unstable management conclusions, the manuscript should be positioned around that narrower contribution.