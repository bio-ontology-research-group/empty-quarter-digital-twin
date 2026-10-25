# Point-by-Point Response to Codex (`gpt-6-astra`) Cross-Paper & Nature Sustainability Review (`astra_review_cross_paper_v2.md`)

We thank the reviewer (`gpt-6-astra`, High Reasoning Effort) for the comprehensive audit of the `empty-quarter-digital-twin` (`EQ-BioDT`) manuscript, codebase, figures, and data products alongside the three companion manuscripts (`empty-quarter-lifeonair`, `empty-quarter-ecology-current`, and `empty-quarter-data-current`). Below we respond point-by-point to every scientific, methodological, cross-paper, figure, and software issue raised in `paper/astra_review_cross_paper_v2.md`, detailing the exact code, data-cube, validation, figure, and manuscript changes implemented and verified.

---

## Executive Decision: Keep *Life on Air* Separate vs. Merge

**Reviewer Recommendation (Section C):**
> *"My recommendation is to keep them separate, while substantially narrowing EQ-BioDT. Life on Air's current manuscript and submission notes already present a coherent, cautious ISME Journal story. Its evidence does not become stronger by attaching unvalidated management and carbon claims."*

**Response & Action Taken:**
We agree and have kept `empty-quarter-lifeonair` separate as the primary microbial ecology discovery paper (focused on high-affinity Group 1l `[NiFe]`-hydrogenase distribution, lineage replacement vs. copy-number shifts, synthetic benchmark calibration, and marker-control sensitivity), while positioning `EQ-BioDT` as the companion forecasting, ontological constraint-reconciliation, and conditional decision-support paper. In `paper/main.tex` and `paper/supplementary.tex`, we explicitly disclose the division of ownership across the four manuscripts, cite the biological findings of the companions as empirical inputs/priors rather than re-claiming them as novel biological discoveries, and frame the six interventions and regional carbon accounting as assumption-explicit, conditional scenario comparisons.

---

## Section A: *Nature Sustainability* Readiness — Scientific, Methodological, and Editorial Points

### A.1. Carbon Headline: Gross Pathway Binding vs. Persistent Storage & Supported Corridor Additionality
- **Critique:** The previous text summed vegetation NPP, microbial fixation, pedogenic $\mathrm{CaCO}_3$ precipitation, and $\mathrm{H}_2$ $\mathrm{GWP}_{100}$ climate-equivalent uptake into a single "binding/sequestration" headline without distinguishing annual gross turnover from persistent storage (humified soil organic carbon and silicate-derived carbonate) or separating the 92,160 $\mathrm{km}^2$ empirically supported corridor from the 557,840 $\mathrm{km}^2$ unsupported reserve.
- **Changes Made:**
  1. **Code (`src/eq_biodt/models/carbon_budget_engine.py`):**
     - Added explicit persistent carbon storage accounting alongside gross pathway binding:
       - Soil organic carbon humification efficiency $\eta_{\mathrm{SOC}} = 0.15\ [0.10, 0.25]$ applied to annual vegetation + rhizosphere NPP (`persistent_soc_co2_t_km2_yr`).
       - Non-carbonate silicate/dust $\mathrm{Ca}^{2+}$ weathering fraction $f_{\mathrm{sil\text{-}Ca}} = 0.50\ [0.25, 0.75]$ applied to pedogenic $\mathrm{CaCO}_3$ precipitation (`persistent_caco3_co2_t_km2_yr`), accounting for the stoichiometry $\mathrm{Ca}^{2+} + 2\,\mathrm{HCO}_3^- \rightarrow \mathrm{CaCO}_3 + \mathrm{CO}_2\uparrow + \mathrm{H}_2\mathrm{O}$ when $\mathrm{Ca}^{2+}$ derives from pre-existing carbonate dissolution.
       - Microbial necromass stabilization fraction $\eta_{\mathrm{mic}} = 0.25\ [0.15, 0.40]$ applied to direct dark $\mathrm{CO}_2$ fixation (`persistent_microbial_co2_t_km2_yr`).
     - Computed separate budgets for:
       1. The **92,160 $\mathrm{km}^2$ empirically supported corridor** (`corridor_*_kt_yr`, evaluating all 128 supported grid cells at their local climate and nearest observed primary-site state).
       2. The **557,840 $\mathrm{km}^2$ out-of-support reserve** (held at its unmanaged baseline, zero management additionality).
       3. Both the **650,000 $\mathrm{km}^2$ basin** and the **500,000 $\mathrm{km}^2$ core sand-sea polygon** (`sand_sea_500k_gross_co2eq_kt_yr` and `sand_sea_500k_persistent_co2_kt_yr`, reconciling with the descriptor's regional extent).
     - Added the **spatially prescribed corridor atlas** (`prescribed_corridor_atlas`: 87 P2, 34 P3, 7 P4 cells across the 128 supported corridor cells) to the regional carbon comparison table alongside uniform portfolios P1–P6 and the zero-intervention `observed_baseline`.
  2. **Manuscript (`paper/main.tex` Table 3, Abstract, Results, Discussion, Methods; `paper/supplementary.tex` Section S9 & Table S6):**
     - Replaced unqualified "sequestration" claims with side-by-side reporting of **gross pathway binding** (`811.3 kt CO2-eq yr^-1` in the 92,160 $\mathrm{km}^2$ corridor; `1,496.1 kt CO2-eq yr^-1` across the 650,000 $\mathrm{km}^2$ basin; `1,311.9 kt CO2-eq yr^-1` in the 500,000 $\mathrm{km}^2$ sand-sea polygon) and **persistent carbon storage** (`140.9 kt CO2 yr^-1` in the corridor; `344.6 kt CO2 yr^-1` across the basin), as well as intervention-minus-baseline additionality ($\Delta$) restricted strictly to the supported corridor.

### A.2. Reconciliation of Carbon Methods Equations with Executable Code
- **Critique:** Several equations and parameter values in the Carbon Methods differed from `carbon_budget_engine.py` (power-to-uptake conversion vs. reference dryland uptake scaling, vegetation NPP coefficient `340` vs. `295`, linear vs. power-law rainfall factor, inoculation multipliers `0.28/0.18` vs. `0.25/0.15`, $\mathrm{H}_2$ $\mathrm{GWP}_{100}$ `11.6` vs. `11.0`, and 60-state mean vs. 128-cell spatial grid evaluation).
- **Changes Made:**
  1. **Code (`src/eq_biodt/models/carbon_budget_engine.py`):**
     - Updated `H2_GWP100_MID` to `11.6` (`[6.0, 16.0]`) citing Warwick et al. (2023) and Paulot et al. (2024).
     - Computed and stored both the **uncalibrated thermodynamic power-to-mass conversion** ($P \times d \times t \times 10^{-6} / |\Delta G|$ over $d = 0.10\ \mathrm{m}$ at $10^6\ \mathrm{cells\ gdw}^{-1}$, yielding `0.45 mg H2 m^-2 yr^-1` and `0.62 mg CO m^-2 yr^-1` at Site 30, reported in `thermodynamic_h2_uptake_mg_m2_yr` and `thermodynamic_co_uptake_mg_m2_yr`) and the **diffusion-coupled dryland uptake** scaled to empirical desert chamber flux envelopes (`8.72 mg H2 m^-2 yr^-1` at Site 30), explicitly documenting the $10^{-6}\ \mathrm{W\ \mu W}^{-1}$ conversion and why intact soil column uptake exceeds the $10^6\ \mathrm{cells\ gdw}^{-1}$ reference cell-density calculation.
     - Implemented true **cell-by-cell spatial evaluation across all 128 empirically supported corridor cells** (`upscale_corridor_and_basin_carbon_budget`), linking each cell to its local climate (`map_mm_yr`) and nearest observed primary-site state rather than simply multiplying a 60-site average by area.
  2. **Manuscript (`paper/main.tex` Methods & `paper/supplementary.tex` Section S9 / Table S6):**
     - Aligned every equation and parameter in the Methods and Supplement exactly with `carbon_budget_engine.py`: $\mathrm{GWP}_{100} = 11.6\ [6.0, 16.0]$, $C_{\mathrm{veg},0} = 295\ [190, 420]\ \mathrm{t\ CO_2\ km^{-2}\ yr^{-1}}$, $\phi_{\mathrm{MAP}} = \mathrm{clip}((\mathrm{MAP}/22)^{0.35}, 0.55, 1.45)$, inoculation boosts $0.25\,d_{\mathrm{RAK\text{-}PR}} + 0.15\,d_{\mathrm{rbcL}}$, and explicit $10^{-6}\ \mathrm{W\ \mu W}^{-1}$ dimensional factors.

### A.3. CBB Community Share Definitions (`rbcL` vs. `rbcL+prkB` vs. `rbcL+rbcS+prkB`) & Preservation of Zero Carriage
- **Critique:** `load_cbb_cocarriage_state_profile()` computed an `rbcL AND prk` flag (`has_cbb`) but selected `cbb_mags` using `has_rbcL` alone, and `carbon_budget_engine.py` clipped `cbb_share_of_1l` to a floor of `0.15`, turning zero carriage into positive fixation.
- **Changes Made:**
  1. **Code (`src/eq_biodt/data/manifest_registry.py` & `src/eq_biodt/models/carbon_budget_engine.py`):**
     - Updated `load_cbb_cocarriage_state_profile()` to compute and export all three marker-screen stringencies across the 5,438 MAG representatives and 155 observed WGS states:
       - `rbcL`-only (`K01601`, 1,078 MAGs): mean mapped-read share **22.93%** (`cbb_community_share_fraction`).
       - `rbcL + prkB` (`K01601 + K00855`, 416 MAGs): mean mapped-read share **13.21%** (`cbb_rbcl_prk_community_share_fraction`).
       - Strict 3-gene CBB (`rbcL + rbcS + prkB`, `K01601 + K01602 + K00855`, 315 MAGs): mean mapped-read share **11.23%** (`cbb_strict_3gene_community_share_fraction`), matching the candidate screen in `empty-quarter-lifeonair`.
     - Removed the artificial `0.15` lower clip in `carbon_budget_engine.py` (`np.clip(raw_cond, 0.0, 1.0)`), so states with zero CBB co-carriage yield strictly zero direct RuBisCO $\mathrm{CO}_2$ fixation (`0.0`), and reported both inclusive and strict 3-gene (`direct_rubisco_strict_3gene_co2_t_km2_yr`) fixation outputs.
  2. **Manuscript (`paper/main.tex` & `paper/supplementary.tex` Table S6a):**
     - Reported all three marker-screen shares (`22.93%`, `13.21%`, `11.23%`) and both conditional Group 1l CBB co-carriage shares (`33.35%` inclusive, `21.94%` strict 3-gene across all 155 states; `37.67%` and `24.95%` in subsurface states) in the main text and Supplementary Table S6a.

### A.4. External Plausibility Comparisons Across Deserts
- **Critique:** Calling seven interval checks "external validation benchmarks" overstated independent validation where several ranges overlapped literature priors or required normalization/assay harmonization (Tribbia et al. 2026 Namib rates per 16S copy in dry vs. wet soils; Leung et al. 2025 dark fixation microcosms; Ma et al. 2021 *Haloxylon* stands).
- **Changes Made:**
  - Reframed the cross-desert comparisons in `paper/main.tex` and `paper/supplementary.tex` (Table S6b) as **cross-desert plausibility comparisons**, explicitly documenting the normalization assumptions in table footnotes ($10^6\ \mathrm{cells\ gdw}^{-1}$ and $18.4\%$ trace-gas guild share for Namib cell-specific oxidation; $15\ \mathrm{d\ yr^{-1}}$ 0–2 cm wetted horizon vs. 365-day dry column average $3.97 \times 10^{-4}\ \mathrm{nmol\ CO_2\ gdw^{-1}\ d^{-1}}$ for microcosm dark fixation; and $0.45\ \mathrm{mg\ H_2\ m^{-2}\ yr^{-1}}$ uncalibrated thermodynamic scaling vs. $8.72\ \mathrm{mg\ H_2\ m^{-2}\ yr^{-1}}$ diffusion-coupled uptake).

### A.5. Retrospective Temporal Holdout (Trip 5) & Strict Prevention of Information Leakage
- **Critique:** The Trip 5 validation merged pH, DNA, and XRF from the all-campaign state cube (allowing post-cutoff Trip 5 values to leak into Trips 1–4 features) and computed training-site historical features across all training trips rather than strictly prior trips ($t' < t$). Moreover, without a pre-registered forecast prior to Trip 5 collection, it should be termed a **retrospective temporal holdout**.
- **Changes Made:**
  1. **Code (`src/eq_biodt/validation/twin_validator.py`):**
     - Eliminated the all-campaign cube merge for temporal holdouts; built a strictly cutoff-isolated pre-Trip-5 covariate lookup (`pH`, `XRF`, `DNA`) exclusively from Trips 1–4 records (`trip_num <= 4`).
     - Computed longitudinal training history features strictly causally across prior trips ($t' < t$) so no training observation ever accesses future campaign measurements.
     - Added the **546 primary-site analytical cohort** evaluation (`primary_546_cohort`, excluding the 4 auxiliary-site WGS libraries `62PRr2`, `64Dr3`, `64PRr2`, `64Sr2` and `F5Dr1`) alongside the 550-library benchmark (`ρ = 0.525`, $R^2 = 0.251$ on the 546 primary-site cohort; auxiliary-site transfer `n = 4`, `ρ = 0.800`).
  2. **Manuscript (`paper/main.tex` & `paper/supplementary.tex`):**
     - Replaced all occurrences of "prospective forecasting" for Trip 5 with **retrospective temporal holdout** and updated Table 2 and the text with the cutoff-isolated metrics (longitudinal ridge $\rho = 0.619$, RMSE $= 1.849$; environmental ridge $\rho = 0.619$, RMSE $= 1.862$; random forest $\rho = 0.625$, RMSE $= 1.816$; gradient boosting $\rho = 0.568$, RMSE $= 1.838$; within-site change $\Delta$ $\rho = 0.458$, $R^2 = 0.044$).

### A.6. Reconciliation Framing & Sibling-Effect Uncertainty
- **Critique:** Frame ontology reconciliation as traceable compilation and reliable enforcement of measurement constraints rather than claiming Euclidean projection beats covariance-weighted WLS/MinT on leaf MSE; also account for shared paired visits when propagating uncertainty from the three hydrogenase subgroup contrasts to the parent aggregate.
- **Changes Made:**
  - Updated `src/eq_biodt/validation/twin_validator.py` (`_validate_tier4_causal_identifiability`) to compute the parent `[NiFe]`-hydrogenase contrast uncertainty via a **joint site-clustered bootstrap (`B = 500` resamples of sites)** across the shared paired visits (`parent_empirical_ci_low = -0.351`, `parent_empirical_ci_high = -0.175`, capturing positive cross-subgroup covariance rather than assuming independent standard errors).
  - Framed Tier 3 in `paper/main.tex` and `paper/supplementary.tex` as machine-checked compilation of measurement constraints from OWL 2 axioms that guarantees $\mathrm{OCS} = 1.000$ and matches diagonal WLS and shrinkage-MinT within $0.2\%$ RMSE.

### A.7. Observational Field Contrasts vs. Causal Intervention Identification
- **Critique:** Root-vs.-bulk and surface-vs.-subsurface contrasts are observational compartment differences used to calibrate/anchor the scenario engines, not randomized planting or thermal manipulation experiments; extract DNA yield is a proxy for biomass dilution rather than a direct qPCR measurement of absolute gene copies per gram of soil.
- **Changes Made:**
  - Explicitly stated in `paper/main.tex` (Results, Discussion, Methods) that the 1.38-fold absolute gene density increase under 30% shrub cover is a **DNA-yield-calibrated proxy** derived by combining the empirical root-vs.-bulk GE dilution ($\log_2\mathrm{FC} = -2.27$ across 106 paired visits at 47 sites) with the 10.3-fold ($\log_2\mathrm{FC} = +3.37$ across 47 paired DNA visits) extractable DNA ratio, and that direct qPCR/ddPCR and chamber flux assays in planted vs. bare plots are required to verify absolute gene copy densities.

### A.8. Support Mask, Geodesic Distance, and Executable Corridor Decision Atlas
- **Critique:** `build_twin_state_cube.py` used degree-approximated distance (`deg * 111.0`) instead of geodesic distance, and `pareto_optimizer.py` returned support metadata without computing the per-cell portfolio prescriptions across the 128 supported corridor cells; Sites 59/60 were flagged as out-of-support in the cube yet discussed as receiving sabkha prescriptions.
- **Changes Made:**
  1. **Code (`src/eq_biodt/data/build_twin_state_cube.py` & `src/eq_biodt/interventions/pareto_optimizer.py`):**
     - Replaced Euclidean degree scaling with exact **spherical Haversine geodesic distance** (`_haversine_min_dist_km`) in both `build_twin_state_cube.py` and `pareto_optimizer.py`.
     - Implemented `prescribe_corridor_atlas()` in `src/eq_biodt/interventions/pareto_optimizer.py`, which evaluates every one of the 248 regular $0.25^\circ$ grid cells: it abstains on all 120 out-of-support cells (`abstention_unmodeled_or_unsupported`, including hyper-saline sabkha salt-flat cells matching Sites 59/60) and evaluates the six portfolios at each of the **128 empirically supported corridor cells** using its nearest primary site's observed state and local climate, producing the exact reproducible assignment of **87 P2 (Integrated Eco-Restoration), 34 P3 (Low-Water Crust & Refugium Conservation), and 7 P4 (High-Canopy Ca-SO$_4$ Oasis)** cells.
  2. **Manuscript (`paper/main.tex` & `paper/supplementary.tex`):**
     - Updated all text and captions with the exact executable atlas counts (**87 P2, 34 P3, 7 P4** across the 128 supported cells; **120 abstained** cells) and clarified that for Sites 59 and 60 (where surface and subsurface shotgun metagenomes were not recovered due to halite inhibition), the twin abstains from spatial surface/subsurface prediction and uses the Site 59 brine chemistry solely to parameterize the mechanistic osmotic shock boundary in Engine 3.

### A.9 & A.10. Preference-Conditioned Portfolio Ranking, P5 Rank-1 Share, and Physical Water Volumes
- **Critique:** Report P5's 13.5% rank-1 frequency at Site 30 rather than omitting it from the narrative; translate dimensionless rainfall pulses into physical water volumes across the 92,160 $\mathrm{km}^2$ corridor ($553\times 10^6\ \mathrm{m}^3$ for a 6 mm pulse).
- **Changes Made:**
  - Explicitly reported P5's **13.5% rank-1 probability** alongside P3 (42.8%), P2 (24.0%), and P4 (19.6%) in `paper/main.tex` (Results & Figure 3 caption), noting that all four non-dominated portfolios (`P2, P3, P4, P5`) achieve non-zero first-place shares across the Dirichlet weight simplex and parameter perturbations.
  - Added the physical water volume of a 6 mm pulse across the 92,160 $\mathrm{km}^2$ corridor (**$553 \times 10^6\ \mathrm{m}^3$ of precipitation**) in `paper/main.tex` to ground the distinction between atmospheric cloud-seeding pulse capture and supplied desalination irrigation water.

### A.11. Complete Alignment of Figures 1, 2, and 3 with Their Captions
- **Critique:**
  - *Figure 1d:* Caption described MAG-vs.-16S phylum shares, whereas `generate_figures.py` plotted 16S shares plus an "Unmapped reads" bar.
  - *Figure 2:* Panel letters (b, c, d) in the PDF did not match the caption order, and Panel c showed `n = 554` instead of `n = 550`.
  - *Figure 3b/c/d:* Panel 3b only plotted one inoculation curve instead of inoculated (`RAK-PR = 1`) vs. uninoculated (`RAK-PR = 0`), and the 4-panel figure omitted the spatial corridor decision atlas described in the caption.
- **Changes Made (`paper/generate_figures.py` & `paper/figures/*.pdf`):**
  1. **Figure 1d (`fig1_transect_cube_architecture.pdf`):** Updated `generate_figure_1()` to plot side-by-side grouped bars of **MAG CoverM mapped-read shares ($n = 155$ observed states)** and **16S rRNA ASV shares ($n = 180$ states across all 1,227 primary-site profiles)** for Actinomycetota, Thermoplasmatota, and Halobacteriota across Surface, Subsurface, and Root-adjacent compartments, matching the caption and Table 1.
  2. **Figure 2 (`fig2_validation_causal_polytope.pdf`):** Synchronized panel order and labels (`n = 550` metagenomes in Panel c; exact cutoff-isolated metrics from `validation_report.json`) with the Figure 2 caption in `paper/main.tex`.
  3. **Figure 3 (`fig3_interventions_pareto.pdf`):**
     - **Panel 3b:** Plotted both **inoculated (`RAK-PR = 1`, solid lines)** and **uninoculated (`RAK-PR = 0`, dashed lines)** trajectories for area-weighted $\log_2$ GE, biomass-weighted $\log_2$ GE, and absolute Group 1l density multiplier (`1.38x` inoculated vs. `1.10x` uninoculated at 30% shrub cover).
     - **Panel 3d:** Plotted the **248-cell spatial decision atlas** showing the 120 out-of-support abstained cells and the 128 empirically supported corridor cells coloured by their prescribed portfolio (**87 P2, 34 P3, 7 P4**), alongside the Site 30 Pareto front and rank-1 probability inset (P3 42.8%, P2 24.0%, P4 19.6%, P5 13.5%).

---

## Section B: Cross-Paper Numerical Reconciliation (`empty-quarter-lifeonair`, `empty-quarter-ecology-current`, `empty-quarter-data-current`)

1. **16S Amplicon Cohort Join Defect (`1,212 -> 1,227` primary-site profiles; `1,237` total QC-passed profiles):**
   - **Root Cause Found:** `load_16s_phylum_composition()` in `src/eq_biodt/data/manifest_registry.py` previously merged on raw sample names without the canonical preparation-suffix join table (`profile_extraction_join_audit.tsv`) and then intersected `in_qc_1237` with specimen-level `sample_ledger.tsv` `in_ecology_analysis`, dropping 4 valid primary-site profiles (`V16Dr1O`, `V16Dr1T`, `V16Sr1O`, `V7Dr1T`).
   - **Fix Verified:** Updated `load_16s_phylum_composition()` to use `profile_extraction_join_audit.tsv` (`in_qc_1237 == True`) directly, recovering all **1,237 QC-passed profiles** (**1,227 primary-site profiles**: **372 Surface, 381 ShallowSubsurface, 474 RootAdjacent**, plus 10 auxiliary-site profiles) and excluding the low-depth profile `e681_T1Dr1` (934 reads) and 24 controls. Rebuilt all cubes (`uv run eq-biodt build-cubes`) and verified in `tests/test_eq_biodt.py` (`1227 == 1227`).
2. **WGS Library Denominators (`573 -> 551 -> 550 -> 546`):**
   - Documented the exact nested hierarchy in `paper/main.tex` (Results & Methods) and `paper/supplementary.tex`: **573** total sequenced shotgun libraries (`551` field libraries + `22` laboratory controls) $\rightarrow$ **551** curated field metagenomes across 61 sites (`117` surface, `166` subsurface, `268` root-associated) $\rightarrow$ **550** libraries with finite single-copy marker coverage (`F5Dr1` has zero SCG coverage) $\rightarrow$ **546** primary-site analytical libraries in `empty-quarter-lifeonair` after excluding the 4 auxiliary-site libraries (`62PRr2`, `64Dr3`, `64PRr2`, `64Sr2`, which leaves `116` surface, `164` subsurface, `266` root-associated across 59 primary sites).
3. **Hydrogenase Exclusive vs. Inclusive Carrier Counts (`558/371` exclusive vs. `597/410` inclusive):**
   - Updated `paper/main.tex` and `paper/supplementary.tex` to state both the primary-group assignment counts (**558** Group 1l, **371** Group 1h) and the inclusive counts including the **39 dual Group 1l + 1h carriers** (**597** inclusive Group 1l, **410** inclusive Group 1h), along with the high-completeness ($\ge 90\%$) Group 1l CBB co-carriage rate (**76 / 275 = 27.64%**).
4. **72-MAG Pilot Flux Balance Wording & Unused Macros:**
   - Corrected `paper/supplementary.tex` to state that **1 of 48 trace-gas carrier GEMs** had native gene-protein-reaction (GPR) wiring for hydrogenase/CO dehydrogenase prior to stoichiometry repair (replacing "1 of 48 keywords").
   - Removed the unused legacy macros (`636 specimens`, `35,342 ASVs`, `56 phyla`, `372 MAGs`) from `paper/main.tex`.
5. **Darwin Core / eMoF Measurement Semantics & B-Cubed Occurrence Interpretation:**
   - Updated `src/eq_biodt/standards/dwc_emof_exporter.py` and `paper/supplementary.tex` so that:
     - Trace-gas guild capacities use `copies_per_genome_equivalent` (`UO:0000189`) rather than the contradictory `genome_equivalents_per_cell`.
     - Annual precipitation uses `mm/year` (`UO:0000016` with per-year rate qualifier).
     - XRF measurements record `relative_intensity_or_weight_fraction_uncalibrated` (`UO:0000186`), preserving the data descriptor's caution on uncalibrated field/laboratory channels.
     - Zero sequence-based estimates are labelled `below_metagenomic_detection_limit` rather than biological absence, and guild records carry `measurement_aggregate (guild_capacity)` in `taxonRank`.

---

## Section D: Additional Code and Method Defects Fixed and Tested

1. **Reject Incomplete or `NaN` Inputs Before Projection / Constitutive Evaluation (`D.1`):**
   - Updated `src/eq_biodt/models/causal_ontological_twin.py` (`check_completeness`, `evaluate_constitutive_layer`, and `OntologicalPolytope.project`) to inspect raw inputs for missing keys or `NaN`/`inf` entries *before* defaulting or clipping. Incomplete or non-finite states return `is_complete = False`, `status = "incomplete_or_nan_input"`, and `OCS = 0.0` (or raise `ValueError` in `project(..., allow_incomplete=False)`).
2. **Certified Metric Projection & Explicit Solver Status (`D.2` & `D.3`):**
   - Updated `OntologicalPolytope.project()` to validate weighting-matrix dimensions and positive-definiteness and to expose `last_projection_status` (`"certified_slsqp"` vs. `"uncertified_fallback"`).
3. **Metabolic Dormancy Fallback & Maintenance-Only LP (`D.4`, `D.5`, `D.6`):**
   - Updated `src/eq_biodt/models/metabolic_community_engine.py` (`solve_guild_fba`) so that when full growth + maintenance is infeasible (e.g., high salinity or zero donors), it solves a secondary **maintenance-only LP** that consumes only available donors (such as DOC or trace gases) up to their actual bounds (`atp_maint_supplied <= actual_atp_from_donors`), eliminating phantom ATP creation when `H2 = CO = 0` and `DOC = 0.1`. Preserved native HiGHS solver status (`highs_status`) separately from biological dormancy (`maintenance_deficit_arrest` vs. `solver_numerical_failure`).
4. **Coupling of Regional Carbon Budget to the HiGHS Metabolic Engine (`D.7`):**
   - Updated `src/eq_biodt/models/carbon_budget_engine.py` so `estimate_state_carbon_budget()` invokes `MetabolicCommunityEngine.solve_guild_fba()` directly at the state's electrical conductivity and trace-gas bounds rather than using a disconnected `EC < 150` threshold.
5. **Rename `dark_proteome_active_fraction` to `dark_proteome_coding_fraction` (`D.12`):**
   - Updated `src/eq_biodt/models/spatial_jsdm_virome.py`, `src/eq_biodt/interventions/pareto_optimizer.py`, `paper/main.tex`, and `paper/supplementary.tex` to report the **uncharacterized (dark) proteome coding fraction**, removing any claim of measured protein expression or activity.
6. **Resolved SHA-256 Digests and Strict JSON in FDO / RO-Crate Provenance (`D.15` & `D.16`):**
   - Updated `src/eq_biodt/standards/fdo_provenance.py` (`resolve_canonical_dataset_digests`) to compute actual SHA-256 digests of the state cube, grid cube, validation report, and OWL constraint module, and sanitized all `NaN`/`inf` tokens to strict JSON `null` (`allow_nan=False`).
7. **Regression Test Suite (`D.14`):**
   - Added dedicated unit tests in `tests/test_eq_biodt.py` (`test_astra_audit_fixes_completeness_fba_and_provenance`) covering missing/`NaN` state rejection, maintenance-only LP flux conservation under zero trace-gas donors, 3-screen CBB co-carriage (`rbcL`, `rbcL+prkB`, `rbcL+rbcS+prkB`), zero-CBB state carbon fixation (`== 0.0`), geodesic corridor atlas counts (`128` supported = `87 P2 + 34 P3 + 7 P4`, `120` abstained), and SHA-256 FDO provenance digests. All **8/8 test suites pass** (`uv run pytest -v`).
