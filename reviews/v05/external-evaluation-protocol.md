# Retrospective external evaluation protocol

Written 2026-10-29 after reading the studies' published conclusions and methods,
before fitting their numerical outcomes. This is not preregistration or a blind
validation. Retain all eligible outcomes, including failures.

## Scientific target

Test which mechanisms in the Rub' al-Khali twin transfer: water limitation of
hydrogen oxidation, competing loss of gas diffusion on wetting, and the need to
separate genomic capacity from measured activity. No external result by itself
calibrates absolute Rub' al-Khali field flux. Independent source studies, rather
than individual bottles, are the main units of geographic replication.

## Studies and roles

* Reji et al., doi:10.1038/s41467-025-67208-3: existing temperate forest,
  meadow and sand bottle data. Development source for the -85 MPa threshold;
  cannot independently validate that threshold. Leave-soil-out response tests
  assess transfer between these environments, separately from desert tests.
* Jordaan et al., doi:10.1128/mSystems.01131-20: Australian, Namib, Gobi and
  Mojave desert soils. Test hydration-response direction where numeric raw
  outcomes are obtainable. Genomic marker comparisons use matching estimands;
  do not substitute marker presence for complete pathways or realized flux.
* Buzzard et al., doi:10.7717/peerj.12966: Sonoran Desert urban soil incubations
  and field measurements. Analyse controlled hydration separately from seasonal
  associations. Four treatment plots at one residence do not constitute four
  independent desert sites. Field calibration/test splits hold out entire
  campaigns; matched baselines receive the same training observations.
* Tribbia et al., doi:10.1128/aem.00265-26: ten Namib transect sites, dry and
  hydrated bottles. Test hydration contrasts and observed dry activity. Historic
  site water potentials are not contemporaneous bottle matric potentials;
  do not treat them as direct threshold validation or convert climate-grid soil
  moisture to bottle water content.
* Leung et al., doi:10.1101/2025.11.29.691269: candidate global-desert comparison;
  eligibility depends on retrieving primary data and compatible measurements.

## Comparisons and reporting

1. Preserve native outcome units, signs, replicate identifiers, intervention
   versions and measurement context. Record every exclusion and missing input.
2. Evaluate the frozen v0.4 activation-diffusion mechanism against diffusion-only
   and no-change/constant baselines. Separate default-parameter predictions,
   input-only physical sensitivity ranges, and calibrated predictions.
3. For controlled wetting, report site-level effect direction and magnitude
   with uncertainty when replicates exist. Keep nondetects/censoring distinct
   from zero. Unpaired bottles must not be analysed as paired replicates.
4. For quantitative held-out predictions, fit amplitudes/other parameters only
   on training observations. Publish row-level predictions and folds; compare
   RMSE/MAE to constant and simple empirical alternatives. No pooling of
   incomparable bottle, biomass-normalized and areal flux values.
5. Study-/site-level uncertainty respects dependence. For one-site Sonoran
   data, report descriptive campaign-held-out errors without claiming global
   uncertainty from many within-plot replicates.
6. Any additional model suggested by these outcomes is exploratory and requires
   a distinct evaluation split. Do not silently overwrite the frozen comparator.
7. Bundle redistributable numeric inputs, retrieval URLs, licenses and hashes;
   make analysis run offline after environment installation.

## Ontological test

A measurement must describe its observed property, unit, specimen, spatial/
temporal support and intervention context. Admit a measurement to a mechanism
only through a declared compatible observation map. Compare with a flat model
that has identical information (expected numerical equivalence) and with
explicit semantic ablations (capacity-as-activity, mixed measurement supports,
unsupported route completion). Benefit means changed admissibility or bounded
interpretation, not automatic superiority of the numerical solver.
