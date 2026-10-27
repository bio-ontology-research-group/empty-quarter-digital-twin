# Manuscript integration and ontology experiment design

Completed 2026-10-27. The manuscript is maintained in the separate Overleaf
repository. The ontology experiment in [PROTOCOL.md](PROTOCOL.md) is planned;
only its descriptive feasibility audit has run.

## Manuscript

Robert requested manuscript integration through `claude-personal` with Opus 5.5.
The local CLI completed successfully with model `claude-opus-5-5`, applying
`write-edit-scientific-paper/SKILL.md` and its required Robert scientific style
reference. The primary agent then checked and refined the scientific claims.
No new laboratory or field results were introduced.

The revised manuscript connects the metagenomic screen, campaign prediction,
external moisture-envelope comparison and parameter sensitivity. A new table
distinguishes ontology, missingness, observations, causal assumptions, prediction,
design and accounting. The generated supplement retains all 31 clauses and
23 parameter defaults unchanged. The main Methods identify the forecast target,
the external threshold and the proposed experiment's withheld condition.

The evidence check retained the negative comparison: compartment means beat
ridge on RMSE in two campaigns, and persistence beats ridge on MAE in campaign 2.
The Discussion links to the planned ontology protocol, explains the fixed-catalogue
backtest and preserves the need for independent local physiological validation.
Missingness remains an observation policy, and the fitted moisture envelope is
distinguished from validation of the full process model. Covariate contributions
are not inferred from aggregate prediction errors alone.

## Soil-label correction

The primary publication identifies `PB_sand` as **Pine Barrens loamy sand** from
the New Jersey oak-pine forest. Its other two reanalysed soils are also temperate;
the Arizona arid soils in the source study are separate. This was verified in
the Methods, “Soil collection, characterization, and experimental setup,” on
page 8 of the [published article](https://repository.library.noaa.gov/view/noaa/72500/noaa_72500_DS1.pdf)
([DOI](https://doi.org/10.1038/s41467-025-67208-3)).

The inherited “Desert sand” figure label and manuscript wording were corrected.
The generator now maps labels by source soil identifier, and the figure was
regenerated from the unchanged stored fits. No observations, estimates or
comparison results changed. This correction narrows the geographical description
of the external test and is relevant to transfer to Rub' al-Khali soils.

## Experiment and validation

The protocol tests same-MAG pathway witnesses, intersecting carrier sets,
measurement compatibility and admission to conditional resource optimization.
It includes equally informed flat computations, reduced-information bounds,
rule ablations, independent source adjudication and explicit failure decisions.
Its feasibility audit preserves input hashes and all sample/MAG denominators.

Checks completed:

- Successful Claude completion with the requested model and no tool denials.
- Manuscript estimates checked against stored v0.3 results; the single-pulse
  uptake ranking checked in all 24 parameter members.
- Repeated feasibility-audit execution produced identical JSON; MAG strata
  sum to 5,438 and usable plus unavailable mapping queries sum to 551.
- Main manuscript and supplement compile with resolved references and no
  overfull boxes; the revised process figure was visually inspected.
- The numerical result JSONs, model code, generated assumptions and generated
  forecast tables remain unchanged. The figure-generator change affects labels.

The independent ontology comparison, its generic rule evaluator and source-review
exercise remain to be implemented and run. Their expected equality to a correct
flat implementation is a prespecified outcome, not an observed result.
