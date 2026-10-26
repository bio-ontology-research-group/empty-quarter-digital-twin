# EQ-BioDT Manuscript Reproducibility Package (`reproducibility/`)

This directory contains a self-contained reproducibility pipeline that regenerates model outputs, tables, figures and the web bundle, and checks **205 selected numerical fixtures** for:

> **A causal-ontological biodiversity digital twin of the Rub' al-Khali desert predicts trace-gas scavenger capacity across sites and campaigns and quantifies chemosynthetic and rhizosphere trade-offs of land management**  
> Robert Hoehndorf and Magnus Rueping

It operates in two modes:
1. **Standalone Repository Mode (`--standalone`, default on external clones)**: Uses the pre-materialized Parquet state cubes in `data/processed/` and the bundled reference tables in `reproducibility/data/` with **zero external file or server dependencies**.
2. **Full Multi-Repository Mode**: When executed on the KAUST workstation with sibling repositories (`empty-quarter`, `empty-quarter-ecology-data`, `empty-quarter-lifeonair`, `empty-quarter-workbench`, `empty-quarter-mobilome`), verifies all 23 upstream canonical files via SHA-256 and rebuilds the state cubes from raw tables.

---

## Quick Start (One-Command Full Reproduction)

From the root of the repository (`empty-quarter-digital-twin/`):

```bash
# Install dependencies via uv
uv sync

# Run the complete end-to-end reproducibility pipeline + 205-assertion numerical audit
uv run python reproducibility/reproduce_all.py

# Or explicitly force standalone mode (no external paths outside this repo)
uv run python reproducibility/reproduce_all.py --standalone

# Run the standalone numerical assertion verifier alone
uv run python reproducibility/verify_all_paper_numbers.py --standalone

# Run the full pytest verification suite
uv run pytest -q
```

---

## What `reproduce_all.py` Executes and Verifies

| Pipeline Stage | Script / Module | Regenerated Outputs | Manuscript Sections, Tables & Figures Verified |
| :--- | :--- | :--- | :--- |
| **Stage 1: State Cube & Five-Tier Validation Suite** | [`src/eq_biodt/validation/twin_validator.py`](../src/eq_biodt/validation/twin_validator.py) | `data/validation/validation_report.json`<br>`data/validation/validation_metrics_summary.tsv` | **Section 2.1, Table 1** (70 sites, 180 states, 155 observed WGS states, 551 WGS samples, 1,227 primary-site 16S profiles, 712 admitted pH measurements, 33 XRF analytes, 248 corridor grid cells)<br>**Section 3.1, Table 2, Fig. 2a–d, Suppl. Section S1** (Tier 1 5-fold spatial block CV; Tier 2 LOCO CV, Space-for-Time transport on $P_{90}^{\mathrm{ERA5}}$ & $P_{365}^{\mathrm{ERA5}}$, and Trip 5 retrospective temporal holdout; Tier 3 OWL 2 DL polytope projection & Nonnegative Shrinkage-MinT; Tier 4 paired root-vs-bulk & surface-vs-subsurface contrasts, 97.4% genus-turnover decomposition, and Group 1l CBB co-carriage; Tier 5 synthetic read-to-guild benchmark & 72-MAG 10,584-run GSMM audit) |
| **Stage 2: Six Biophysical Interventions & DESDEO Pareto Optimization** | [`src/eq_biodt/interventions/`](../src/eq_biodt/interventions/)<br>[`reproducibility/reproduce_all.py`](reproduce_all.py) | `data/validation/all_scenarios_summary.json`<br>`paper/tables/table_s3_scores.tex` | **Section 3.2, Table 3, Fig. 3a–c, Suppl. Tables S2–S4** (Cloud seeding 6 mm vs 25 mm at S30; *Haloxylon* revegetation + RAK-PR / `rbcL-hyaB` inoculation at S30; $4 \times 2 \times 2$ factorial desalination brine audit at S60/S59; Subsurface thermal refugium & choline/betaine at S30; Oryx/gazelle rewilding vs camel overgrazing at S30; Nodal seismic + `alkB`/`ladA` bioremediation vs graded Vibroseis at S30)<br>**Section 3.3, Fig. 3d, Suppl. Table S5** (Portfolios $P_1$–$P_6$, $B=200$ Monte Carlo dominance & rank-1 probabilities, and 128-cell supported corridor allocation: 122 $P_3$ + 6 $P_4$, with 120 abstained cells) |
| **Stage 3: Four-Pathway Dimensional Carbon Budget & Cross-Desert Context** | [`src/eq_biodt/models/carbon_budget_engine.py`](../src/eq_biodt/models/carbon_budget_engine.py) | `data/validation/validation_report.json` (`carbon_budget_cross_desert_validation`) | **Section 3.4, Table 4, Suppl. Section S4, Suppl. Table S6** (Microbial dark RuBisCO $\mathrm{CO_2}$ fixation, indirect $\mathrm{H_2}$ soil sink $\mathrm{CO_2}$-eq via $\mathrm{GWP}_{100}$, photosynthetic shrub NPP, and pedogenic carbonate precipitation across the $92{,}160\text{ km}^2$ supported corridor, $500{,}000\text{ km}^2$ sand-sea polygon, and $650{,}000\text{ km}^2$ basin; 7 configured literature-context intervals, partly shared with calibration) |
| **Stage 4: 12 Focal Bacterial Genera & Full Arabian Desert WebGL Grid** | [`reproducibility/reproduce_all.py`](reproduce_all.py) (`refit_12_genus_ridge_models`) | `reproducibility/outputs/taxon_12genera_refit_summary.json`<br>`data/processed/taxon_responses.json`<br>`data/processed/arabian_desert_grid_summary.json` | **Section 2.6, Suppl. Table S9** (7 Hallmark desert genera: *Domibacillus*, *Massilia*, *Bacillus*, *Microvirga*, *Flavisolibacter*, *Noviherbaspirillum*, *CSP1-2*; 5 root-adjacent indicator genera: *TM7a*, *Pelagibacterium*, *Arsenicitalea*, *Ohtaekwangia*, *Neorhizobium*; site-mean CLR ridge regression on $\mathrm{MAT}$, $\mathrm{MAP}/12$, and XRF $\mathrm{S}$ with 200 site bootstraps and 5 contiguous longitude CV blocks; Arabian Desert $14^\circ\text{--}32^\circ\text{N}$, $36^\circ\text{--}60^\circ\text{E}$ multi-resolution grids at $0.10^\circ$, $0.05^\circ$, and $0.04^\circ$) |
| **Stage 5: Interactive Web Digital Twin Bundle** | [`reproducibility/export_web_twin_bundle.py`](export_web_twin_bundle.py) | `reproducibility/outputs/eq_biodt_bundle.json` | **Section 2.6, Section 4.5** (Complete JSON data payload powering the interactive 2D/3D WebGL Digital Twin at `https://rubalkhali.science/twin` and `https://rubalkhali.science/model`) |
| **Stage 6: Vector PDF Figures & 205-Assertion Numerical Audit** | [`scripts/generate_figures.py`](../scripts/generate_figures.py)<br>[`reproducibility/verify_all_paper_numbers.py`](verify_all_paper_numbers.py) | `paper/figures/fig1_transect_cube_architecture.pdf`<br>`paper/figures/fig2_validation_causal_polytope.pdf`<br>`paper/figures/fig3_interventions_pareto.pdf`<br>`reproducibility/outputs/paper_numbers_audit.json`<br>`reproducibility/outputs/reproducibility_summary.tsv` | **Figures 1–3** and **all 205 exact numerical assertions** across `paper/main.tex` and `paper/supplementary.tex` |

---

## Directory Contents

- `reproduce_all.py` — Master end-to-end reproduction script.
- `verify_all_paper_numbers.py` — Automated assertion verifier checking 205 exact manuscript values.
- `export_web_twin_bundle.py` — Exporter for the `eq_biodt_bundle.json` web twin payload.
- `data/` — Bundled standalone reference tables:
  - `xrf_chemical_mapping.yml`, `biome_codes.yml`, `site_iri_registry.tsv`, `site_aliases.tsv`, `taxonomy_domain_phylum_iris.tsv` (programmatic ontology and site registries; Supplementary Table S8).
  - `aquatic_water_chemistry_s2.tsv` (measured hydrochemistry for Sites 59 and 60; Supplementary Table S4).
  - `taxon_12genera_site_clr_panel.tsv` (180-row site $\times$ compartment CLR + $\mathrm{MAT}$ + $\mathrm{MAP}/12$ + $\mathrm{S}$ panel for the 12 focal bacterial genera; Supplementary Table S9).
  - `arabian_desert_grid_summary.json` (Arabian Desert $14^\circ\text{--}32^\circ\text{N}$, $36^\circ\text{--}60^\circ\text{E}$ raster resolution and mask statistics).
  - `h123_v2/` (paired visit $\log_2$ fold-change tables, lineage-turnover decomposition tables, CBB co-carriage tables, and within/between-site temporal slope tables used by Tiers 2 and 4).
  - `benchmark_r2/` (synthetic metagenome read-to-guild benchmark tables used by Tier 5).
  - `canonical_sources_manifest.json` (byte sizes and SHA-256 checksums of the 23 upstream tables).
- `outputs/` — Generated audit reports (`paper_numbers_audit.json`, `reproducibility_summary.tsv`, `taxon_12genera_refit_summary.json`, `eq_biodt_bundle.json`).

## Scientific scope of the checks

The numerical fixture suite compares executable outputs with constants; it does not parse the manuscript or establish intervention effects. Read the [26 October 2026 audit](../reviews/validation-2026-10-26/report.md) alongside these results. Ontological consistency concerns compiled algebraic constraints, causal scenarios depend on prescribed response functions, and carbon low/high values are sensitivity envelopes. The source-level [coefficient inventory](outputs/model_coefficient_inventory.json) records implementation constants without assigning empirical support.
