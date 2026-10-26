# Empty Quarter Biodiversity Digital Twin (`EQ-BioDT`)

[![Live Digital Twin](https://img.shields.io/badge/Live_Digital_Twin-rubalkhali.science%2Ftwin-0ea5e9)](https://rubalkhali.science/twin)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](https://opensource.org/licenses/MIT)

**EQ-BioDT** is a causal-ontological biodiversity digital twin of the Rub' al-Khali (Empty Quarter) and the wider Arabian Desert ecoregion. It accompanies the manuscript:

> **A causal-ontological biodiversity digital twin of the Rub' al-Khali desert predicts trace-gas scavenger capacity across sites and campaigns and quantifies chemosynthetic and rhizosphere trade-offs of land management**  
> Robert Hoehndorf & Magnus Rueping

The live interactive 2D/3D WebGL Digital Twin is deployed at **[https://rubalkhali.science/twin](https://rubalkhali.science/twin)** (and **[https://rubalkhali.science/model](https://rubalkhali.science/model)**).

---

## Architecture & Scientific Overview

1. **Transect State Cube (`180` site-compartment states across `70` sites)**:
   - Harmonizes `551` shotgun metagenomes (`550` with non-missing trace-gas guild estimates, `549` with mobilome assays, `5,438` dereplicated species-level MAGs), `1,227` quality-controlled primary-site 16S rRNA amplicon profiles (`1,237` across all sites; `1,271` total sequenced profiles), `712` admitted $\mathrm{CaCl}_2$ pH records (`178` observed states), X-ray fluorescence (XRF) geochemistry (`11` core elements + vanadium), and ERA5 / Open-Meteo climate products across `60` primary sites (`S01`–`S60`) and `10` auxiliary sites (`RAK_1000061`–`RAK_1000070`) along a `1,043 km` transect of the Rub' al-Khali across five campaigns (`Trips 1–5`, March 2023 – October 2025).
   - Indexed on Uber H3 hexagonal cells (`h3_5`, `h3_6`) and `0.25°` ERA5 corridor grid cells (`248` cells total: `128` inside the $\le 45\text{ km}$ empirical domain of support covering `92,160 km²`, and `120` abstained cells).
2. **Three-Layer Causal-Ontological Engine**:
   - **Predictive Layer**: Out-of-sample ridge regressions for trace-gas guilds (`hydB_1l`, `hydB_1h`, `hydB_2a`, `coxL_formI`), 16S Shannon diversity ($H'$), viral and plasmid richness (`vOTU`, `pOTU`), and `12` focal bacterial genera (`7` hallmark desert genera: *Domibacillus*, *Massilia*, *Bacillus*, *Microvirga*, *Flavisolibacter*, *Noviherbaspirillum*, *CSP1-2*; and `5` root-adjacent indicator genera: *TM7a*, *Pelagibacterium*, *Arsenicitalea*, *Ohtaekwangia*, *Neorhizobium*).
   - **Constitutive Layer ($\Phi_{\mathrm{const}}$)**: Compiles `13` subsumption (`SubClassOf`), `5` exhaustive partition (`DisjointUnion`), and `8` unit upper-bound axioms over `17` state variables from the OWL 2 DL module [`ontology/rubalkhali_twin_constraints.owl`](ontology/rubalkhali_twin_constraints.owl) into a closed convex polyhedron $\mathcal{C}(\mathcal{T}_{\mathrm{ont}})$ and enforces ontological consistency ($\mathrm{OCS} = 1.0$, $\mathrm{CVM} = 0.0$) via Euclidean projection $\Pi_{\mathcal{C}(\mathcal{T}_{\mathrm{ont}})}$.
   - **Mechanistic Scenario Layer ($\mathcal{M}_{\mathrm{caus}}$)**: Couples Millington–Quirk porous-media gas diffusion, distributed-lag moisture-pulse resuscitation, a curated 10-flux guild-level linear program solved with HiGHS, and spatial virome/mobilome response models across six desert management engines.
3. **Multi-Objective Portfolio Optimization & Dimensional Carbon Accounting**:
   - Evaluates the current hyperarid baseline (`P0`) and six management portfolios (`P1`–`P6`) across four objectives ($J_1$–$J_4$) using DESDEO achievement scalarizing functions and `B = 200` Monte Carlo robustness draws.
   - Upscales four physical carbon and greenhouse-gas pathways (direct microbial RuBisCO dark $\mathrm{CO}_2$ fixation, the indirect tropospheric $\mathrm{OH}^{\bullet}$ / $\mathrm{CH}_4$ sink of soil $\mathrm{H}_2$ oxidation at $\mathrm{GWP}_{100} = 11.6$, *Haloxylon persicum* and rhizosphere net primary production, and pedogenic $\mathrm{CaCO}_3$ precipitation) across the `92,160 km²` empirically supported corridor (`122` `P3` and `6` `P4` cells in the prescribed corridor atlas), the `500,000 km²` core sand-sea polygon, and the `650,000 km²` Rub' al-Khali basin, validated against seven published cross-desert benchmarks.
4. **European Biodiversity Digital Twin Interoperability**:
   - Exports Darwin Core Archives (`event.txt`, `occurrence.txt`, `emof.txt`, `meta.xml`, `eml.xml`), B-Cubed occurrence cubes (`datapackage.json`), and content-addressed FAIR Digital Objects (`FDO`) / Workflow Run RO-Crate v1.1 (`ro-crate-metadata.json`) provenance records.

---

## Repository Structure

```text
empty-quarter-digital-twin/
├── data/
│   ├── processed/                 # Self-contained Parquet state cubes, 551-WGS panel & 12-genus JSON
│   ├── dwca/                      # Darwin Core Archive (Event core + Occurrence + OBIS eMoF)
│   ├── b3_cubes/                  # B-Cubed occurrence cube & H3 indicator tables + Frictionless spec
│   ├── ro_crate/                  # Workflow Run RO-Crate v1.1 & FDO kernel metadata
│   └── validation/                # Five-tier validation report JSON & metrics summary TSV
├── ontology/
│   └── rubalkhali_twin_constraints.owl  # OWL 2 DL constraint module (17 classes, 5 partitions, 13 subsumptions)
├── paper/
│   ├── main.tex                   # Main manuscript LaTeX source
│   ├── supplementary.tex          # Supplementary Information LaTeX source (Notes 1-3, Tables S1-S9)
│   ├── references.bib             # BibTeX bibliography
│   ├── figures/                   # Multi-panel vector PDF figures (fig1, fig2, fig3)
│   └── tables/                    # Auto-generated LaTeX table fragments (table_s3_scores.tex)
├── reproducibility/               # Full standalone reproducibility package for every number in the paper
│   ├── README.md                  # Detailed table-by-table & section-by-section reproduction guide
│   ├── reproduce_all.py           # Master script regenerating all tables, figures, bundles & audits
│   ├── verify_all_paper_numbers.py# Automated assertion suite checking 205 exact paper numbers
│   ├── export_web_twin_bundle.py  # Exporter for the interactive WebGL Digital Twin JSON bundle
│   ├── data/                      # Bundled standalone empirical contrasts, benchmarks & registries
│   └── outputs/                   # Generated reproducibility reports, CSVs & audit logs
├── scripts/
│   ├── generate_figures.py        # Generates vector PDF Figures 1, 2, and 3
│   └── verify_manuscript_numbers.py # Inspects cube denominators and scenario outputs
├── src/eq_biodt/
│   ├── api/                       # FastAPI REST backend & interactive WebGL console (/twin)
│   ├── cli.py                     # Command-line interface (eq-biodt)
│   ├── data/                      # State cube builder & programmatic ontology/manifest registries
│   ├── interventions/             # 6 biophysical scenario engines & DESDEO Pareto portfolio optimizer
│   ├── models/                    # OWL 2 polytope projection, Millington-Quirk, HiGHS FBA, CarbonBudget
│   ├── standards/                 # DwC-A/eMoF, B-Cubed & RO-Crate/FDO exporters
│   └── validation/                # Five-tier quantitative validation suite
└── tests/
    └── test_eq_biodt.py           # Comprehensive pytest suite (unit, API, 5-tier validation & standalone audit)
```

---

## Installation & Quick Start

All Python execution uses [`uv`](https://docs.astral.sh/uv/) with locked dependencies in [`uv.lock`](uv.lock):

```bash
# 1. Install locked environment
uv sync

# 2. Run the full test suite (including standalone reproducibility audit)
uv run pytest -q
```

### Reproducing Every Number, Table, and Figure in the Paper

The [`reproducibility/`](reproducibility/) directory is completely self-contained and runs without external mounts or companion repositories:

```bash
# Run the complete end-to-end reproduction pipeline in standalone mode:
uv run python reproducibility/reproduce_all.py --standalone

# Or run the 205-assertion numerical audit directly:
uv run python reproducibility/verify_all_paper_numbers.py --standalone
```

See [`reproducibility/README.md`](reproducibility/README.md) for the complete mapping between manuscript sections, tables (`Tables 1–4`, `Suppl. Tables S1–S9`), and figures (`Figs. 1–3`) and their generating functions.

### Command-Line Interface (`eq-biodt`)

```bash
# Verify all 23 canonical data and ontology sources
uv run eq-biodt verify-sources

# Build all Parquet state cubes (70 sites, 180 site-compartments, 414 visit-compartments, 248 grid cells)
uv run eq-biodt build-cubes

# Export Darwin Core Archive, B-Cubed occurrence cube, and RO-Crate metadata
uv run eq-biodt export-standards

# Simulate any of the 6 desert interventions at a specific site and compartment
uv run eq-biodt simulate --site 30 --compartment Surface --intervention cloud_seeding --param rain_pulse_mm=6.0
uv run eq-biodt simulate --site 60 --compartment Surface --intervention desalination_brine --param discharge_salinity_ms_cm=45.0 --param ca_so4_to_nacl_ratio=2.4

# Run multi-objective Pareto portfolio optimization (P1-P6) with 200 Monte Carlo draws
uv run eq-biodt optimize --site 30 --compartment Surface

# Run the Five-Tier Quantitative Validation Suite
uv run eq-biodt validate

# Launch the FastAPI REST server and interactive web console at http://localhost:8090/twin
uv run eq-biodt serve --port 8090
```

### Exporting the Interactive WebGL Digital Twin Bundle

To export the self-contained JSON bundle (`eq_biodt_bundle.json`) that powers the interactive 2D/3D WebGL Digital Twin at [https://rubalkhali.science/twin](https://rubalkhali.science/twin):

```bash
uv run python reproducibility/export_web_twin_bundle.py
```

---

## License & Citation

Released under the **MIT License**. When citing **EQ-BioDT** or the Rub' al-Khali Knowledge Graph (`v3.0.1`), please cite:

- Hoehndorf, R. & Rueping, M. *A causal-ontological biodiversity digital twin of the Rub' al-Khali desert predicts trace-gas scavenger capacity across sites and campaigns and quantifies chemosynthetic and rhizosphere trade-offs of land management* (2026).
- Rub' al-Khali Knowledge Graph v3.0.1, Zenodo: [https://doi.org/10.5281/zenodo.23134168](https://doi.org/10.5281/zenodo.23134168).
