# Empty Quarter Biodiversity Digital Twin (`EQ-BioDT`)

A standards-compliant, causally grounded **Biodiversity Digital Twin for the Rub' al-Khali (Empty Quarter)** integrating:
1. **Empirical Multi-Omic & Geochemical Atlas**: 60 primary + 10 auxiliary sites across a 1,043 km transect, 5 field campaigns (2023–2025), and 3 soil compartments (`Surface` 0 cm, `ShallowSubsurface` 5–10 cm, `RootAdjacent` 5–10 cm), combining 1,237 16S rRNA profiles, 551 shotgun soil metagenomes (5,438 species-level MAGs, 237 subsurface `Thermoplasmatota` species, 106,042 recurrent dark protein families), trace-gas scavenging guilds (`group 1l`/`1h`/`2a` [NiFe]-hydrogenases, form-I `CoxL`, `RuBisCO`), 19,300 viral OTUs, 322,755 plasmid OTUs, 722 long-read cultured genomes (`RAK-PR`), 210 aquatic + 316 soil + 61 *Bacillus* genome-scale metabolic models (GSMMs), 712 CaCl2 pH records, and 79 ChEBI-mapped XRF analytes.
2. **Two-Layer Causal-Ontological Engine (`OCM`)**: Bridges the Rub' al-Khali OWL Knowledge Graph (`rubalkhali.owl`, 45.7M triples) with the measure-theoretic Two-Layer Causal-Ontological framework (`causal-knowledge-graphs` / `onto-causal`), separating constitutive ontological relations ($\Phi_{\mathrm{const}}$) from causal mechanisms ($\mathcal{M}_{\mathrm{caus}}$) with regime-universal ontological consistency ($\text{OCS} = 1.0$).
3. **European Biodiversity Digital Twin Standards**: Built on locally archived specifications from **BioDT** (FAIR Digital Objects, RO-Crate v1.3, ODMAP model cards, Workflow Run RO-Crate v0.4), **B-Cubed (`b3verse`)** (Uber H3 + 0.25° Open-Meteo/ERA5 occurrence & functional cubes), **LTER-LIFE (`FDFDT`)** (Darwin Core Archive + OBIS `eMoF`), and **DESDEO**-style multi-objective Pareto optimization.

## Environmental Management Simulators
1. **Cloud Seeding & Artificial Precipitation Pulses** (`cloud_seeding`)
2. **Revegetation & Native `RAK-PR` Bio-Inoculation** (`revegetation`)
3. **Desalination, Brine Disposal & Oasis Hydrology** (`desalination_brine`)
4. **Subsurface Refugium & Thermal/UV Shielding** (`refugium_thermal`)
5. **Megafauna Reintroduction & Grazing Management (Arabian Oryx, Sand Gazelle & Camel Stocking)** (`megafauna_grazing`)
6. **Seismic Exploration (Vibroseis) & Hydrocarbon Bioremediation** (`seismic_hydrocarbon`)

## Quick Start

```bash
uv run eq-biodt build-cubes
uv run eq-biodt export-dwca
uv run eq-biodt simulate --site 60 --intervention desalination_brine --ec-mscm 220
uv run eq-biodt optimize --site 25
uv run pytest
```
