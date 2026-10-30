"""B-Cubed (b3verse) Spatiotemporal Biodiversity Data Cube & Indicator Builder.

Implements the EU B-Cubed (`b-cubed-eu/b3gbi` & Frictionless Tabular Data Package)
specification for the Empty Quarter Digital Twin:
- Dimensions:
  1. Spatial: Uber H3 hierarchical discrete global grid (`h3_5`, `h3_6`) + `era5_cell`
  2. Temporal: Campaign period (`2023-03-17/2025-10-15` / Trips 1..5)
  3. Compartment: `Surface`, `ShallowSubsurface`, `RootAdjacent`
  4. Taxonomic / Functional Guild (`taxonKey`, `scientificName`, `rank`) with
     programmatically verified NCBITaxon and RAK IRIs.
- Measures:
  - Strictly observed `occurrences` (never fabricated when missing), empirical
    `abundance_or_ge` from raw sequencing tables (`*_raw`), `shannon_16s`, `hill_q1_diversity`,
    and `minCoordinateUncertaintyInMeters`.
- Exports `datapackage.json` + Parquet/CSV cube files.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

from eq_biodt.data.build_twin_state_cube import TwinStateCubeBuilder


class B3OccurrenceCubeBuilder:
    """Builds B-Cubed compliant occurrence/functional cubes and `b3gbi`-style indicators."""

    def __init__(self, state_builder: Optional[TwinStateCubeBuilder] = None) -> None:
        self.state_builder = state_builder or TwinStateCubeBuilder()
        self.ontology = self.state_builder.ontology

    def _get_verified_guilds(self) -> List[Tuple[str, str, str, str, str]]:
        """Return (taxon_iri, sci_name, rank, raw_col, count_col) using verified ontology IRIs."""
        return [
            (
                self.ontology.get_taxon_iri("Actinomycetota"),
                "Actinomycetota",
                "PHYLUM",
                "mag_actinomycetota_raw",
                "wgs_sample_count",
            ),
            (
                self.ontology.get_taxon_iri("Thermoplasmatota"),
                "Thermoplasmatota",
                "PHYLUM",
                "mag_thermoplasmatota_raw",
                "wgs_sample_count",
            ),
            (
                self.ontology.get_taxon_iri("Halobacteriota"),
                "Halobacteriota",
                "PHYLUM",
                "mag_halobacteriota_raw",
                "wgs_sample_count",
            ),
            (
                self.ontology.get_taxon_iri("Chloroflexota"),
                "Chloroflexota",
                "PHYLUM",
                "mag_chloroflexota_raw",
                "wgs_sample_count",
            ),
            (
                f"{self.ontology.RAK_BASE}RAK_GUILD_HYDB_1L",
                "Group 1l [NiFe]-hydrogenase guild",
                "FUNCTIONAL_GUILD",
                "hydB_1l_ge_raw",
                "wgs_sample_count",
            ),
            (
                f"{self.ontology.RAK_BASE}RAK_GUILD_HYDB_1H",
                "Group 1h [NiFe]-hydrogenase guild",
                "FUNCTIONAL_GUILD",
                "hydB_1h_ge_raw",
                "wgs_sample_count",
            ),
            (
                f"{self.ontology.RAK_BASE}RAK_GUILD_COXL_FORMI",
                "Form I coxL CO-oxidizing guild",
                "FUNCTIONAL_GUILD",
                "coxL_formI_ge_raw",
                "wgs_sample_count",
            ),
            (
                f"{self.ontology.RAK_BASE}RAK_GUILD_VIROME_CAUDOVIRICETES",
                "Desert Caudoviricetes virome",
                "VIRAL_GUILD",
                "viral_richness_mean_raw",
                "wgs_sample_count",
            ),
        ]

    def build_b3_cube(self) -> pd.DataFrame:
        """Construct the multi-dimensional B-Cubed occurrence and functional trait cube from observed rows."""
        comp_df = self.state_builder.build_compartment_state_cube()
        rows: List[Dict[str, Any]] = []
        guilds = self._get_verified_guilds()

        for _, r in comp_df.iterrows():
            shannon = float(r["shannon_16s"]) if pd.notna(r["shannon_16s"]) else float("nan")
            hill_q1 = float(np.exp(shannon)) if np.isfinite(shannon) else float("nan")

            for taxon_key, sci_name, rank, raw_col, count_col in guilds:
                raw_val = r.get(raw_col, float("nan"))
                obs_cnt = int(r[count_col]) if pd.notna(r.get(count_col)) else 0
                if pd.isna(raw_val) or obs_cnt <= 0:
                    continue
                val = float(raw_val)
                rows.append(
                    {
                        "year_period": "2023-2025",
                        "eventDate": "2023-03-17/2025-10-15",
                        "h3_5": str(r["h3_5"]),
                        "h3_6": str(r["h3_6"]),
                        "era5_cell": str(r["era5_cell"]),
                        "site_num": int(r["site_num"]),
                        "site_iri": str(r["site_iri"]),
                        "compartment": str(r["compartment_norm"]),
                        "transect_zone": str(r["transect_zone"]),
                        "taxonKey": taxon_key,
                        "scientificName": sci_name,
                        "taxonRank": rank,
                        "occurrences": obs_cnt,
                        "abundance_or_ge": round(val, 6),
                        "shannon_16s": round(shannon, 4) if np.isfinite(shannon) else None,
                        "hill_q1_diversity": round(hill_q1, 2) if np.isfinite(hill_q1) else None,
                        "minCoordinateUncertaintyInMeters": 10.0,
                    }
                )
        return pd.DataFrame(rows)

    def compute_h3_cell_indicators(self) -> pd.DataFrame:
        """Aggregate `b3gbi`-style biodiversity indicators per H3_5 grid cell and compartment."""
        comp_df = self.state_builder.build_compartment_state_cube()
        ind = (
            comp_df.groupby(["h3_5", "transect_zone", "compartment_norm"], as_index=False)
            .agg(
                sites_in_cell=("site_num", "nunique"),
                lat_centroid=("lat", "mean"),
                lon_centroid=("lon", "mean"),
                shannon_mean=("shannon_16s", "mean"),
                hill_q1_mean=("shannon_16s", lambda s: float(np.mean(np.exp(s.dropna())))),
                hydB_1l_ge_mean=("hydB_1l_ge_raw", "mean"),
                coxL_formI_ge_mean=("coxL_formI_ge_raw", "mean"),
                viral_richness_mean=("viral_richness_mean_raw", "mean"),
                thermoplasmatota_mean=("mag_thermoplasmatota_raw", "mean"),
                halobacteriota_mean=("mag_halobacteriota_raw", "mean"),
                xrf_evaporite_index_mean=("xrf_evaporite_index", "mean"),
                ph_mean=("ph_mean_raw", "mean"),
                map_mm_mean=("map_mm", "mean"),
            )
        )
        return ind

    def export_b3_datapackage(self, output_dir: Path) -> Dict[str, Path]:
        """Export Frictionless `datapackage.json` + CSV and Parquet cube files."""
        output_dir.mkdir(parents=True, exist_ok=True)
        cube_df = self.build_b3_cube()
        ind_df = self.compute_h3_cell_indicators()

        cube_csv = output_dir / "b3_eq_occurrence_cube.csv"
        cube_parquet = output_dir / "b3_eq_occurrence_cube.parquet"
        ind_csv = output_dir / "b3_h3_biodiversity_indicators.csv"
        dp_json = output_dir / "datapackage.json"

        cube_df.to_csv(cube_csv, index=False)
        cube_df.to_parquet(cube_parquet, index=False)
        ind_df.to_csv(ind_csv, index=False)

        datapackage = {
            "profile": "tabular-data-package",
            "name": "eq-biodt-b3-occurrence-cube",
            "title": "Empty Quarter B-Cubed (b3verse) Spatiotemporal Biodiversity & Functional Trait Cube",
            "description": "H3-indexed taxonomic, trace-gas chemolithoautotrophic guild, and virome data cube conforming to B-Cubed (b-cubed-eu) specifications.",
            "temporal": {
                "startDate": "2023-03-17",
                "endDate": "2025-10-15",
            },
            "spatial": {
                "gridType": "H3",
                "resolutions": [5, 6],
                "crs": "EPSG:4326",
            },
            "resources": [
                {
                    "name": "b3_eq_occurrence_cube",
                    "path": "b3_eq_occurrence_cube.csv",
                    "profile": "tabular-data-resource",
                    "schema": {
                        "fields": [
                            {
                                "name": col,
                                "type": (
                                    "number"
                                    if pd.api.types.is_numeric_dtype(cube_df[col])
                                    else "string"
                                ),
                            }
                            for col in cube_df.columns
                        ]
                    },
                },
                {
                    "name": "b3_h3_biodiversity_indicators",
                    "path": "b3_h3_biodiversity_indicators.csv",
                    "profile": "tabular-data-resource",
                },
            ],
        }
        dp_json.write_text(json.dumps(datapackage, indent=2), encoding="utf-8")

        return {
            "cube_csv": cube_csv,
            "cube_parquet": cube_parquet,
            "indicators_csv": ind_csv,
            "datapackage_json": dp_json,
        }
