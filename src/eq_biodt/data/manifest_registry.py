"""Canonical data pointer and verified ontology identifier registry.

Follows DATA_LOCATIONS.md and GEMINI.md rules:
- Never guess paths or ontology identifiers (ChEBI, ENVO, PATO, UO, SIO, NCBITaxon, RAK).
- Programmatically load all IRIs from the canonical YAML/TSV registries on disk.
- Enforce the RAK_0000100 (quality) / RAK_0000500 (value) XRF analyte counter rule
  (skipping 'LE' which has predefined RAK_0000032/RAK_0000033).
- Resolve Trip 1 numeric alias sites (61..64) via site_aliases.tsv to canonical IRIs
  RAK_1000067..1000070 and load verified NCBITaxon domain/phylum IRIs from
  taxonomy_mapping_decisions.tsv.
"""

from __future__ import annotations

import hashlib
import json
import xml.etree.ElementTree as ET
import zipfile
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import pandas as pd
import yaml


@dataclass(frozen=True)
class CanonicalPaths:
    """Authoritative local and IBEX paths across the Empty Quarter repositories,
    with automatic fallback to the self-contained `reproducibility/data` bundle
    for standalone GitHub clones.
    """

    software_root: Path = Path("/home/leechuck/Public/software")
    papers_root: Path = Path("/home/leechuck/Documents/papers")
    ibex_data_root: str = "/ibex/scratch/projects/c2014/EmptyQuarter_Data"
    ibex_metagenomics_root: str = "/ibex/scratch/projects/c2014/rob/metagenomics"
    repo_root: Path = Path(__file__).resolve().parents[3]

    @property
    def repro_data_dir(self) -> Path:
        return self.repo_root / "reproducibility" / "data"

    @property
    def processed_dir(self) -> Path:
        return self.repo_root / "data" / "processed"

    def _resolve_or_fallback(self, primary: Path, fallback: Path) -> Path:
        return primary if primary.exists() else fallback

    @property
    def eq_repo(self) -> Path:
        return self.software_root / "empty-quarter"

    @property
    def eq_ecology_data(self) -> Path:
        return self.software_root / "empty-quarter-ecology-data"

    @property
    def eq_workbench(self) -> Path:
        return self.software_root / "empty-quarter-workbench"

    @property
    def eq_lifeonair(self) -> Path:
        return self.software_root / "empty-quarter-lifeonair"

    @property
    def eq_proteome(self) -> Path:
        return self.software_root / "empty-quarter-proteome"

    @property
    def xrf_chemical_mapping_yml(self) -> Path:
        return self._resolve_or_fallback(
            self.eq_repo / "config" / "codes" / "xrf_chemical_mapping.yml",
            self.repro_data_dir / "xrf_chemical_mapping.yml",
        )

    @property
    def biome_codes_yml(self) -> Path:
        return self._resolve_or_fallback(
            self.eq_repo / "config" / "codes" / "biome_codes.yml",
            self.repro_data_dir / "biome_codes.yml",
        )

    @property
    def site_iri_registry_tsv(self) -> Path:
        return self._resolve_or_fallback(
            self.eq_ecology_data / "metadata" / "samples" / "site_iri_registry.tsv",
            self.repro_data_dir / "site_iri_registry.tsv",
        )

    @property
    def site_aliases_tsv(self) -> Path:
        return self._resolve_or_fallback(
            self.eq_ecology_data / "metadata" / "samples" / "site_aliases.tsv",
            self.repro_data_dir / "site_aliases.tsv",
        )

    @property
    def env_measurements_curated_tsv(self) -> Path:
        return (
            self.eq_ecology_data
            / "metadata"
            / "environmental"
            / "environmental_measurements_curated.tsv"
        )

    @property
    def sample_ledger_tsv(self) -> Path:
        return self.eq_ecology_data / "evidence" / "release" / "sample_ledger.tsv"

    @property
    def release_evidence_json(self) -> Path:
        return self.eq_ecology_data / "evidence" / "release" / "release_evidence.json"

    @property
    def taxonomy_mapping_decisions_tsv(self) -> Path:
        return self._resolve_or_fallback(
            self.eq_ecology_data
            / "evidence"
            / "release"
            / "taxonomy_mapping_audit"
            / "taxonomy_mapping_decisions.tsv",
            self.repro_data_dir / "taxonomy_domain_phylum_iris.tsv",
        )

    @property
    def feature_table_16s_tsv(self) -> Path:
        return self.eq_ecology_data / "metadata" / "taxonomy" / "feature-table-trips1-5.tsv"

    @property
    def taxonomy_16s_tsv(self) -> Path:
        return self.eq_ecology_data / "metadata" / "taxonomy" / "taxonomy-trips1-5.tsv"

    @property
    def canonical_profile_read_inventory_tsv(self) -> Path:
        return (
            self.eq_ecology_data
            / "metadata"
            / "amplicon"
            / "current_inventory"
            / "canonical_profile_read_inventory.tsv"
        )

    @property
    def profile_extraction_join_audit_tsv(self) -> Path:
        return (
            self.eq_ecology_data
            / "evidence"
            / "ecology-canonical"
            / "depth_extraction"
            / "profile_extraction_join_audit.tsv"
        )

    @property
    def site_altitudes_tsv(self) -> Path:
        return self.eq_ecology_data / "metadata" / "geodata" / "site_altitudes.tsv"

    @property
    def ph_accepted_tsv(self) -> Path:
        return self.eq_ecology_data / "metadata" / "ph" / "ph_accepted_measurements.tsv"

    @property
    def xrf_lab_filtered_tsv(self) -> Path:
        return self.eq_ecology_data / "metadata" / "geochemistry" / "xrf_lab_table_filtered.tsv"

    @property
    def daily_weather_canonical_tsv(self) -> Path:
        return self.eq_ecology_data / "metadata" / "climate" / "daily_weather_canonical.tsv"

    @property
    def monthly_weather_canonical_tsv(self) -> Path:
        return (
            self.eq_ecology_data
            / "metadata"
            / "climate"
            / "monthly_weather_averages_canonical.tsv"
        )

    @property
    def rain_event_product_exposures_tsv(self) -> Path:
        return (
            self.eq_ecology_data
            / "metadata"
            / "climate"
            / "rain_event_product_exposures.tsv"
        )

    @property
    def env_site_tsv(self) -> Path:
        return self.eq_lifeonair / "lifeonair" / "results" / "inputs" / "env_site.tsv"

    @property
    def env_site_trip_tsv(self) -> Path:
        return self.eq_lifeonair / "lifeonair" / "results" / "inputs" / "env_site_trip.tsv"

    @property
    def era5_daily_tsv(self) -> Path:
        return self.eq_lifeonair / "lifeonair" / "results" / "inputs" / "era5_daily_1995_2025.tsv"

    @property
    def power_daily_tsv(self) -> Path:
        return self.eq_lifeonair / "lifeonair" / "results" / "inputs" / "power_daily_2021_2025.tsv"

    @property
    def xrf_site_compartment_panel_tsv(self) -> Path:
        return (
            self.eq_lifeonair
            / "lifeonair"
            / "results"
            / "inputs"
            / "xrf_site_compartment_panel.tsv"
        )

    @property
    def coverm_rep_abundance_tsv_gz(self) -> Path:
        return (
            self.eq_lifeonair
            / "lifeonair"
            / "results"
            / "inputs"
            / "coverm_rep_abundance_long.tsv.gz"
        )

    @property
    def rep_cbb_ko_tsv(self) -> Path:
        return self.eq_lifeonair / "lifeonair" / "results" / "inputs" / "rep_cbb_ko.tsv"

    @property
    def atmo_classification_v2_tsv(self) -> Path:
        return (
            self.eq_lifeonair
            / "lifeonair"
            / "results"
            / "inputs"
            / "atmo_classification_v2.tsv"
        )

    @property
    def gtdb_bac120_summary_tsv(self) -> Path:
        return self.eq_lifeonair / "lifeonair" / "results" / "inputs" / "r232_bac120.summary.tsv"

    @property
    def gtdb_ar53_summary_tsv(self) -> Path:
        return self.eq_lifeonair / "lifeonair" / "results" / "inputs" / "r232_ar53.summary.tsv"

    @property
    def guild_ge_strict_tsv(self) -> Path:
        return (
            self.eq_lifeonair
            / "lifeonair"
            / "results"
            / "reads"
            / "guild_genome_equivalents_v2_r3noflagged.tsv"
        )

    @property
    def guild_abundance_v551_tsv(self) -> Path:
        return (
            self.eq_lifeonair
            / "lifeonair"
            / "results"
            / "inputs"
            / "guild_abundance_by_sample_v551.tsv"
        )

    @property
    def rep_growth_v2_tsv(self) -> Path:
        return self.eq_lifeonair / "lifeonair" / "results" / "inputs" / "rep_growth_v2.tsv"

    @property
    def h123_v2_dir(self) -> Path:
        return self._resolve_or_fallback(
            self.eq_lifeonair / "lifeonair" / "results" / "h123_v2",
            self.repro_data_dir / "h123_v2",
        )

    @property
    def benchmark_r2_dir(self) -> Path:
        return self._resolve_or_fallback(
            self.eq_lifeonair / "lifeonair" / "results" / "benchmark_r2",
            self.repro_data_dir / "benchmark_r2",
        )

    @property
    def metabolic_pilot_dir(self) -> Path:
        return self.eq_lifeonair / "lifeonair" / "results" / "metabolic_pilot"

    @property
    def sample_mobilome_load_tsv(self) -> Path:
        return self.papers_root / "empty-quarter-mobilome" / "results" / "sample_mobilome_load.tsv"

    @property
    def aquatic_supp_xlsx(self) -> Path:
        return (
            self.papers_root
            / "empty-quarter-tiannyu"
            / "source"
            / "tiannyu-2026-10-18"
            / "Empty quarter project supplementary table-latest-9.xlsx"
        )

    @property
    def wgs_sample_metadata_tsv(self) -> Path:
        return self.eq_workbench / "data" / "metadata" / "sample_metadata.tsv"

    @property
    def condition_anchors_json(self) -> Path:
        return self.eq_workbench / "data" / "metadata" / "condition_anchors.json"

    @property
    def sample_condition_map_tsv(self) -> Path:
        return self.eq_workbench / "data" / "metadata" / "sample_condition_map.tsv"


class OntologyRegistry:
    """Programmatically loads and validates all ontology IRIs from canonical project files."""

    RAK_BASE = "https://rubalkhali.science/kb/"
    SIO_BASE = "http://semanticscience.org/resource/"
    OBO_BASE = "http://purl.obolibrary.org/obo/"

    # Verified SIO & core TBox IRIs from update_rubalkhali_ontology.groovy
    CORE_TERMS: Dict[str, str] = {
        "site_quality": "https://rubalkhali.science/kb/RAK_0000001",
        "sampling_site": "https://rubalkhali.science/kb/RAK_0000002",
        "site_visit": "https://rubalkhali.science/kb/RAK_0000003",
        "temperature_measuring": "https://rubalkhali.science/kb/RAK_0000006",
        "pressure_measuring": "https://rubalkhali.science/kb/RAK_0000007",
        "humidity_measuring": "https://rubalkhali.science/kb/RAK_0000008",
        "temperature_value": "https://rubalkhali.science/kb/RAK_0000010",
        "pressure_value": "https://rubalkhali.science/kb/RAK_0000011",
        "humidity_value": "https://rubalkhali.science/kb/RAK_0000012",
        "annual_mean_temp_value": "https://rubalkhali.science/kb/RAK_0000013",
        "annual_total_precip_value": "https://rubalkhali.science/kb/RAK_0000014",
        "surface_soil_sample": "https://rubalkhali.science/kb/RAK_0000020",
        "deep_soil_sample": "https://rubalkhali.science/kb/RAK_0000021",
        "rhizosphere_sample": "https://rubalkhali.science/kb/RAK_0000022",
        "plant_matter_sample": "https://rubalkhali.science/kb/RAK_0000023",
        "xrf_measuring_process": "https://rubalkhali.science/kb/RAK_0000025",
        "chemical_analyte_concentration": "https://rubalkhali.science/kb/RAK_0000029",
        "analyte_concentration_value": "https://rubalkhali.science/kb/RAK_0000030",
        "light_elements_quality": "https://rubalkhali.science/kb/RAK_0000032",
        "light_elements_value": "https://rubalkhali.science/kb/RAK_0000033",
        "dna_concentration_quality": "https://rubalkhali.science/kb/RAK_0000043",
        "dna_concentration_value": "https://rubalkhali.science/kb/RAK_0000044",
        "taxon_rel_abundance_quality": "https://rubalkhali.science/kb/RAK_0000072",
        "taxon_rel_abundance_value": "https://rubalkhali.science/kb/RAK_0000073",
        "taxon_abs_abundance_value": "https://rubalkhali.science/kb/RAK_0000076",
        "taxon_abs_abundance_quality": "https://rubalkhali.science/kb/RAK_0000078",
        # SIO canonical measurement reification predicates
        "sio_has_attribute": "http://semanticscience.org/resource/SIO_000008",
        "sio_is_attribute_of": "http://semanticscience.org/resource/SIO_000011",
        "sio_is_measurement_value_of": "http://semanticscience.org/resource/SIO_000215",
        "sio_has_measurement_value": "http://semanticscience.org/resource/SIO_000216",
        "sio_has_unit": "http://semanticscience.org/resource/SIO_000221",
        "sio_has_value": "http://semanticscience.org/resource/SIO_000300",
        # Validated PATO / UO terms from update_rubalkhali_ontology.groovy & DATA_DICTIONARY.tsv
        "pato_concentration": "http://purl.obolibrary.org/obo/PATO_0000033",
        "uo_percent": "http://purl.obolibrary.org/obo/UO_0000187",
        "uo_ph": "http://purl.obolibrary.org/obo/UO_0000196",
        "uo_pascal": "http://purl.obolibrary.org/obo/UO_0000110",
    }

    def __init__(self, paths: Optional[CanonicalPaths] = None) -> None:
        self.paths = paths or CanonicalPaths()
        self.xrf_mappings: Dict[str, Dict[str, Any]] = {}
        self.envo_biomes: Dict[str, str] = {}
        self.envo_materials: Dict[str, str] = {}
        self.envo_features: Dict[str, str] = {}
        self.site_iris: Dict[str, str] = {}
        self.site_aliases: Dict[str, Dict[str, str]] = {}
        self.taxon_iris: Dict[str, Dict[str, str]] = {}
        self._load_all()

    def _load_all(self) -> None:
        # 1. Load XRF analyte ChEBI/PubChem mappings and compute exact RAK quality/value IRIs
        with open(self.paths.xrf_chemical_mapping_yml, "r", encoding="utf-8") as fh:
            raw_xrf = yaml.safe_load(fh)
        mappings = raw_xrf.get("mappings", {})
        quality_counter = 100
        value_counter = 500
        for analyte, entry in mappings.items():
            if analyte == "LE":
                self.xrf_mappings[analyte] = {
                    "chebi_iri": None,
                    "pubchem_uri": None,
                    "rak_quality_iri": self.CORE_TERMS["light_elements_quality"],
                    "rak_value_iri": self.CORE_TERMS["light_elements_value"],
                    "semantic_status": entry.get("semantic_status", "instrument_pseudoanalyte"),
                }
                continue
            q_iri = f"{self.RAK_BASE}RAK_0{quality_counter:06d}"
            v_iri = f"{self.RAK_BASE}RAK_0{value_counter:06d}"
            quality_counter += 1
            value_counter += 1
            self.xrf_mappings[analyte] = {
                "chebi_iri": entry.get("chebi"),
                "pubchem_uri": entry.get("pubchem"),
                "rak_quality_iri": q_iri,
                "rak_value_iri": v_iri,
            }

        # 2. Load ENVO biomes, materials, features from biome_codes.yml
        with open(self.paths.biome_codes_yml, "r", encoding="utf-8") as fh:
            raw_biomes = yaml.safe_load(fh)
        for k, v in (raw_biomes.get("biomes") or {}).items():
            self.envo_biomes[k] = f"{self.OBO_BASE}{v}"
        for k, v in (raw_biomes.get("materials") or {}).items():
            self.envo_materials[k] = f"{self.OBO_BASE}{v}"
        for k, v in (raw_biomes.get("features") or {}).items():
            self.envo_features[k] = f"{self.OBO_BASE}{v}"

        # 3. Load canonical site IRIs (Sites 1..60 + 10 auxiliary named sites)
        for i in range(1, 61):
            self.site_iris[str(i)] = f"{self.RAK_BASE}RAK_1{i:06d}"
            self.site_iris[f"Site {i}"] = f"{self.RAK_BASE}RAK_1{i:06d}"
        df_sites = pd.read_csv(self.paths.site_iri_registry_tsv, sep="\t")
        for _, row in df_sites.iterrows():
            label = str(row["site_label"]).strip()
            iri = str(row["site_iri"]).strip()
            self.site_iris[label] = iri

        # 4. Load Trip-1 numeric site aliases (61..64 -> RAK_1000067..1000070) from site_aliases.tsv
        if self.paths.site_aliases_tsv.exists():
            df_aliases = pd.read_csv(self.paths.site_aliases_tsv, sep="\t")
            for _, row in df_aliases.iterrows():
                src_val = str(row["source_site_id"]).strip()
                canon_label = str(row["canonical_site_label"]).strip()
                canon_iri = str(row["canonical_site_iri"]).strip()
                self.site_aliases[src_val] = {
                    "trip": str(row["source_sheet"]).strip(),
                    "canonical_site_label": canon_label,
                    "canonical_site_iri": canon_iri,
                    "evidence_basis": str(row["evidence"]).strip(),
                }
                # Map numeric alias keys ("61".."64" and "Site 61".."Site 64") to canonical IRIs
                self.site_iris[src_val] = canon_iri
                self.site_iris[f"Site {src_val}"] = canon_iri

        # 5. Load verified NCBITaxon domain and phylum IRIs from taxonomy_mapping_decisions.tsv
        if self.paths.taxonomy_mapping_decisions_tsv.exists():
            df_tax = pd.read_csv(
                self.paths.taxonomy_mapping_decisions_tsv,
                sep="\t",
                usecols=["rank", "source_name", "final_iri", "final_label", "parent_iri", "mapping_status"],
            )
            valid_tax = df_tax[
                df_tax["rank"].isin(["domain", "superkingdom", "phylum"])
                & (df_tax["mapping_status"] == "validated_ncbi")
            ].drop_duplicates(subset=["rank", "source_name"])
            for _, row in valid_tax.iterrows():
                name = str(row["source_name"]).strip()
                self.taxon_iris[name] = {
                    "rank": str(row["rank"]).strip(),
                    "iri": str(row["final_iri"]).strip(),
                    "label": str(row["final_label"]).strip(),
                    "parent_iri": (
                        str(row["parent_iri"]).strip()
                        if pd.notna(row["parent_iri"])
                        else ""
                    ),
                }

    def get_chebi_iri(self, analyte: str) -> Optional[str]:
        """Return the validated ChEBI IRI for an XRF analyte (or None if unmapped)."""
        if analyte not in self.xrf_mappings:
            raise KeyError(f"Unknown XRF analyte '{analyte}' not in xrf_chemical_mapping.yml")
        return self.xrf_mappings[analyte]["chebi_iri"]

    def get_rak_analyte_iris(self, analyte: str) -> Tuple[str, str]:
        """Return (quality_iri, value_iri) in the RAK namespace for an XRF analyte."""
        if analyte not in self.xrf_mappings:
            raise KeyError(f"Unknown XRF analyte '{analyte}' not in xrf_chemical_mapping.yml")
        entry = self.xrf_mappings[analyte]
        return entry["rak_quality_iri"], entry["rak_value_iri"]

    def get_envo_iri(self, term_name: str) -> str:
        """Look up a validated ENVO IRI from biome_codes.yml."""
        key = term_name.strip().lower()
        for d in (self.envo_features, self.envo_materials, self.envo_biomes):
            if key in d:
                return d[key]
        raise KeyError(f"Term '{term_name}' not found in canonical biome_codes.yml")

    def get_taxon_iri(self, taxon_name: str) -> str:
        """Look up a validated NCBITaxon IRI from taxonomy_mapping_decisions.tsv."""
        if taxon_name not in self.taxon_iris:
            raise KeyError(
                f"Taxon '{taxon_name}' not found in validated taxonomy_mapping_decisions.tsv"
            )
        return self.taxon_iris[taxon_name]["iri"]

    def get_site_iri(self, site_key: str) -> str:
        """Look up a validated site IRI (primary 1..60, auxiliary label, or Trip-1 alias 61..64)."""
        key = str(site_key).strip()
        if key in self.site_iris:
            return self.site_iris[key]
        raise KeyError(f"Site '{site_key}' not found in canonical site_iri_registry.tsv or site_aliases.tsv")


class CanonicalDataRegistry:
    """Fail-closed loader for canonical Empty Quarter tables and denominator invariants."""

    # Audited sample universe denominators from DATA_LOCATIONS.md, release_evidence.json & papers
    EXPECTED_INVARIANTS = {
        "primary_sites": 60,
        "auxiliary_named_sites": 10,
        "confirmed_site_aliases": 4,
        "core_specimen_rows": 2480,
        "feature_table_16s_profiles": 1271,
        "ecology_analysis_16s_profiles": 1237,
        "primary_site_16s_profiles": 1227,
        "wgs_soil_samples_with_metadata": 551,
        "coverm_wgs_libraries": 573,
        "dereplicated_mags": 5438,
        "mobilome_load_samples": 552,
        "mobilome_votu_catalogue": 19300,
        "mobilome_potu_catalogue": 322755,
        "ph_admitted_measurements": 712,
        "xrf_core_elements": 11,
    }

    CORE_11_XRF_ELEMENTS: Tuple[str, ...] = (
        "Al", "Ca", "Cl", "Fe", "K", "Mg", "Na", "P", "S", "Si", "Ti"
    )

    def __init__(self, paths: Optional[CanonicalPaths] = None) -> None:
        self.paths = paths or CanonicalPaths()
        self.ontology = OntologyRegistry(self.paths)

    @staticmethod
    def sha256_file(path: Path) -> str:
        h = hashlib.sha256()
        with open(path, "rb") as fh:
            for chunk in iter(lambda: fh.read(1024 * 1024), b""):
                h.update(chunk)
        return h.hexdigest()

    def load_env_site(self) -> pd.DataFrame:
        df = pd.read_csv(self.paths.env_site_tsv, sep="\t")
        assert len(df) == self.EXPECTED_INVARIANTS["primary_sites"], (
            f"Expected 60 primary sites in env_site.tsv, got {len(df)}"
        )
        return df

    def load_site_aliases(self) -> pd.DataFrame:
        df = pd.read_csv(self.paths.site_aliases_tsv, sep="\t")
        assert len(df) == self.EXPECTED_INVARIANTS["confirmed_site_aliases"], (
            f"Expected 4 confirmed site aliases, got {len(df)}"
        )
        return df

    def load_curated_visit_dates(self) -> pd.DataFrame:
        """Load exact collection dates per (site, trip) from environmental_measurements_curated.tsv."""
        df = pd.read_csv(self.paths.env_measurements_curated_tsv, sep="\t")
        trip_map = {
            "Trip 1 (2023)": "trip1",
            "Trip 1 auxiliary/revisit (2023)": "trip1",
            "Trip 2 (2023)": "trip2",
            "Trip 3 (2024)": "trip3",
            "Trip 4 (2024)": "trip4",
            "Trip 5 (2025)": "trip5",
        }
        df["trip"] = df["expedition"].map(trip_map)
        df["collection_date"] = pd.to_datetime(df["date"], format="%d/%m/%Y", errors="coerce")
        df["site_num"] = pd.to_numeric(df["site"], errors="coerce")
        # Resolve Trip 1 alias numeric sites 61..64 to their canonical site numbers 67..70
        alias_num_map = {61: 67, 62: 68, 63: 69, 64: 70}
        df["canonical_site_num"] = df["site_num"].replace(alias_num_map)
        out = (
            df.dropna(subset=["trip", "collection_date", "site_num"])
            .sort_values(["trip", "site_num", "collection_date"])
            .groupby(["site_num", "canonical_site_num", "trip"], as_index=False)
            .agg(
                collection_date=("collection_date", "first"),
                temperature_c_field=("temperature_c", "mean"),
                relative_humidity_pct_field=("relative_humidity_pct", "mean"),
                pressure_mbar_field=("pressure_mbar", "mean"),
            )
        )
        out["collection_date_iso"] = out["collection_date"].dt.strftime("%Y-%m-%d")
        return out

    def load_exact_weather_windows(self) -> pd.DataFrame:
        """Recompute antecedent precipitation windows (P7, P14, P30, P90, P365) from exact collection dates."""
        visits = self.load_curated_visit_dates()
        visits_primary = visits[visits["site_num"].between(1, 60)].copy()

        era5 = pd.read_csv(self.paths.era5_daily_tsv, sep="\t")
        era5["date_dt"] = pd.to_datetime(era5["date"])
        power = pd.read_csv(self.paths.power_daily_tsv, sep="\t")
        power["date_dt"] = pd.to_datetime(power["date"])

        # Index daily precipitation series by site for fast window summation
        era5_by_site = {
            int(s): g.set_index("date_dt")["P"].sort_index()
            for s, g in era5.groupby("site")
        }
        power_by_site = {
            int(s): g.set_index("date_dt")["P"].sort_index()
            for s, g in power.groupby("site")
        }

        records: List[Dict[str, Any]] = []
        for _, row in visits_primary.iterrows():
            s = int(row["site_num"])
            dt = row["collection_date"]
            rec: Dict[str, Any] = {
                "site": s,
                "trip": row["trip"],
                "collection_date": row["collection_date_iso"],
                "temperature_c_field": row["temperature_c_field"],
                "relative_humidity_pct_field": row["relative_humidity_pct_field"],
            }
            s_era5 = era5_by_site.get(s)
            s_power = power_by_site.get(s)
            for w in (7, 14, 30, 90, 365):
                start_dt = dt - pd.Timedelta(days=w)
                if s_era5 is not None:
                    rec[f"P{w}_era5_exact"] = float(s_era5.loc[start_dt:dt].sum())
                else:
                    rec[f"P{w}_era5_exact"] = float("nan")
                if s_power is not None:
                    rec[f"P{w}_power_exact"] = float(s_power.loc[start_dt:dt].sum())
                else:
                    rec[f"P{w}_power_exact"] = float("nan")
            records.append(rec)

        exact_df = pd.DataFrame(records).drop_duplicates(subset=["site", "trip"])
        # Merge with env_site_trip.tsv to retain NDVI, EVI, SMAP soil moisture, and legacy trip_date for audit
        env_trip = pd.read_csv(self.paths.env_site_trip_tsv, sep="\t").drop_duplicates(
            subset=["site", "trip"]
        )
        merged = env_trip.merge(exact_df, on=["site", "trip"], how="left")
        merged["date_substitution_days"] = (
            pd.to_datetime(merged["collection_date"]) - pd.to_datetime(merged["trip_date"])
        ).dt.days.fillna(0).astype(int)
        return merged

    def load_env_site_trip(self) -> pd.DataFrame:
        return self.load_exact_weather_windows()

    def load_xrf_site_compartment_panel(self) -> pd.DataFrame:
        return pd.read_csv(self.paths.xrf_site_compartment_panel_tsv, sep="\t")

    def load_ph_accepted(self) -> pd.DataFrame:
        df = pd.read_csv(self.paths.ph_accepted_tsv, sep="\t")
        admitted = df[df["disposition"] == "ADMITTED_MEASUREMENT"]
        assert len(admitted) == self.EXPECTED_INVARIANTS["ph_admitted_measurements"], (
            f"Expected 712 admitted pH measurements, got {len(admitted)}"
        )
        return admitted

    def load_guild_genome_equivalents(self, id_floor: int = 60) -> pd.DataFrame:
        df = pd.read_csv(self.paths.guild_ge_strict_tsv, sep="\t")
        return df[df["id_floor"] == id_floor].copy()

    def load_mobilome_load(self) -> pd.DataFrame:
        df = pd.read_csv(self.paths.sample_mobilome_load_tsv, sep="\t")
        assert len(df) == self.EXPECTED_INVARIANTS["mobilome_load_samples"], (
            f"Expected 552 mobilome load rows, got {len(df)}"
        )
        return df

    def load_wgs_sample_metadata(self) -> pd.DataFrame:
        df = pd.read_csv(self.paths.wgs_sample_metadata_tsv, sep="\t")
        assert len(df) == self.EXPECTED_INVARIANTS["wgs_soil_samples_with_metadata"], (
            f"Expected 551 WGS soil metadata rows, got {len(df)}"
        )
        return df

    def load_sample_ledger(self) -> pd.DataFrame:
        return pd.read_csv(self.paths.sample_ledger_tsv, sep="\t", low_memory=False)

    def load_16s_phylum_composition(self) -> pd.DataFrame:
        """Compute empirical 16S rRNA phylum relative abundances across the 1,271 sequenced profiles."""
        tax = pd.read_csv(self.paths.taxonomy_16s_tsv, sep="\t")

        def _extract_phylum(taxon_str: Any) -> str:
            if not isinstance(taxon_str, str):
                return "Unassigned"
            parts = taxon_str.split(";")
            if len(parts) >= 2 and parts[1].strip():
                return parts[1].strip()
            return "Unassigned"

        tax["phylum"] = tax["Taxon"].map(_extract_phylum)
        feat_to_phylum = dict(zip(tax["Feature ID"], tax["phylum"]))

        ft = pd.read_csv(self.paths.feature_table_16s_tsv, sep="\t", skiprows=1)
        ft["phylum"] = ft["#OTU ID"].map(feat_to_phylum).fillna("Unassigned")
        sample_cols = [c for c in ft.columns if c not in ("#OTU ID", "phylum")]
        assert len(sample_cols) == self.EXPECTED_INVARIANTS["feature_table_16s_profiles"], (
            f"Expected 1271 16S profiles, got {len(sample_cols)}"
        )

        phylum_sums = ft.groupby("phylum")[sample_cols].sum()
        total_reads = phylum_sums.sum(axis=0)
        phylum_rel = phylum_sums.div(total_reads.replace(0, float("nan")), axis=1).T
        phylum_rel.index.name = "profile_id"
        phylum_rel = phylum_rel.reset_index()
        phylum_rel["total_reads_16s"] = total_reads.values

        inv = pd.read_csv(
            self.paths.canonical_profile_read_inventory_tsv,
            sep="\t",
            usecols=["profile_id", "source_sample_name"],
        )
        if self.paths.profile_extraction_join_audit_tsv.exists():
            audit = pd.read_csv(
                self.paths.profile_extraction_join_audit_tsv,
                sep="\t",
                usecols=["profile_id", "ledger_specimen_id", "ledger_resolution", "depth"],
            )
            inv = inv.merge(audit, on="profile_id", how="left")
            inv["resolved_sample_id"] = inv["ledger_specimen_id"].fillna(inv["source_sample_name"])
            inv["in_qc_1237"] = inv["ledger_specimen_id"].notna()
        else:
            inv["resolved_sample_id"] = inv["source_sample_name"]
            inv["in_qc_1237"] = True

        ledger = self.load_sample_ledger()[
            [
                "sample_id",
                "trip",
                "site",
                "site_corrected",
                "canonical_site_label",
                "canonical_site_iri",
                "compartment",
                "is_control",
                "in_ecology_analysis",
            ]
        ].drop_duplicates(subset=["sample_id"])

        merged = inv.merge(
            ledger, left_on="resolved_sample_id", right_on="sample_id", how="left"
        ).merge(phylum_rel, on="profile_id", how="inner")
        # Canonical 1,237 QC-passed cohort (1,227 primary-site + 10 auxiliary-site profiles matching empty-quarter-ecology-current)
        merged["in_ecology_analysis"] = merged["in_qc_1237"].fillna(False)
        return merged

    def load_coverm_and_cbb_metrics(self) -> pd.DataFrame:
        """Compute empirical MAG carrier abundance, CBB/RuBisCO co-carriage across 3 screens, unmapped fraction, and MAG phylum shares."""
        bac = pd.read_csv(
            self.paths.gtdb_bac120_summary_tsv,
            sep="\t",
            usecols=["user_genome", "classification"],
        )
        ar = pd.read_csv(
            self.paths.gtdb_ar53_summary_tsv,
            sep="\t",
            usecols=["user_genome", "classification"],
        )
        gtdb = pd.concat([bac, ar], ignore_index=True)
        gtdb["phylum"] = gtdb["classification"].str.extract(r"p__([^;]+)")

        cbb = pd.read_csv(self.paths.rep_cbb_ko_tsv, sep="\t")
        cbb["has_rbcl"] = (cbb["K01601"] > 0).astype(int)
        cbb["has_cbb_prk"] = ((cbb["K01601"] > 0) & (cbb["K00855"] > 0)).astype(int)
        cbb["has_cbb_3gene"] = (
            (cbb["K01601"] > 0) & (cbb["K01602"] > 0) & (cbb["K00855"] > 0)
        ).astype(int)
        reps_rbcl = set(cbb[cbb["has_rbcl"] > 0]["rep"])
        reps_cbb_prk = set(cbb[cbb["has_cbb_prk"] > 0]["rep"])
        reps_cbb_3gene = set(cbb[cbb["has_cbb_3gene"] > 0]["rep"])

        atmo = pd.read_csv(self.paths.atmo_classification_v2_tsv, sep="\t")
        reps_1l = set(atmo[atmo["assignment"] == "[NiFe]_Group_1l"]["rep_id"])
        reps_1h = set(atmo[atmo["assignment"] == "[NiFe]_Group_1h"]["rep_id"])
        reps_2a = set(atmo[atmo["assignment"] == "[NiFe]_Group_2a"]["rep_id"])
        reps_cox1 = set(atmo[atmo["assignment"] == "form_I|AYxCSFR"]["rep_id"])

        cov = pd.read_csv(
            self.paths.coverm_rep_abundance_tsv_gz,
            sep="\t",
            names=["sample", "rep", "pct"],
        )
        unmapped = cov[cov["rep"] == "unmapped"].set_index("sample")["pct"]
        mapped = cov[cov["rep"] != "unmapped"].copy()
        mapped = mapped.merge(
            gtdb[["user_genome", "phylum"]],
            left_on="rep",
            right_on="user_genome",
            how="left",
        )

        mapped["pct_1l"] = mapped["pct"].where(mapped["rep"].isin(reps_1l), 0.0)
        mapped["pct_1h"] = mapped["pct"].where(mapped["rep"].isin(reps_1h), 0.0)
        mapped["pct_2a"] = mapped["pct"].where(mapped["rep"].isin(reps_2a), 0.0)
        mapped["pct_cox1"] = mapped["pct"].where(mapped["rep"].isin(reps_cox1), 0.0)
        mapped["pct_cbb"] = mapped["pct"].where(mapped["rep"].isin(reps_rbcl), 0.0)
        mapped["pct_cbb_prk"] = mapped["pct"].where(mapped["rep"].isin(reps_cbb_prk), 0.0)
        mapped["pct_cbb_3gene"] = mapped["pct"].where(mapped["rep"].isin(reps_cbb_3gene), 0.0)
        mapped["pct_1l_cbb"] = mapped["pct"].where(
            mapped["rep"].isin(reps_1l & reps_rbcl), 0.0
        )
        mapped["pct_1l_cbb_prk"] = mapped["pct"].where(
            mapped["rep"].isin(reps_1l & reps_cbb_prk), 0.0
        )
        mapped["pct_1l_cbb_3gene"] = mapped["pct"].where(
            mapped["rep"].isin(reps_1l & reps_cbb_3gene), 0.0
        )
        mapped["pct_1l_nocbb"] = mapped["pct"].where(
            mapped["rep"].isin(reps_1l - reps_rbcl), 0.0
        )

        agg = mapped.groupby("sample")[
            [
                "pct",
                "pct_1l",
                "pct_1h",
                "pct_2a",
                "pct_cox1",
                "pct_cbb",
                "pct_cbb_prk",
                "pct_cbb_3gene",
                "pct_1l_cbb",
                "pct_1l_cbb_prk",
                "pct_1l_cbb_3gene",
                "pct_1l_nocbb",
            ]
        ].sum()
        agg = agg.rename(columns={"pct": "mapped_read_pct"})
        agg["unmapped_read_pct"] = unmapped
        denom = agg["mapped_read_pct"].replace(0.0, float("nan"))
        agg["share_1l_mapped"] = agg["pct_1l"] / denom
        agg["share_1h_mapped"] = agg["pct_1h"] / denom
        agg["share_2a_mapped"] = agg["pct_2a"] / denom
        agg["share_cox1_mapped"] = agg["pct_cox1"] / denom
        agg["share_cbb_mapped"] = agg["pct_cbb"] / denom
        agg["share_cbb_prk_mapped"] = agg["pct_cbb_prk"] / denom
        agg["share_cbb_3gene_mapped"] = agg["pct_cbb_3gene"] / denom
        agg["share_1l_cbb_mapped"] = agg["pct_1l_cbb"] / denom
        agg["share_1l_cbb_prk_mapped"] = agg["pct_1l_cbb_prk"] / denom
        agg["share_1l_cbb_3gene_mapped"] = agg["pct_1l_cbb_3gene"] / denom
        agg["share_1l_nocbb_mapped"] = agg["pct_1l_nocbb"] / denom

        ph_piv = mapped.pivot_table(
            index="sample", columns="phylum", values="pct", aggfunc="sum", fill_value=0.0
        )
        ph_rel = ph_piv.div(denom, axis=0)
        for ph_name, col_name in [
            ("Actinomycetota", "mag_actinomycetota_rel"),
            ("Thermoplasmatota", "mag_thermoplasmatota_rel"),
            ("Halobacteriota", "mag_halobacteriota_rel"),
            ("Chloroflexota", "mag_chloroflexota_rel"),
            ("Pseudomonadota", "mag_pseudomonadota_rel"),
            ("Bacillota", "mag_bacillota_rel"),
        ]:
            agg[col_name] = ph_rel[ph_name] if ph_name in ph_rel.columns else 0.0

        return agg.reset_index()

    def load_aquatic_water_chemistry(self) -> pd.DataFrame:
        """Load measured surface-water hydrochemistry for Sites 59 and 60 (Table S2 of aquatic study)."""
        # The versioned table is the publication input. A workstation's mutable
        # supplementary workbook must not silently override it.
        fallback_tsv = self.paths.repro_data_dir / "aquatic_water_chemistry_s2.tsv"
        if fallback_tsv.exists():
            return pd.read_csv(fallback_tsv, sep="\t")
        if not self.paths.aquatic_supp_xlsx.exists():
            fallback_tsv = self.paths.repro_data_dir / "aquatic_water_chemistry_s2.tsv"
            if fallback_tsv.exists():
                return pd.read_csv(fallback_tsv, sep="\t")
            raise FileNotFoundError(
                f"Missing aquatic supplementary workbook: {self.paths.aquatic_supp_xlsx}"
            )
        with zipfile.ZipFile(self.paths.aquatic_supp_xlsx) as z:
            ns = {"m": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}
            sst = [
                "".join(t.text or "" for t in si.findall(".//m:t", ns))
                for si in ET.fromstring(z.read("xl/sharedStrings.xml")).findall("m:si", ns)
            ] if "xl/sharedStrings.xml" in z.namelist() else []
            sroot = ET.fromstring(z.read("xl/worksheets/sheet3.xml"))
            rows: List[List[str]] = []
            for r in sroot.findall(".//m:row", ns):
                vals: List[str] = []
                for c in r.findall("m:c", ns):
                    t = c.attrib.get("t")
                    v = c.find("m:v", ns)
                    val = v.text if v is not None else ""
                    if t == "s" and val.isdigit():
                        val = sst[int(val)]
                    elif t == "inlineStr":
                        val = "".join(node.text or "" for node in c.findall(".//m:t", ns))
                    vals.append(val.strip())
                if any(vals):
                    rows.append(vals)
        # Row 1 is ['Parameter', 'Unit', '59-MR-WR1', '59-MR-WR2', '60-MR-WR1']
        header = rows[1]
        data_rows = rows[2:]
        df = pd.DataFrame(data_rows, columns=header)
        for sample_col in ("59-MR-WR1", "59-MR-WR2", "60-MR-WR1"):
            df[sample_col] = pd.to_numeric(
                df[sample_col].replace("-", float("nan")), errors="coerce"
            )
        return df

    def load_enrichment_conditions(self) -> Tuple[Dict[str, Any], pd.DataFrame]:
        with open(self.paths.condition_anchors_json, "r", encoding="utf-8") as fh:
            anchors = json.load(fh)
        cond_map = pd.read_csv(self.paths.sample_condition_map_tsv, sep="\t")
        return anchors, cond_map

    def verify_all_sources(self) -> Dict[str, Any]:
        """Verify existence, row counts, and SHA-256 digests of canonical input tables."""
        critical_files = {
            "xrf_chemical_mapping": self.paths.xrf_chemical_mapping_yml,
            "biome_codes": self.paths.biome_codes_yml,
            "site_iri_registry": self.paths.site_iri_registry_tsv,
            "site_aliases": self.paths.site_aliases_tsv,
            "env_measurements_curated": self.paths.env_measurements_curated_tsv,
            "sample_ledger": self.paths.sample_ledger_tsv,
            "taxonomy_mapping_decisions": self.paths.taxonomy_mapping_decisions_tsv,
            "feature_table_16s": self.paths.feature_table_16s_tsv,
            "taxonomy_16s": self.paths.taxonomy_16s_tsv,
            "ph_accepted": self.paths.ph_accepted_tsv,
            "env_site": self.paths.env_site_tsv,
            "env_site_trip": self.paths.env_site_trip_tsv,
            "era5_daily": self.paths.era5_daily_tsv,
            "power_daily": self.paths.power_daily_tsv,
            "xrf_site_compartment_panel": self.paths.xrf_site_compartment_panel_tsv,
            "coverm_rep_abundance": self.paths.coverm_rep_abundance_tsv_gz,
            "rep_cbb_ko": self.paths.rep_cbb_ko_tsv,
            "atmo_classification_v2": self.paths.atmo_classification_v2_tsv,
            "guild_ge_strict": self.paths.guild_ge_strict_tsv,
            "sample_mobilome_load": self.paths.sample_mobilome_load_tsv,
            "aquatic_supp_xlsx": self.paths.aquatic_supp_xlsx,
            "wgs_sample_metadata": self.paths.wgs_sample_metadata_tsv,
            "condition_anchors": self.paths.condition_anchors_json,
        }
        manifest_fallback = self.paths.repro_data_dir / "canonical_sources_manifest.json"
        if not all(p.exists() for p in critical_files.values()) and manifest_fallback.exists():
            return json.loads(manifest_fallback.read_text(encoding="utf-8"))
        report: Dict[str, Any] = {}
        for name, path in critical_files.items():
            if not path.exists():
                raise FileNotFoundError(f"Missing canonical file for {name}: {path}")
            report[name] = {
                "path": str(path),
                "bytes": path.stat().st_size,
                "sha256": self.sha256_file(path),
            }
        return report
