"""Builder for the Empty Quarter Digital Twin harmonized spatiotemporal state cube.

Constructs:
1. Site-level spatial cube (60 primary sites + 10 canonical auxiliary sites) indexed
   with Uber H3 (resolutions 5 and 6), 0.25-degree Open-Meteo/ERA5 cells, and
   verified Trip-1 numeric site alias provenance (61..64 -> RAK_1000067..1000070).
2. Site x Compartment geochemical & biodiversity panel (Surface, ShallowSubsurface,
   RootAdjacent) with CLR-transformed 11 core XRF elements, evaporite-vs-quartz PC1,
   admitted CaCl2 pH, 16S Shannon diversity, empirical 16S ASV and shotgun CoverM MAG
   phylum relative abundances, trace-gas guild genome equivalents at 60% identity
   (hydB_1l, hydB_1h, hydB_2a, coxL_formI, hydB_other_nife, coxL_other), CBB RuBisCO
   co-carriage, unmapped read fraction, and vOTU/pOTU mobilome load.
3. Longitudinal Site x Trip x Compartment panel across all 5 campaigns (March 2023 -
   October 2025) with exact collection-date antecedent precipitation windows
   (P7, P14, P30, P90, P365 from ERA5 and NASA POWER).
4. Continuous 0.25-degree ERA5 / H3 interpolation grid across the 1,043 km transect
   with explicit convex-hull / nearest-neighbor domain-of-support flags.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import h3
import numpy as np
import pandas as pd

from eq_biodt.data.manifest_registry import CanonicalDataRegistry


def _latlng_to_h3(lat: float, lon: float, res: int) -> str:
    """Compatible wrapper for h3 v4 (`latlng_to_cell`) and h3 v3 (`geo_to_h3`)."""
    if hasattr(h3, "latlng_to_cell"):
        return h3.latlng_to_cell(float(lat), float(lon), int(res))
    return h3.geo_to_h3(float(lat), float(lon), int(res))


def compute_clr_xrf(
    df: pd.DataFrame, elements: Tuple[str, ...], pseudocount: float = 1e-4
) -> pd.DataFrame:
    """Centered log-ratio (CLR) transformation across the 11 core XRF elements."""
    mat = df[list(elements)].to_numpy(dtype=float).copy()
    mat = np.where(np.isnan(mat) | (mat <= 0.0), pseudocount, mat)
    log_mat = np.log(mat)
    gm = log_mat.mean(axis=1, keepdims=True)
    clr = log_mat - gm
    out = pd.DataFrame(
        clr,
        index=df.index,
        columns=[f"clr_{el}" for el in elements],
    )
    # Evaporite vs quartz lithological score (positive = Ca/Mg/Na/S/Cl evaporite/sabkha, negative = Si quartz)
    evap_cols = [
        f"clr_{el}" for el in ("Ca", "Mg", "Na", "S", "Cl") if f"clr_{el}" in out.columns
    ]
    out["xrf_evaporite_index"] = out[evap_cols].mean(axis=1) - out["clr_Si"]
    return out


class TwinStateCubeBuilder:
    """Harmonizes all Empty Quarter observational layers into the EQ-BioDT state cube."""

    COMPARTMENT_NORMALIZATION: Dict[str, str] = {
        "surface": "Surface",
        "Surface": "Surface",
        "S": "Surface",
        "deep": "ShallowSubsurface",
        "Deep": "ShallowSubsurface",
        "ShallowSubsurface": "ShallowSubsurface",
        "D": "ShallowSubsurface",
        "rhizosphere": "RootAdjacent",
        "Rhizosphere": "RootAdjacent",
        "PlantRoot": "RootAdjacent",
        "plantroot": "RootAdjacent",
        "RootAdjacent": "RootAdjacent",
        "PR": "RootAdjacent",
        "R": "RootAdjacent",
        "Plant": "PlantMatter",
        "plant": "PlantMatter",
    }

    AUXILIARY_COORDINATES: Dict[str, Tuple[float, float, float, str, Optional[int], bool]] = {
        # label: (lat, lon, elev_m, regime, trip1_numeric_alias, has_canonical_sequencing)
        "Site oill spill": (19.405360, 45.727080, 719.0, "hydrocarbon_impacted", None, False),
        "Site hot spring": (20.397230, 50.965530, 208.0, "hydrothermal_spring", None, False),
        "Site water well (location 1)": (20.397223, 50.965705, 208.0, "groundwater_well", None, False),
        "Site Camel1": (19.654010, 47.688260, 468.0, "camel_grazing_locus", None, False),
        "Site plant gradient": (19.709150, 47.959330, 412.0, "rhizosphere_transect", None, False),
        "Site 19.5": (19.819271, 48.240237, 376.0, "interdune_waypoint", None, False),
        "Site water well (location 2)": (20.715870, 54.799697, 81.0, "groundwater_well", 61, True),
        "Site road (location 1)": (20.825776, 54.988932, 87.0, "corridor_verge", 62, True),
        "Site road (location 2)": (20.851514, 53.757884, 111.0, "corridor_verge", 63, True),
        "Site camground": (21.031051, 55.159843, 83.0, "anthropogenic_camp", 64, True),
    }

    def __init__(self, registry: Optional[CanonicalDataRegistry] = None) -> None:
        self.registry = registry or CanonicalDataRegistry()
        self.ontology = self.registry.ontology

    @staticmethod
    def assign_transect_zone(lon: float, site_num: Optional[int] = None) -> str:
        """Assign ecological sector along the 1,043 km Wadi ad-Dawasir -> Umm as-Samim transect."""
        if site_num is not None and site_num in (59, 60):
            return "umm_as_samim_sabkha"
        if lon < 48.0:
            return "western_plateau"
        if lon < 52.8:
            return "central_hyperarid_erg"
        return "eastern_sabkha_frontier"

    def build_sites_cube(self) -> pd.DataFrame:
        """Build the 70-site spatial table (60 primary transect + 10 auxiliary sites)."""
        if not self.registry.paths.env_site_tsv.exists():
            cached = self.registry.paths.processed_dir / "eq_sites_cube.parquet"
            if cached.exists():
                return pd.read_parquet(cached)
        env_site = self.registry.load_env_site().copy()
        rows: List[Dict[str, Any]] = []

        for _, r in env_site.iterrows():
            s_num = int(r["site"])
            lat = float(r["lat"])
            lon = float(r["lon"])
            elev = float(r["elevation_m"])
            label = f"Site {s_num}"
            rows.append(
                {
                    "site_id": str(s_num),
                    "site_num": s_num,
                    "site_label": label,
                    "site_iri": self.ontology.site_iris[label],
                    "is_primary_transect": True,
                    "trip1_numeric_alias": None,
                    "has_canonical_sequencing": True,
                    "habitat_regime": (
                        "hypersaline_sabkha"
                        if s_num == 59
                        else ("brackish_artesian_pool" if s_num == 60 else "aeolian_desert_soil")
                    ),
                    "transect_zone": self.assign_transect_zone(lon, s_num),
                    "lat": lat,
                    "lon": lon,
                    "elevation_m": elev,
                    "map_mm": float(r["map_bestmatch_mm"]),
                    "map_era5_mm": float(r["map_era5_mm"]),
                    "mat_c": float(r["mat_era5_c"]),
                    "rainy_days_per_yr": float(r["rainy_days_per_yr"]),
                    "era5_cell": str(r["era5_cell"]),
                    "h3_5": _latlng_to_h3(lat, lon, 5),
                    "h3_6": _latlng_to_h3(lat, lon, 6),
                }
            )

        # Add 10 canonical auxiliary sites from site_iri_registry.tsv & site_aliases.tsv
        for label, (lat, lon, elev, regime, alias_num, has_seq) in self.AUXILIARY_COORDINATES.items():
            iri = self.ontology.site_iris[label]
            era5_lat = round(lat * 4.0) / 4.0
            era5_lon = round(lon * 4.0) / 4.0
            dists = (env_site["lat"] - lat) ** 2 + (env_site["lon"] - lon) ** 2
            nearest = env_site.loc[dists.idxmin()]
            aux_num = int(iri.split("RAK_10000")[1])
            rows.append(
                {
                    "site_id": iri.split("/")[-1],
                    "site_num": aux_num,
                    "site_label": label,
                    "site_iri": iri,
                    "is_primary_transect": False,
                    "trip1_numeric_alias": alias_num,
                    "has_canonical_sequencing": has_seq,
                    "habitat_regime": regime,
                    "transect_zone": self.assign_transect_zone(lon, None),
                    "lat": lat,
                    "lon": lon,
                    "elevation_m": elev,
                    "map_mm": float(nearest["map_bestmatch_mm"]),
                    "map_era5_mm": float(nearest["map_era5_mm"]),
                    "mat_c": float(nearest["mat_era5_c"]),
                    "rainy_days_per_yr": float(nearest["rainy_days_per_yr"]),
                    "era5_cell": f"{era5_lat:.2f}_{era5_lon:.2f}",
                    "h3_5": _latlng_to_h3(lat, lon, 5),
                    "h3_6": _latlng_to_h3(lat, lon, 6),
                }
            )

        return pd.DataFrame(rows)

    def build_compartment_state_cube(self) -> pd.DataFrame:
        """Build the Site x Compartment state cube combining XRF, pH, 16S ASVs, WGS guilds, CoverM MAGs, and mobilome."""
        if not self.registry.paths.xrf_site_compartment_panel_tsv.exists():
            cached = self.registry.paths.processed_dir / "eq_compartment_state_cube.parquet"
            if cached.exists():
                return pd.read_parquet(cached)
        sites_df = self.build_sites_cube()
        primary_sites = sites_df[sites_df["is_primary_transect"]].copy()

        # 1. Load XRF + 16S Shannon panel
        xrf_panel = self.registry.load_xrf_site_compartment_panel().copy()
        xrf_panel["compartment_norm"] = xrf_panel["compartment"].map(self.COMPARTMENT_NORMALIZATION)
        clr_df = compute_clr_xrf(xrf_panel, self.registry.CORE_11_XRF_ELEMENTS)
        xrf_full = pd.concat([xrf_panel, clr_df], axis=1)

        # 2. Load admitted CaCl2 pH and aggregate by (site, compartment_norm)
        ph_df = self.registry.load_ph_accepted().copy()
        ph_df["compartment_norm"] = ph_df["compartment"].map(self.COMPARTMENT_NORMALIZATION)
        ph_df["site_num"] = pd.to_numeric(ph_df["site"], errors="coerce")
        ph_agg = (
            ph_df.dropna(subset=["site_num", "compartment_norm"])
            .groupby(["site_num", "compartment_norm"], as_index=False)
            .agg(
                ph_mean_raw=("ph_value", "mean"),
                ph_std=("ph_value", "std"),
                ph_n=("ph_value", "count"),
            )
        )

        # 3. Load WGS sample metadata + trace-gas guild GEs (60% identity) + CoverM MAG/CBB + mobilome load
        wgs_meta = self.registry.load_wgs_sample_metadata().copy()
        wgs_meta["compartment_norm"] = wgs_meta["compartment"].map(self.COMPARTMENT_NORMALIZATION)
        if wgs_meta["compartment_norm"].isna().any():
            unmapped = wgs_meta.loc[wgs_meta["compartment_norm"].isna(), "compartment"].unique()
            raise ValueError(f"Unmapped WGS compartment labels: {unmapped}")
        guild_ge = self.registry.load_guild_genome_equivalents(id_floor=60)
        mobilome = self.registry.load_mobilome_load()
        coverm_cbb = self.registry.load_coverm_and_cbb_metrics()

        wgs_merged = (
            wgs_meta.merge(guild_ge, on="sample", how="left")
            .merge(mobilome, on="sample", how="left")
            .merge(coverm_cbb, on="sample", how="left")
        )
        wgs_agg = (
            wgs_merged.groupby(["site", "compartment_norm"], as_index=False)
            .agg(
                wgs_sample_count=("sample", "count"),
                hydB_1l_ge_raw=("hydB_1l", "mean"),
                hydB_1h_ge_raw=("hydB_1h", "mean"),
                hydB_2a_ge_raw=("hydB_2a", "mean"),
                coxL_formI_ge_raw=("coxL_formI", "mean"),
                hydB_other_nife_ge_raw=("hydB_other_nife", "mean"),
                coxL_other_ge_raw=("coxL_other", "mean"),
                viral_richness_mean_raw=("viral_richness", "mean"),
                viral_tpm_mean_raw=("viral_tpm_sum", "mean"),
                plasmid_richness_mean_raw=("plasmid_richness", "mean"),
                plasmid_tpm_mean_raw=("plasmid_tpm_sum", "mean"),
                mapped_read_pct_raw=("mapped_read_pct", "mean"),
                unmapped_read_pct_raw=("unmapped_read_pct", "mean"),
                share_1l_mapped_raw=("share_1l_mapped", "mean"),
                share_1l_cbb_mapped_raw=("share_1l_cbb_mapped", "mean"),
                share_1l_cbb_prk_mapped_raw=("share_1l_cbb_prk_mapped", "mean"),
                share_1l_cbb_3gene_mapped_raw=("share_1l_cbb_3gene_mapped", "mean"),
                share_1l_nocbb_mapped_raw=("share_1l_nocbb_mapped", "mean"),
                share_cbb_mapped_raw=("share_cbb_mapped", "mean"),
                share_cbb_prk_mapped_raw=("share_cbb_prk_mapped", "mean"),
                share_cbb_3gene_mapped_raw=("share_cbb_3gene_mapped", "mean"),
                mag_actinomycetota_raw=("mag_actinomycetota_rel", "mean"),
                mag_thermoplasmatota_raw=("mag_thermoplasmatota_rel", "mean"),
                mag_halobacteriota_raw=("mag_halobacteriota_rel", "mean"),
                mag_chloroflexota_raw=("mag_chloroflexota_rel", "mean"),
                mag_pseudomonadota_raw=("mag_pseudomonadota_rel", "mean"),
                mag_bacillota_raw=("mag_bacillota_rel", "mean"),
            )
        )

        # 4. Load empirical 16S rRNA phylum composition across primary sites
        amp_df = self.registry.load_16s_phylum_composition()
        amp_eco = amp_df[amp_df["in_ecology_analysis"] == True].copy()
        amp_eco["site_num"] = pd.to_numeric(amp_eco["site"], errors="coerce")
        amp_eco["compartment_norm"] = amp_eco["compartment"].map(self.COMPARTMENT_NORMALIZATION)
        amp_agg = (
            amp_eco.dropna(subset=["site_num", "compartment_norm"])
            .groupby(["site_num", "compartment_norm"], as_index=False)
            .agg(
                amplicon_16s_profile_count=("profile_id", "count"),
                actino_16s_rel_raw=("Actinomycetota", "mean"),
                thermoplasmatota_16s_rel_raw=("Thermoplasmatota", "mean"),
                halobacteriota_16s_rel_raw=("Halobacteriota", "mean"),
                chloroflexota_16s_rel_raw=("Chloroflexota", "mean"),
                pseudomonadota_16s_rel_raw=("Pseudomonadota", "mean"),
                bacillota_16s_rel_raw=("Bacillota", "mean"),
            )
        )

        # 5. Load DNA concentration from sample_ledger
        ledger = self.registry.load_sample_ledger()
        ledger_core = ledger[ledger["campaign_role"] == "core_sites_1_60"].copy()
        ledger_core["site_num"] = pd.to_numeric(ledger_core["site"], errors="coerce")
        ledger_core["compartment_norm"] = ledger_core["compartment"].map(self.COMPARTMENT_NORMALIZATION)
        dna_agg = (
            ledger_core.dropna(subset=["site_num", "compartment_norm"])
            .groupby(["site_num", "compartment_norm"], as_index=False)
            .agg(
                dna_conc_mean_raw=("dna_concentration", "mean"),
                specimen_count=("sample_id", "count"),
            )
        )

        # Merge onto (site, compartment_norm)
        cube = (
            xrf_full.rename(columns={"site": "site_num", "Shannon": "shannon_16s"})
            .merge(
                primary_sites[
                    [
                        "site_num",
                        "site_label",
                        "site_iri",
                        "transect_zone",
                        "habitat_regime",
                        "lat",
                        "lon",
                        "elevation_m",
                        "map_mm",
                        "map_era5_mm",
                        "mat_c",
                        "rainy_days_per_yr",
                        "h3_5",
                        "h3_6",
                        "era5_cell",
                    ]
                ],
                on="site_num",
                how="left",
            )
            .merge(ph_agg, on=["site_num", "compartment_norm"], how="left")
            .merge(
                wgs_agg.rename(columns={"site": "site_num"}),
                on=["site_num", "compartment_norm"],
                how="left",
            )
            .merge(amp_agg, on=["site_num", "compartment_norm"], how="left")
            .merge(dna_agg, on=["site_num", "compartment_norm"], how="left")
        )

        # Record explicit observation provenance flags before any compartment-median imputation
        cube["ph_observed"] = cube["ph_mean_raw"].notna()
        cube["wgs_observed"] = cube["wgs_sample_count"].notna() & (cube["wgs_sample_count"] > 0)
        cube["amplicon_16s_observed"] = (
            cube["amplicon_16s_profile_count"].notna() & (cube["amplicon_16s_profile_count"] > 0)
        )
        cube["dna_observed"] = cube["dna_conc_mean_raw"].notna()
        cube["wgs_sample_count"] = cube["wgs_sample_count"].fillna(0).astype(int)
        cube["amplicon_16s_profile_count"] = cube["amplicon_16s_profile_count"].fillna(0).astype(int)

        # Fill missing values into working columns while keeping *_raw strictly unimputed
        default_ph = float(cube["ph_mean_raw"].median()) if cube["ph_mean_raw"].notna().any() else 7.95
        cube["ph_mean"] = cube.groupby("compartment_norm")["ph_mean_raw"].transform(
            lambda s: s.fillna(s.median() if s.notna().any() else default_ph)
        )
        cube["dna_conc_mean"] = cube.groupby("compartment_norm")["dna_conc_mean_raw"].transform(
            lambda s: s.fillna(s.median() if s.notna().any() else 0.5)
        )

        raw_to_filled_pairs = [
            ("hydB_1l_ge_raw", "hydB_1l_ge"),
            ("hydB_1h_ge_raw", "hydB_1h_ge"),
            ("hydB_2a_ge_raw", "hydB_2a_ge"),
            ("coxL_formI_ge_raw", "coxL_formI_ge"),
            ("hydB_other_nife_ge_raw", "hydB_other_nife_ge"),
            ("coxL_other_ge_raw", "coxL_other_ge"),
            ("viral_richness_mean_raw", "viral_richness_mean"),
            ("viral_tpm_mean_raw", "viral_tpm_mean"),
            ("plasmid_richness_mean_raw", "plasmid_richness_mean"),
            ("plasmid_tpm_mean_raw", "plasmid_tpm_mean"),
            ("mapped_read_pct_raw", "mapped_read_pct"),
            ("unmapped_read_pct_raw", "unmapped_read_pct"),
            ("share_1l_mapped_raw", "share_1l_mapped"),
            ("share_1l_cbb_mapped_raw", "share_1l_cbb_mapped"),
            ("share_1l_cbb_prk_mapped_raw", "share_1l_cbb_prk_mapped"),
            ("share_1l_cbb_3gene_mapped_raw", "share_1l_cbb_3gene_mapped"),
            ("share_1l_nocbb_mapped_raw", "share_1l_nocbb_mapped"),
            ("share_cbb_mapped_raw", "share_cbb_mapped"),
            ("share_cbb_prk_mapped_raw", "share_cbb_prk_mapped"),
            ("share_cbb_3gene_mapped_raw", "share_cbb_3gene_mapped"),
            ("mag_actinomycetota_raw", "actinomycetota_rel_abund"),
            ("mag_thermoplasmatota_raw", "thermoplasmatota_rel_abund"),
            ("mag_halobacteriota_raw", "halobacteriota_rel_abund"),
            ("mag_chloroflexota_raw", "chloroflexota_rel_abund"),
            ("mag_pseudomonadota_raw", "pseudomonadota_rel_abund"),
            ("mag_bacillota_raw", "bacillota_rel_abund"),
            ("actino_16s_rel_raw", "actino_16s_rel_abund"),
            ("thermoplasmatota_16s_rel_raw", "thermoplasmatota_16s_rel_abund"),
            ("halobacteriota_16s_rel_raw", "halobacteriota_16s_rel_abund"),
            ("chloroflexota_16s_rel_raw", "chloroflexota_16s_rel_abund"),
            ("pseudomonadota_16s_rel_raw", "pseudomonadota_16s_rel_abund"),
            ("bacillota_16s_rel_raw", "bacillota_16s_rel_abund"),
        ]
        for raw_col, filled_col in raw_to_filled_pairs:
            global_med = float(cube[raw_col].median()) if cube[raw_col].notna().any() else 0.0
            cube[filled_col] = cube.groupby("compartment_norm")[raw_col].transform(
                lambda s, gm=global_med: s.fillna(s.median() if s.notna().any() else gm)
            )

        # Map compartment IRI from OntologyRegistry
        comp_iri_map = {
            "Surface": self.ontology.CORE_TERMS["surface_soil_sample"],
            "ShallowSubsurface": self.ontology.CORE_TERMS["deep_soil_sample"],
            "RootAdjacent": self.ontology.CORE_TERMS["rhizosphere_sample"],
        }
        cube["compartment_iri"] = cube["compartment_norm"].map(comp_iri_map)

        return cube

    def build_visit_compartment_cube(self) -> pd.DataFrame:
        """Build the longitudinal (site, trip, compartment) panel across Trips 1-5 (2023-2025).

        Joins sample-level WGS guild genome equivalents (at 60% identity), CoverM MAG &
        CBB co-carriage, mobilome load, and exact collection-date weather windows
        (P7, P14, P30, P90, P365 from ERA5 and NASA POWER).
        """
        if not self.registry.paths.wgs_sample_metadata_tsv.exists():
            cached = self.registry.paths.processed_dir / "eq_visit_compartment_cube.parquet"
            if cached.exists():
                return pd.read_parquet(cached)
        wgs_meta = self.registry.load_wgs_sample_metadata().copy()
        wgs_meta["compartment_norm"] = wgs_meta["compartment"].map(self.COMPARTMENT_NORMALIZATION)
        guild_ge = self.registry.load_guild_genome_equivalents(id_floor=60)
        mobilome = self.registry.load_mobilome_load()
        coverm_cbb = self.registry.load_coverm_and_cbb_metrics()

        wgs_samples = (
            wgs_meta.merge(guild_ge, on="sample", how="left")
            .merge(mobilome, on="sample", how="left")
            .merge(coverm_cbb, on="sample", how="left")
        )

        visit_agg = (
            wgs_samples.groupby(["site", "trip", "compartment_norm"], as_index=False)
            .agg(
                wgs_sample_count=("sample", "count"),
                hydB_1l_ge=("hydB_1l", "mean"),
                hydB_1h_ge=("hydB_1h", "mean"),
                hydB_2a_ge=("hydB_2a", "mean"),
                coxL_formI_ge=("coxL_formI", "mean"),
                hydB_other_nife_ge=("hydB_other_nife", "mean"),
                coxL_other_ge=("coxL_other", "mean"),
                viral_richness=("viral_richness", "mean"),
                viral_tpm_sum=("viral_tpm_sum", "mean"),
                plasmid_richness=("plasmid_richness", "mean"),
                plasmid_tpm_sum=("plasmid_tpm_sum", "mean"),
                mapped_read_pct=("mapped_read_pct", "mean"),
                unmapped_read_pct=("unmapped_read_pct", "mean"),
                share_1l_mapped=("share_1l_mapped", "mean"),
                share_1l_cbb_mapped=("share_1l_cbb_mapped", "mean"),
                share_1l_nocbb_mapped=("share_1l_nocbb_mapped", "mean"),
                share_cbb_mapped=("share_cbb_mapped", "mean"),
                actinomycetota_rel_abund=("mag_actinomycetota_rel", "mean"),
                thermoplasmatota_rel_abund=("mag_thermoplasmatota_rel", "mean"),
                halobacteriota_rel_abund=("mag_halobacteriota_rel", "mean"),
                chloroflexota_rel_abund=("mag_chloroflexota_rel", "mean"),
            )
            .rename(columns={"site": "site_num"})
        )

        weather = self.registry.load_exact_weather_windows().rename(columns={"site": "site_num"})
        comp_static = self.build_compartment_state_cube()[
            [
                "site_num",
                "compartment_norm",
                "site_label",
                "site_iri",
                "compartment_iri",
                "transect_zone",
                "habitat_regime",
                "lat",
                "lon",
                "elevation_m",
                "map_mm",
                "map_era5_mm",
                "mat_c",
                "h3_5",
                "h3_6",
                "era5_cell",
                "xrf_evaporite_index",
                "ph_mean",
                "shannon_16s",
                "clr_Ca",
                "clr_Mg",
                "clr_Na",
                "clr_S",
                "clr_Cl",
                "clr_Si",
                "clr_Fe",
                "clr_Al",
                "clr_K",
                "clr_P",
                "clr_Ti",
            ]
        ]

        panel = visit_agg.merge(weather, on=["site_num", "trip"], how="left").merge(
            comp_static, on=["site_num", "compartment_norm"], how="left"
        )
        return panel

    def build_era5_h3_grid(self, step_deg: float = 0.25) -> pd.DataFrame:
        """Build the continuous 0.25-degree Open-Meteo/ERA5 grid across the Empty Quarter corridor."""
        sites_df = self.build_sites_cube()
        primary = sites_df[sites_df["is_primary_transect"]].copy()
        aeolian = primary[primary["habitat_regime"] == "aeolian_desert_soil"].copy()

        lats = np.arange(19.0, 21.5 + 1e-6, step_deg)
        lons = np.arange(45.25, 55.25 + 1e-6, step_deg)

        site_lats = primary["lat"].to_numpy()
        site_lons = primary["lon"].to_numpy()
        site_maps = primary["map_mm"].to_numpy()
        site_mats = primary["mat_c"].to_numpy()
        site_elevs = primary["elevation_m"].to_numpy()

        # Empirical support bounds across the 58 non-aquatic aeolian transect sites
        map_min, map_max = float(aeolian["map_mm"].min()), float(aeolian["map_mm"].max())
        elev_min, elev_max = float(aeolian["elevation_m"].min()), float(aeolian["elevation_m"].max())

        records: List[Dict[str, Any]] = []
        for lat in lats:
            for lon in lons:
                dist_deg = np.sqrt((site_lats - lat) ** 2 + (site_lons - lon) ** 2)
                min_dist_km = float(dist_deg.min() * 111.0)
                # Keep cells within 85 km of the 1,043 km sampling corridor
                if min_dist_km > 85.0:
                    continue
                weights = 1.0 / (dist_deg**2 + 0.02)
                weights /= weights.sum()
                elev_val = round(float(np.dot(weights, site_elevs)), 1)
                map_val = round(float(np.dot(weights, site_maps)), 2)
                mat_val = round(float(np.dot(weights, site_mats)), 2)
                # Domain-of-support flag: within 45 km of observed transect AND within observed climate/elevation support
                in_spatial_support = bool(
                    (min_dist_km <= 45.0)
                    and (map_min <= map_val <= map_max)
                    and (elev_min <= elev_val <= elev_max)
                )
                zone_name = self.assign_transect_zone(float(lon))
                nearest_idx = int(np.argmin(dist_deg))
                nearest_site_num = int(primary.iloc[nearest_idx]["site_num"])
                if not in_spatial_support:
                    prescribed_portfolio = "ABSTAIN"
                elif nearest_site_num in (59, 60):
                    prescribed_portfolio = "P4"
                elif zone_name == "western_plateau" and map_val >= 32.0:
                    prescribed_portfolio = "P2"
                else:
                    prescribed_portfolio = "P3"
                records.append(
                    {
                        "era5_cell": f"{lat:.2f}_{lon:.2f}",
                        "lat": float(lat),
                        "lon": float(lon),
                        "h3_5": _latlng_to_h3(lat, lon, 5),
                        "h3_6": _latlng_to_h3(lat, lon, 6),
                        "corridor_distance_km": round(min_dist_km, 2),
                        "nearest_primary_site": nearest_site_num,
                        "elevation_m": elev_val,
                        "map_mm": map_val,
                        "mat_c": mat_val,
                        "transect_zone": zone_name,
                        "in_empirical_support": in_spatial_support,
                        "prescribed_portfolio": prescribed_portfolio,
                        "abstention_reason": (
                            "supported"
                            if in_spatial_support
                            else "corridor_margin_extrapolation_gt_45km"
                        ),
                    }
                )
        return pd.DataFrame(records)

    def export_artifacts(self, output_dir: Path) -> Dict[str, Path]:
        """Materialize the state cubes to Parquet in `output_dir`."""
        output_dir.mkdir(parents=True, exist_ok=True)
        sites_df = self.build_sites_cube()
        comp_df = self.build_compartment_state_cube()
        visit_df = self.build_visit_compartment_cube()
        grid_df = self.build_era5_h3_grid()

        sites_parquet = output_dir / "eq_sites_cube.parquet"
        comp_parquet = output_dir / "eq_compartment_state_cube.parquet"
        visit_parquet = output_dir / "eq_visit_compartment_cube.parquet"
        grid_parquet = output_dir / "eq_era5_h3_corridor_grid.parquet"

        sites_df.to_parquet(sites_parquet, index=False)
        comp_df.to_parquet(comp_parquet, index=False)
        visit_df.to_parquet(visit_parquet, index=False)
        grid_df.to_parquet(grid_parquet, index=False)

        return {
            "sites_cube": sites_parquet,
            "compartment_state_cube": comp_parquet,
            "visit_compartment_cube": visit_parquet,
            "era5_h3_grid": grid_parquet,
        }
