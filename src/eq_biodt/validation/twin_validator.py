"""Five-Tier Quantitative Validation & Benchmarking Suite for EQ-BioDT.

Computes out-of-sample predictive accuracy, site-clustered bootstrap confidence intervals,
Shrinkage-MinT / Nonnegative-MinT hierarchical reconciliation, space-for-time transport tests,
strictly time-ordered Trip 5 prospective forecasting, and natural-experiment calibration
across the canonical Empty Quarter datasets:

- Tier 1: 5-Fold Spatially-Blocked Out-of-Sample Cross-Validation (contiguous longitude blocks,
  zero site overlap across folds, train-only median imputation, observed-only target scoring,
  site-clustered bootstrap 95% CIs, and compartment/climate/depth baseline comparisons).
- Tier 2: Longitudinal Validation across Trips 1..5 (March 2023 - October 2025):
  (a) Leave-One-Campaign-Out (LOCO) CV across Trips 1..5,
  (b) Space-for-Time Transport Test comparing between-site spatial aridity slopes against
      within-site 90-day and 365-day precipitation anomaly slopes (ERA5 vs NASA POWER), and
  (c) Strictly Time-Ordered Prospective Forecasting (Train on Trips 1..4 -> Forecast Trip 5
      for both absolute levels and within-site temporal anomalies Delta_Trip5 against
      site-history persistence, environmental ridge, random forest, and gradient boosting).
- Tier 3: OWL-Compiled Polyhedron & Hierarchical Reconciliation Benchmark:
  Compares Unconstrained Base Predictions, Ordinary Euclidean Projection (W = I),
  Diagonal WLS (W = diag(sigma_i^2)), Bottom-Up Nonnegative Leaf Sum, Shrinkage-MinT,
  Nonnegative Shrinkage-MinT (Ledoit-Wolf covariance W_LW estimated strictly from inner
  training-fold base residuals), and Matched-Feature Flat Ridge.
- Tier 4: Paired Natural-Experiment Calibration & Mechanistic Decomposition:
  - n=106 paired Root vs Bulk visits across 47 sites (`h1_visit_log2fc_root_vs_bulk.tsv`)
  - Genus-turnover (97.4%) vs within-genus prevalence decomposition (`A2_root_bulk_decomposition_summary.tsv`)
  - Group 1l H2-oxidation vs CBB RuBisCO co-carriage (`101/558` MAGs, `A4_cbb_cocarriage_by_group.tsv`)
  - n=68 paired Surface vs Subsurface visits across 36 sites (`h1_visit_log2fc_surface_vs_bulk.tsv`)
  - n=51 observed bulk-soil site aridity gradient (`H3_spatial_by_marker_position.tsv`)
- Tier 5: Synthetic Metagenome Read-to-Guild Benchmark (`mock_ge_accuracy.tsv`),
  72-MAG 10,584-run GSMM Audit (`metabolic_pilot/REPORT.txt`), Hard Maintenance Feasibility
  Invariants, and Common-Baseline 4x2x2 Factorial Brine Analysis.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Dict, List, Optional

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.ensemble import HistGradientBoostingRegressor, RandomForestRegressor

from eq_biodt.data.build_twin_state_cube import TwinStateCubeBuilder
from eq_biodt.data.manifest_registry import CanonicalDataRegistry
from eq_biodt.interventions.desalination_brine import DesalinationBrineIntervention
from eq_biodt.interventions.pareto_optimizer import ParetoPortfolioOptimizer
from eq_biodt.interventions.refugium_thermal import RefugiumThermalIntervention
from eq_biodt.interventions.revegetation import RevegetationIntervention
from eq_biodt.models.carbon_budget_engine import CarbonBudgetEngine
from eq_biodt.models.causal_ontological_twin import (
    CausalOntologicalTwinEngine,
    MinTReconciler,
    OntologicalPolytope,
    PartitionSumConstraint,
    SubsumptionRelation,
)
from eq_biodt.models.metabolic_community_engine import (
    CommunityMetabolicInputs,
    MetabolicCommunityEngine,
)


def _impute_train_test(
    X_train: np.ndarray, X_test: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """Impute missing values using strictly training-fold column medians."""
    X_tr = np.asarray(X_train, dtype=float).copy()
    X_te = np.asarray(X_test, dtype=float).copy()
    col_med = np.nanmedian(X_tr, axis=0)
    col_med = np.where(np.isfinite(col_med), col_med, 0.0)
    for j in range(X_tr.shape[1]):
        tr_nan = ~np.isfinite(X_tr[:, j])
        te_nan = ~np.isfinite(X_te[:, j])
        if np.any(tr_nan):
            X_tr[tr_nan, j] = col_med[j]
        if np.any(te_nan):
            X_te[te_nan, j] = col_med[j]
    return X_tr, X_te


def _ridge_fit_predict(
    X_train: np.ndarray,
    y_train: np.ndarray,
    X_test: np.ndarray,
    alpha: float = 5.0,
) -> np.ndarray:
    """Closed-form standardized ridge regression with strictly training-fold median imputation."""
    X_tr, X_te = _impute_train_test(X_train, X_test)
    y_tr = np.asarray(y_train, dtype=float)

    mu = X_tr.mean(axis=0, keepdims=True)
    sd = X_tr.std(axis=0, keepdims=True) + 1e-8
    Z_tr = (X_tr - mu) / sd
    Z_te = (X_te - mu) / sd

    Z_tr_b = np.column_stack([np.ones(Z_tr.shape[0]), Z_tr])
    Z_te_b = np.column_stack([np.ones(Z_te.shape[0]), Z_te])

    p = Z_tr_b.shape[1]
    reg = alpha * np.eye(p)
    reg[0, 0] = 0.0  # do not penalize intercept
    w = np.linalg.solve(Z_tr_b.T @ Z_tr_b + reg, Z_tr_b.T @ y_tr)
    return Z_te_b @ w


def _regression_metrics(
    y_true: np.ndarray,
    y_pred: np.ndarray,
    site_ids: Optional[np.ndarray] = None,
    n_boot: int = 300,
    seed: int = 42,
) -> Dict[str, Any]:
    """Compute RMSE, MAE, R^2, Pearson r, Spearman rho, and site-clustered bootstrap 95% CIs."""
    yt = np.asarray(y_true, dtype=float)
    yp = np.asarray(y_pred, dtype=float)
    if not np.all(np.isfinite(yp)):
        raise ValueError("Out-of-fold predictions contain NaN or infinite values.")
    mask = np.isfinite(yt) & np.isfinite(yp)
    yt = yt[mask]
    yp = yp[mask]

    err = yp - yt
    mse = float(np.mean(err**2))
    rmse = float(np.sqrt(mse))
    mae = float(np.mean(np.abs(err)))
    ss_tot = float(np.sum((yt - np.mean(yt)) ** 2))
    r2 = float(1.0 - np.sum(err**2) / max(1e-12, ss_tot))

    if len(yt) > 2 and np.std(yt) > 1e-12 and np.std(yp) > 1e-12:
        pr, pp = stats.pearsonr(yt, yp)
        sr, sp = stats.spearmanr(yt, yp)
    else:
        pr, pp, sr, sp = 0.0, 1.0, 0.0, 1.0

    out: Dict[str, Any] = {
        "n": int(len(yt)),
        "rmse": round(rmse, 4),
        "mae": round(mae, 4),
        "r2": round(r2, 4),
        "pearson_r": round(float(pr), 4),
        "pearson_p": float(pp),
        "spearman_rho": round(float(sr), 4),
        "spearman_p": float(sp),
    }

    if site_ids is not None:
        s_arr = np.asarray(site_ids)[mask]
        uniq_s = np.unique(s_arr)
        out["n_clusters_sites"] = int(len(uniq_s))
        if len(uniq_s) >= 5 and n_boot > 0:
            rng = np.random.default_rng(seed)
            idx_by_site = {s: np.where(s_arr == s)[0] for s in uniq_s}
            boot_rho: List[float] = []
            boot_r2: List[float] = []
            for _ in range(n_boot):
                sampled_sites = rng.choice(uniq_s, size=len(uniq_s), replace=True)
                b_idx = np.concatenate([idx_by_site[s] for s in sampled_sites])
                yt_b = yt[b_idx]
                yp_b = yp[b_idx]
                if np.std(yt_b) > 1e-12 and np.std(yp_b) > 1e-12:
                    boot_rho.append(float(stats.spearmanr(yt_b, yp_b)[0]))
                    ss_b = float(np.sum((yt_b - np.mean(yt_b)) ** 2))
                    boot_r2.append(
                        float(1.0 - np.sum((yp_b - yt_b) ** 2) / max(1e-12, ss_b))
                    )
            if boot_rho:
                out["spearman_rho_cluster_ci95"] = [
                    round(float(np.percentile(boot_rho, 2.5)), 4),
                    round(float(np.percentile(boot_rho, 97.5)), 4),
                ]
                out["r2_cluster_ci95"] = [
                    round(float(np.percentile(boot_r2, 2.5)), 4),
                    round(float(np.percentile(boot_r2, 97.5)), 4),
                ]
    return out


class TwinValidationSuite:
    """Executes all 5 quantitative validation tiers for the Empty Quarter Digital Twin."""

    def __init__(self, builder: Optional[TwinStateCubeBuilder] = None) -> None:
        self.builder = builder or TwinStateCubeBuilder()
        self.registry: CanonicalDataRegistry = self.builder.registry
        self.paths = self.registry.paths
        self.causal_engine = CausalOntologicalTwinEngine()

    def build_sample_level_validation_panel(self, id_floor: int = 60) -> pd.DataFrame:
        """Merge the 551 WGS soil samples with site coordinates, exact weather windows, XRF, pH, guild GEs, and mobilome."""
        wgs = self.registry.load_wgs_sample_metadata().copy()
        wgs["compartment_norm"] = wgs["compartment"].map(self.builder.COMPARTMENT_NORMALIZATION)
        if wgs["compartment_norm"].isna().any():
            raise ValueError("Unmapped compartment labels in WGS metadata.")

        guild = self.registry.load_guild_genome_equivalents(id_floor=id_floor)
        mobilome = self.registry.load_mobilome_load()
        sites_cube = self.builder.build_sites_cube()
        env_trip = self.registry.load_env_site_trip()
        comp_cube = self.builder.build_compartment_state_cube()

        site_climate = sites_cube[
            [
                "site_num",
                "lat",
                "lon",
                "elevation_m",
                "map_mm",
                "map_era5_mm",
                "mat_c",
                "rainy_days_per_yr",
                "h3_5",
                "transect_zone",
            ]
        ].rename(
            columns={
                "site_num": "site",
                "map_mm": "map_bestmatch_mm",
                "mat_c": "mat_era5_c",
            }
        )

        trip_cols = [
            "site",
            "trip",
            "collection_date",
            "date_substitution_days",
            "ndvi_W60d",
            "evi_W60d",
            "smap_sm_W30d",
            "P7_era5_exact",
            "P14_era5_exact",
            "P30_era5_exact",
            "P90_era5_exact",
            "P365_era5_exact",
            "P90_power_exact",
            "P365_power_exact",
            "P90_era5",
            "P365_era5",
        ]
        avail_trip_cols = [c for c in trip_cols if c in env_trip.columns]

        df = (
            wgs.merge(guild, on="sample", how="inner")
            .merge(mobilome, on="sample", how="left")
            .merge(site_climate, on="site", how="left")
            .merge(env_trip[avail_trip_cols], on=["site", "trip"], how="left")
            .merge(
                comp_cube[
                    [
                        "site_num",
                        "compartment_norm",
                        "xrf_evaporite_index",
                        "ph_mean",
                        "shannon_16s",
                        "dna_conc_mean",
                        "clr_Si",
                        "clr_Ca",
                        "clr_Mg",
                        "clr_Fe",
                        "clr_Al",
                        "clr_K",
                        "clr_Na",
                        "clr_S",
                        "clr_Cl",
                        "clr_P",
                        "clr_Ti",
                    ]
                ].rename(columns={"site_num": "site"}),
                on=["site", "compartment_norm"],
                how="left",
            )
        )

        if df["lon"].isna().any():
            raise ValueError("Missing longitude coordinates after merging 70-site spatial table.")

        for gcol in (
            "hydB_1l",
            "hydB_1h",
            "hydB_2a",
            "coxL_formI",
            "viral_richness",
            "plasmid_richness",
        ):
            df[gcol] = pd.to_numeric(df[gcol], errors="coerce")

        df["guild_observed"] = (
            df["hydB_1l"].notna()
            & df["hydB_1h"].notna()
            & df["hydB_2a"].notna()
            & df["coxL_formI"].notna()
        )
        df["mobilome_observed"] = df["viral_richness"].notna() & df["plasmid_richness"].notna()

        df["is_surface"] = (df["compartment_norm"] == "Surface").astype(float)
        df["is_deep"] = (df["compartment_norm"] == "ShallowSubsurface").astype(float)
        df["is_root"] = (df["compartment_norm"] == "RootAdjacent").astype(float)

        pc = 0.01
        for gcol in ("hydB_1l", "hydB_1h", "hydB_2a", "coxL_formI"):
            vals = df[gcol].to_numpy(dtype=float)
            df[f"log2_{gcol}"] = np.where(
                np.isfinite(vals), np.log2(np.maximum(0.0, vals) + pc), np.nan
            )

        df["total_h2ase_ge"] = np.where(
            df["guild_observed"],
            np.maximum(0.0, df["hydB_1l"])
            + np.maximum(0.0, df["hydB_1h"])
            + np.maximum(0.0, df["hydB_2a"]),
            np.nan,
        )
        df["total_trace_gas_scavenger_ge"] = np.where(
            df["guild_observed"],
            df["total_h2ase_ge"] + np.maximum(0.0, df["coxL_formI"]),
            np.nan,
        )
        df["log2_total_h2ase_ge"] = np.where(
            df["guild_observed"],
            np.log2(df["total_h2ase_ge"] + pc),
            np.nan,
        )
        df["log1p_dna_conc"] = np.log1p(
            np.maximum(0.0, pd.to_numeric(df["dna_conc_mean"], errors="coerce"))
        )
        df["log1p_scg_rpk"] = np.log1p(
            np.maximum(0.0, pd.to_numeric(df["scg_median_rpk"], errors="coerce"))
        )
        df["scg_genes_detected"] = pd.to_numeric(df["scg_genes_detected"], errors="coerce")
        df["log1p_viral_richness"] = np.where(
            df["mobilome_observed"],
            np.log1p(np.maximum(0.0, df["viral_richness"])),
            np.nan,
        )
        df["log1p_plasmid_richness"] = np.where(
            df["mobilome_observed"],
            np.log1p(np.maximum(0.0, df["plasmid_richness"])),
            np.nan,
        )
        return df

    def run_tier1_spatial_block_cv(self, n_folds: int = 5) -> Dict[str, Any]:
        """Tier 1: 5-Fold Spatially-Blocked Out-of-Sample Cross-Validation (zero site leakage)."""
        df = self.build_sample_level_validation_panel(id_floor=60)
        comp_cube = self.builder.build_compartment_state_cube().copy()

        unique_sites = df[["site", "lon"]].drop_duplicates().sort_values(["lon", "site"])
        unique_sites["spatial_fold"] = np.floor(
            np.linspace(0, n_folds - 1e-6, len(unique_sites))
        ).astype(int)
        fold_map = dict(zip(unique_sites["site"], unique_sites["spatial_fold"]))
        df["spatial_fold"] = df["site"].map(fold_map)
        if df["spatial_fold"].isna().any():
            raise ValueError("Unassigned spatial folds in WGS panel.")

        feature_cols = [
            "is_surface",
            "is_deep",
            "is_root",
            "map_bestmatch_mm",
            "mat_era5_c",
            "ndvi_W60d",
            "xrf_evaporite_index",
            "ph_mean",
            "log1p_dna_conc",
            "clr_Si",
            "clr_Ca",
            "clr_Fe",
        ]
        comp_only_cols = ["is_surface", "is_deep", "is_root"]
        clim_only_cols = [
            "map_bestmatch_mm",
            "map_era5_mm",
            "mat_era5_c",
            "elevation_m",
            "P365_era5_exact",
        ]

        guild_targets = [
            "log2_hydB_1l",
            "log2_hydB_1h",
            "log2_hydB_2a",
            "log2_coxL_formI",
            "log2_total_h2ase_ge",
        ]
        df_guild = df[df["guild_observed"]].reset_index(drop=True)
        oof_preds: Dict[str, np.ndarray] = {
            t: np.full(len(df_guild), np.nan) for t in guild_targets
        }
        oof_comp_base = np.full(len(df_guild), np.nan)
        oof_clim_base = np.full(len(df_guild), np.nan)

        for f in range(n_folds):
            tr_mask = (df_guild["spatial_fold"] != f).to_numpy()
            te_mask = (df_guild["spatial_fold"] == f).to_numpy()
            X_tr = df_guild.loc[tr_mask, feature_cols].to_numpy(dtype=float)
            X_te = df_guild.loc[te_mask, feature_cols].to_numpy(dtype=float)
            for t in guild_targets:
                y_tr = df_guild.loc[tr_mask, t].to_numpy(dtype=float)
                oof_preds[t][te_mask] = _ridge_fit_predict(X_tr, y_tr, X_te, alpha=100.0)
            y_tr_1l = df_guild.loc[tr_mask, "log2_hydB_1l"].to_numpy(dtype=float)
            oof_comp_base[te_mask] = _ridge_fit_predict(
                df_guild.loc[tr_mask, comp_only_cols].to_numpy(dtype=float),
                y_tr_1l,
                df_guild.loc[te_mask, comp_only_cols].to_numpy(dtype=float),
                alpha=1.0,
            )
            oof_clim_base[te_mask] = _ridge_fit_predict(
                df_guild.loc[tr_mask, clim_only_cols].to_numpy(dtype=float),
                y_tr_1l,
                df_guild.loc[te_mask, clim_only_cols].to_numpy(dtype=float),
                alpha=25.0,
            )

        results: Dict[str, Any] = {}
        guild_sites = df_guild["site"].to_numpy()
        for t in guild_targets:
            results[t] = _regression_metrics(
                df_guild[t].to_numpy(dtype=float),
                oof_preds[t],
                site_ids=guild_sites,
            )

        df_mob = df[df["mobilome_observed"] & df["guild_observed"]].reset_index(drop=True)
        mob_sites = df_mob["site"].to_numpy()
        mob_features = feature_cols + [
            "log2_hydB_1l",
            "log2_hydB_1h",
            "log2_coxL_formI",
            "shannon_16s",
            "log1p_scg_rpk",
            "scg_genes_detected",
        ]
        depth_only_cols = ["log1p_scg_rpk", "scg_genes_detected"]
        oof_depth_viral = np.full(len(df_mob), np.nan)

        for mt in ("log1p_viral_richness", "log1p_plasmid_richness"):
            pred_arr = np.full(len(df_mob), np.nan)
            for f in range(n_folds):
                tr_mask = (df_mob["spatial_fold"] != f).to_numpy()
                te_mask = (df_mob["spatial_fold"] == f).to_numpy()
                X_tr = df_mob.loc[tr_mask, mob_features].to_numpy(dtype=float)
                X_te = df_mob.loc[te_mask, mob_features].to_numpy(dtype=float)
                y_tr = df_mob.loc[tr_mask, mt].to_numpy(dtype=float)
                pred_arr[te_mask] = _ridge_fit_predict(X_tr, y_tr, X_te, alpha=30.0)
                if mt == "log1p_viral_richness":
                    oof_depth_viral[te_mask] = _ridge_fit_predict(
                        df_mob.loc[tr_mask, depth_only_cols].to_numpy(dtype=float),
                        y_tr,
                        df_mob.loc[te_mask, depth_only_cols].to_numpy(dtype=float),
                        alpha=5.0,
                    )
            results[mt] = _regression_metrics(
                df_mob[mt].to_numpy(dtype=float),
                pred_arr,
                site_ids=mob_sites,
            )

        comp_Cube_clean = comp_cube.dropna(subset=["shannon_16s"]).copy().reset_index(drop=True)
        cube_sites = (
            comp_Cube_clean[["site_num", "lon"]]
            .drop_duplicates()
            .sort_values(["lon", "site_num"])
        )
        cube_sites["spatial_fold"] = np.floor(
            np.linspace(0, n_folds - 1e-6, len(cube_sites))
        ).astype(int)
        cube_fold_map = dict(zip(cube_sites["site_num"], cube_sites["spatial_fold"]))
        comp_Cube_clean["spatial_fold"] = comp_Cube_clean["site_num"].map(cube_fold_map)

        comp_Cube_clean["is_surface"] = (
            comp_Cube_clean["compartment_norm"] == "Surface"
        ).astype(float)
        comp_Cube_clean["is_deep"] = (
            comp_Cube_clean["compartment_norm"] == "ShallowSubsurface"
        ).astype(float)
        comp_Cube_clean["is_root"] = (
            comp_Cube_clean["compartment_norm"] == "RootAdjacent"
        ).astype(float)
        comp_Cube_clean["log1p_dna_conc"] = np.log1p(
            np.maximum(0.0, pd.to_numeric(comp_Cube_clean["dna_conc_mean"], errors="coerce"))
        )
        shannon_feats = [
            "is_surface",
            "is_deep",
            "is_root",
            "map_mm",
            "mat_c",
            "ph_mean",
            "xrf_evaporite_index",
            "log1p_dna_conc",
            "clr_Si",
            "clr_Ca",
            "clr_Mg",
            "clr_Fe",
            "clr_K",
        ]
        sh_oof = np.full(len(comp_Cube_clean), np.nan)
        for f in range(n_folds):
            tr = (comp_Cube_clean["spatial_fold"] != f).to_numpy()
            te = (comp_Cube_clean["spatial_fold"] == f).to_numpy()
            X_tr = comp_Cube_clean.loc[tr, shannon_feats].to_numpy(dtype=float)
            X_te = comp_Cube_clean.loc[te, shannon_feats].to_numpy(dtype=float)
            y_tr = comp_Cube_clean.loc[tr, "shannon_16s"].to_numpy(dtype=float)
            sh_oof[te] = _ridge_fit_predict(X_tr, y_tr, X_te, alpha=15.0)
        results["shannon_16s_site_compartment"] = _regression_metrics(
            comp_Cube_clean["shannon_16s"].to_numpy(dtype=float),
            sh_oof,
            site_ids=comp_Cube_clean["site_num"].to_numpy(),
        )

        baselines = {
            "log2_hydB_1l_compartment_only": _regression_metrics(
                df_guild["log2_hydB_1l"].to_numpy(dtype=float),
                oof_comp_base,
                site_ids=guild_sites,
            ),
            "log2_hydB_1l_climate_only": _regression_metrics(
                df_guild["log2_hydB_1l"].to_numpy(dtype=float),
                oof_clim_base,
                site_ids=guild_sites,
            ),
            "log1p_viral_richness_sequencing_depth_only": _regression_metrics(
                df_mob["log1p_viral_richness"].to_numpy(dtype=float),
                oof_depth_viral,
                site_ids=mob_sites,
            ),
        }

        return {
            "validation_scheme": f"{n_folds}-Fold Spatial Block Cross-Validation (zero site overlap across folds, train-only median imputation)",
            "total_wgs_samples": int(len(df)),
            "observed_guild_wgs_samples": int(len(df_guild)),
            "observed_mobilome_wgs_samples": int(len(df_mob)),
            "total_site_compartments": int(len(comp_Cube_clean)),
            "targets": results,
            "baseline_comparisons": baselines,
        }

    def _run_space_for_time_and_trip5_forecast(
        self, df_guild: pd.DataFrame
    ) -> Dict[str, Any]:
        """Compute Space-for-Time between- vs within-site slopes and strictly time-ordered Trip 5 forecasts."""
        h123_dir = self.paths.h123_v2_dir
        wb_df = pd.read_csv(h123_dir / "H3_temporal_within_between.tsv", sep="\t")
        ag_df = pd.read_csv(h123_dir / "H3_temporal_product_agreement.tsv", sep="\t")

        # Extract headline between-site vs within-site slopes for hydB_1l across Deep, Surface, PlantRoot
        wb_1l = wb_df[
            (wb_df["marker"] == "hydB_1l")
            & (wb_df["exposure"].isin(["P90_era5", "P365_era5", "P90_power", "P365_power"]))
        ][
            [
                "compartment",
                "exposure",
                "n_visits",
                "n_sites",
                "between_slope",
                "within_slope",
                "within_se",
                "within_LRT_p",
                "within_boot_lo",
                "within_boot_hi",
                "q_BH_family",
            ]
        ].to_dict(orient="records")

        # Recompute exact collection-date within-site vs between-site slope on Deep/Surface/RootAdjacent
        visit_panel = self.builder.build_visit_compartment_cube()
        visit_panel["log2_hydB_1l"] = np.log2(np.maximum(0.0, visit_panel["hydB_1l_ge"]) + 0.01)
        visit_panel["log1p_P90_era5_exact"] = np.log1p(
            np.maximum(0.0, visit_panel["P90_era5_exact"])
        )
        visit_panel["log1p_P365_era5_exact"] = np.log1p(
            np.maximum(0.0, visit_panel["P365_era5_exact"])
        )

        exact_date_wb: List[Dict[str, Any]] = []
        for comp_name in ("ShallowSubsurface", "Surface", "RootAdjacent"):
            sub = visit_panel[visit_panel["compartment_norm"] == comp_name].dropna(
                subset=["log2_hydB_1l", "log1p_P90_era5_exact", "log1p_P365_era5_exact"]
            )
            site_counts = sub.groupby("site_num")["trip"].nunique()
            multi_sites = site_counts[site_counts >= 2].index
            sub_m = sub[sub["site_num"].isin(multi_sites)].copy()
            if len(sub_m) < 10:
                continue
            for exp_col in ("log1p_P90_era5_exact", "log1p_P365_era5_exact"):
                site_mean_x = sub_m.groupby("site_num")[exp_col].transform("mean")
                site_mean_y = sub_m.groupby("site_num")["log2_hydB_1l"].transform("mean")
                dx = sub_m[exp_col] - site_mean_x
                dy = sub_m["log2_hydB_1l"] - site_mean_y
                within_beta = float(np.dot(dx, dy) / max(1e-9, np.dot(dx, dx)))
                s_agg = sub_m.groupby("site_num", as_index=False).agg(
                    xm=(exp_col, "mean"), ym=("log2_hydB_1l", "mean")
                )
                bx = s_agg["xm"] - s_agg["xm"].mean()
                by = s_agg["ym"] - s_agg["ym"].mean()
                between_beta = float(np.dot(bx, by) / max(1e-9, np.dot(bx, bx)))
                exact_date_wb.append(
                    {
                        "compartment": comp_name,
                        "exposure": exp_col,
                        "n_visits": int(len(sub_m)),
                        "n_sites": int(sub_m["site_num"].nunique()),
                        "between_site_slope": round(between_beta, 4),
                        "within_site_slope_exact_dates": round(within_beta, 4),
                    }
                )

        # Strictly Time-Ordered Prospective Forecasting Benchmark: Train on Trips 1..4 -> Test on Trip 5
        # Aggregate to (site, trip, compartment_norm) visit-compartment observations
        guild_cols = ["hydB_1l", "hydB_1h", "hydB_2a", "coxL_formI"]
        env_cols = [
            "is_surface",
            "is_deep",
            "is_root",
            "map_bestmatch_mm",
            "mat_era5_c",
            "ndvi_W60d",
            "smap_sm_W30d",
            "P90_era5_exact",
            "P365_era5_exact",
            "xrf_evaporite_index",
            "ph_mean",
            "log1p_dna_conc",
            "clr_Si",
            "clr_Ca",
            "clr_Fe",
        ]
        vc_df = (
            df_guild.groupby(["site", "trip", "compartment_norm"], as_index=False)
            .agg(
                {
                    "log2_hydB_1l": "mean",
                    "log2_hydB_1h": "mean",
                    "log2_hydB_2a": "mean",
                    "log2_coxL_formI": "mean",
                    "log2_total_h2ase_ge": "mean",
                    "hydB_1l": "mean",
                    "hydB_1h": "mean",
                    "hydB_2a": "mean",
                    "total_h2ase_ge": "mean",
                    "coxL_formI": "mean",
                    "total_trace_gas_scavenger_ge": "mean",
                    **{c: "mean" for c in env_cols},
                }
            )
        )

        train_vc = vc_df[vc_df["trip"].isin(["trip1", "trip2", "trip3", "trip4"])].copy().reset_index(drop=True)
        test_vc = vc_df[vc_df["trip"] == "trip5"].copy().reset_index(drop=True)

        # Compute frozen Trips 1..4 site-compartment historical means (strictly no Trip 5 leakage)
        hist_means = (
            train_vc.groupby(["site", "compartment_norm"], as_index=False)
            .agg(
                hist_log2_1l=("log2_hydB_1l", "mean"),
                hist_log2_h2tot=("log2_total_h2ase_ge", "mean"),
                hist_P90=("P90_era5_exact", "mean"),
                hist_ndvi=("ndvi_W60d", "mean"),
            )
        )
        comp_fallback = train_vc.groupby("compartment_norm")["log2_hydB_1l"].mean().to_dict()
        comp_fallback_h2 = train_vc.groupby("compartment_norm")["log2_total_h2ase_ge"].mean().to_dict()

        # Leave-one-trip-out historical mean inside Trips 1..4 for training the longitudinal twin
        train_hist_1l = np.zeros(len(train_vc), dtype=float)
        train_dP90 = np.zeros(len(train_vc), dtype=float)
        for idx, row in train_vc.iterrows():
            other = train_vc[
                (train_vc["site"] == row["site"])
                & (train_vc["compartment_norm"] == row["compartment_norm"])
                & (train_vc["trip"] != row["trip"])
            ]
            if len(other) > 0:
                train_hist_1l[idx] = float(other["log2_hydB_1l"].mean())
                train_dP90[idx] = float(row["P90_era5_exact"] - other["P90_era5_exact"].mean())
            else:
                train_hist_1l[idx] = float(comp_fallback.get(row["compartment_norm"], -1.0))
                train_dP90[idx] = 0.0
        train_vc["hist_log2_1l"] = train_hist_1l
        train_vc["delta_P90_era5"] = train_dP90

        test_vc = test_vc.merge(hist_means, on=["site", "compartment_norm"], how="left")
        has_prior_history = test_vc["hist_log2_1l"].notna().to_numpy()
        test_vc["hist_log2_1l"] = test_vc.apply(
            lambda r: r["hist_log2_1l"]
            if pd.notna(r["hist_log2_1l"])
            else comp_fallback.get(r["compartment_norm"], -1.0),
            axis=1,
        )
        test_vc["delta_P90_era5"] = (
            test_vc["P90_era5_exact"] - test_vc["hist_P90"].fillna(test_vc["P90_era5_exact"])
        )

        # 1) Persistence baseline (Frozen Trips 1..4 site-compartment mean)
        pred_persist_1l = test_vc["hist_log2_1l"].to_numpy(dtype=float)

        # 2) Environmental Ridge (spatial + weather features only, no site history)
        X_tr_env = train_vc[env_cols].to_numpy(dtype=float)
        X_te_env = test_vc[env_cols].to_numpy(dtype=float)
        y_tr_1l = train_vc["log2_hydB_1l"].to_numpy(dtype=float)
        y_te_1l = test_vc["log2_hydB_1l"].to_numpy(dtype=float)
        pred_env_ridge_1l = _ridge_fit_predict(X_tr_env, y_tr_1l, X_te_env, alpha=60.0)

        # 3) Random Forest & HistGradientBoosting baselines (environmental + site history)
        long_cols = env_cols + ["hist_log2_1l", "delta_P90_era5"]
        X_tr_long, X_te_long = _impute_train_test(
            train_vc[long_cols].to_numpy(dtype=float),
            test_vc[long_cols].to_numpy(dtype=float),
        )
        rf = RandomForestRegressor(n_estimators=200, max_depth=6, random_state=42)
        rf.fit(X_tr_long, y_tr_1l)
        pred_rf_1l = rf.predict(X_te_long)

        hgb = HistGradientBoostingRegressor(max_iter=100, max_depth=4, l2_regularization=5.0, random_state=42)
        hgb.fit(X_tr_long, y_tr_1l)
        pred_hgb_1l = hgb.predict(X_te_long)

        # 4) Longitudinal Causal-Ontological Twin (anchored site history + compartment-specific rain anomaly + ridge)
        pred_twin_1l = _ridge_fit_predict(X_tr_long, y_tr_1l, X_te_long, alpha=45.0)

        # Evaluate both Trip 5 absolute levels (n=89) and within-site temporal change Delta_Trip5 (n=85 with prior history)
        te_sites = test_vc["site"].to_numpy()
        y_te_delta = (y_te_1l - pred_persist_1l)[has_prior_history]
        s_te_delta = te_sites[has_prior_history]

        return {
            "space_for_time_transport": {
                "h123_within_vs_between_hydB_1l": wb_1l,
                "exact_collection_date_within_vs_between": exact_date_wb,
                "era5_vs_nasa_power_agreement": ag_df.where(pd.notna(ag_df), None).to_dict(
                    orient="records"
                ),
            },
            "trip5_prospective_forecast_benchmark": {
                "train_campaigns": ["trip1", "trip2", "trip3", "trip4"],
                "test_campaign": "trip5 (12-15 October 2025)",
                "n_train_visit_compartments": int(len(train_vc)),
                "n_test_visit_compartments": int(len(test_vc)),
                "n_test_with_prior_site_compartment_history": int(has_prior_history.sum()),
                "level_forecasts_log2_hydB_1l": {
                    "longitudinal_twin": _regression_metrics(
                        y_te_1l, pred_twin_1l, site_ids=te_sites
                    ),
                    "site_compartment_persistence": _regression_metrics(
                        y_te_1l, pred_persist_1l, site_ids=te_sites
                    ),
                    "environmental_ridge_spatial_only": _regression_metrics(
                        y_te_1l, pred_env_ridge_1l, site_ids=te_sites
                    ),
                    "random_forest_baseline": _regression_metrics(
                        y_te_1l, pred_rf_1l, site_ids=te_sites
                    ),
                    "gradient_boosting_baseline": _regression_metrics(
                        y_te_1l, pred_hgb_1l, site_ids=te_sites
                    ),
                },
                "temporal_anomaly_forecasts_delta_trip5_minus_history": {
                    "longitudinal_twin": _regression_metrics(
                        y_te_delta,
                        (pred_twin_1l - pred_persist_1l)[has_prior_history],
                        site_ids=s_te_delta,
                    ),
                    "environmental_ridge_spatial_only": _regression_metrics(
                        y_te_delta,
                        (pred_env_ridge_1l - pred_persist_1l)[has_prior_history],
                        site_ids=s_te_delta,
                    ),
                    "random_forest_baseline": _regression_metrics(
                        y_te_delta,
                        (pred_rf_1l - pred_persist_1l)[has_prior_history],
                        site_ids=s_te_delta,
                    ),
                },
            },
        }

    def run_tier2_temporal_loco_cv(self) -> Dict[str, Any]:
        """Tier 2: Leave-One-Campaign-Out (LOCO) CV, Space-for-Time Transport, and Trip 5 Prospective Forecast."""
        df = self.build_sample_level_validation_panel(id_floor=60)
        df_guild = df[df["guild_observed"]].reset_index(drop=True)
        df_mob = df[df["mobilome_observed"] & df["guild_observed"]].reset_index(drop=True)

        feature_cols = [
            "is_surface",
            "is_deep",
            "is_root",
            "map_bestmatch_mm",
            "mat_era5_c",
            "ndvi_W60d",
            "xrf_evaporite_index",
            "ph_mean",
            "log1p_dna_conc",
            "clr_Si",
            "clr_Ca",
            "clr_Fe",
        ]
        trips = sorted(df_guild["trip"].unique())
        oof_1l = np.full(len(df_guild), np.nan)
        oof_h2tot = np.full(len(df_guild), np.nan)
        oof_viral = np.full(len(df_mob), np.nan)

        per_trip_metrics: Dict[str, Any] = {}
        for trip in trips:
            tr_mask = (df_guild["trip"] != trip).to_numpy()
            te_mask = (df_guild["trip"] == trip).to_numpy()
            X_tr = df_guild.loc[tr_mask, feature_cols].to_numpy(dtype=float)
            X_te = df_guild.loc[te_mask, feature_cols].to_numpy(dtype=float)

            p_1l = _ridge_fit_predict(
                X_tr,
                df_guild.loc[tr_mask, "log2_hydB_1l"].to_numpy(dtype=float),
                X_te,
                alpha=50.0,
            )
            p_h2 = _ridge_fit_predict(
                X_tr,
                df_guild.loc[tr_mask, "log2_total_h2ase_ge"].to_numpy(dtype=float),
                X_te,
                alpha=50.0,
            )
            oof_1l[te_mask] = p_1l
            oof_h2tot[te_mask] = p_h2

            tr_mob = (df_mob["trip"] != trip).to_numpy()
            te_mob = (df_mob["trip"] == trip).to_numpy()
            X_tr_mob = df_mob.loc[
                tr_mob, feature_cols + ["log1p_scg_rpk", "scg_genes_detected"]
            ].to_numpy(dtype=float)
            X_te_mob = df_mob.loc[
                te_mob, feature_cols + ["log1p_scg_rpk", "scg_genes_detected"]
            ].to_numpy(dtype=float)
            oof_viral[te_mob] = _ridge_fit_predict(
                X_tr_mob,
                df_mob.loc[tr_mob, "log1p_viral_richness"].to_numpy(dtype=float),
                X_te_mob,
                alpha=25.0,
            )

            te_sites = df_guild.loc[te_mask, "site"].to_numpy()
            per_trip_metrics[str(trip)] = {
                "n_samples": int(te_mask.sum()),
                "log2_hydB_1l": _regression_metrics(
                    df_guild.loc[te_mask, "log2_hydB_1l"].to_numpy(dtype=float),
                    p_1l,
                    site_ids=te_sites,
                ),
                "log2_total_h2ase_ge": _regression_metrics(
                    df_guild.loc[te_mask, "log2_total_h2ase_ge"].to_numpy(dtype=float),
                    p_h2,
                    site_ids=te_sites,
                ),
            }

        sft_and_t5 = self._run_space_for_time_and_trip5_forecast(df_guild)

        return {
            "validation_scheme": "Leave-One-Campaign-Out (LOCO), Space-for-Time Transport, and Prospective Trip 5 Forecast",
            "pooled_metrics": {
                "log2_hydB_1l": _regression_metrics(
                    df_guild["log2_hydB_1l"].to_numpy(dtype=float),
                    oof_1l,
                    site_ids=df_guild["site"].to_numpy(),
                ),
                "log2_total_h2ase_ge": _regression_metrics(
                    df_guild["log2_total_h2ase_ge"].to_numpy(dtype=float),
                    oof_h2tot,
                    site_ids=df_guild["site"].to_numpy(),
                ),
                "log1p_viral_richness": _regression_metrics(
                    df_mob["log1p_viral_richness"].to_numpy(dtype=float),
                    oof_viral,
                    site_ids=df_mob["site"].to_numpy(),
                ),
            },
            "per_campaign_breakdown": per_trip_metrics,
            **sft_and_t5,
        }

    def run_tier3_causal_ontological_polytope_benchmark(self) -> Dict[str, Any]:
        """Tier 3: Benchmark Unconstrained ML vs Euclidean, WLS, Bottom-Up, Shrinkage-MinT, and Nonnegative-MinT."""
        df_all = self.build_sample_level_validation_panel(id_floor=60)
        df = df_all[df_all["guild_observed"]].reset_index(drop=True)

        guild_vars = [
            "hydB_1l",
            "hydB_1h",
            "hydB_2a",
            "total_h2ase_ge",
            "coxL_formI",
            "total_trace_gas_scavenger_ge",
        ]
        leaf_indices = [0, 1, 2, 4]
        agg_indices = [3, 5]

        # Summing matrix S (6 x 4) mapping bottom leaves [hydB_1l, hydB_1h, hydB_2a, coxL_formI] to all 6 nodes
        S_mat = np.array(
            [
                [1.0, 0.0, 0.0, 0.0],  # hydB_1l
                [0.0, 1.0, 0.0, 0.0],  # hydB_1h
                [0.0, 0.0, 1.0, 0.0],  # hydB_2a
                [1.0, 1.0, 1.0, 0.0],  # total_h2ase_ge = 1l + 1h + 2a
                [0.0, 0.0, 0.0, 1.0],  # coxL_formI
                [1.0, 1.0, 1.0, 1.0],  # total_trace_gas_scavenger_ge = total_h2ase + coxL_formI
            ],
            dtype=float,
        )

        # Compile the 6-variable trace-gas sub-polytope directly from the formal OWL 2 DL specification
        owl_poly = OntologicalPolytope.from_owl_file()
        name_map = {
            "hydB_1l_ge": "hydB_1l",
            "hydB_1h_ge": "hydB_1h",
            "hydB_2a_ge": "hydB_2a",
            "total_h2ase_ge": "total_h2ase_ge",
            "coxL_formI_ge": "coxL_formI",
            "total_trace_gas_scavenger_ge": "total_trace_gas_scavenger_ge",
        }
        polytope = OntologicalPolytope(
            variable_names=guild_vars,
            subsumptions=[
                SubsumptionRelation(
                    name_map[s.child_var], name_map[s.parent_var], s.ontology_axiom
                )
                for s in owl_poly.subsumptions
                if s.child_var in name_map and s.parent_var in name_map
            ],
            partitions=[
                PartitionSumConstraint(
                    children_vars=tuple(name_map[c] for c in p.children_vars),
                    parent_var=name_map[p.parent_var],
                    exhaustive=p.exhaustive,
                    ontology_axiom=p.ontology_axiom,
                )
                for p in owl_poly.partitions
                if p.parent_var in name_map and all(c in name_map for c in p.children_vars)
            ],
        )

        unique_sites = df[["site", "lon"]].drop_duplicates().sort_values(["lon", "site"])
        unique_sites["spatial_fold"] = np.floor(
            np.linspace(0, 5 - 1e-6, len(unique_sites))
        ).astype(int)
        df["spatial_fold"] = df["site"].map(
            dict(zip(unique_sites["site"], unique_sites["spatial_fold"]))
        )

        leaf_feats = [
            "is_surface",
            "is_deep",
            "is_root",
            "map_era5_mm",
            "mat_era5_c",
            "elevation_m",
            "xrf_evaporite_index",
            "ph_mean",
            "log1p_dna_conc",
        ]
        agg_feats = [
            "is_surface",
            "is_deep",
            "is_root",
            "map_bestmatch_mm",
            "ndvi_W60d",
            "smap_sm_W30d",
            "P365_era5_exact",
            "clr_Si",
            "clr_Ca",
        ]

        n = len(df)
        Y_true = df[guild_vars].to_numpy(dtype=float)
        Y_unconstrained = np.full_like(Y_true, np.nan)
        Y_matched_ridge = np.full_like(Y_true, np.nan)
        Y_wls_nn = np.full_like(Y_true, np.nan)
        Y_mint_shrunk = np.full_like(Y_true, np.nan)
        Y_mint_nn = np.full_like(Y_true, np.nan)
        fold_shrinkage_lambdas: List[float] = []

        for f in range(5):
            tr = (df["spatial_fold"] != f).to_numpy()
            te = (df["spatial_fold"] == f).to_numpy()
            X_tr_leaf = df.loc[tr, leaf_feats].to_numpy(dtype=float)
            X_te_leaf = df.loc[te, leaf_feats].to_numpy(dtype=float)
            X_tr_agg = df.loc[tr, agg_feats].to_numpy(dtype=float)
            X_te_agg = df.loc[te, agg_feats].to_numpy(dtype=float)

            Y_tr_hat = np.zeros_like(Y_true[tr])
            for idx, vname in enumerate(guild_vars):
                y_tr = Y_true[tr, idx]
                Y_matched_ridge[te, idx] = _ridge_fit_predict(
                    X_tr_leaf, y_tr, X_te_leaf, alpha=8.0
                )
                if vname in ("total_h2ase_ge", "total_trace_gas_scavenger_ge"):
                    Y_unconstrained[te, idx] = _ridge_fit_predict(
                        X_tr_agg, y_tr, X_te_agg, alpha=1.0
                    )
                    Y_tr_hat[:, idx] = _ridge_fit_predict(
                        X_tr_agg, y_tr, X_tr_agg, alpha=1.0
                    )
                else:
                    Y_unconstrained[te, idx] = _ridge_fit_predict(
                        X_tr_leaf, y_tr, X_te_leaf, alpha=8.0
                    )
                    Y_tr_hat[:, idx] = _ridge_fit_predict(
                        X_tr_leaf, y_tr, X_tr_leaf, alpha=8.0
                    )

            # Estimate residual covariance strictly inside the training fold
            E_tr = Y_tr_hat - Y_true[tr]
            W_lw, shrink_lam = MinTReconciler.estimate_ledoit_wolf_covariance(E_tr)
            fold_shrinkage_lambdas.append(shrink_lam)
            W_diag = np.diag(np.maximum(1e-6, np.var(E_tr, axis=0)))

            Y_wls_nn[te] = MinTReconciler.reconcile_matrix(
                Y_unconstrained[te], S_mat, W_diag, nonnegative=True
            )
            Y_mint_shrunk[te] = MinTReconciler.reconcile_matrix(
                Y_unconstrained[te], S_mat, W_lw, nonnegative=False
            )
            Y_mint_nn[te] = MinTReconciler.reconcile_matrix(
                Y_unconstrained[te], S_mat, W_lw, nonnegative=True
            )

        Y_projected = np.zeros_like(Y_unconstrained)
        Y_bottom_up = np.zeros_like(Y_unconstrained)
        unconstrained_states: List[Dict[str, float]] = []
        projected_states: List[Dict[str, float]] = []
        mint_nn_states: List[Dict[str, float]] = []
        cvm_unconstrained_list: List[float] = []
        cvm_projected_list: List[float] = []

        for i in range(n):
            raw_dict = {
                guild_vars[j]: float(Y_unconstrained[i, j]) for j in range(len(guild_vars))
            }
            unconstrained_states.append(raw_dict)
            v_neg = sum(max(0.0, -raw_dict[k]) for k in guild_vars)
            v_part1 = abs(
                (raw_dict["hydB_1l"] + raw_dict["hydB_1h"] + raw_dict["hydB_2a"])
                - raw_dict["total_h2ase_ge"]
            )
            v_part2 = abs(
                (raw_dict["total_h2ase_ge"] + raw_dict["coxL_formI"])
                - raw_dict["total_trace_gas_scavenger_ge"]
            )
            cvm_unconstrained_list.append(v_neg + v_part1 + v_part2)

            proj_dict = polytope.project(raw_dict)
            projected_states.append(proj_dict)
            for j, k in enumerate(guild_vars):
                Y_projected[i, j] = proj_dict[k]

            mint_nn_states.append(
                {guild_vars[j]: float(Y_mint_nn[i, j]) for j in range(len(guild_vars))}
            )

            b1l = max(0.0, raw_dict["hydB_1l"])
            b1h = max(0.0, raw_dict["hydB_1h"])
            b2a = max(0.0, raw_dict["hydB_2a"])
            bcox = max(0.0, raw_dict["coxL_formI"])
            btot_h2 = b1l + b1h + b2a
            Y_bottom_up[i, :] = [b1l, b1h, b2a, btot_h2, bcox, btot_h2 + bcox]

            v_neg_p = sum(max(0.0, -proj_dict[k]) for k in guild_vars)
            v_part1_p = abs(
                (proj_dict["hydB_1l"] + proj_dict["hydB_1h"] + proj_dict["hydB_2a"])
                - proj_dict["total_h2ase_ge"]
            )
            v_part2_p = abs(
                (proj_dict["total_h2ase_ge"] + proj_dict["coxL_formI"])
                - proj_dict["total_trace_gas_scavenger_ge"]
            )
            cvm_projected_list.append(v_neg_p + v_part1_p + v_part2_p)

        ocs_before = polytope.compute_ocs(unconstrained_states, tol=1e-4)
        ocs_after = polytope.compute_ocs(projected_states, tol=1e-4)
        ocs_mint_nn = polytope.compute_ocs(mint_nn_states, tol=1e-4)

        mse_before = float(np.mean((Y_unconstrained - Y_true) ** 2))
        mse_after = float(np.mean((Y_projected - Y_true) ** 2))
        mse_mint_nn = float(np.mean((Y_mint_nn - Y_true) ** 2))
        rmse_before = float(np.sqrt(mse_before))
        rmse_after = float(np.sqrt(mse_after))
        rmse_mint_nn = float(np.sqrt(mse_mint_nn))
        rmse_mint_shrunk = float(np.sqrt(np.mean((Y_mint_shrunk - Y_true) ** 2)))
        rmse_wls_nn = float(np.sqrt(np.mean((Y_wls_nn - Y_true) ** 2)))
        pct_mse_reduction = 100.0 * (mse_before - mse_after) / max(1e-12, mse_before)
        pct_mse_reduction_mint_nn = (
            100.0 * (mse_before - mse_mint_nn) / max(1e-12, mse_before)
        )

        sites_arr = df["site"].to_numpy()
        uniq_s = np.unique(sites_arr)
        rng = np.random.default_rng(42)
        idx_by_s = {s: np.where(sites_arr == s)[0] for s in uniq_s}
        boot_red: List[float] = []
        for _ in range(300):
            s_draw = rng.choice(uniq_s, size=len(uniq_s), replace=True)
            b_idx = np.concatenate([idx_by_s[s] for s in s_draw])
            m_b0 = float(np.mean((Y_unconstrained[b_idx] - Y_true[b_idx]) ** 2))
            m_b1 = float(np.mean((Y_projected[b_idx] - Y_true[b_idx]) ** 2))
            boot_red.append(100.0 * (m_b0 - m_b1) / max(1e-12, m_b0))
        mse_red_ci95 = [
            round(float(np.percentile(boot_red, 2.5)), 2),
            round(float(np.percentile(boot_red, 97.5)), 2),
        ]

        leaf_rmse_before = float(
            np.sqrt(np.mean((Y_unconstrained[:, leaf_indices] - Y_true[:, leaf_indices]) ** 2))
        )
        leaf_rmse_after = float(
            np.sqrt(np.mean((Y_projected[:, leaf_indices] - Y_true[:, leaf_indices]) ** 2))
        )
        leaf_rmse_mint_nn = float(
            np.sqrt(np.mean((Y_mint_nn[:, leaf_indices] - Y_true[:, leaf_indices]) ** 2))
        )
        agg_rmse_before = float(
            np.sqrt(np.mean((Y_unconstrained[:, agg_indices] - Y_true[:, agg_indices]) ** 2))
        )
        agg_rmse_after = float(
            np.sqrt(np.mean((Y_projected[:, agg_indices] - Y_true[:, agg_indices]) ** 2))
        )
        agg_rmse_mint_nn = float(
            np.sqrt(np.mean((Y_mint_nn[:, agg_indices] - Y_true[:, agg_indices]) ** 2))
        )

        rmse_bottom_up = float(np.sqrt(np.mean((Y_bottom_up - Y_true) ** 2)))
        rmse_matched_ridge = float(np.sqrt(np.mean((Y_matched_ridge - Y_true) ** 2)))

        sample_l2_before = np.sum((Y_unconstrained - Y_true) ** 2, axis=1)
        sample_l2_after = np.sum((Y_projected - Y_true) ** 2, axis=1)
        theorem2_pass_rate = float(np.mean(sample_l2_after <= sample_l2_before + 1e-7))
        strict_improvement_rate = float(np.mean(sample_l2_after < sample_l2_before - 1e-9))

        sub_df = df[df["compartment_norm"] == "ShallowSubsurface"]
        baseline_abund = {
            "hydB_1l": float(sub_df["hydB_1l"].mean()),
            "hydB_1h": float(sub_df["hydB_1h"].mean()),
            "hydB_2a": float(sub_df["hydB_2a"].mean()),
        }
        sibling_slopes = {
            "hydB_1l": (-2.2747, 0.33),
            "hydB_1h": (-0.3277, 0.16),
            "hydB_2a": (-0.3432, 0.13),
        }
        prop2_res = self.causal_engine.agreement_weighted_class_transport(
            sibling_slopes,
            baseline_abundances=baseline_abund,
        )

        return {
            "n_samples_evaluated": n,
            "ocs_unconstrained": round(ocs_before, 4),
            "ocs_projected_two_layer": round(ocs_after, 4),
            "ocs_nonnegative_mint": round(ocs_mint_nn, 4),
            "mean_cvm_unconstrained": round(float(np.mean(cvm_unconstrained_list)), 4),
            "mean_cvm_projected_two_layer": round(float(np.mean(cvm_projected_list)), 6),
            "rmse_unconstrained": round(rmse_before, 4),
            "rmse_projected_two_layer": round(rmse_after, 4),
            "rmse_nonnegative_mint": round(rmse_mint_nn, 4),
            "rmse_shrinkage_mint_unconstrained": round(rmse_mint_shrunk, 4),
            "rmse_diagonal_wls_nonnegative": round(rmse_wls_nn, 4),
            "leaf_rmse_unconstrained": round(leaf_rmse_before, 4),
            "leaf_rmse_projected_two_layer": round(leaf_rmse_after, 4),
            "leaf_rmse_nonnegative_mint": round(leaf_rmse_mint_nn, 4),
            "aggregate_rmse_unconstrained": round(agg_rmse_before, 4),
            "aggregate_rmse_projected_two_layer": round(agg_rmse_after, 4),
            "aggregate_rmse_nonnegative_mint": round(agg_rmse_mint_nn, 4),
            "mse_reduction_pct": round(pct_mse_reduction, 2),
            "mse_reduction_pct_nonnegative_mint": round(pct_mse_reduction_mint_nn, 2),
            "mse_reduction_pct_cluster_ci95": mse_red_ci95,
            "mean_ledoit_wolf_shrinkage_intensity": round(float(np.mean(fold_shrinkage_lambdas)), 4),
            "theorem2_samplewise_error_reduction_rate": round(theorem2_pass_rate, 4),
            "strict_samplewise_error_reduction_rate": round(strict_improvement_rate, 4),
            "reconciliation_baselines": {
                "bottom_up_nonneg_leaf_sum_rmse": round(rmse_bottom_up, 4),
                "diagonal_wls_nonneg_rmse": round(rmse_wls_nn, 4),
                "shrinkage_mint_unconstrained_rmse": round(rmse_mint_shrunk, 4),
                "nonnegative_shrinkage_mint_rmse": round(rmse_mint_nn, 4),
                "matched_feature_linear_ridge_rmse": round(rmse_matched_ridge, 4),
            },
            "proposition2_sibling_class_transport": prop2_res,
        }

    def run_tier4_natural_experiment_counterfactual_calibration(self) -> Dict[str, Any]:
        """Tier 4: Consistency checks against paired-visit field contrasts, turnover decomposition, and CBB co-carriage."""
        h123_dir = self.paths.h123_v2_dir
        contrasts_df = pd.read_csv(h123_dir / "contrasts_by_position.tsv", sep="\t")
        rvb_visits_df = pd.read_csv(h123_dir / "h1_visit_log2fc_root_vs_bulk.tsv", sep="\t")
        svb_visits_df = pd.read_csv(h123_dir / "h1_visit_log2fc_surface_vs_bulk.tsv", sep="\t")
        a2_rvb_df = pd.read_csv(h123_dir / "A2_root_bulk_decomposition_summary.tsv", sep="\t")
        a2_arid_df = pd.read_csv(h123_dir / "A2_aridity_decomposition_summary.tsv", sep="\t")
        a4_cbb_df = pd.read_csv(h123_dir / "A4_cbb_cocarriage_by_group.tsv", sep="\t")
        a4_sub_df = pd.read_csv(h123_dir / "A4_cbb_subsets_position_and_aridity.tsv", sep="\t")

        comp_cube = self.builder.build_compartment_state_cube()
        reveg = RevegetationIntervention()
        refugium = RefugiumThermalIntervention()

        # 1. Root vs Bulk consistency check across all 47 paired sites (n=106 visits)
        paired_rvb_sites = sorted(rvb_visits_df["site"].unique())
        twin_rvb_log2fc_1l: List[float] = []
        twin_rvb_log2fc_dna: List[float] = []

        for s in paired_rvb_sites:
            sub = comp_cube[
                (comp_cube["site_num"] == int(s))
                & (comp_cube["compartment_norm"] == "ShallowSubsurface")
            ]
            if sub.empty:
                continue
            row = sub.iloc[0].to_dict()
            sim_root = reveg.simulate(
                row,
                shrub_cover_fraction=0.85,
                rak_pr_pgpr_dose=0.5,
                site60_rbcl_hyab_plasmid_dose=0.0,
            )
            base_1l = float(sim_root["summary"]["baseline_hydB_1l_ge"])
            root_1l = float(sim_root["summary"]["root_zone_hydB_1l_ge"])
            twin_rvb_log2fc_1l.append(float(np.log2((root_1l + 0.01) / (base_1l + 0.01))))
            twin_rvb_log2fc_dna.append(
                float(np.log2(sim_root["summary"]["biomass_dna_multiplier"]))
            )

        emp_rvb_1l_row = contrasts_df[
            (contrasts_df["contrast"] == "root/bulk")
            & (contrasts_df["quantity"] == "hydB_1l")
        ].iloc[0]
        emp_rvb_median = float(emp_rvb_1l_row["median_log2"])
        emp_rvb_ci_lo = float(emp_rvb_1l_row["ci_lo"])
        emp_rvb_ci_hi = float(emp_rvb_1l_row["ci_hi"])
        twin_rvb_median = float(np.median(twin_rvb_log2fc_1l))

        site_emp_medians = rvb_visits_df.groupby("site")["hydB_1l"].median()
        sign_matches = sum(1 for v in site_emp_medians.to_numpy() if v < 0.0) / len(
            site_emp_medians
        )

        # 2. Surface vs Subsurface (5-10 cm refugium) across 36 paired sites (n=68 visits)
        paired_svb_sites = sorted(svb_visits_df["site"].unique())
        twin_svb_log2fc_1l: List[float] = []
        for s in paired_svb_sites:
            sub = comp_cube[
                (comp_cube["site_num"] == int(s))
                & (comp_cube["compartment_norm"] == "ShallowSubsurface")
            ]
            if sub.empty:
                continue
            row = sub.iloc[0].to_dict()
            sim_intact = refugium.simulate(
                row, subsurface_refugium_intact_fraction=1.0, surface_temp_c=35.0
            )
            sim_exposed = refugium.simulate(
                row, subsurface_refugium_intact_fraction=0.0, surface_temp_c=58.0
            )
            h1l_intact = float(sim_intact["projected_state"]["hydB_1l_ge"])
            h1l_exposed = float(sim_exposed["projected_state"]["hydB_1l_ge"])
            twin_svb_log2fc_1l.append(
                float(np.log2((h1l_exposed + 0.01) / (h1l_intact + 0.01)))
            )

        emp_svb_1l_row = contrasts_df[
            (contrasts_df["contrast"] == "surface/bulk")
            & (contrasts_df["quantity"] == "hydB_1l")
        ].iloc[0]
        emp_svb_median = float(emp_svb_1l_row["median_log2"])
        emp_svb_ci_lo = float(emp_svb_1l_row["ci_lo"])
        emp_svb_ci_hi = float(emp_svb_1l_row["ci_hi"])
        twin_svb_median = float(np.median(twin_svb_log2fc_1l))

        # 3. Transect Aridity Gradient (observed ShallowSubsurface sites, n=51)
        deep_obs = comp_cube[
            (comp_cube["compartment_norm"] == "ShallowSubsurface")
            & comp_cube["wgs_observed"]
        ].dropna(subset=["hydB_1l_ge_raw", "map_mm"])
        rho_aridity, p_aridity = stats.spearmanr(
            deep_obs["hydB_1l_ge_raw"], deep_obs["map_mm"]
        )

        a2_rvb_row = a2_rvb_df.iloc[0].to_dict()
        a2_arid_row = a2_arid_df.iloc[0].to_dict()

        return {
            "root_vs_bulk_rhizosphere_paradox": {
                "n_paired_visits": int(emp_rvb_1l_row["n_visits"]),
                "n_paired_sites": int(emp_rvb_1l_row["n_sites"]),
                "empirical_median_log2fc_hydB_1l": round(emp_rvb_median, 4),
                "empirical_95ci": [round(emp_rvb_ci_lo, 4), round(emp_rvb_ci_hi, 4)],
                "twin_predicted_median_log2fc_hydB_1l": round(twin_rvb_median, 4),
                "abs_calibration_error_log2": round(abs(twin_rvb_median - emp_rvb_median), 4),
                "within_empirical_95ci": bool(
                    emp_rvb_ci_lo <= twin_rvb_median <= emp_rvb_ci_hi
                ),
                "empirical_root_bulk_log2_dna_ratio": 3.37,
                "twin_predicted_log2_dna_ratio": round(
                    float(np.median(twin_rvb_log2fc_dna)), 4
                ),
                "empirical_site_depletion_fraction": round(float(sign_matches), 4),
            },
            "lineage_turnover_vs_copy_number_decomposition": {
                "root_vs_bulk_decomposition": a2_rvb_row,
                "aridity_gradient_decomposition": a2_arid_row,
            },
            "h2_oxidation_vs_cbb_autotrophy_decoupling": {
                "cbb_cocarriage_by_group": a4_cbb_df.to_dict(orient="records"),
                "cbb_subsets_position_and_aridity": a4_sub_df.to_dict(orient="records"),
            },
            "surface_vs_subsurface_refugium": {
                "n_paired_visits": int(emp_svb_1l_row["n_visits"]),
                "n_paired_sites": int(emp_svb_1l_row["n_sites"]),
                "empirical_median_log2fc_hydB_1l": round(emp_svb_median, 4),
                "empirical_95ci": [round(emp_svb_ci_lo, 4), round(emp_svb_ci_hi, 4)],
                "twin_predicted_median_log2fc_hydB_1l": round(twin_svb_median, 4),
                "abs_calibration_error_log2": round(abs(twin_svb_median - emp_svb_median), 4),
                "within_empirical_95ci": bool(
                    emp_svb_ci_lo <= twin_svb_median <= emp_svb_ci_hi
                ),
            },
            "aridity_gradient_hydB_1l_vs_MAP": {
                "n_sites": int(len(deep_obs)),
                "cube_spearman_rho": round(float(rho_aridity), 4),
                "cube_spearman_p": float(p_aridity),
                "paper_reference_rho_bestmatch": -0.5119,
                "sign_and_significance_match": bool(
                    rho_aridity < -0.35 and p_aridity < 0.01
                ),
            },
        }

    def run_tier5_synthetic_and_metabolic_benchmark(self) -> Dict[str, Any]:
        """Tier 5: Synthetic Metagenome Read-to-Guild Accuracy, 72-MAG GSMM Pilot Audit & Factorial Brine Analysis."""
        bench_dir = self.paths.eq_lifeonair / "lifeonair" / "results" / "benchmark_r2"
        mock_acc = pd.read_csv(bench_dir / "mock_ge_accuracy.tsv", sep="\t")
        mock_rb = pd.read_csv(bench_dir / "mock_ge_root_bulk.tsv", sep="\t")

        sub_acc = mock_acc[
            (mock_acc["db"] == "full")
            & (mock_acc["id_floor"] == 60)
            & (mock_acc["group"].isin(["hydB_1l", "hydB_1h", "hydB_2a", "coxL_formI"]))
        ]
        sub_rb = mock_rb[
            (mock_rb["db"] == "full")
            & (mock_rb["id_floor"] == 60)
            & (mock_rb["group"] == "hydB_1l")
        ].dropna(subset=["true_log2_root_bulk", "est_log2_root_bulk"])

        rb_metrics = _regression_metrics(
            sub_rb["true_log2_root_bulk"].to_numpy(dtype=float),
            sub_rb["est_log2_root_bulk"].to_numpy(dtype=float),
        )

        fba = MetabolicCommunityEngine()
        closed_res = fba.solve_community_fba(
            CommunityMetabolicInputs(
                h2_uptake_bound=0.0,
                co_uptake_bound=0.0,
                root_exudate_doc_bound=0.0,
                exogenous_nh4_no3_bound=0.0,
                sulfate_so4_bound=0.0,
                choline_betaine_bound=0.0,
            )
        )
        closed_exchange_pass = bool(
            closed_res["total_biomass_flux"] <= 1e-7
            and not closed_res["hard_maintenance_feasible"]
        )

        so4_only_res = fba.solve_community_fba(
            CommunityMetabolicInputs(
                h2_uptake_bound=0.0,
                co_uptake_bound=0.0,
                root_exudate_doc_bound=0.0,
                exogenous_nh4_no3_bound=0.0,
                sulfate_so4_bound=1.0,
                choline_betaine_bound=0.0,
            )
        )
        sulfate_donor_closure_pass = bool(
            so4_only_res["total_biomass_flux"] <= 1e-7
            and so4_only_res["sulfur_redox_flux"] <= 1e-7
        )

        h2_grid = np.linspace(0.0, 2.5, 25)
        h2_biomass = [
            fba.solve_community_fba(
                CommunityMetabolicInputs(
                    h2_uptake_bound=float(h), exogenous_nh4_no3_bound=0.5
                )
            )["total_biomass_flux"]
            for h in h2_grid
        ]
        h2_monotonic_pass = bool(np.all(np.diff(h2_biomass) >= -1e-7))

        ec_grid = np.linspace(2.0, 220.0, 25)
        ec_atpm = [
            fba.solve_community_fba(CommunityMetabolicInputs(salinity_ms_cm=float(ec)))[
                "atpm_required_unmitigated"
            ]
            for ec in ec_grid
        ]
        salinity_atpm_monotonic_pass = bool(np.all(np.diff(ec_atpm) >= -1e-7))

        # Common-baseline 4x2x2 factorial brine analysis at Site 60 reference baseline
        comp_cube = self.builder.build_compartment_state_cube()
        s60_row = (
            comp_cube[
                (comp_cube["site_num"] == 60)
                & (comp_cube["compartment_norm"] == "Surface")
            ]
            .iloc[0]
            .to_dict()
        )
        desal = DesalinationBrineIntervention()
        factorial_brine = desal.evaluate_factorial_brine_grid(s60_row)

        return {
            "synthetic_metagenome_ge_benchmark": {
                "markers_evaluated": sub_acc[
                    ["group", "n", "spearman_log", "pearson_log", "median_log2_bias"]
                ].to_dict(orient="records"),
                "mean_pearson_r_across_guilds": round(
                    float(sub_acc["pearson_log"].mean()), 4
                ),
                "mean_spearman_rho_across_guilds": round(
                    float(sub_acc["spearman_log"].mean()), 4
                ),
                "root_vs_bulk_synthetic_contrast_hydB_1l": rb_metrics,
            },
            "metabolic_pilot_72mag_10584run_audit": fba.METABOLIC_PILOT_AUDIT_SUMMARY,
            "measured_aquatic_hydrochemistry": fba.MEASURED_AQUATIC_HYDROCHEMISTRY,
            "factorial_brine_common_baseline_grid": factorial_brine,
            "highs_fba_thermodynamic_audit": {
                "closed_exchange_zero_growth_pass": closed_exchange_pass,
                "sulfate_electron_donor_closure_pass": sulfate_donor_closure_pass,
                "h2_uptake_monotonicity_25_steps_pass": h2_monotonic_pass,
                "salinity_atpm_monotonicity_25_steps_pass": salinity_atpm_monotonic_pass,
            },
        }

    def run_carbon_budget_cross_desert_validation(self) -> Dict[str, Any]:
        """Validate the 4-pathway CarbonBudgetEngine internally and against 7 published cross-desert benchmarks."""
        comp_cube = self.builder.build_compartment_state_cube()
        grid_df = self.builder.build_era5_h3_grid()
        engine = CarbonBudgetEngine()
        optimizer = ParetoPortfolioOptimizer()

        # 1. Empirical CBB / RuBisCO and H2-scavenger read shares across the 155 observed WGS site-compartments
        obs_comp = comp_cube[comp_cube["wgs_observed"]].copy()
        cbb_by_comp: Dict[str, Dict[str, float]] = {}
        for comp_name, grp in obs_comp.groupby("compartment_norm"):
            s_cbb = grp["share_cbb_mapped_raw"].dropna()
            s_cbb_prk = grp["share_cbb_prk_mapped_raw"].dropna()
            s_cbb_3g = grp["share_cbb_3gene_mapped_raw"].dropna()
            s_1l = grp["share_1l_mapped_raw"].dropna()
            s_1l_cbb = grp["share_1l_cbb_mapped_raw"].dropna()
            s_1l_3g = grp["share_1l_cbb_3gene_mapped_raw"].dropna()
            valid_1l = s_1l > 1e-9
            cond_ratio = (
                (s_1l_cbb[valid_1l] / s_1l[valid_1l]).mean() if valid_1l.any() else 0.0
            )
            cond_3g_ratio = (
                (s_1l_3g[valid_1l] / s_1l[valid_1l]).mean() if valid_1l.any() else 0.0
            )
            cbb_by_comp[str(comp_name)] = {
                "n_observed_sites": int(len(s_cbb)),
                "share_cbb_mapped_mean": round(float(s_cbb.mean()), 4),
                "share_cbb_prkb_mapped_mean": round(float(s_cbb_prk.mean()), 4),
                "share_cbb_3gene_mapped_mean": round(float(s_cbb_3g.mean()), 4),
                "share_1l_mapped_mean": round(float(s_1l.mean()), 4),
                "share_1l_cbb_mapped_mean": round(float(s_1l_cbb.mean()), 4),
                "share_1l_cbb_3gene_mapped_mean": round(float(s_1l_3g.mean()), 4),
                "conditional_cbb_given_1l_mean": round(float(cond_ratio), 4),
                "conditional_cbb_3gene_given_1l_mean": round(float(cond_3g_ratio), 4),
            }

        overall_cbb_mean = float(obs_comp["share_cbb_mapped_raw"].dropna().mean())
        overall_cbb_prkb_mean = float(obs_comp["share_cbb_prk_mapped_raw"].dropna().mean())
        overall_cbb_3gene_mean = float(obs_comp["share_cbb_3gene_mapped_raw"].dropna().mean())

        # MAG-level CBB co-carriage from A4_cbb_cocarriage_by_group.tsv
        a4_df = pd.read_csv(
            self.paths.h123_v2_dir / "A4_cbb_cocarriage_by_group.tsv", sep="\t"
        )
        g1l_row = a4_df[a4_df["group"] == "1l"].iloc[0]
        mag_1l_cbb_frac = float(g1l_row["frac_cbb"])
        mag_1l_hc90_cbb_frac = float(g1l_row["hc90_frac_1l"])

        # 2. Internal thermodynamic, stoichiometric, and physical invariant checks
        s35 = (
            comp_cube[
                (comp_cube["site_num"] == 35) & (comp_cube["compartment_norm"] == "Surface")
            ]
            .iloc[0]
            .to_dict()
        )
        b_infeasible = engine.compute_state_carbon_budget(
            s35, hard_maintenance_feasible=False
        )
        maintenance_zero_pass = (
            b_infeasible["direct_microbial_rubisco_t_co2_km2_yr"]["point"] == 0.0
            and b_infeasible["direct_microbial_rubisco_t_co2_km2_yr"]["high"] == 0.0
        )

        s35_zero_gas = dict(s35)
        for gk in ("hydB_1l_ge", "hydB_1h_ge", "hydB_2a_ge", "coxL_formI_ge"):
            s35_zero_gas[gk] = 0.0
        b_zero_gas = engine.compute_state_carbon_budget(s35_zero_gas)
        zero_gas_pass = (
            b_zero_gas["direct_microbial_rubisco_t_co2_km2_yr"]["point"] == 0.0
            and b_zero_gas["indirect_h2_sink_t_co2eq_km2_yr"]["point"] == 0.0
        )

        plasmid_grid = np.linspace(0.0, 1.0, 11)
        rubisco_vs_plasmid = [
            engine.compute_state_carbon_budget(
                s35, shrub_cover_fraction=0.25, rbcl_hyab_plasmid_dose=float(p_dose)
            )["direct_microbial_rubisco_t_co2_km2_yr"]["point"]
            for p_dose in plasmid_grid
        ]
        plasmid_monotonic_pass = bool(np.all(np.diff(rubisco_vs_plasmid) >= -1e-9))

        cover_grid = np.linspace(0.0, 0.60, 13)
        veg_vs_cover = [
            engine.compute_state_carbon_budget(s35, shrub_cover_fraction=float(c))[
                "vegetation_npp_t_co2_km2_yr"
            ]["point"]
            for c in cover_grid
        ]
        sic_vs_cover = [
            engine.compute_state_carbon_budget(s35, shrub_cover_fraction=float(c))[
                "pedogenic_carbonate_t_co2_km2_yr"
            ]["point"]
            for c in cover_grid
        ]
        cover_monotonic_pass = bool(
            np.all(np.diff(veg_vs_cover) >= -1e-9) and np.all(np.diff(sic_vs_cover) >= -1e-9)
        )

        porosity_grid = np.linspace(0.22, 0.42, 11)
        h2_sink_vs_porosity = [
            engine.compute_state_carbon_budget(s35, effective_porosity=float(eps))[
                "indirect_h2_sink_t_co2eq_km2_yr"
            ]["point"]
            for eps in porosity_grid
        ]
        porosity_monotonic_pass = bool(np.all(np.diff(h2_sink_vs_porosity) > 0.0))

        # 3. Upscaled corridor (92,160 km^2) and basin (650,000 km^2) budgets across Baseline and P1..P6
        basin_eval = optimizer.evaluate_basin_carbon_budgets(comp_cube, grid_df)
        base_b = basin_eval["baseline_observed_state"]
        port_b = {p["policy_name"]: p for p in basin_eval["candidate_portfolio_budgets"]}

        domain_abstention_pass = (
            base_b["supported_corridor_area_km2"] == 92160.0
            and base_b["unmanaged_reserve_area_km2"] == 557840.0
            and base_b["total_basin_area_km2"] == 650000.0
        )

        # 4. Cross-desert benchmark comparisons across all 60 primary surface sites
        surf_df = comp_cube[comp_cube["compartment_norm"] == "Surface"].copy()
        site_base_budgets = [
            engine.compute_state_carbon_budget(r.to_dict()) for _, r in surf_df.iterrows()
        ]
        site_bare_budgets = [
            engine.compute_state_carbon_budget(r.to_dict(), shrub_cover_fraction=0.002)
            for _, r in surf_df.iterrows()
        ]
        eq_dark_co2_annual_nmol_gdw_d = float(
            np.mean([b["dark_co2_fixation_annual_avg_nmol_gdw_d"] for b in site_base_budgets])
        )
        eq_dark_co2_wetted_nmol_gdw_d = float(
            np.mean([b["dark_co2_fixation_wetted_active_nmol_gdw_d"] for b in site_base_budgets])
        )
        eq_h2_mg_m2_yr = float(
            np.mean([b["h2_uptake_mg_m2_yr"]["point"] for b in site_base_budgets])
        )
        eq_bare_sic_t_km2_yr = float(
            np.mean([b["pedogenic_carbonate_t_co2_km2_yr"]["point"] for b in site_bare_budgets])
        )
        eq_cell_h2_nmol_h_cell = float(
            np.mean([b["guild_cell_specific_h2_oxidation_nmol_h_cell"] for b in site_base_budgets])
        )
        p2_corridor_t_km2 = port_b[
            "Integrated Eco-Restoration & Oryx Rewilding (Balanced Optimum)"
        ]["total_direct_co2_bound"]["corridor_t_km2_yr_point"]

        eq_values_by_metric = {
            "cbb_rubisco_community_share": round(overall_cbb_mean, 4),
            "dark_and_h2_stimulated_co2_fixation_nmol_gdw_d": round(
                eq_dark_co2_wetted_nmol_gdw_d, 5
            ),
            "cell_specific_h2_oxidation_nmol_h_cell": float(f"{eq_cell_h2_nmol_h_cell:.3e}"),
            "areal_h2_uptake_mg_m2_yr": round(eq_h2_mg_m2_yr, 3),
            "actinomycetota_rubisco_cocarriage_fraction": round(mag_1l_cbb_frac, 4),
            "bare_pedogenic_carbonate_t_co2_km2_yr": round(eq_bare_sic_t_km2_yr, 4),
            "revegetated_corridor_total_t_co2_km2_yr": round(p2_corridor_t_km2, 3),
        }

        cross_desert_checks: List[Dict[str, Any]] = []
        for bench in engine.CROSS_DESERT_BENCHMARKS:
            m_name = bench["metric"]
            eq_val = eq_values_by_metric[m_name]
            within_range = bool(bench["literature_low"] <= eq_val <= bench["literature_high"])
            cross_desert_checks.append(
                {
                    **bench,
                    "eq_biodt_value": eq_val,
                    "within_published_desert_range": within_range,
                }
            )

        return {
            "empirical_cbb_rubisco_by_compartment": cbb_by_comp,
            "overall_wgs_cbb_community_share_mean": round(overall_cbb_mean, 4),
            "overall_wgs_cbb_prkb_community_share_mean": round(overall_cbb_prkb_mean, 4),
            "overall_wgs_cbb_3gene_community_share_mean": round(overall_cbb_3gene_mean, 4),
            "mag_hydB_1l_cbb_cocarriage_fraction": round(mag_1l_cbb_frac, 4),
            "mag_hydB_1l_hc90_cbb_cocarriage_fraction": round(mag_1l_hc90_cbb_frac, 4),
            "dark_co2_fixation_annual_dry_column_nmol_gdw_d": round(
                eq_dark_co2_annual_nmol_gdw_d, 6
            ),
            "dark_co2_fixation_wetted_active_horizon_nmol_gdw_d": round(
                eq_dark_co2_wetted_nmol_gdw_d, 5
            ),
            "internal_invariants_audit": {
                "maintenance_infeasible_zero_rubisco_pass": maintenance_zero_pass,
                "zero_trace_gas_zero_rubisco_pass": zero_gas_pass,
                "plasmid_cbb_monotonicity_pass": plasmid_monotonic_pass,
                "shrub_cover_npp_and_sic_monotonicity_pass": cover_monotonic_pass,
                "compaction_diffusivity_h2_sink_monotonicity_pass": porosity_monotonic_pass,
                "domain_of_support_abstention_preserved_pass": domain_abstention_pass,
            },
            "cross_desert_benchmark_validation": cross_desert_checks,
            "all_cross_desert_benchmarks_pass": all(
                c["within_published_desert_range"] for c in cross_desert_checks
            ),
            "corridor_and_basin_budgets": {
                "baseline_observed": base_b,
                "prescribed_corridor_atlas": basin_eval["prescribed_corridor_atlas_budget"],
                "candidate_portfolios": basin_eval["candidate_portfolio_budgets"],
            },
        }

    def run_all_validations(self, output_dir: Optional[Path] = None) -> Dict[str, Any]:
        """Run all 5 validation tiers + carbon budget cross-desert validation and optionally write reports."""
        t1 = self.run_tier1_spatial_block_cv()
        t2 = self.run_tier2_temporal_loco_cv()
        t3 = self.run_tier3_causal_ontological_polytope_benchmark()
        t4 = self.run_tier4_natural_experiment_counterfactual_calibration()
        t5 = self.run_tier5_synthetic_and_metabolic_benchmark()
        t_carbon = self.run_carbon_budget_cross_desert_validation()

        report = {
            "twin_version": "0.2.0",
            "tier1_spatial_block_cv": t1,
            "tier2_temporal_loco_cv": t2,
            "tier3_causal_ontological_projection": t3,
            "tier4_natural_experiment_calibration": t4,
            "tier5_synthetic_and_metabolic_benchmark": t5,
            "carbon_budget_cross_desert_validation": t_carbon,
        }

        if output_dir is not None:
            output_dir.mkdir(parents=True, exist_ok=True)
            json_path = output_dir / "validation_report.json"

            def _clean_nans(obj: Any) -> Any:
                if isinstance(obj, float) and not np.isfinite(obj):
                    return None
                if isinstance(obj, dict):
                    return {k: _clean_nans(v) for k, v in obj.items()}
                if isinstance(obj, list):
                    return [_clean_nans(v) for v in obj]
                return obj

            clean_report = _clean_nans(report)
            json_path.write_text(
                json.dumps(clean_report, indent=2, allow_nan=False), encoding="utf-8"
            )

            t5_bench = t2["trip5_prospective_forecast_benchmark"]["level_forecasts_log2_hydB_1l"]["longitudinal_twin"]
            summary_rows = [
                {
                    "tier": "Tier 1 (5-Fold Spatial Block CV)",
                    "metric_target": f"log2(hydB_1l_ge) [n={t1['targets']['log2_hydB_1l']['n']} WGS]",
                    "rmse": t1["targets"]["log2_hydB_1l"]["rmse"],
                    "mae": t1["targets"]["log2_hydB_1l"]["mae"],
                    "r2": t1["targets"]["log2_hydB_1l"]["r2"],
                    "pearson_r": t1["targets"]["log2_hydB_1l"]["pearson_r"],
                    "spearman_rho": t1["targets"]["log2_hydB_1l"]["spearman_rho"],
                    "p_value": t1["targets"]["log2_hydB_1l"]["spearman_p"],
                },
                {
                    "tier": "Tier 1 (5-Fold Spatial Block CV)",
                    "metric_target": f"log2(total_h2ase_ge) [n={t1['targets']['log2_total_h2ase_ge']['n']} WGS]",
                    "rmse": t1["targets"]["log2_total_h2ase_ge"]["rmse"],
                    "mae": t1["targets"]["log2_total_h2ase_ge"]["mae"],
                    "r2": t1["targets"]["log2_total_h2ase_ge"]["r2"],
                    "pearson_r": t1["targets"]["log2_total_h2ase_ge"]["pearson_r"],
                    "spearman_rho": t1["targets"]["log2_total_h2ase_ge"]["spearman_rho"],
                    "p_value": t1["targets"]["log2_total_h2ase_ge"]["spearman_p"],
                },
                {
                    "tier": "Tier 1 (5-Fold Spatial Block CV)",
                    "metric_target": f"log1p(viral_richness) [n={t1['targets']['log1p_viral_richness']['n']} WGS]",
                    "rmse": t1["targets"]["log1p_viral_richness"]["rmse"],
                    "mae": t1["targets"]["log1p_viral_richness"]["mae"],
                    "r2": t1["targets"]["log1p_viral_richness"]["r2"],
                    "pearson_r": t1["targets"]["log1p_viral_richness"]["pearson_r"],
                    "spearman_rho": t1["targets"]["log1p_viral_richness"]["spearman_rho"],
                    "p_value": t1["targets"]["log1p_viral_richness"]["spearman_p"],
                },
                {
                    "tier": "Tier 1 (5-Fold Spatial Block CV)",
                    "metric_target": f"log1p(plasmid_richness) [n={t1['targets']['log1p_plasmid_richness']['n']} WGS]",
                    "rmse": t1["targets"]["log1p_plasmid_richness"]["rmse"],
                    "mae": t1["targets"]["log1p_plasmid_richness"]["mae"],
                    "r2": t1["targets"]["log1p_plasmid_richness"]["r2"],
                    "pearson_r": t1["targets"]["log1p_plasmid_richness"]["pearson_r"],
                    "spearman_rho": t1["targets"]["log1p_plasmid_richness"]["spearman_rho"],
                    "p_value": t1["targets"]["log1p_plasmid_richness"]["spearman_p"],
                },
                {
                    "tier": "Tier 1 (5-Fold Spatial Block CV)",
                    "metric_target": "shannon_16s [n=180 site-compartments]",
                    "rmse": t1["targets"]["shannon_16s_site_compartment"]["rmse"],
                    "mae": t1["targets"]["shannon_16s_site_compartment"]["mae"],
                    "r2": t1["targets"]["shannon_16s_site_compartment"]["r2"],
                    "pearson_r": t1["targets"]["shannon_16s_site_compartment"]["pearson_r"],
                    "spearman_rho": t1["targets"]["shannon_16s_site_compartment"]["spearman_rho"],
                    "p_value": t1["targets"]["shannon_16s_site_compartment"]["spearman_p"],
                },
                {
                    "tier": "Tier 2A (Leave-One-Campaign-Out CV)",
                    "metric_target": f"log2(hydB_1l_ge) across Trips 1..5 [n={t2['pooled_metrics']['log2_hydB_1l']['n']}]",
                    "rmse": t2["pooled_metrics"]["log2_hydB_1l"]["rmse"],
                    "mae": t2["pooled_metrics"]["log2_hydB_1l"]["mae"],
                    "r2": t2["pooled_metrics"]["log2_hydB_1l"]["r2"],
                    "pearson_r": t2["pooled_metrics"]["log2_hydB_1l"]["pearson_r"],
                    "spearman_rho": t2["pooled_metrics"]["log2_hydB_1l"]["spearman_rho"],
                    "p_value": t2["pooled_metrics"]["log2_hydB_1l"]["spearman_p"],
                },
                {
                    "tier": "Tier 2B (Prospective Trips 1..4 -> Trip 5 Forecast)",
                    "metric_target": f"log2(hydB_1l_ge) Trip 5 [n={t5_bench['n']} site-compartments]",
                    "rmse": t5_bench["rmse"],
                    "mae": t5_bench["mae"],
                    "r2": t5_bench["r2"],
                    "pearson_r": t5_bench["pearson_r"],
                    "spearman_rho": t5_bench["spearman_rho"],
                    "p_value": t5_bench["spearman_p"],
                },
            ]
            pd.DataFrame(summary_rows).to_csv(
                output_dir / "validation_metrics_summary.tsv", sep="\t", index=False
            )

        return report
