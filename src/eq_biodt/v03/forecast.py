"""Campaign-blocked backtests and immutable, scoreable prospective forecasts."""
from __future__ import annotations

from datetime import date, datetime, timezone
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.compose import ColumnTransformer
from sklearn.linear_model import Ridge
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import OneHotEncoder, StandardScaler

from .evidence import sha256

TARGET = "total_h2_marker_ge"


def visits(panel: pd.DataFrame) -> pd.DataFrame:
    columns = ["site", "compartment", "trip_num", "collection_date", "date_substitution_days", "lat", "lon", TARGET]
    frame = panel[columns].dropna(subset=["site", "compartment", "trip_num", "collection_date", "lat", "lon", TARGET]).copy()
    # date_substitution_days is the correction from a nominal trip date to the
    # curated actual collection date, not uncertainty in that actual date.
    frame["site"] = frame.site.astype(str)
    frame["collection_date"] = pd.to_datetime(frame.collection_date)
    if (frame[TARGET] < 0).any() or not np.isfinite(frame[TARGET]).all():
        raise ValueError("Invalid observed genome equivalents")
    # Biological replicates are aggregated before testing; visits are the unit.
    frame = frame.groupby(["site", "compartment", "trip_num"], as_index=False).agg(
        collection_date=("collection_date", "max"), lat=("lat", "mean"), lon=("lon", "mean"),
        total_h2_marker_ge=(TARGET, "mean"))
    frame["y"] = np.log1p(frame[TARGET])
    return frame.sort_values(["collection_date", "site", "compartment"]).reset_index(drop=True)


def features(frame: pd.DataFrame) -> pd.DataFrame:
    x = frame[["lat", "lon", "compartment"]].copy()
    day = pd.to_datetime(frame.collection_date).dt.day_of_year
    x["season_sin"] = np.sin(2*np.pi*day/365.25)
    x["season_cos"] = np.cos(2*np.pi*day/365.25)
    return x


def predict(train: pd.DataFrame, target: pd.DataFrame) -> dict[str, np.ndarray]:
    if train.empty:
        raise ValueError("Training visits required")
    transformer = ColumnTransformer([
        ("numeric", StandardScaler(), ["lat", "lon"]),
        # Cyclic features already have a fixed [-1,1] scale. Dividing by the
        # tiny variance of a three-day first campaign causes unstable extrapolation.
        ("season", "passthrough", ["season_sin", "season_cos"]),
        ("compartment", OneHotEncoder(handle_unknown="ignore"), ["compartment"])])
    model = make_pipeline(transformer, Ridge(alpha=10.))
    model.fit(features(train), train.y)
    compartment = train.groupby("compartment").y.mean()
    fallback = target.compartment.map(compartment).fillna(train.y.mean()).to_numpy()
    latest = train.sort_values("collection_date").drop_duplicates(["site", "compartment"], keep="last")
    history = dict(zip(zip(latest.site, latest.compartment), latest.y))
    persistence = np.array([history.get((r.site, r.compartment), fallback[i]) for i,r in enumerate(target.itertuples())])
    return {"geography_season_ridge": np.maximum(model.predict(features(target)), 0),
            "persistence": persistence, "compartment_mean": fallback,
            "global_mean": np.repeat(train.y.mean(), len(target))}


def backtest(panel: pd.DataFrame) -> dict:
    frame = visits(panel)
    predictions = []
    for campaign in sorted(frame.trip_num.unique())[1:]:
        test = frame.loc[frame.trip_num.eq(campaign)]
        # All data must precede the first sample of the held-out campaign.
        train = frame.loc[(frame.trip_num < campaign) & (frame.collection_date < test.collection_date.min())]
        for name, values in predict(train, test).items():
            for row, value in zip(test.itertuples(), values):
                predictions.append({"campaign": int(campaign), "site": row.site, "compartment": row.compartment,
                                    "model": name, "observed_log1p_ge": row.y, "predicted_log1p_ge": float(value),
                                    "training_visits": len(train), "training_end": train.collection_date.max().date().isoformat()})
    metrics = []
    pred = pd.DataFrame(predictions)
    if not pred.empty:
        for (campaign, model), group in pred.groupby(["campaign", "model"]):
            residual = group.predicted_log1p_ge-group.observed_log1p_ge
            metrics.append({"campaign": int(campaign), "model": model, "n_visits": len(group),
                            "mae_log1p_ge": float(abs(residual).mean()), "rmse_log1p_ge": float(np.sqrt((residual**2).mean()))})
    return {"target": TARGET, "unit": "log(1+marker genome equivalents)",
            "split": "whole_future_campaign", "metrics": metrics, "predictions": predictions,
            "reference_scope": "Retrospective fixed 2026 MAG/marker catalogue; catalogue discovery is not time-held-out.",
            "interpretation": "Genomic-capacity prediction; causal flux accuracy requires independent experiments."}


def digest(payload: dict) -> str:
    return hashlib.sha256(json.dumps(payload, sort_keys=True, allow_nan=False).encode()).hexdigest()


def freeze(panel_path: Path, target_date: str, output: Path, *, issued_on: str | None = None,
           window_days: int = 14) -> dict:
    issued = date.fromisoformat(issued_on) if issued_on else datetime.now(timezone.utc).date()
    target = date.fromisoformat(target_date)
    if target <= issued or window_days < 0:
        raise ValueError("Forecast target must be after issue date; window must be nonnegative")
    panel = pd.read_csv(panel_path, sep="\t")
    frame = visits(panel)
    train = frame.loc[frame.collection_date.dt.date <= issued]
    latest = train.drop_duplicates(["site", "compartment"], keep="last").copy()
    latest["last_observed"] = latest.collection_date.dt.date.astype(str)
    latest["collection_date"] = pd.Timestamp(target)
    forecasts = predict(train, latest)
    rows = []
    for i, row in enumerate(latest.itertuples()):
        rows.append({"site": row.site, "compartment": row.compartment, "last_observed": row.last_observed,
                     "days_since_observation": (target-date.fromisoformat(row.last_observed)).days,
                     "predictions_ge": {name: float(np.expm1(v[i])) for name,v in forecasts.items()}})
    payload = {"version": "0.3.0", "issued_on": issued.isoformat(), "target_date": target_date,
               "window_days": window_days, "source_sha256": sha256(panel_path),
               "target": TARGET, "training_visits": len(train), "predictions": rows,
               "method": "ridge alpha=10; fixed geography, season, compartment; train-only preprocessing",
               "uncertainty": "No calibrated predictive intervals; compare frozen baseline forecasts.",
               "reference_catalogue": "fixed corrected v2 annotations (2026)",
               "status": "test_fixture" if issued_on is not None else "frozen_prospective_forecast"}
    payload["sha256"] = digest(payload)
    # Never overwrite a registered forecast.
    with output.open("x") as f:
        json.dump(payload, f, indent=2, allow_nan=False)
        f.write("\n")
    return payload


def score(frozen: dict, observed: pd.DataFrame) -> dict:
    payload = {k:v for k,v in frozen.items() if k != "sha256"}
    if digest(payload) != frozen.get("sha256"):
        raise ValueError("Frozen forecast checksum mismatch")
    target = pd.Timestamp(frozen["target_date"])
    issued = pd.Timestamp(frozen["issued_on"])
    data = observed.copy()
    data["collection_date"] = pd.to_datetime(data.collection_date)
    data = data.loc[(data.collection_date > issued) &
                    ((data.collection_date-target).dt.days.abs() <= frozen["window_days"])]
    data["site"] = data.site.astype(str)
    data = data.groupby(["site", "compartment"], as_index=False)[TARGET].mean()
    errors = []
    for row in frozen["predictions"]:
        match = data.loc[data.site.eq(row["site"]) & data.compartment.eq(row["compartment"])]
        if match.empty:
            continue
        value = float(match.iloc[0][TARGET])
        if not np.isfinite(value) or value < 0:
            raise ValueError("Invalid validation observation")
        for name, predicted in row["predictions_ge"].items():
            errors.append({"model": name, "error": float(np.log1p(predicted)-np.log1p(value))})
    results = []
    if errors:
        for name, group in pd.DataFrame(errors).groupby("model"):
            results.append({"model": name, "n_visits": len(group),
                            "mae_log1p_ge": float(group.error.abs().mean()),
                            "rmse_log1p_ge": float(np.sqrt((group.error**2).mean()))})
    return {"forecast_sha256": frozen["sha256"], "metrics": results,
            "status": "scored" if results else "awaiting_matching_future_observations"}
