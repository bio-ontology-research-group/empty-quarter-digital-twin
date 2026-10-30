"""Export the external-transfer summary served by /v3/transfer.

Reads the reviewed v0.5 outputs (reviews/v05) and writes data/v03/transfer_summary.json, which is
packaged with the wheel. No number is recomputed here; values are copied from the reviewed tables.
Run from the repository root:  uv run python scripts/web/export_transfer_summary.py
"""
import json
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
REVIEW = ROOT / "reviews" / "v05"
OUT = ROOT / "data" / "v03" / "transfer_summary.json"


def records(frame: pd.DataFrame) -> list[dict]:
    return json.loads(frame.to_json(orient="records"))


def main() -> None:
    ext = json.loads((REVIEW / "external_transfer.json").read_text())
    namib = pd.read_csv(REVIEW / "namib_hydration_contrasts.tsv", sep="\t")
    pred = pd.read_csv(REVIEW / "namib_conditional_predictions.tsv", sep="\t")
    sonoran_lab = pd.read_csv(REVIEW / "sonoran_hydration_contrasts.tsv", sep="\t")
    quantified = namib.dropna(subset=["fold"])
    summary = {
        "source": "reviews/v05 (external_transfer.json and companion tables); protocol in reviews/v05/external-evaluation-protocol.md",
        "seed": ext["seed"],
        "namib": {
            "sites": ext["namib"]["sites"],
            "dry_sites_with_quantified_atmospheric_rates": ext["namib"]["dry_sites_with_quantified_atmospheric_rates"],
            "contrasts_quantified": records(quantified),
            "other_sites": records(namib[namib.fold.isna()][["site", "n_dry_rates", "n_wet_rates", "wet_mean", "status"]]),
            "reference_predictions": records(pred[["site", "assumed_bulk_density_g_cm3", "dry_reference_flux",
                                                    "wet_reference_flux", "predicted_direction"]]),
            "prediction_scope": ext["namib"]["prediction_scope"],
            "reference_failure": ext["namib"]["reference_failure"],
            "measurement_exclusion": ext["namib"]["measurement_exclusion"],
            "rate_unit": "bottle rate in the native units of the source study (not equated with areal flux)",
        },
        "sonoran": {
            "lab_bottles": ext["sonoran"]["lab_bottles"],
            "field_collars": ext["sonoran"]["field_collars"],
            "field_plot_campaign_means": ext["sonoran"]["field_plot_campaign_means"],
            "field_validation": ext["sonoran"]["field_validation"],
            "field_metrics": ext["sonoran"]["field_metrics"],
            "field_metric_unit": "nmol H2 m-2 s-1",
            "lab_contrasts": records(sonoran_lab[["site", "dry_gwc", "wet_gwc", "difference", "lcl", "ucl", "unit"]]),
            "lab_scope": ext["sonoran"]["lab_scope"],
            "lab_default_predictions": ext["sonoran"]["lab_default_predictions"],
        },
        "temperate": {
            "design": "leave-one-soil-out; common response amplitude fitted on two soils, third predicted",
            "rmse": ext["temperate_leave_soil_out"],
            "unit": "bottle first-order rate (min-1)",
        },
    }
    OUT.write_text(json.dumps(summary, indent=1) + "\n")
    print(f"wrote {OUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
