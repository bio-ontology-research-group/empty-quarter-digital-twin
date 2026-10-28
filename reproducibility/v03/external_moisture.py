"""Reanalyse published bottle experiments; do not transfer absolute rates to soils.

Primary source: Reji et al., doi:10.1038/s41467-025-67208-3.
Fixed-C0 exponential fits follow the authors' corrected analysis. Rates can be
negative; air-dried observations are fitted here rather than assigned zero.
"""
from pathlib import Path
import json
import urllib.request

import numpy as np
import pandas as pd
from scipy.optimize import least_squares

from eq_biodt.v03.evidence import sha256, source_directory
from eq_biodt.v03.process import Parameters, activation

URL = "https://raw.githubusercontent.com/Linta-Reji/Reji2025_H2_SoilMoisture/1f860701b8572d3bef7f0657bd5e40862feca477/data/H2_experiments_allData_finalCorrected.txt"
SHA = "44bd240facea7db00348fe49046560784652b72a0b9cf43c84c863f80d48dcd0"


def analyse(output: Path, source: Path | None = None) -> dict:
    output.mkdir(parents=True, exist_ok=True)
    if source is None:
        source = source_directory()/"reji-bottles.txt"
    if sha256(source) != SHA:
        raise ValueError("External data digest mismatch")
    raw = pd.read_csv(source, sep="\t")
    rows = []
    for key,g in raw.groupby(["soil_type", "sample_id", "replicate"]):
        g = g.sort_values("time_min")
        t = g.time_min.to_numpy(float)
        y = g.h2_ppmv_dil_corr.to_numpy(float)
        if len(g) < 3 or not np.isfinite(y).all() or y[0] <= 0:
            raise ValueError(f"Invalid bottle series {key}")
        t = t-t[0]
        fit = least_squares(lambda x: y[0]*np.exp(np.clip(-x[0]*t, -700, 50))-y,
                            [.001], bounds=([-1/max(t)], [1.]), xtol=1e-12, gtol=1e-12, ftol=1e-12)
        if not fit.success:
            raise RuntimeError(f"Bottle fit failed: {key}")
        row = g.iloc[0]
        porosity = float(row.porosity)
        theta = float(row.volumetric_moisture_perDrySoil_rhoW)
        psi = float(row.WP)
        # Supersaturated bottles retain their observations; the gas-pore model
        # assigns zero air space and exposes its missing aqueous transport.
        diffusion = max(0., porosity-theta)**(10/3)/porosity**2
        p = Parameters(k_h2_s=.001, k_co_s=.001)
        # Directional moisture-envelope test only, not a chamber transport model.
        rows.append({"soil":key[0], "sample_id":key[1], "replicate":int(key[2]),
                     "psi_mpa":psi, "theta":theta, "porosity":porosity,
                     "bottle_k_min":float(fit.x[0]), "n_timepoints":len(g),
                     "diffusion_only":diffusion, "activation_diffusion":diffusion*activation(psi,p)})
    rates = pd.DataFrame(rows)
    rates.to_csv(output/"external_bottle_rates.tsv",sep="\t",index=False)
    fits = []
    for soil,g in rates.groupby("soil"):
        for rep in sorted(g.replicate.unique()):
            train, test = g.loc[g.replicate.ne(rep)], g.loc[g.replicate.eq(rep)]
            for model in ["diffusion_only", "activation_diffusion"]:
                x = train[model].to_numpy()
                coefficient = max(0., float(x@train.bottle_k_min/(x@x)))
                residual = test[model]*coefficient-test.bottle_k_min
                fits.append({"soil":soil,"heldout_replicate":int(rep),"model":model,
                             "amplitude":coefficient,"rmse_per_min":float(np.sqrt((residual**2).mean())),
                             "test_bottles":len(test)})
    report = {"source_url":URL,"source_sha256":SHA,"bottles":len(rates),"timepoints":len(raw),
              "fixed_threshold_mpa":-85,"fits":fits,
              "supersaturated_bottles":int((rates.theta > rates.porosity).sum()),
              "scope":"Directional moisture-envelope test with leave-replicate-out amplitude fits. Threshold informed by the same published study; not independent threshold validation or local kinetic calibration.",
              "model":"Observed bottle k = fitted nonnegative soil amplitude times D_ratio, with or without the fixed water-potential activation ramp.",
              "air_dried":"Fitted from measured concentrations; never assigned a zero rate."}
    (output/"external_moisture.json").write_text(json.dumps(report,indent=2)+"\n")
    return report


if __name__ == "__main__":
    import argparse
    p=argparse.ArgumentParser();p.add_argument("--source",type=Path);p.add_argument("--output",type=Path,required=True)
    args=p.parse_args();print(json.dumps(analyse(args.output,args.source),indent=2))
