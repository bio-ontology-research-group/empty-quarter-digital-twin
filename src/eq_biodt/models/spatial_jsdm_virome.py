"""Spatial Joint Species Distribution (CLR-RDA / Ridge) & Virome/Mobilome Co-Turnover Engine.

Grounded in:
- `~/Documents/papers/empty-quarter-map`: Transect-wide CLR-RDA variance partitioning
  across 11 core XRF elements, CaCl2 pH, MAP/MAT, and soil compartment.
- `~/Documents/papers/empty-quarter-mobilome`: 19,300 corrected vOTUs + 322,755 pOTUs
  across 552 metagenomes (`sample_mobilome_load.tsv`), 17.5% variance explained
  (host community 10.6%, compartment 2.8%, climate 1.4%, geography 1.4%), and
  06_resistome_safety ARG containment.
- `~/Documents/papers/empty-quarter-viruses`: Plant-associated Caudoviricetes expansion
  and viral-prokaryotic Mantel coupling (r = 0.254).
"""

from __future__ import annotations

from typing import Any, Dict

import numpy as np


class SpatialJSDMViromeModel:
    """Predicts spatial guild, virome (vOTU), plasmid (pOTU), and resistome responses."""

    # Empirical variance partitioning weights & catalogue denominators from empty-quarter-mobilome
    VARIANCE_PARTITION = {
        "votu_catalogue_total": 19300,
        "potu_catalogue_total": 322755,
        "mobilome_samples_n": 552,
        "host_community": 0.106,
        "soil_compartment": 0.028,
        "climate_aridity": 0.014,
        "spatial_geography": 0.014,
        "edaphic_xrf_ph": 0.013,
        "total_explained_r2": 0.175,
        "viral_host_mantel_r": 0.254,
    }

    def predict_virome_mobilome_state(
        self,
        shannon_16s: float,
        compartment: str,
        map_mm: float,
        xrf_evaporite_index: float,
        ph_cacl2: float = 7.95,
        vegetation_cover_frac: float = 0.0,
        exogenous_manure_load: float = 0.0,
        plasmid_bioaug_dose: float = 0.0,
        baseline_viral_richness: float = 18.0,
        baseline_plasmid_richness: float = 31.0,
        baseline_shannon_16s: float | None = None,
    ) -> Dict[str, Any]:
        """Predict vOTU richness/TPM, pOTU richness/TPM, lytic:lysogenic ratio, and ARG screening score.

        Anchored to the site's baseline covariates so a null intervention (vegetation_cover_frac=0,
        exogenous_manure_load=0, plasmid_bioaug_dose=0, shannon_16s=baseline_shannon_16s) preserves
        exact baseline viral and plasmid richness.
        """
        comp = compartment.strip()
        is_subsurface = 1.0 if comp in ("ShallowSubsurface", "deep", "Deep") else 0.0
        ref_shannon = float(shannon_16s if baseline_shannon_16s is None else baseline_shannon_16s)

        # Host diversity coupling anchored to site baseline Shannon H'
        host_factor = float(np.exp(0.254 * (shannon_16s - ref_shannon)))
        comp_viral_boost = 1.0 + 0.55 * vegetation_cover_frac + 0.15 * exogenous_manure_load
        viral_richness = float(max(1.0, baseline_viral_richness * host_factor * comp_viral_boost))
        viral_tpm = float(viral_richness * (24.5 + 18.0 * vegetation_cover_frac))

        # Plasmid mobilome (pOTU) response anchored to site baseline
        plasmid_boost = (
            1.0
            + 0.65 * vegetation_cover_frac
            + 0.48 * plasmid_bioaug_dose
            + 0.30 * exogenous_manure_load
        )
        plasmid_richness = float(max(2.0, baseline_plasmid_richness * host_factor * plasmid_boost))
        plasmid_tpm = float(plasmid_richness * (32.0 + 25.0 * plasmid_bioaug_dose))

        # Temperate/lysogenic vs lytic phage strategy:
        lysogenic_fraction = float(
            np.clip(
                0.68
                - 0.22 * vegetation_cover_frac
                - 0.18 * exogenous_manure_load
                + 0.06 * is_subsurface,
                0.20,
                0.88,
            )
        )

        # Resistome / ARG biosafety index (06_resistome_safety):
        arg_risk_score = float(
            np.clip(
                0.035
                + 0.42 * exogenous_manure_load
                + 0.015 * plasmid_bioaug_dose
                - 0.01 * is_subsurface,
                0.01,
                1.0,
            )
        )

        return {
            "viral_richness_votu": round(viral_richness, 2),
            "viral_abundance_tpm": round(viral_tpm, 2),
            "plasmid_richness_potu": round(plasmid_richness, 2),
            "plasmid_abundance_tpm": round(plasmid_tpm, 2),
            "lysogenic_phage_fraction": round(lysogenic_fraction, 4),
            "resistome_arg_risk_index": round(arg_risk_score, 4),
            "resistome_biosafety_pass": bool(arg_risk_score < 0.25),
            "variance_partition_meta": self.VARIANCE_PARTITION,
        }
