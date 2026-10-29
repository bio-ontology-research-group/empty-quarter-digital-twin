"""Build metagenomic evidence from corrected MAG annotations and read mappings."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

from .semantics import Observation, additive_total

# KEGG M00165 module-block alternatives, snapshot 2026-10-28.
# K01602 is optional in the KEGG definition (minus operator), not a universal
# prerequisite for RuBisCO. The additional three-gene screen remains available.
CBB_STEPS = {
    "carboxylation": ("K01601",),
    "phosphoglycerate_kinase": ("K00927",),
    "gap_dehydrogenase": ("K05298", "K00150", "K00134"),
    "triose_isomerase": ("K01803",),
    "aldolase": ("K01623", "K01624"),
    "fructose_bisphosphatase": ("K03841", "K02446", "K11532", "K01086"),
    "transketolase": ("K00615",),
    # Alternative regeneration routes; transaldolase is not an SBPase enzyme.
    "sedoheptulose_regeneration": ("K01100", "K11532", "K01086", "K00616", "K13810"),
    "ribose_isomerase": ("K01807", "K01808"),
    "ribulose_epimerase": ("K01783",),
    "phosphoribulokinase": ("K00855", "K28462"),
}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024*1024), b""):
            h.update(block)
    return h.hexdigest()


def source_directory() -> Path:
    packaged = Path(__file__).parent / "sources"
    return packaged if packaged.exists() else Path(__file__).resolve().parents[3] / "data/sources"


def build_bundled_evidence(output: Path, sources: Path | None = None) -> dict:
    """Rebuild every evidence column from checksum-verified, portable source tables."""
    sources = sources or source_directory()
    manifest = json.loads((sources / "manifest.json").read_text())["sources"]
    paths = {}
    for key, info in manifest.items():
        path = sources / info["file"]
        if sha256(path) != info["sha256"]:
            raise ValueError(f"Source digest mismatch: {key}")
        paths[key] = path
    return build_evidence(full_ko=paths["corrected_full_MAG_KOs"], quality=paths["MAG_quality"],
                          metadata_panel=paths["metadata_dates"], read_mapping=paths["read_mapping"],
                          marker_classification=paths["gas_marker_classification"],
                          marker_ge=paths["marker_genome_equivalents"], output=output)


def build_evidence(*, full_ko: Path, quality: Path, metadata_panel: Path,
                   read_mapping: Path, marker_classification: Path, marker_ge: Path,
                   output: Path) -> dict:
    output.mkdir(parents=True, exist_ok=True)
    kos = pd.read_csv(full_ko, sep="\t", names=["rep", "kos"], keep_default_na=False)
    if kos.rep.duplicated().any():
        raise ValueError("Duplicate reference MAG annotation")
    ko_sets = kos.kos.map(lambda x: set(x.split(",")) if x else set())
    mag = kos[["rep"]].copy()
    for step, alternatives in CBB_STEPS.items():
        mag[step] = ko_sets.map(lambda s: bool(s.intersection(alternatives)))
    mag["annotation_available"] = kos.kos.ne("")
    mag["cbb_steps_detected"] = mag[list(CBB_STEPS)].sum(axis=1)
    mag["cbb_all_steps_detected"] = mag.cbb_steps_detected.eq(len(CBB_STEPS)) & mag.annotation_available
    mag["sbpase_detected"] = ko_sets.map(lambda s: bool(s.intersection({"K01100", "K11532", "K01086"})))
    mag["transaldolase_detected"] = ko_sets.map(lambda s: bool(s.intersection({"K00616", "K13810"})))
    mag["cbb_canonical_detected"] = mag.cbb_all_steps_detected & mag.sbpase_detected
    mag["cbb_transaldolase_detected"] = mag.cbb_all_steps_detected & mag.transaldolase_detected
    mag["rbcl_detected"] = mag.carboxylation
    mag["rbcl_prk_detected"] = mag.carboxylation & mag.phosphoribulokinase
    mag["three_gene_screen"] = mag.rbcl_prk_detected & ko_sets.map(lambda s: "K01602" in s)
    target_kos = {ko for alts in CBB_STEPS.values() for ko in alts} | {"K01602"}
    mag["cbb_kos"] = ko_sets.map(lambda s: ";".join(sorted(s.intersection(target_kos))))
    mag["missing_steps"] = mag.apply(lambda r: ";".join(k for k in CBB_STEPS if not r[k]), axis=1)
    q = pd.read_csv(quality, sep="\t", usecols=["mag_id", "completeness", "contamination", "quality_tier", "gtdb_classification"])
    for rank, prefix in [("domain", "d__"), ("phylum", "p__"), ("order", "o__"), ("family", "f__"), ("genus", "g__")]:
        q[rank] = q["gtdb_classification"].str.extract(rf"{prefix}([^;]+)").fillna("")
    q = q.drop(columns=["gtdb_classification"])
    mag = mag.merge(q, left_on="rep", right_on="mag_id", how="left", validate="one_to_one")
    if mag.completeness.isna().any():
        raise ValueError("Missing MAG quality; identity join must be exact")
    mag["quality_90_5"] = (mag.completeness >= 90) & (mag.contamination < 5)
    atmo = pd.read_csv(marker_classification, sep="\t")
    for label, assignment in {"h2_1l": "[NiFe]_Group_1l", "h2_1h": "[NiFe]_Group_1h",
                              "h2_2a": "[NiFe]_Group_2a", "co_form1": "form_I|AYxCSFR"}.items():
        mag[label] = mag.rep.isin(set(atmo.loc[atmo.assignment.eq(assignment), "rep_id"]))
    mag["h2_carrier_union"] = mag[["h2_1l", "h2_1h", "h2_2a"]].any(axis=1)
    mag["h2_cbb_all_steps"] = mag.h2_carrier_union & mag.cbb_all_steps_detected
    mag["h2_cbb_quality_90_5"] = mag.h2_cbb_all_steps & mag.quality_90_5
    mag["dual_h2_co_cbb_all_steps"] = mag.h2_cbb_all_steps & mag.co_form1
    mag.to_csv(output/"mag_evidence.tsv", sep="\t", index=False)
    cover = pd.read_csv(read_mapping, sep="\t", names=["sample", "rep", "pct"])
    if cover.duplicated(["sample", "rep"]).any() or (cover.pct < 0).any() or np.isinf(cover.pct).any():
        raise ValueError("Invalid/duplicate read mapping weights")
    missing_mapping_samples = sorted(cover.loc[cover.pct.isna(), "sample"].unique().tolist())
    # Entire observation excluded: NaN mapping weights are not zero abundance.
    cover = cover.loc[~cover["sample"].isin(missing_mapping_samples)].copy()
    mapped = cover.loc[cover.rep.ne("unmapped")].copy()
    if not mapped.rep.isin(mag.rep).all():
        raise ValueError("Mapped reference outside annotation catalogue")
    mapped = mapped.merge(mag, on="rep", how="left", validate="many_to_one")
    total = mapped.groupby("sample").pct.sum()
    panel = total.rename("mapped_read_pct").to_frame()
    panel["unmapped_read_pct"] = cover.loc[cover.rep.eq("unmapped")].set_index("sample").pct
    if panel.unmapped_read_pct.isna().any() or not np.allclose(panel.sum(axis=1), 100, atol=.02, rtol=0):
        raise ValueError("Mapping fractions do not close to 100 percent")
    screens = ["rbcl_detected", "rbcl_prk_detected", "three_gene_screen", "cbb_all_steps_detected",
               "h2_carrier_union", "h2_cbb_all_steps", "h2_cbb_quality_90_5", "quality_90_5",
               "cbb_canonical_detected", "cbb_transaldolase_detected"]
    for label in screens:
        value = mapped.pct.where(mapped[label], 0).groupby(mapped["sample"]).sum()
        panel[label+"_share_mapped"] = value/total.replace(0, np.nan)
        panel[label+"_pct_all_reads"] = value
    # Shannon describes mapped MAG read proportions, with no 16S substitution.
    w = mapped.pct/ mapped["sample"].map(total)
    entropy = -w*np.log(w.where(w > 0, 1))
    panel["mapped_mag_shannon"] = entropy.groupby(mapped["sample"]).sum()
    ge = pd.read_csv(marker_ge, sep="\t")
    ge = ge.loc[ge.id_floor.eq(60), ["sample", "hydB_1l", "hydB_1h", "hydB_2a", "coxL_formI",
                                  "reads_screened", "reads_scg", "scg_genes_detected"]]
    # Metadata fields only: pooled 16S, pH, DNA and filled environmental values are excluded.
    columns = ["sample", "site", "compartment", "trip_num", "collection_date", "date_substitution_days", "lat", "lon"]
    meta = (pd.read_parquet(metadata_panel) if metadata_panel.suffix == ".parquet"
            else pd.read_csv(metadata_panel, sep="\t"))[columns]
    panel = meta.merge(panel.reset_index(), on="sample", how="left", validate="one_to_one").merge(ge, on="sample", how="left", validate="one_to_one")
    def total_marker(row):
        parts = [Observation(None if pd.isna(row[k]) else float(row[k]), "GE",
                             "SCG:"+str(row["sample"]), marker_ge.name,
                             "gene_copy_capacity", k) for k in ["hydB_1l", "hydB_1h", "hydB_2a"]]
        return additive_total(parts, disjoint=True).value
    panel["total_h2_marker_ge"] = panel.apply(total_marker, axis=1)
    panel.to_csv(output/"sample_evidence.tsv", sep="\t", index=False, na_rep="NA")
    inputs = {"corrected_full_MAG_KOs": full_ko, "MAG_quality": quality,
              "read_mapping": read_mapping,
              "gas_marker_classification": marker_classification,
              "marker_genome_equivalents": marker_ge,
              "metadata_dates": metadata_panel}
    manifest = {"version": "0.5.0", "sources": {k: {"file": p.name, "sha256": sha256(p)} for k,p in inputs.items()},
                "reference_MAGs": len(mag), "samples": len(panel),
                "all_CBB_steps_MAGs": int(mag.cbb_all_steps_detected.sum()),
                "h2_CBB_all_steps_MAGs": int(mag.h2_cbb_all_steps.sum()),
                "h2_CBB_quality_90_5_MAGs": int(mag.h2_cbb_quality_90_5.sum()),
                "annotated_MAGs": int(mag.annotation_available.sum()),
                "marker_GE_observed_samples": int(panel.total_h2_marker_ge.notna().sum()),
                "dated_samples": int(panel.collection_date.notna().sum()),
                "missing_mapping_samples": panel.loc[panel.mapped_read_pct.isna(), "sample"].tolist(),
                "source_missing_mapping_samples": missing_mapping_samples,
                "CBB_step_alternatives": CBB_STEPS,
                "endpoint": "detected_genomic_capacity",
                "negative_screen_interpretation": "incomplete_genome_or_undetected_annotation_or_true_absence",
                "activity_absolute_biomass": "unobserved",
                "artifacts": {p.name: sha256(p) for p in (output/"mag_evidence.tsv", output/"sample_evidence.tsv")}}
    (output/"manifest.json").write_text(json.dumps(manifest, indent=2)+"\n")
    return manifest
