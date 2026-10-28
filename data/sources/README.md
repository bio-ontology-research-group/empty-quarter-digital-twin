# Analysis input snapshot

These are the inputs required to reproduce EQ-BioDT v0.4.0 from the annotated
MAG catalogue and read summaries. `manifest.json` records each file's SHA-256,
the original uncompressed digest where applicable, and external source URLs.
The evidence builder validates all snapshot digests before running.

- `rep_ko_v2.tsv.gz`: corrected, post-refinement MAG KO annotations.
- `refined_mag_manifest.tsv.gz`: MAG quality and GTDB taxonomy.
- `atmo_classification_v2.tsv.gz`: trace-gas marker assignments.
- `coverm_rep_abundance_long.tsv.gz`: per-sample MAG and unmapped read percentages.
- `guild_genome_equivalents_v2_r3noflagged.tsv.gz`: read-screen marker genome equivalents.
- `eq_wgs_validation_panel.parquet`: frozen metadata panel; evidence building uses sample, site, compartment, campaign, location and curated collection dates. Other legacy columns are not prediction inputs.
- `M00165.txt`: KEGG M00165 module snapshot, retrieved 2026-10-28 from the official REST endpoint. KEGG is the source of the pathway definition and retains its rights. Standard KO identifiers and module alternatives are used for annotation interpretation.
- `reji-bottles.txt`: dilution-corrected observations from Reji et al., *Global implications of a low soil moisture threshold for microbial hydrogen uptake*, Nature Communications 17, 515 (2026), https://doi.org/10.1038/s41467-025-67208-3. Upstream commit `1f860701b8572d3bef7f0657bd5e40862feca477` in [Reji2025_H2_SoilMoisture](https://github.com/Linta-Reji/Reji2025_H2_SoilMoisture). Original authors retain ownership and applicable data terms.

The first six inputs are the project’s released analysis inputs. This package
reproduces the Twin analyses starting from those annotations and summaries;
it does not rerun sequencing, assembly, bin refinement or primary gene annotation.
No sibling checkout, institutional storage or private metadata lookup is needed.
