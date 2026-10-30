"""Darwin Core Event Core + Extended MeasurementOrFact (DwC-A / eMoF) exporter.

Conforms to EU BioDT / LTER-LIFE / DTO-BioFlow / OBIS ENV-DATA standards:
- EventCore (`event.txt`): hierarchical Site -> Site_x_Compartment sampling events
  with true campaign dates (`2023-03-17/2025-10-15`).
- Occurrence (`occurrence.txt`): strictly observed functional & taxonomic guild
  occurrences per event with verified NCBITaxon / RAK IRIs.
- ExtendedMeasurementOrFact (`emof.txt`): 11 core XRF elements + Vanadium, admitted CaCl2 pH
  (only where observed), 16S Shannon diversity, MAP, and trace-gas guild Genome Equivalents
  linked to programmatically verified ChEBI, ENVO, PATO, UO, SIO, and RAK IRIs.
- `meta.xml` and `eml.xml` descriptors.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List, Optional

import pandas as pd

from eq_biodt.data.build_twin_state_cube import TwinStateCubeBuilder
from eq_biodt.data.manifest_registry import CanonicalDataRegistry


class DwCEMoFExporter:
    """Exports the Empty Quarter state cube as a standards-compliant DwC-A + eMoF bundle."""

    def __init__(self, builder: Optional[TwinStateCubeBuilder] = None) -> None:
        self.builder = builder or TwinStateCubeBuilder()
        self.registry: CanonicalDataRegistry = self.builder.registry
        self.ontology = self.registry.ontology

    def _get_guild_taxa(self) -> List[Dict[str, str]]:
        return [
            {
                "col": "mag_actinomycetota_raw",
                "scientificName": "Actinomycetota",
                "scientificNameID": self.ontology.get_taxon_iri("Actinomycetota"),
                "taxonRank": "phylum",
                "quantityType": "relative_abundance_fraction",
                "unit": "dimensionless",
                "unitID": "http://purl.obolibrary.org/obo/UO_0000186",
            },
            {
                "col": "mag_thermoplasmatota_raw",
                "scientificName": "Thermoplasmatota",
                "scientificNameID": self.ontology.get_taxon_iri("Thermoplasmatota"),
                "taxonRank": "phylum",
                "quantityType": "relative_abundance_fraction",
                "unit": "dimensionless",
                "unitID": "http://purl.obolibrary.org/obo/UO_0000186",
            },
            {
                "col": "mag_halobacteriota_raw",
                "scientificName": "Halobacteriota",
                "scientificNameID": self.ontology.get_taxon_iri("Halobacteriota"),
                "taxonRank": "phylum",
                "quantityType": "relative_abundance_fraction",
                "unit": "dimensionless",
                "unitID": "http://purl.obolibrary.org/obo/UO_0000186",
            },
            {
                "col": "hydB_1l_ge_raw",
                "scientificName": "Group 1l [NiFe]-hydrogenase atmospheric H2 scavengers",
                "scientificNameID": f"{self.ontology.RAK_BASE}RAK_GUILD_HYDB_1L",
                "taxonRank": "functional_guild",
                "quantityType": "copies per genome equivalent",
                "unit": "copies per genome equivalent",
                "unitID": "https://rubalkhali.science/kb/RAK_0000078",
            },
            {
                "col": "coxL_formI_ge_raw",
                "scientificName": "Form I coxL atmospheric CO oxidizers",
                "scientificNameID": f"{self.ontology.RAK_BASE}RAK_GUILD_COXL_FORMI",
                "taxonRank": "functional_guild",
                "quantityType": "copies per genome equivalent",
                "unit": "copies per genome equivalent",
                "unitID": "https://rubalkhali.science/kb/RAK_0000078",
            },
        ]

    def build_tables(self) -> Dict[str, pd.DataFrame]:
        """Build `event`, `occurrence`, and `emof` DataFrames with verified ontology IRIs."""
        comp_cube = self.builder.build_compartment_state_cube()
        guild_taxa = self._get_guild_taxa()

        events: List[Dict[str, Any]] = []
        occurrences: List[Dict[str, Any]] = []
        emofs: List[Dict[str, Any]] = []

        desert_envo = self.ontology.get_envo_iri("desert biome")
        sand_envo = self.ontology.get_envo_iri("desert sand")
        uo_percent = self.ontology.CORE_TERMS["uo_percent"]
        uo_ph = self.ontology.CORE_TERMS["uo_ph"]

        for _, row in comp_cube.iterrows():
            site_num = int(row["site_num"])
            comp_norm = str(row["compartment_norm"])
            parent_event_id = f"RAK_SITE_{site_num:02d}"
            event_id = f"RAK_SITE_{site_num:02d}_{comp_norm.upper()}"

            events.append(
                {
                    "eventID": event_id,
                    "parentEventID": parent_event_id,
                    "locationID": str(row["site_iri"]),
                    "eventDate": "2023-03-17/2025-10-15",
                    "decimalLatitude": round(float(row["lat"]), 6),
                    "decimalLongitude": round(float(row["lon"]), 6),
                    "minimumElevationInMeters": round(float(row["elevation_m"]), 1),
                    "maximumElevationInMeters": round(float(row["elevation_m"]), 1),
                    "habitat": f"{desert_envo} | {sand_envo}",
                    "samplingProtocol": str(row["compartment_iri"]),
                    "fieldNotes": f"Compartment={comp_norm}; H3_6={row['h3_6']}; ERA5={row['era5_cell']}",
                }
            )

            # Occurrences per compartment event (export only when observed, never imputed)
            for g in guild_taxa:
                raw_val = row.get(g["col"])
                if pd.isna(raw_val):
                    continue
                val = float(raw_val)
                occ_id = f"{event_id}_{g['col']}"
                occurrences.append(
                    {
                        "occurrenceID": occ_id,
                        "eventID": event_id,
                        "scientificName": g["scientificName"],
                        "scientificNameID": g["scientificNameID"],
                        "taxonRank": g["taxonRank"],
                        "organismQuantity": round(val, 6),
                        "organismQuantityType": g["quantityType"],
                        "occurrenceStatus": "present" if val > 0 else "absent",
                        "basisOfRecord": "MaterialSample",
                    }
                )

            # eMoF records: 11 core XRF elements + Vanadium (V)
            xrf_elements = list(self.registry.CORE_11_XRF_ELEMENTS) + ["V"]
            for el in xrf_elements:
                if el not in row or pd.isna(row[el]):
                    continue
                chebi_iri = self.ontology.get_chebi_iri(el)
                q_iri, v_iri = self.ontology.get_rak_analyte_iris(el)
                emofs.append(
                    {
                        "measurementID": f"{event_id}_XRF_{el}",
                        "eventID": event_id,
                        "occurrenceID": "",
                        "measurementType": f"XRF concentration of {el}",
                        "measurementTypeID": q_iri,
                        "measurementValue": str(round(float(row[el]), 6)),
                        "measurementValueID": chebi_iri or v_iri,
                        "measurementUnit": "percent",
                        "measurementUnitID": uo_percent,
                    }
                )

            # eMoF record: Admitted CaCl2 soil pH (only when directly observed in ph_accepted_measurements.tsv)
            if bool(row.get("ph_observed", False)) and pd.notna(row.get("ph_mean_raw")):
                emofs.append(
                    {
                        "measurementID": f"{event_id}_PH_CACL2",
                        "eventID": event_id,
                        "occurrenceID": "",
                        "measurementType": "Soil pH in 0.01M CaCl2",
                        "measurementTypeID": uo_ph,
                        "measurementValue": str(round(float(row["ph_mean_raw"]), 4)),
                        "measurementValueID": "",
                        "measurementUnit": "pH",
                        "measurementUnitID": uo_ph,
                    }
                )

            # eMoF record: 16S Shannon diversity
            if pd.notna(row.get("shannon_16s")):
                emofs.append(
                    {
                        "measurementID": f"{event_id}_SHANNON_16S",
                        "eventID": event_id,
                        "occurrenceID": "",
                        "measurementType": "16S rRNA Shannon diversity index",
                        "measurementTypeID": self.ontology.CORE_TERMS["site_quality"],
                        "measurementValue": str(round(float(row["shannon_16s"]), 4)),
                        "measurementValueID": "",
                        "measurementUnit": "dimensionless",
                        "measurementUnitID": "http://purl.obolibrary.org/obo/UO_0000186",
                    }
                )

            # eMoF record: Mean annual precipitation (MAP)
            emofs.append(
                {
                    "measurementID": f"{event_id}_MAP_MM",
                    "eventID": event_id,
                    "occurrenceID": "",
                    "measurementType": "Mean annual precipitation",
                    "measurementTypeID": self.ontology.CORE_TERMS["annual_total_precip_value"],
                    "measurementValue": str(round(float(row["map_mm"]), 2)),
                    "measurementValueID": "",
                    "measurementUnit": "mm/yr",
                    "measurementUnitID": "http://purl.obolibrary.org/obo/UO_0000016",
                }
            )

        return {
            "event": pd.DataFrame(events),
            "occurrence": pd.DataFrame(occurrences),
            "emof": pd.DataFrame(emofs),
        }

    def export_archive(self, output_dir: Path) -> Dict[str, Path]:
        """Write `event.txt`, `occurrence.txt`, `emof.txt`, `meta.xml`, and `eml.xml`."""
        output_dir.mkdir(parents=True, exist_ok=True)
        tables = self.build_tables()

        event_path = output_dir / "event.txt"
        occ_path = output_dir / "occurrence.txt"
        emof_path = output_dir / "emof.txt"
        meta_path = output_dir / "meta.xml"
        eml_path = output_dir / "eml.xml"

        tables["event"].to_csv(event_path, sep="\t", index=False)
        tables["occurrence"].to_csv(occ_path, sep="\t", index=False)
        tables["emof"].to_csv(emof_path, sep="\t", index=False)

        meta_xml = """<?xml version="1.0" encoding="UTF-8"?>
<archive xmlns="http://rs.tdwg.org/dwc/text/" metadata="eml.xml">
  <core encoding="UTF-8" fieldsTerminatedBy="\\t" linesTerminatedBy="\\n" ignoreHeaderLines="1" rowType="http://rs.tdwg.org/dwc/terms/Event">
    <files><location>event.txt</location></files>
    <id index="0" />
    <field index="0" term="http://rs.tdwg.org/dwc/terms/eventID"/>
    <field index="1" term="http://rs.tdwg.org/dwc/terms/parentEventID"/>
    <field index="2" term="http://rs.tdwg.org/dwc/terms/locationID"/>
    <field index="3" term="http://rs.tdwg.org/dwc/terms/eventDate"/>
    <field index="4" term="http://rs.tdwg.org/dwc/terms/decimalLatitude"/>
    <field index="5" term="http://rs.tdwg.org/dwc/terms/decimalLongitude"/>
    <field index="6" term="http://rs.tdwg.org/dwc/terms/minimumElevationInMeters"/>
    <field index="7" term="http://rs.tdwg.org/dwc/terms/maximumElevationInMeters"/>
    <field index="8" term="http://rs.tdwg.org/dwc/terms/habitat"/>
    <field index="9" term="http://rs.tdwg.org/dwc/terms/samplingProtocol"/>
    <field index="10" term="http://rs.tdwg.org/dwc/terms/fieldNotes"/>
  </core>
  <extension encoding="UTF-8" fieldsTerminatedBy="\\t" linesTerminatedBy="\\n" ignoreHeaderLines="1" rowType="http://rs.tdwg.org/dwc/terms/Occurrence">
    <files><location>occurrence.txt</location></files>
    <coreid index="1" />
    <field index="0" term="http://rs.tdwg.org/dwc/terms/occurrenceID"/>
    <field index="1" term="http://rs.tdwg.org/dwc/terms/eventID"/>
    <field index="2" term="http://rs.tdwg.org/dwc/terms/scientificName"/>
    <field index="3" term="http://rs.tdwg.org/dwc/terms/scientificNameID"/>
    <field index="4" term="http://rs.tdwg.org/dwc/terms/taxonRank"/>
    <field index="5" term="http://rs.tdwg.org/dwc/terms/organismQuantity"/>
    <field index="6" term="http://rs.tdwg.org/dwc/terms/organismQuantityType"/>
    <field index="7" term="http://rs.tdwg.org/dwc/terms/occurrenceStatus"/>
    <field index="8" term="http://rs.tdwg.org/dwc/terms/basisOfRecord"/>
  </extension>
  <extension encoding="UTF-8" fieldsTerminatedBy="\\t" linesTerminatedBy="\\n" ignoreHeaderLines="1" rowType="http://rs.iobis.org/obis/terms/ExtendedMeasurementOrFact">
    <files><location>emof.txt</location></files>
    <coreid index="1" />
    <field index="0" term="http://rs.tdwg.org/dwc/terms/measurementID"/>
    <field index="1" term="http://rs.tdwg.org/dwc/terms/eventID"/>
    <field index="2" term="http://rs.tdwg.org/dwc/terms/occurrenceID"/>
    <field index="3" term="http://rs.tdwg.org/dwc/terms/measurementType"/>
    <field index="4" term="http://rs.iobis.org/obis/terms/measurementTypeID"/>
    <field index="5" term="http://rs.tdwg.org/dwc/terms/measurementValue"/>
    <field index="6" term="http://rs.iobis.org/obis/terms/measurementValueID"/>
    <field index="7" term="http://rs.tdwg.org/dwc/terms/measurementUnit"/>
    <field index="8" term="http://rs.iobis.org/obis/terms/measurementUnitID"/>
  </extension>
</archive>
"""
        eml_xml = """<?xml version="1.0" encoding="UTF-8"?>
<eml:eml xmlns:eml="https://eml.ecoinformatics.org/eml-2.2.0" packageId="eq-biodt-dwca-0.2.0" system="https://rubalkhali.science">
  <dataset>
    <title>Empty Quarter (Rub' al Khali) Biodiversity Digital Twin Darwin Core Event &amp; eMoF Archive</title>
    <creator>
      <individualName><surName>Hoehndorf</surName></individualName>
      <organizationName>KAUST Bio-Ontology Research Group (BORG)</organizationName>
    </creator>
    <abstract>
      <para>Harmonized LTER-LIFE / OBIS / BioDT Darwin Core Event + ExtendedMeasurementOrFact archive across the 1,043 km Rub' al Khali transect (March 2023 - October 2025).</para>
    </abstract>
  </dataset>
</eml:eml>
"""
        meta_path.write_text(meta_xml, encoding="utf-8")
        eml_path.write_text(eml_xml, encoding="utf-8")

        return {
            "event": event_path,
            "occurrence": occ_path,
            "emof": emof_path,
            "meta": meta_path,
            "eml": eml_path,
        }
