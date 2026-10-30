"""OWL 2 ontology generation, Kobayashi-MaRust (km) classification and justification extraction."""
from __future__ import annotations

import json
import os
import shutil
import subprocess
import time
import hashlib
import tempfile
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd

from .evidence import CBB_STEPS, sha256

RAK = "https://rubalkhali.science/kb/"
KEGG = "https://www.kegg.jp/entry/"
OBO = "http://purl.obolibrary.org/obo/"
SIO = "http://semanticscience.org/resource/"
SKOS = "http://www.w3.org/2004/02/skos/core#"

# Pinned reasoner (must equal VERSION/SHA256 in scripts/install_reasoner.py; checked by tests).
KM_VERSION = "1.5.3"
KM_SHA256 = "2470b0415f1a78533024e15b5c317f9e743a98916b597739db93b9452a8dd8d4"
OFN_NAME = "eq_biodt_v03.ofn"

ALL_CBB_KOS = tuple(sorted({ko for alts in CBB_STEPS.values() for ko in alts} | {"K01602"}))

ENDPOINT_CLASSES = (
    "GenomicCapacityEndpoint",
    "MappedReadCompositionEndpoint",
    "DNARecoveryEndpoint",
    "ViableDNAEndpoint",
    "BiomassEndpoint",
    "PhysiologicalActivityEndpoint",
)


# ---------------------------------------------------------------------------
# Measured-quantity (endpoint) partition, units and observation maps.
# Name-based rak: IRIs (the twin style) never collide with the KG's numeric RAK_nnnnnnn ranges.
# ---------------------------------------------------------------------------
# name -> (label, definition, UO unit IRIs, UCUM, denominator, CHEBI measured entities, PATO quality type)
NEW_ENDPOINTS = {
    "ArealH2FluxEndpoint": (
        "areal H2 flux endpoint",
        "Net exchange of dihydrogen between a soil surface and the atmosphere per unit ground area and time, in nmol H2 m-2 s-1 (negative = uptake in the source sign convention).",
        ("UO_0000041", "UO_0000080", "UO_0000010"), "nmol.m-2.s-1", "ground surface area and time", ("CHEBI_18276",), "PATO_0000161"),
    "MassNormalizedBottleUptakeRateEndpoint": (
        "mass-normalized bottle uptake rate endpoint",
        "Gas uptake rate measured in a closed bottle incubation and normalised by soil mass, in nmol H2 h-1 g-1 soil.",
        ("UO_0000041", "UO_0000032", "UO_0000021"), "nmol.h-1.g-1", "soil mass and time", ("CHEBI_18276",), "PATO_0000161"),
    "FirstOrderRateConstantEndpoint": (
        "first-order rate constant endpoint",
        "Rate constant k of a first-order decline of headspace gas amount, in s-1 (convertible from h-1 by dividing by 3600).",
        ("UO_0000010",), "s-1", "time (headspace amount is not a denominator)", ("CHEBI_18276",), "PATO_0000161"),
    "GravimetricWaterContentEndpoint": (
        "gravimetric water content endpoint",
        "Mass of water per mass of dry soil (g g-1).",
        ("UO_0000191",), "g.g-1", "dry soil mass", ("CHEBI_15377",), "PATO_0001470"),
    "VolumetricWaterContentEndpoint": (
        "volumetric water content endpoint",
        "Volume of water per bulk soil volume (cm3 cm-3).",
        ("UO_0000191",), "cm3.cm-3", "bulk soil volume", ("CHEBI_15377",), "PATO_0001470"),
    "MatricWaterPotentialEndpoint": (
        "matric water potential endpoint",
        "Matric component of soil water potential (pressure units; reported in MPa, negative when dry).",
        ("UO_0000110",), "Pa", "none (energy per unit volume of water)", ("CHEBI_15377",), None),
    "TotalWaterPotentialEndpoint": (
        "total water potential endpoint",
        "Total soil water potential including at least matric and osmotic components (pressure units).",
        ("UO_0000110",), "Pa", "none (energy per unit volume of water)", ("CHEBI_15377",), None),
    "GenomeEquivalentsEndpoint": (
        "genome equivalents endpoint",
        "Estimated number of genome copies of a marker-defined guild from single-copy-gene-normalised read counts.",
        ("UO_0000189",), "{genome_equivalent}", "sample as normalised in the source table", (), None),
    "MappedReadShareEndpoint": (
        "mapped-read share endpoint",
        "Fraction of reads mapped to the reference MAG set that belong to a stated MAG class (denominator: mapped reads). A kind of mapped-read composition.",
        ("UO_0000191",), "1", "mapped reads", (), "PATO_0001470"),
}
# Partition members (pairwise disjoint). MappedReadShareEndpoint is a subclass of MappedReadCompositionEndpoint.
SUBCLASS_ENDPOINTS = {"MappedReadShareEndpoint": "MappedReadCompositionEndpoint"}
PARTITION_ENDPOINTS = ENDPOINT_CLASSES + tuple(k for k in NEW_ENDPOINTS if k not in SUBCLASS_ENDPOINTS)
ALL_ENDPOINTS = ENDPOINT_CLASSES + tuple(NEW_ENDPOINTS)

ENDPOINT_LABELS = {
    "GenomicCapacityEndpoint": "genomic capacity endpoint",
    "MappedReadCompositionEndpoint": "mapped-read composition endpoint",
    "DNARecoveryEndpoint": "DNA recovery endpoint",
    "ViableDNAEndpoint": "viable DNA endpoint",
    "BiomassEndpoint": "biomass endpoint",
    "PhysiologicalActivityEndpoint": "physiological activity endpoint",
}
ENDPOINT_DEFS = {
    "GenomicCapacityEndpoint": "Presence or share of genomes encoding a capacity (gene content); says nothing about expression or activity.",
    "MappedReadCompositionEndpoint": "Composition of the mapped read pool (shares of reads mapped to reference genomes).",
    "DNARecoveryEndpoint": "Total extractable DNA recovered per unit sample, including relic DNA.",
    "ViableDNAEndpoint": "DNA from cells with intact membranes (viable fraction).",
    "BiomassEndpoint": "Living biomass of a stated organism group per unit sample.",
    "PhysiologicalActivityEndpoint": "Measured physiological rate (for example gas uptake) of a community.",
}

# map name -> (label, definition, from, to, invertible, assumption/input classes)
OBSERVATION_MAPS = {
    "GravimetricToVolumetricMap": (
        "gravimetric to volumetric water content map",
        "Converts gravimetric to volumetric water content: theta_v = theta_g * rho_b / rho_w; requires a bulk-density assumption.",
        "GravimetricWaterContentEndpoint", "VolumetricWaterContentEndpoint", True, ("BulkDensityAssumption",)),
    "FirstOrderToBottleUptakeMap": (
        "first-order constant to bottle uptake map",
        "Converts a first-order rate constant to a mass-normalised bottle uptake rate: rate = k * n_gas / m_soil; requires headspace gas amount and soil mass.",
        "FirstOrderRateConstantEndpoint", "MassNormalizedBottleUptakeRateEndpoint", True,
        ("HeadspaceGasAmountInput", "SoilMassInput")),
    "VolumetricToMatricPotentialMap": (
        "volumetric water content to matric potential map",
        "Converts volumetric water content to matric potential through a water-retention curve; requires retention-curve parameters (including porosity) as assumptions.",
        "VolumetricWaterContentEndpoint", "MatricWaterPotentialEndpoint", True, ("WaterRetentionCurveAssumption",)),
    "MatricToTotalPotentialMap": (
        "matric to total water potential map",
        "Adds the osmotic (and other) components to the matric potential; requires a solute (here NaCl) osmotic-potential assumption.",
        "MatricWaterPotentialEndpoint", "TotalWaterPotentialEndpoint", True, ("OsmoticPotentialAssumption",)),
}
MAP_ASSUMPTIONS = {
    "BulkDensityAssumption": ("bulk density assumption", "Assumed (not measured) soil dry bulk density used to convert between mass- and volume-based quantities."),
    "HeadspaceGasAmountInput": ("headspace gas amount input", "Amount of gas (mol) in the incubation headspace needed to scale a rate constant to an uptake rate."),
    "SoilMassInput": ("soil mass input", "Mass of soil in the incubation bottle used as denominator."),
    "WaterRetentionCurveAssumption": ("water retention curve assumption", "Assumed soil water-retention relation (for example Campbell/Clapp-Hornberger parameters, air-entry potential, porosity)."),
    "OsmoticPotentialAssumption": ("osmotic potential assumption", "Assumed solute composition (NaCl) and concentration for the osmotic component of water potential."),
}
# Maps deliberately NOT declared: no direct observation map exists between these endpoints.
UNDECLARED_MAPS = {
    ("MassNormalizedBottleUptakeRateEndpoint", "ArealH2FluxEndpoint"): (
        "BottleRateToArealFlux",
        "Not declared: a bottle (mass-normalised) rate does not map to an areal flux without soil geometry (depth, bulk density), gas transport (diffusion, boundary layer) and moisture-dependent reaction-layer assumptions. This is a process-model prediction, not an observation map.",
        ("TransportGeometryAssumption",)),
}

COMPARTMENT_CLASSES = {
    "SurfaceSoilCompartment": ("surface soil sampling compartment", "Surface soil sampled at the transect ('Surface' compartment in sample_evidence.tsv)."),
    "DeepSoilCompartment": ("deep soil sampling compartment", "Deeper soil sampled at the transect ('Deep' compartment in sample_evidence.tsv)."),
    "PlantRootCompartment": ("plant root sampling compartment", "Plant root and adhering rhizosphere material sampled at the transect ('PlantRoot' compartment in sample_evidence.tsv)."),
}

# (rak term, relation, target local name (obo/sio), expected OLS label, namespace)
#   relation: skos:closeMatch | skos:relatedMatch | rak:measuredEntity | rak:unit | SubClassOf
_MAPPINGS_RAW = [
    ("CBBCompleteMAG", "skos:relatedMatch", "GO_0019253", "reductive pentose-phosphate cycle"),
    ("H2AutotrophicCandidateMAG", "skos:relatedMatch", "GO_0015977", "carbon fixation"),
    ("RuBisCOMAG", "skos:relatedMatch", "GO_0016984", "ribulose-bisphosphate carboxylase activity"),
    ("RuBisCOPRKMAG", "skos:relatedMatch", "GO_0016984", "ribulose-bisphosphate carboxylase activity"),
    ("RuBisCOPRKMAG", "skos:relatedMatch", "GO_0008974", "phosphoribulokinase activity"),
    ("CBBStep_carboxylation", "skos:relatedMatch", "GO_0016984", "ribulose-bisphosphate carboxylase activity"),
    ("CBBStep_phosphoribulokinase", "skos:relatedMatch", "GO_0008974", "phosphoribulokinase activity"),
    ("H2CarrierMAG", "skos:relatedMatch", "GO_0033748", "hydrogenase (acceptor) activity"),
    ("H2CarrierMAG", "skos:relatedMatch", "GO_1902421", "hydrogen metabolic process"),
    ("RAK_GUILD_TOTAL_H2ASE", "skos:relatedMatch", "GO_0033748", "hydrogenase (acceptor) activity"),
    ("RAK_FLUX_H2_OXIDATION_POWER", "skos:relatedMatch", "GO_1902421", "hydrogen metabolic process"),
    ("COCarrierMAG", "skos:relatedMatch", "GO_0008805", "aerobic carbon monoxide dehydrogenase activity"),
    ("RAK_GUILD_COXL_FORMI", "skos:relatedMatch", "GO_0008805", "aerobic carbon monoxide dehydrogenase activity"),
    ("RAK_FLUX_CO_OXIDATION_POWER", "skos:relatedMatch", "GO_0019410", "aerobic respiration, using carbon monoxide as electron donor"),
    ("RAK_FLUX_H2_OXIDATION_POWER", "rak:measuredEntity", "CHEBI_18276", "dihydrogen"),
    ("RAK_FLUX_H2_OXIDATION_POWER", "rak:measuredEntity", "CHEBI_15379", "dioxygen"),
    ("RAK_FLUX_H2_OXIDATION_POWER", "rak:measuredEntity", "CHEBI_15377", "water"),
    ("RAK_FLUX_CO_OXIDATION_POWER", "rak:measuredEntity", "CHEBI_17245", "carbon monoxide"),
    ("RAK_FLUX_CO_OXIDATION_POWER", "rak:measuredEntity", "CHEBI_15379", "dioxygen"),
    ("RAK_FLUX_CO_OXIDATION_POWER", "rak:measuredEntity", "CHEBI_16526", "carbon dioxide"),
    ("OsmoticPotentialAssumption", "skos:relatedMatch", "CHEBI_26710", "sodium chloride"),
    ("SurfaceSoilCompartment", "skos:relatedMatch", "ENVO_00001998", "soil"),
    ("SurfaceSoilCompartment", "skos:relatedMatch", "ENVO_01000017", "sand"),
    ("DeepSoilCompartment", "skos:relatedMatch", "ENVO_00001998", "soil"),
    ("PlantRootCompartment", "skos:relatedMatch", "ENVO_00005801", "rhizosphere"),
    ("BulkDensityAssumption", "skos:relatedMatch", "ENVO_00010005", "soil bulk density"),
    ("VolumetricWaterContentEndpoint", "skos:relatedMatch", "ENVO_09000026", "concentration of liquid water in soil"),
    ("MeasuredQuantity", "SubClassOf", "SIO_000070", "measurement value"),
    ("RAK_TAXON_TOTAL_PROKARYOTA", "(existing) child", "NCBITaxon_201174", "Actinomycetota"),
    ("RAK_TAXON_TOTAL_PROKARYOTA", "(existing) child", "NCBITaxon_2283796", "Thermoplasmatota"),
    ("RAK_TAXON_TOTAL_PROKARYOTA", "(existing) child", "NCBITaxon_3100166", "Halobacteriota"),
]
for _n, (_l, _d, _uo, _u, _den, _che, _pato) in NEW_ENDPOINTS.items():
    for _x in _uo:
        _MAPPINGS_RAW.append((_n, "rak:unit", _x, None))
    for _x in _che:
        _MAPPINGS_RAW.append((_n, "rak:measuredEntity", _x, None))
    if _pato:
        _MAPPINGS_RAW.append((_n, "skos:relatedMatch", _pato, None))
_EXPECTED_LABELS = {
    "UO_0000041": "nanomole", "UO_0000080": "square meter", "UO_0000010": "second", "UO_0000032": "hour",
    "UO_0000021": "gram", "UO_0000191": "fraction", "UO_0000110": "pascal", "UO_0000189": "count unit",
    "CHEBI_18276": "dihydrogen", "CHEBI_15377": "water", "PATO_0000161": "rate", "PATO_0001470": "ratio",
}
MAPPINGS = [(t, r, i, lab if lab else _EXPECTED_LABELS[i]) for t, r, i, lab in _MAPPINGS_RAW]
# Vocabulary terms used as annotation/object properties or superclasses: (IRI, expected label)
VOCABULARY = [
    (SKOS + "relatedMatch", "has related match"), (OBO + "IAO_0000115", "definition"),
    (SIO + "SIO_000070", "measurement value"),
]


def _mapping_iri(local: str) -> str:
    return (SIO if local.startswith("SIO_") else OBO) + local


def _q(text: str) -> str:
    return '"' + text.replace("\\", "\\\\").replace('"', '\\"') + '"'


def _step_defs() -> dict:
    out = {}
    for step, alts in CBB_STEPS.items():
        out[f"StepKO_{step}"] = (f"KEGG orthology groups for CBB step {step}",
                                 f"Union of the KEGG orthology groups ({', '.join(alts)}) accepted as evidence for the Calvin-Benson-Bassham step '{step}'.")
        out[f"CBBStep_{step}"] = (f"genome encoding CBB step {step}",
                                 f"A genome that encodes at least one gene assigned to a KEGG orthology group of the CBB step '{step}'.")
    return out


def class_annotations() -> dict:
    """rdfs:label and IAO:0000115 definition for every rak: class and property declared in the TBox."""
    d = {
        "hasKO": ("has KEGG ortholog", "Relates a genome (MAG) to a KEGG orthology group (KO) when the genome encodes at least one gene assigned to that KO by the annotation pipeline. Absence of the assertion means not detected, not absent."),
        "hasMarkerGuild": ("has marker guild", "Relates a genome to a marker-defined functional guild when it carries at least one marker gene of that guild (for example a [NiFe] hydrogenase group or form I CO dehydrogenase)."),
        "hasEndpoint": ("has endpoint", "Relates a quantity observation to the endpoint (measured quantity kind) it reports."),
        "mapsFrom": ("maps from", "Source endpoint of a declared observation map."),
        "mapsTo": ("maps to", "Target endpoint of a declared observation map."),
        "requiresAssumption": ("requires assumption", "Assumption or auxiliary input that a declared observation map needs before it can convert values."),
        "ReferenceMAG": ("reference metagenome-assembled genome", "A dereplicated representative MAG of the Rub' al Khali reference set."),
        "BacterialMAG": ("bacterial MAG", "A reference MAG classified in the domain Bacteria."),
        "ArchaealMAG": ("archaeal MAG", "A reference MAG classified in the domain Archaea."),
        "Quality90_5MAG": ("MAG with completeness >= 90% and contamination < 5%", "A reference MAG with CheckM completeness at least 90 and contamination below 5 percent."),
        "QuantityObservation": ("quantity observation", "An observation recording a value of a measured quantity with endpoint, unit and denominator."),
        "MeasuredQuantity": ("measured quantity", "Superclass of endpoint classes: kinds of numeric quantity that a measurement or derived estimate reports. Each is a measurement value in the SIO sense; subclasses in the partition are pairwise disjoint."),
        "ObservationMap": ("observation map", "A declared, assumption-bearing conversion between two disjoint endpoints. Endpoints without a declared map must not be compared or conflated."),
        "MapAssumption": ("map assumption", "An assumption or auxiliary input a declared observation map requires."),
        "RAK_TAXON_TOTAL_PROKARYOTA": ("total prokaryota (taxon partition)", "All prokaryotic reference MAGs; partitioned into Actinomycetota, Thermoplasmatota, Halobacteriota and other."),
        "RAK_TAXON_OTHER_PROKARYOTA": ("other prokaryota (taxon partition)", "Prokaryotic MAGs outside the three named phyla of the partition."),
        "RAK_DNA_TOTAL_SOIL_FRACTION": ("total soil DNA fraction", "Total extractable soil DNA; partitioned into viable and relic DNA."),
        "RAK_DNA_VIABLE_FRACTION": ("viable DNA fraction", "DNA from cells with intact membranes."),
        "RAK_DNA_RELIC_FRACTION": ("relic DNA fraction", "Extracellular or compromised-cell DNA."),
        "RAK_FLUX_TOTAL_TRACE_GAS_POWER": ("total trace-gas oxidation power", "Combined potential trace-gas (H2 plus CO) oxidation capacity."),
        "RAK_FLUX_H2_OXIDATION_POWER": ("H2 oxidation power", "Potential H2 oxidation capacity; component of total trace-gas power."),
        "RAK_FLUX_CO_OXIDATION_POWER": ("CO oxidation power", "Potential CO oxidation capacity; component of total trace-gas power."),
        "RAK_GUILD_TOTAL_TRACE_GAS": ("trace-gas oxidiser guild (total)", "All marker-defined trace-gas oxidisers: hydrogenase plus CO dehydrogenase guilds."),
        "RAK_GUILD_TOTAL_H2ASE": ("hydrogenase guild (total)", "Marker-defined [NiFe] hydrogenase guild: groups 1l, 1h and 2a."),
        "RAK_GUILD_HYDB_1L": ("[NiFe] hydrogenase group 1l guild", "Guild defined by the [NiFe] group 1l hydrogenase large-subunit (HydB) marker."),
        "RAK_GUILD_HYDB_1H": ("[NiFe] hydrogenase group 1h guild", "Guild defined by the [NiFe] group 1h hydrogenase large-subunit (HydB) marker."),
        "RAK_GUILD_HYDB_2A": ("[NiFe] hydrogenase group 2a guild", "Guild defined by the [NiFe] group 2a hydrogenase large-subunit (HydB) marker."),
        "RAK_GUILD_COXL_FORMI": ("form I CO dehydrogenase guild", "Guild defined by the form I carbon monoxide dehydrogenase large-subunit (CoxL) marker."),
        "H2CarrierMAG": ("H2 carrier MAG", "A MAG carrying at least one [NiFe] hydrogenase marker of group 1l, 1h or 2a (genomic capacity only)."),
        "COCarrierMAG": ("CO carrier MAG", "A MAG carrying a form I CO dehydrogenase (CoxL) marker (genomic capacity only)."),
        "RuBisCOMAG": ("RuBisCO-encoding MAG", "A MAG with at least one gene assigned to the RuBisCO large-subunit KO (carboxylation step); a single-marker screen, not a complete CBB cycle."),
        "ArchaealRuBisCOMAG": ("archaeal RuBisCO-encoding MAG", "A RuBisCO-encoding archaeal MAG."),
        "RuBisCOPRKMAG": ("RuBisCO plus phosphoribulokinase MAG", "A MAG encoding both RuBisCO and phosphoribulokinase (two-gene screen)."),
        "ThreeGeneCBBMAG": ("three-gene CBB screen MAG", "A MAG encoding RuBisCO, phosphoribulokinase and K01602 (RuBisCO small subunit)."),
        "CBBCompleteMAG": ("complete CBB (11-step) MAG", "A MAG encoding at least one KO for each of the 11 CBB steps; a genomic-capacity screen, not evidence of operation."),
        "H2AutotrophicCandidateMAG": ("H2-autotrophic candidate MAG", "A MAG that is an H2 carrier and encodes all 11 CBB steps; a candidate defined by gene content, not demonstrated autotrophy."),
        "HighQualityH2AutotrophicCandidateMAG": ("high-quality H2-autotrophic candidate MAG", "An H2-autotrophic candidate MAG with completeness >= 90% and contamination < 5%."),
        "DualH2COAutotrophicCandidateMAG": ("dual H2/CO autotrophic candidate MAG", "An H2-autotrophic candidate MAG that also carries a form I CO dehydrogenase marker."),
        "ConflatedObservationProbe": ("conflated observation probe", "Test class used to show that two disjoint endpoints cannot be conflated (unsatisfiable by design)."),
    }
    for n, (lab, definition, *_r) in NEW_ENDPOINTS.items():
        d[n] = (lab, definition)
    for n, definition in ENDPOINT_DEFS.items():
        d[n] = (ENDPOINT_LABELS[n], definition)
    for n, (lab, definition, *_r) in OBSERVATION_MAPS.items():
        d[n] = (lab, definition)
    for n, (lab, definition) in MAP_ASSUMPTIONS.items():
        d[n] = (lab, definition)
    for n, (lab, definition) in COMPARTMENT_CLASSES.items():
        d[n] = (lab, definition)
    d.update(_step_defs())
    return d


def _annotation_lines() -> list[str]:
    lines = [f"Declaration(AnnotationProperty({p}))" for p in
             ("obo:IAO_0000115", "skos:closeMatch", "skos:relatedMatch", "rak:unit", "rak:unitUCUM",
              "rak:denominator", "rak:measuredEntity", "rak:invertible")]
    for name, (lab, definition) in sorted(class_annotations().items()):
        lines.append(f"AnnotationAssertion(rdfs:label rak:{name} {_q(lab)}@en)")
        lines.append(f"AnnotationAssertion(obo:IAO_0000115 rak:{name} {_q(definition)}@en)")
    for n, (lab, definition, uo, ucum, den, che, pato) in NEW_ENDPOINTS.items():
        lines.append(f"AnnotationAssertion(rak:unitUCUM rak:{n} {_q(ucum)})")
        lines.append(f"AnnotationAssertion(rak:denominator rak:{n} {_q(den)})")
    for n, (lab, definition, a, b, inv, ass) in OBSERVATION_MAPS.items():
        lines.append(f"AnnotationAssertion(rak:invertible rak:{n} {_q(str(inv).lower())})")
    for term, rel, local, _lab in MAPPINGS:
        if rel in ("skos:closeMatch", "skos:relatedMatch", "rak:unit", "rak:measuredEntity"):
            lines.append(f"AnnotationAssertion({rel} rak:{term} <{_mapping_iri(local)}>)")
    return lines


def _measurement_lines() -> list[str]:
    lines = [
        "Declaration(Class(rak:MeasuredQuantity))",
        f"Declaration(Class(sio:SIO_000070))",
        "SubClassOf(rak:MeasuredQuantity sio:SIO_000070)",
        "Declaration(ObjectProperty(rak:mapsFrom))", "Declaration(ObjectProperty(rak:mapsTo))",
        "Declaration(ObjectProperty(rak:requiresAssumption))",
        "Declaration(Class(rak:ObservationMap))", "Declaration(Class(rak:MapAssumption))",
        "Declaration(Class(rak:SurfaceSoilCompartment))", "Declaration(Class(rak:DeepSoilCompartment))",
        "Declaration(Class(rak:PlantRootCompartment))",
        "DisjointClasses(rak:SurfaceSoilCompartment rak:DeepSoilCompartment rak:PlantRootCompartment)",
    ]
    for ep in ALL_ENDPOINTS:
        if ep in NEW_ENDPOINTS:
            lines.append(f"Declaration(Class(rak:{ep}))")
        lines.append(f"SubClassOf(rak:{ep} rak:MeasuredQuantity)")
    for sub, sup in SUBCLASS_ENDPOINTS.items():
        lines.append(f"SubClassOf(rak:{sub} rak:{sup})")
    lines.append("DisjointClasses(" + " ".join(f"rak:{ep}" for ep in PARTITION_ENDPOINTS) + ")")
    for a in MAP_ASSUMPTIONS:
        lines.append(f"Declaration(Class(rak:{a}))")
        lines.append(f"SubClassOf(rak:{a} rak:MapAssumption)")
    for m, (_l, _d, src, dst, _inv, assumptions) in OBSERVATION_MAPS.items():
        lines.append(f"Declaration(Class(rak:{m}))")
        lines.append(f"SubClassOf(rak:{m} rak:ObservationMap)")
        lines.append(f"SubClassOf(rak:{m} ObjectSomeValuesFrom(rak:mapsFrom rak:{src}))")
        lines.append(f"SubClassOf(rak:{m} ObjectSomeValuesFrom(rak:mapsTo rak:{dst}))")
        for a in assumptions:
            lines.append(f"SubClassOf(rak:{m} ObjectSomeValuesFrom(rak:requiresAssumption rak:{a}))")
    return lines


def find_km() -> Path:
    """Locate the Kobayashi-MaRust (km) OWL 2 reasoner executable."""
    env_bin = os.environ.get("KM_BIN")
    if env_bin and Path(env_bin).exists():
        return Path(env_bin)
    pinned = Path(__file__).resolve().parents[3] / ".tools" / "km"
    if pinned.exists():
        return pinned
    which_bin = shutil.which("km")
    if which_bin:
        return Path(which_bin)
    local_bin = Path.home() / ".local" / "bin" / "km"
    if local_bin.exists():
        return local_bin
    raise FileNotFoundError("Kobayashi-MaRust executable 'km' not found in KM_BIN, PATH or ~/.local/bin/km")


def _tbox_lines(annotations: bool = True) -> list[str]:
    lines = [
        f"Prefix(:=<{RAK}>)",
        f"Prefix(rak:=<{RAK}>)",
        f"Prefix(kegg:=<{KEGG}>)",
        f"Prefix(obo:=<{OBO}>)",
        f"Prefix(sio:=<{SIO}>)",
        f"Prefix(skos:=<{SKOS}>)",
        f"Ontology(<{RAK}eq_biodt_v03.owl>",
        "Declaration(ObjectProperty(rak:hasKO))",
        "Declaration(ObjectProperty(rak:hasMarkerGuild))",
        "Declaration(ObjectProperty(rak:hasEndpoint))",
        "Declaration(Class(rak:ReferenceMAG))",
        "Declaration(Class(rak:BacterialMAG))",
        "Declaration(Class(rak:ArchaealMAG))",
        "DisjointClasses(rak:BacterialMAG rak:ArchaealMAG)",
        "Declaration(Class(rak:Quality90_5MAG))",
        "Declaration(Class(rak:QuantityObservation))",
    ]
    for ep in ENDPOINT_CLASSES:
        lines.append(f"Declaration(Class(rak:{ep}))")
    lines.append("DisjointClasses(" + " ".join(f"rak:{ep}" for ep in ENDPOINT_CLASSES) + ")")

    # Verified taxonomic, DNA-partition, power and trace-gas guild hierarchy from rubalkhali_twin_constraints.owl
    lines.extend([
        "Declaration(Class(rak:RAK_TAXON_TOTAL_PROKARYOTA))",
        "Declaration(Class(obo:NCBITaxon_201174))",
        "Declaration(Class(obo:NCBITaxon_2283796))",
        "Declaration(Class(obo:NCBITaxon_3100166))",
        "Declaration(Class(rak:RAK_TAXON_OTHER_PROKARYOTA))",
        "SubClassOf(obo:NCBITaxon_201174 rak:RAK_TAXON_TOTAL_PROKARYOTA)",
        "SubClassOf(obo:NCBITaxon_2283796 rak:RAK_TAXON_TOTAL_PROKARYOTA)",
        "SubClassOf(obo:NCBITaxon_3100166 rak:RAK_TAXON_TOTAL_PROKARYOTA)",
        "SubClassOf(rak:RAK_TAXON_OTHER_PROKARYOTA rak:RAK_TAXON_TOTAL_PROKARYOTA)",
        "DisjointUnion(rak:RAK_TAXON_TOTAL_PROKARYOTA obo:NCBITaxon_201174 obo:NCBITaxon_2283796 obo:NCBITaxon_3100166 rak:RAK_TAXON_OTHER_PROKARYOTA)",
        "Declaration(Class(rak:RAK_DNA_TOTAL_SOIL_FRACTION))",
        "Declaration(Class(rak:RAK_DNA_VIABLE_FRACTION))",
        "Declaration(Class(rak:RAK_DNA_RELIC_FRACTION))",
        "SubClassOf(rak:RAK_DNA_VIABLE_FRACTION rak:RAK_DNA_TOTAL_SOIL_FRACTION)",
        "SubClassOf(rak:RAK_DNA_RELIC_FRACTION rak:RAK_DNA_TOTAL_SOIL_FRACTION)",
        "DisjointUnion(rak:RAK_DNA_TOTAL_SOIL_FRACTION rak:RAK_DNA_VIABLE_FRACTION rak:RAK_DNA_RELIC_FRACTION)",
        "Declaration(Class(rak:RAK_FLUX_TOTAL_TRACE_GAS_POWER))",
        "Declaration(Class(rak:RAK_FLUX_H2_OXIDATION_POWER))",
        "Declaration(Class(rak:RAK_FLUX_CO_OXIDATION_POWER))",
        "SubClassOf(rak:RAK_FLUX_H2_OXIDATION_POWER rak:RAK_FLUX_TOTAL_TRACE_GAS_POWER)",
        "SubClassOf(rak:RAK_FLUX_CO_OXIDATION_POWER rak:RAK_FLUX_TOTAL_TRACE_GAS_POWER)",
        "DisjointUnion(rak:RAK_FLUX_TOTAL_TRACE_GAS_POWER rak:RAK_FLUX_H2_OXIDATION_POWER rak:RAK_FLUX_CO_OXIDATION_POWER)",
        "Declaration(Class(rak:RAK_GUILD_TOTAL_TRACE_GAS))",
        "Declaration(Class(rak:RAK_GUILD_TOTAL_H2ASE))",
        "Declaration(Class(rak:RAK_GUILD_HYDB_1L))",
        "Declaration(Class(rak:RAK_GUILD_HYDB_1H))",
        "Declaration(Class(rak:RAK_GUILD_HYDB_2A))",
        "Declaration(Class(rak:RAK_GUILD_COXL_FORMI))",
        "SubClassOf(rak:RAK_GUILD_HYDB_1L rak:RAK_GUILD_TOTAL_H2ASE)",
        "SubClassOf(rak:RAK_GUILD_HYDB_1H rak:RAK_GUILD_TOTAL_H2ASE)",
        "SubClassOf(rak:RAK_GUILD_HYDB_2A rak:RAK_GUILD_TOTAL_H2ASE)",
        "SubClassOf(rak:RAK_GUILD_TOTAL_H2ASE rak:RAK_GUILD_TOTAL_TRACE_GAS)",
        "SubClassOf(rak:RAK_GUILD_COXL_FORMI rak:RAK_GUILD_TOTAL_TRACE_GAS)",
        "DisjointUnion(rak:RAK_GUILD_TOTAL_H2ASE rak:RAK_GUILD_HYDB_1L rak:RAK_GUILD_HYDB_1H rak:RAK_GUILD_HYDB_2A)",
        "DisjointUnion(rak:RAK_GUILD_TOTAL_TRACE_GAS rak:RAK_GUILD_TOTAL_H2ASE rak:RAK_GUILD_COXL_FORMI)",
        "Declaration(Class(rak:H2CarrierMAG))",
        "EquivalentClasses(rak:H2CarrierMAG ObjectSomeValuesFrom(rak:hasMarkerGuild rak:RAK_GUILD_TOTAL_H2ASE))",
        "Declaration(Class(rak:COCarrierMAG))",
        "EquivalentClasses(rak:COCarrierMAG ObjectSomeValuesFrom(rak:hasMarkerGuild rak:RAK_GUILD_COXL_FORMI))",
    ])
    for ko in ALL_CBB_KOS:
        lines.append(f"Declaration(Class(kegg:{ko}))")
    for step, alts in CBB_STEPS.items():
        lines.append(f"Declaration(Class(rak:StepKO_{step}))")
        lines.append(f"Declaration(Class(rak:CBBStep_{step}))")
        for ko in alts:
            lines.append(f"SubClassOf(kegg:{ko} rak:StepKO_{step})")
        lines.append(f"EquivalentClasses(rak:CBBStep_{step} ObjectSomeValuesFrom(rak:hasKO rak:StepKO_{step}))")

    step_classes = " ".join(f"rak:CBBStep_{s}" for s in CBB_STEPS)
    lines.extend([
        "Declaration(Class(rak:RuBisCOMAG))",
        "EquivalentClasses(rak:RuBisCOMAG rak:CBBStep_carboxylation)",
        "Declaration(Class(rak:ArchaealRuBisCOMAG))",
        "EquivalentClasses(rak:ArchaealRuBisCOMAG ObjectIntersectionOf(rak:ArchaealMAG rak:RuBisCOMAG))",
        "Declaration(Class(rak:RuBisCOPRKMAG))",
        "EquivalentClasses(rak:RuBisCOPRKMAG ObjectIntersectionOf(rak:CBBStep_carboxylation rak:CBBStep_phosphoribulokinase))",
        "Declaration(Class(rak:ThreeGeneCBBMAG))",
        "EquivalentClasses(rak:ThreeGeneCBBMAG ObjectIntersectionOf(rak:CBBStep_carboxylation rak:CBBStep_phosphoribulokinase ObjectSomeValuesFrom(rak:hasKO kegg:K01602)))",
        "Declaration(Class(rak:CBBCompleteMAG))",
        f"EquivalentClasses(rak:CBBCompleteMAG ObjectIntersectionOf({step_classes}))",
        "Declaration(Class(rak:H2AutotrophicCandidateMAG))",
        "EquivalentClasses(rak:H2AutotrophicCandidateMAG ObjectIntersectionOf(rak:H2CarrierMAG rak:CBBCompleteMAG))",
        "Declaration(Class(rak:HighQualityH2AutotrophicCandidateMAG))",
        "EquivalentClasses(rak:HighQualityH2AutotrophicCandidateMAG ObjectIntersectionOf(rak:H2AutotrophicCandidateMAG rak:Quality90_5MAG))",
        "Declaration(Class(rak:DualH2COAutotrophicCandidateMAG))",
        "EquivalentClasses(rak:DualH2COAutotrophicCandidateMAG ObjectIntersectionOf(rak:H2AutotrophicCandidateMAG ObjectSomeValuesFrom(rak:hasMarkerGuild rak:RAK_GUILD_COXL_FORMI)))",
    ])
    lines.extend(_measurement_lines())
    if annotations:
        lines.extend(_annotation_lines())
    return lines


def _mag_axiom_lines(row) -> list[str]:
    m_iri = f"rak:MAG_{row.rep}"
    lines = [
        f"Declaration(Class({m_iri}))",
        f"SubClassOf({m_iri} rak:ReferenceMAG)",
    ]
    domain = getattr(row, "domain", "")
    if domain == "Bacteria":
        lines.append(f"SubClassOf({m_iri} rak:BacterialMAG)")
    elif domain == "Archaea":
        lines.append(f"SubClassOf({m_iri} rak:ArchaealMAG)")
    phylum = getattr(row, "phylum", "")
    if phylum == "Actinomycetota":
        lines.append(f"SubClassOf({m_iri} obo:NCBITaxon_201174)")
    elif phylum == "Thermoplasmatota":
        lines.append(f"SubClassOf({m_iri} obo:NCBITaxon_2283796)")
    elif phylum == "Halobacteriota":
        lines.append(f"SubClassOf({m_iri} obo:NCBITaxon_3100166)")
    elif phylum:
        lines.append(f"SubClassOf({m_iri} rak:RAK_TAXON_OTHER_PROKARYOTA)")
    cbb_kos_str = getattr(row, "cbb_kos", "")
    if isinstance(cbb_kos_str, str) and cbb_kos_str:
        for ko in cbb_kos_str.split(";"):
            lines.append(f"SubClassOf({m_iri} ObjectSomeValuesFrom(rak:hasKO kegg:{ko}))")
    if row.h2_1l:
        lines.append(f"SubClassOf({m_iri} ObjectSomeValuesFrom(rak:hasMarkerGuild rak:RAK_GUILD_HYDB_1L))")
    if row.h2_1h:
        lines.append(f"SubClassOf({m_iri} ObjectSomeValuesFrom(rak:hasMarkerGuild rak:RAK_GUILD_HYDB_1H))")
    if row.h2_2a:
        lines.append(f"SubClassOf({m_iri} ObjectSomeValuesFrom(rak:hasMarkerGuild rak:RAK_GUILD_HYDB_2A))")
    if row.co_form1:
        lines.append(f"SubClassOf({m_iri} ObjectSomeValuesFrom(rak:hasMarkerGuild rak:RAK_GUILD_COXL_FORMI))")
    if row.quality_90_5:
        lines.append(f"SubClassOf({m_iri} rak:Quality90_5MAG)")
    return lines


def build_owl_ontology(mag: pd.DataFrame, output_path: Path) -> Path:
    """Write the self-contained OWL 2 Functional Syntax ontology for all reference MAGs."""
    lines = _tbox_lines()
    for row in mag.itertuples():
        lines.extend(_mag_axiom_lines(row))
    lines.append(")")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text("\n".join(lines) + "\n")
    return output_path


def classify_ontology(ofn_path: Path) -> dict:
    """Classify an OWL 2 ontology with Kobayashi-MaRust (km)."""
    km = find_km()
    t0 = time.perf_counter()
    res = subprocess.run([str(km), "classify", str(ofn_path)], capture_output=True, text=True, check=True)
    wall_s = time.perf_counter() - t0
    payload = json.loads(res.stdout)
    if not payload.get("consistent", False) or payload.get("unsatisfiable") or payload.get("dropped", 0) != 0:
        raise RuntimeError(f"OWL ontology classification failed consistency/completeness gate: {payload}")
    prof_res = subprocess.run([str(km), "profile", str(ofn_path)], capture_output=True, text=True, check=True)
    profile = json.loads(prof_res.stdout)
    return {"classification": payload, "profile": profile, "wall_seconds": wall_s, "reasoner": str(km)}


def explain_mag_entailment(mag: pd.DataFrame, mag_id: str, target_class: str, tmp_dir: Path) -> dict:
    """Extract a verified, subset-minimal OWL source-axiom justification with km explain."""
    km = find_km()
    row = mag.loc[mag.rep.eq(mag_id)]
    if row.empty:
        raise ValueError(f"Unknown reference MAG: {mag_id}")
    lines = _tbox_lines(annotations=False)
    lines.extend(_mag_axiom_lines(row.iloc[0]))
    lines.append(")")
    focal_path = tmp_dir / f"explain_{mag_id}_{target_class}.ofn"
    focal_path.write_text("\n".join(lines) + "\n")
    res = subprocess.run(
        [
            str(km), "explain", "--max-checks", "2048",
            str(focal_path), "subclass",
            f"{RAK}MAG_{mag_id}", f"{RAK}{target_class}",
        ],
        capture_output=True, text=True, check=True,
    )
    return json.loads(res.stdout)


def explain_endpoint_conflict(endpoint_a: str, endpoint_b: str, tmp_dir: Path) -> dict:
    """Verify via km explain that conflating two distinct observation endpoints is unsatisfiable."""
    km = find_km()
    lines = _tbox_lines(annotations=False)
    lines.extend([
        "Declaration(Class(rak:ConflatedObservationProbe))",
        f"SubClassOf(rak:ConflatedObservationProbe rak:{endpoint_a})",
        f"SubClassOf(rak:ConflatedObservationProbe rak:{endpoint_b})",
        ")",
    ])
    probe_path = tmp_dir / f"conflict_{endpoint_a}_{endpoint_b}.ofn"
    probe_path.write_text("\n".join(lines) + "\n")
    res = subprocess.run(
        [
            str(km), "explain", "--max-checks", "1024",
            str(probe_path), "unsatisfiable",
            f"{RAK}ConflatedObservationProbe",
        ],
        capture_output=True, text=True, check=True,
    )
    return json.loads(res.stdout)



def _km_json(args: list[str]) -> dict:
    res = subprocess.run([str(find_km()), *args], capture_output=True, text=True, check=True)
    return json.loads(res.stdout)


def build_observation_table(tmp_dir: Path) -> dict:
    """Classify every endpoint pair with km and join the result with the declared observation maps.

    Disjointness is established by one probe class per unordered pair (subclass of both endpoints):
    the pair is disjoint exactly when km reports the probe unsatisfiable. Subsumption comes from the
    classification of the same ontology. Declared maps come from the TBox (OBSERVATION_MAPS).
    """
    names = list(ALL_ENDPOINTS)
    pairs = [(a, b) for i, a in enumerate(names) for b in names[i + 1:]]
    lines = _tbox_lines(annotations=False)
    for k, (a, b) in enumerate(pairs):
        lines += [f"Declaration(Class(rak:ObsProbe_{k}))", f"SubClassOf(rak:ObsProbe_{k} rak:{a})",
                  f"SubClassOf(rak:ObsProbe_{k} rak:{b})"]
    lines.append(")")
    path = tmp_dir / "observation_probes.ofn"
    path.write_text("\n".join(lines) + "\n")
    out = _km_json(["classify", str(path)])
    if not out.get("consistent", False) or out.get("dropped", 0) != 0:
        raise RuntimeError(f"Observation probe ontology failed gate: {out}")
    unsat = set(out["unsatisfiable"])
    sub = {(x.replace("rak:", "").replace(RAK, ""), y.replace("rak:", "").replace(RAK, "")) for x, y in out["subsumptions"]}
    unsat_names = {u.replace("rak:", "").replace(RAK, "") for u in unsat}
    declared = {}
    for m, (lab, definition, src, dst, inv, ass) in OBSERVATION_MAPS.items():
        declared[(src, dst)] = (m, "forward", ass, inv)
        if inv:
            declared[(dst, src)] = (m, "inverse", ass, inv)

    def justify(k, a, b):
        one = tmp_dir / f"observation_probe_{k}.ofn"  # km explain caps source axioms: one probe per file
        one.write_text("\n".join(_tbox_lines(annotations=False) + [
            f"Declaration(Class(rak:ObsProbe_{k}))", f"SubClassOf(rak:ObsProbe_{k} rak:{a})",
            f"SubClassOf(rak:ObsProbe_{k} rak:{b})", ")"]) + "\n")
        res = _km_json(["explain", "--max-checks", "1024", str(one), "unsatisfiable", f"{RAK}ObsProbe_{k}"])
        if res.get("status") != "entailed" or not res.get("justifications"):
            raise RuntimeError(f"No disjointness justification for {a}/{b}")
        j = res["justifications"][0]
        if not (j["verified"] and j["subsetMinimal"]):
            raise RuntimeError(f"Unverified disjointness justification for {a}/{b}")
        return [ax["functionalSyntax"].replace(f"rak:ObsProbe_{k}", "rak:Probe") for ax in j["axioms"]]

    disjoint = [(k, a, b) for k, (a, b) in enumerate(pairs) if f"ObsProbe_{k}" in unsat_names]
    with ThreadPoolExecutor(max_workers=8) as pool:
        justs = dict(zip([k for k, _, _ in disjoint], pool.map(lambda t: justify(*t), disjoint)))
    table = {}
    for k, (a, b) in enumerate(pairs):
        entry = {"a": a, "b": b}
        if (a, b) in sub or (b, a) in sub:
            entry["relation"] = "subsumed"
            entry["subsumption"] = f"{a} SubClassOf {b}" if (a, b) in sub else f"{b} SubClassOf {a}"
        elif k in justs:
            entry["relation"] = "disjoint"
            entry["disjointness_justification"] = justs[k]
        else:
            entry["relation"] = "overlapping"
        entry["declared_maps"] = [{"map": declared[d][0], "assumptions": list(declared[d][2])}
                                  for d in ((a, b), (b, a)) if d in declared]
        table[f"{a}|{b}"] = entry
    return {"pairs": table, "endpoints": names, "probe_count": len(pairs), "disjoint_pairs": len(disjoint)}


def export_observation_maps(report_obs: dict, data: Path, ontology_sha256: str) -> None:
    undeclared = [{"from": a, "to": b, "name": n, "note": note, "would_require": list(req)}
                  for (a, b), (n, note, req) in UNDECLARED_MAPS.items()]
    endpoints = {n: {"label": class_annotations()[n][0], "definition": class_annotations()[n][1],
                     "unit_ucum": NEW_ENDPOINTS[n][3] if n in NEW_ENDPOINTS else None,
                     "unit_uo": [_mapping_iri(u) for u in NEW_ENDPOINTS[n][2]] if n in NEW_ENDPOINTS else [],
                     "denominator": NEW_ENDPOINTS[n][4] if n in NEW_ENDPOINTS else None}
                 for n in ALL_ENDPOINTS}
    maps = {m: {"label": lab, "definition": d, "from": src, "to": dst, "invertible": inv,
                "assumptions": list(ass)} for m, (lab, d, src, dst, inv, ass) in OBSERVATION_MAPS.items()}
    art = {"schema_version": 1, "rules_sha256": rules_digest(), "ontology_sha256": ontology_sha256,
           "endpoints": endpoints, "declared_maps": maps, "undeclared_maps": undeclared,
           "probe_count": report_obs["probe_count"], "disjoint_pairs": report_obs["disjoint_pairs"],
           "pairs": report_obs["pairs"]}
    (data / "observation_maps.json").write_text(json.dumps(art, indent=1, sort_keys=True) + "\n")


def _load_observation_table(data: Path) -> dict:
    manifest = json.loads((data / "manifest.json").read_text())
    path = data / "observation_maps.json"
    if sha256(path) != manifest["artifacts"].get("observation_maps.json"):
        raise ValueError("OWL artifact checksum mismatch: observation_maps.json")
    table = json.loads(path.read_text())
    if table["rules_sha256"] != rules_digest():
        raise ValueError("Stale observation map table: rebuild the ontology")
    return table


def default_data_dir() -> Path:
    """Bundled data in a wheel, or data/v03 in a repository checkout."""
    package = Path(__file__).resolve().parent
    return package/"data" if (package/"data").exists() else package.parents[2]/"data"/"v03"


def observation_compatible(endpoint_a: str, endpoint_b: str, data: Path | None = None,
                           table: dict | None = None) -> dict:
    """Decide whether values of two endpoints may be compared, using the reasoner-derived table.

    status: 'compatible' (same endpoint or one subsumes the other), 'compatible_via_map' (disjoint
    endpoints with a declared observation map; the map and its required assumptions are returned) or
    'incompatible' (disjoint, no declared map; justification returned).
    """
    if table is None:
        table = _load_observation_table(data if data is not None else default_data_dir())
    known = table["endpoints"]
    for e in (endpoint_a, endpoint_b):
        if e not in known:
            raise ValueError(f"Unknown endpoint class: {e}")
    if endpoint_a == endpoint_b:
        return {"status": "compatible", "relation": "same", "a": endpoint_a, "b": endpoint_b, "map": None,
                "assumptions": [], "justification": "Same endpoint class."}
    key = "|".join(sorted((endpoint_a, endpoint_b), key=ALL_ENDPOINTS.index))
    entry = table["pairs"][key]
    base = {"a": endpoint_a, "b": endpoint_b, "relation": entry["relation"], "map": None, "assumptions": []}
    if entry["relation"] == "subsumed":
        return {**base, "status": "compatible", "justification": f"Reasoner subsumption: {entry['subsumption']}."}
    if entry["relation"] == "overlapping":
        return {**base, "status": "incompatible",
                "justification": "Endpoints are neither identical, subsumed nor declared disjoint; no compatibility is entailed."}
    for dm in entry["declared_maps"]:
        spec = table["declared_maps"][dm["map"]]
        forward = (spec["from"], spec["to"]) == (endpoint_a, endpoint_b)
        if forward or spec["invertible"]:
            return {**base, "status": "compatible_via_map", "map": dm["map"],
                    "direction": "forward" if forward else "inverse", "assumptions": list(dm["assumptions"]),
                    "justification": f"Disjoint endpoints; declared observation map {dm['map']} converts them given the listed assumptions/inputs."}
    undeclared = [u for u in table["undeclared_maps"] if {u["from"], u["to"]} == {endpoint_a, endpoint_b}]
    why = undeclared[0]["note"] if undeclared else "No observation map is declared between these disjoint endpoints."
    return {**base, "status": "incompatible", "would_require": undeclared[0]["would_require"] if undeclared else [],
            "justification": f"Disjoint endpoints ({'; '.join(entry['disjointness_justification'][:1])}). {why}"}


def extract_justifications(mag: pd.DataFrame, entailed: dict, work: Path, workers: int = 8) -> dict:
    """Subset-minimal, verified km justification for every entailed member of the two target classes."""
    t0 = time.perf_counter()
    targets = (("H2AutotrophicCandidateMAG", "H2AutotrophicCandidateMAG"), ("COCarrierMAG", "COCarrierMAG"))
    jobs = [(cls, rep) for cls, key in targets for rep in sorted(entailed[key])]

    def one(job):
        cls, rep = job
        res = explain_mag_entailment(mag, rep, cls, work)
        just = res["justifications"][0] if res.get("justifications") else None
        if res.get("status") != "entailed" or just is None or not (just["verified"] and just["subsetMinimal"]):
            raise RuntimeError(f"No verified subset-minimal justification for {rep} / {cls}")
        axioms = [a["functionalSyntax"].replace(f"rak:MAG_{rep}", "?mag") for a in just["axioms"]]
        return cls, rep, axioms, res.get("enumerationComplete"), len(res["justifications"])

    with ThreadPoolExecutor(max_workers=workers) as pool:
        results = list(pool.map(one, jobs))
    tbox_index: dict[str, int] = {}
    out: dict = {c: {} for c, _ in targets}
    for cls, rep, axioms, complete, n in results:
        ids = [tbox_index.setdefault(a, len(tbox_index)) for a in axioms if "?mag" not in a]
        out[cls][rep] = {"tbox_axiom_ids": sorted(ids), "mag_axioms": [a for a in axioms if "?mag" in a],
                         "justifications_returned": n, "enumeration_complete": complete}
    return {"classes": out, "tbox_axioms": [a for a, _ in sorted(tbox_index.items(), key=lambda kv: kv[1])],
            "counts": {c: len(v) for c, v in out.items()}, "wall_seconds": round(time.perf_counter() - t0, 1),
            "workers": workers, "all_verified_subset_minimal": True}


def evaluate_ontology_reasoning(mag: pd.DataFrame, panel: pd.DataFrame, ofn_path: Path) -> dict:
    """Run Kobayashi-MaRust OWL classification, justification checks and semantic ablations."""
    build_owl_ontology(mag, ofn_path)
    run = classify_ontology(ofn_path)
    subsumptions = run["classification"]["subsumptions"]
    profile = run["profile"]

    def entailed_mags(cls_name: str) -> set[str]:
        target_curie = f"rak:{cls_name}"
        target_iri = f"{RAK}{cls_name}"
        prefix_curie = "rak:MAG_"
        prefix_iri = f"{RAK}MAG_"
        out = set()
        for sub, sup in subsumptions:
            if sup in (target_curie, target_iri):
                if sub.startswith(prefix_curie):
                    out.add(sub[len(prefix_curie):])
                elif sub.startswith(prefix_iri):
                    out.add(sub[len(prefix_iri):])
        return out

    entailed = {
        "RuBisCOMAG": entailed_mags("RuBisCOMAG"),
        "ArchaealRuBisCOMAG": entailed_mags("ArchaealRuBisCOMAG"),
        "RuBisCOPRKMAG": entailed_mags("RuBisCOPRKMAG"),
        "ThreeGeneCBBMAG": entailed_mags("ThreeGeneCBBMAG"),
        "CBBCompleteMAG": entailed_mags("CBBCompleteMAG"),
        "H2CarrierMAG": entailed_mags("H2CarrierMAG"),
        "COCarrierMAG": entailed_mags("COCarrierMAG"),
        "H2AutotrophicCandidateMAG": entailed_mags("H2AutotrophicCandidateMAG"),
        "HighQualityH2AutotrophicCandidateMAG": entailed_mags("HighQualityH2AutotrophicCandidateMAG"),
        "DualH2COAutotrophicCandidateMAG": entailed_mags("DualH2COAutotrophicCandidateMAG"),
    }

    expected = {
        "RuBisCOMAG": set(mag.loc[mag.rbcl_detected, "rep"]),
        "ArchaealRuBisCOMAG": set(mag.loc[mag.rbcl_detected & mag.domain.eq("Archaea"), "rep"]),
        "RuBisCOPRKMAG": set(mag.loc[mag.rbcl_prk_detected, "rep"]),
        "ThreeGeneCBBMAG": set(mag.loc[mag.three_gene_screen, "rep"]),
        "CBBCompleteMAG": set(mag.loc[mag.cbb_all_steps_detected, "rep"]),
        "H2CarrierMAG": set(mag.loc[mag.h2_carrier_union, "rep"]),
        "COCarrierMAG": set(mag.loc[mag.co_form1, "rep"]),
        "H2AutotrophicCandidateMAG": set(mag.loc[mag.h2_cbb_all_steps, "rep"]),
        "HighQualityH2AutotrophicCandidateMAG": set(mag.loc[mag.h2_cbb_quality_90_5, "rep"]),
        "DualH2COAutotrophicCandidateMAG": set(mag.loc[mag.h2_cbb_all_steps & mag.co_form1, "rep"]),
    }
    for k in expected:
        if entailed[k] != expected[k]:
            raise RuntimeError(f"OWL reasoner mismatch on {k}: {len(entailed[k])} vs {len(expected[k])}")

    # Extract subset-minimal explanation on the first candidate MAG in lexical order
    first_candidate = sorted(entailed["H2AutotrophicCandidateMAG"])[0]
    with tempfile.TemporaryDirectory(prefix="reason-", dir=ofn_path.parent) as work:
        mag_exp = explain_mag_entailment(mag, first_candidate, "H2AutotrophicCandidateMAG", Path(work))
        conflict_exp = explain_endpoint_conflict("GenomicCapacityEndpoint", "MappedReadCompositionEndpoint", Path(work))

    if mag_exp.get("status") != "entailed" or not mag_exp.get("justifications"):
        raise RuntimeError("Failed to extract subset-minimal OWL justification for candidate MAG")
    if conflict_exp.get("status") != "entailed" or not conflict_exp.get("justifications"):
        raise RuntimeError("Failed to extract unsatisfiability justification for endpoint conflict")

    just = mag_exp["justifications"][0]
    conflict_just = conflict_exp["justifications"][0]
    with tempfile.TemporaryDirectory(prefix="reason-", dir=ofn_path.parent) as work:
        observation = build_observation_table(Path(work))
        # Exhaustive justification extraction (~3 min with 8 workers) only for the shipped data/ontology layout.
        justifications = extract_justifications(mag, entailed, Path(work)) if ofn_path.parent.name == "ontology" else None

    # Semantic ablations & genome-linkage vs marginal interval analysis on the 551 samples
    rbcl_h2_mags = entailed["RuBisCOMAG"] & entailed["H2CarrierMAG"]
    prk_h2_mags = entailed["RuBisCOPRKMAG"] & entailed["H2CarrierMAG"]
    multi_h2_mags = int((mag[["h2_1l", "h2_1h", "h2_2a"]].sum(axis=1) > 1).sum())

    valid = panel.loc[panel.mapped_read_pct.notna()].copy()
    h = valid["h2_carrier_union_share_mapped"].to_numpy(float)
    c = valid["cbb_all_steps_detected_share_mapped"].to_numpy(float)
    linked = valid["h2_cbb_all_steps_share_mapped"].to_numpy(float)
    lower = np.maximum(0.0, h + c - 1.0)
    upper = np.minimum(h, c)
    width = upper - lower

    cbb_df = mag.loc[mag.rep.isin(entailed["CBBCompleteMAG"])]
    cand_df = mag.loc[mag.rep.isin(entailed["H2AutotrophicCandidateMAG"])]
    arch_rbcl_df = mag.loc[mag.rep.isin(entailed["ArchaealRuBisCOMAG"])]

    compartment_summary = {}
    for comp, grp in panel.groupby("compartment"):
        ge_s = grp["total_h2_marker_ge"].dropna()
        rbcl_s = grp["rbcl_detected_share_mapped"].dropna() * 100.0
        cbb_s = grp["cbb_all_steps_detected_share_mapped"].dropna() * 100.0
        h2_s = grp["h2_carrier_union_share_mapped"].dropna() * 100.0
        h2cbb_s = grp["h2_cbb_all_steps_share_mapped"].dropna() * 100.0
        compartment_summary[str(comp)] = {
            "n_samples": int(len(grp)),
            "total_h2_marker_ge_median": round(float(ge_s.median()), 4),
            "total_h2_marker_ge_iqr": [round(float(ge_s.quantile(0.25)), 4), round(float(ge_s.quantile(0.75)), 4)],
            "h2_carrier_mapped_pct_median": round(float(h2_s.median()), 2),
            "rbcl_single_marker_mapped_pct_median": round(float(rbcl_s.median()), 2),
            "cbb_11_step_mapped_pct_median": round(float(cbb_s.median()), 2),
            "h2_cbb_candidate_mapped_pct_median": round(float(h2cbb_s.median()), 3),
            "h2_cbb_candidate_mapped_pct_max": round(float(h2cbb_s.max()), 2),
            "samples_exceeding_1pct_h2_cbb": int((h2cbb_s > 1.0).sum()),
        }

    return {
        "reasoner": "Kobayashi-MaRust (km)",
        "reasoner_version": mag_exp.get("reasonerVersion"),
        "reasoner_binary_sha256": sha256(find_km()),
        "ontology_sha256": sha256(ofn_path),
        "selected_route": profile.get("selected_route"),
        "expressivity": profile.get("profile", {}).get("expressivity", {}).get("code"),
        "ontology_file": str(ofn_path.name),
        "logical_axioms": profile.get("profile", {}).get("source", {}).get("logical_axioms"),
        "declared_classes": profile.get("profile", {}).get("source", {}).get("declared_classes"),
        "consistent": True,
        "unsatisfiable_classes": 0,
        "dropped_clauses": 0,
        "entailed_class_counts": {k: len(v) for k, v in entailed.items()},
        "exact_agreement_all_5438_mags": True,
        "taxonomic_architecture": {
            "archaeal_rubisco_detected_mags": len(entailed["ArchaealRuBisCOMAG"]),
            "archaeal_rubisco_by_phylum": {k: int(v) for k, v in arch_rbcl_df["phylum"].value_counts().items()},
            "cbb_complete_by_phylum": {k: int(v) for k, v in cbb_df["phylum"].value_counts().items()},
            "h2_autotrophic_candidates_by_phylum": {k: int(v) for k, v in cand_df["phylum"].value_counts().items()},
            "h2_autotrophic_candidates_by_order": {k: int(v) for k, v in cand_df["order"].value_counts().items()},
            "h2_autotrophic_candidates_by_genus": {k: int(v) for k, v in cand_df["genus"].replace("", "unclassified").value_counts().items()},
            "h2_autotrophic_hydrogenase_subtypes": {
                "group_1l": int(cand_df["h2_1l"].sum()),
                "group_1h": int(cand_df["h2_1h"].sum()),
                "group_2a": int(cand_df["h2_2a"].sum()),
                "dual_h2_co_form1": len(entailed["DualH2COAutotrophicCandidateMAG"]),
            },
        },
        "compartment_niche_partitioning": compartment_summary,
        "example_candidate_justification": {
            "mag_id": first_candidate,
            "superClass": f"{RAK}H2AutotrophicCandidateMAG",
            "verified": bool(just["verified"]),
            "subset_minimal": bool(just["subsetMinimal"]),
            "axiom_count": len(just["axioms"]),
            "classification_checks": mag_exp["classificationChecks"],
            "axioms": [a["functionalSyntax"] for a in just["axioms"]],
        },
        "endpoint_disjointness_justification": {
            "probe_class": f"{RAK}ConflatedObservationProbe",
            "endpoints": ["GenomicCapacityEndpoint", "MappedReadCompositionEndpoint"],
            "status": conflict_exp["status"],
            "verified": bool(conflict_just["verified"]),
            "subset_minimal": bool(conflict_just["subsetMinimal"]),
            "axiom_count": len(conflict_just["axioms"]),
            "axioms": [a["functionalSyntax"] for a in conflict_just["axioms"]],
        },
        "ablations": {
            "single_marker_rubisco_without_11_step_conjunction": {
                "entailed_cbb_mags": len(entailed["RuBisCOMAG"]),
                "additional_incomplete_screen_mags": len(entailed["RuBisCOMAG"] - entailed["CBBCompleteMAG"]),
                "entailed_h2_cbb_mags": len(rbcl_h2_mags),
                "additional_h2_incomplete_screen_mags": len(rbcl_h2_mags - entailed["H2AutotrophicCandidateMAG"]),
            },
            "two_gene_rubisco_prk_without_11_step_conjunction": {
                "entailed_cbb_mags": len(entailed["RuBisCOPRKMAG"]),
                "additional_incomplete_screen_mags": len(entailed["RuBisCOPRKMAG"] - entailed["CBBCompleteMAG"]),
                "entailed_h2_cbb_mags": len(prk_h2_mags),
                "additional_h2_incomplete_screen_mags": len(prk_h2_mags - entailed["H2AutotrophicCandidateMAG"]),
            },
            "additive_hydrogenase_guilds_without_carrier_union": {
                "multi_hydrogenase_mags_double_counted": multi_h2_mags,
            },
            "sample_marginals_without_genome_linkage": {
                "usable_mapped_samples": int(len(valid)),
                "missing_mapped_samples": int(panel.mapped_read_pct.isna().sum()),
                "nonzero_ambiguity_samples": int((width > 1e-12).sum()),
                "frechet_width_median": float(np.median(width)),
                "frechet_width_max": float(np.max(width)),
                "linked_intersection_median": float(np.median(linked)),
                "linked_intersection_max": float(np.max(linked)),
                "linked_within_frechet_bounds_all_samples": bool(
                    np.all((linked >= lower - 1e-12) & (linked <= upper + 1e-12))
                ),
            },
        },
        "observation_maps_summary": {"endpoints": len(ALL_ENDPOINTS), "disjoint_pairs": observation["disjoint_pairs"],
                                     "declared_maps": sorted(OBSERVATION_MAPS)},
        "justification_summary": ({"counts": justifications["counts"], "wall_seconds": justifications["wall_seconds"],
                                   "workers": justifications["workers"]} if justifications else None),
        "_observation_table": observation,
        "_justifications": justifications,
        "entailed_candidate_mags": sorted(entailed["H2AutotrophicCandidateMAG"]),
        "entailed_high_quality_candidate_mags": sorted(entailed["HighQualityH2AutotrophicCandidateMAG"]),
        "entailed_memberships": {str(rep): sorted(cls for cls, ids in entailed.items() if rep in ids)
                                  for rep in mag.rep},
    }


def rules_digest() -> str:
    return hashlib.sha256("\n".join(_tbox_lines()).encode()).hexdigest()


def export_classification(report: dict, data: Path) -> None:
    """Publish actual reasoner outputs tied to their exact evidence and TBox."""
    report = dict(report)
    observation = report.pop("_observation_table", None)
    justifications = report.pop("_justifications", None)
    artifact = {"schema_version": 1, "rules_sha256": rules_digest(),
                "evidence_sha256": sha256(data / "mag_evidence.tsv"),
                "ontology_sha256": report["ontology_sha256"],
                "reasoner_binary_sha256": report["reasoner_binary_sha256"],
                "classes_by_mag": report["entailed_memberships"]}
    names = ["entailed_memberships.json", "ontology_report.json"]
    if observation is not None:
        export_observation_maps(observation, data, report["ontology_sha256"])
        names.append("observation_maps.json")
    if justifications is not None:
        jart = {"schema_version": 1, "rules_sha256": rules_digest(), "ontology_sha256": report["ontology_sha256"],
                "reasoner_binary_sha256": report["reasoner_binary_sha256"], **justifications}
        (data / "justifications.json").write_text(json.dumps(jart, indent=1, sort_keys=True) + "\n")
        names.append("justifications.json")
    (data / "entailed_memberships.json").write_text(json.dumps(artifact, indent=2, sort_keys=True) + "\n")
    (data / "ontology_report.json").write_text(json.dumps(report, indent=2) + "\n")
    manifest_path = data / "manifest.json"
    manifest = json.loads(manifest_path.read_text())
    if (data / "ontology_mappings.tsv").exists():
        names.append("ontology_mappings.tsv")
    for name in names:
        manifest["artifacts"][name] = sha256(data / name)
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")


def ontology_path(data: Path) -> Path:
    """Shipped OWL file: data/../ontology (repository and wheel layouts share this relative layout)."""
    for cand in (data.parent / "ontology" / OFN_NAME, data / OFN_NAME):
        if cand.exists():
            return cand
    raise FileNotFoundError(f"Shipped ontology {OFN_NAME} not found next to {data}")


def load_classification(data: Path) -> dict:
    manifest = json.loads((data / "manifest.json").read_text())
    for name in ("entailed_memberships.json", "ontology_report.json"):
        if sha256(data / name) != manifest["artifacts"].get(name):
            raise ValueError("OWL artifact checksum mismatch")
    artifact = json.loads((data / "entailed_memberships.json").read_text())
    if artifact["rules_sha256"] != rules_digest() or artifact["evidence_sha256"] != sha256(data / "mag_evidence.tsv"):
        raise ValueError("Stale OWL classification: rebuild evidence and run the reasoner")
    report = json.loads((data / "ontology_report.json").read_text())
    ofn_digest = sha256(ontology_path(data))
    if not (ofn_digest == artifact["ontology_sha256"] == report["ontology_sha256"]):
        raise ValueError("Ontology digest mismatch: shipped OWL file differs from the classified ontology")
    if not (artifact["reasoner_binary_sha256"] == report["reasoner_binary_sha256"] == KM_SHA256):
        raise ValueError("Reasoner digest mismatch: recorded reasoner binary is not the pinned km release")
    return artifact
