"""v0.6 ontology layer: standard-ontology mappings, observation semantics, digest checks, justifications."""
import csv
import json
import re
import shutil
from pathlib import Path

import pytest

from eq_biodt.v03 import ontology as ont
from eq_biodt.v03.ontology import load_classification, observation_compatible

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data/v03"
OFN = ROOT / "data/ontology" / ont.OFN_NAME

COUNTS = {"H2CarrierMAG": 1094, "RuBisCOMAG": 1078, "RuBisCOPRKMAG": 416, "ThreeGeneCBBMAG": 315,
          "CBBCompleteMAG": 232, "COCarrierMAG": 198, "H2AutotrophicCandidateMAG": 136,
          "HighQualityH2AutotrophicCandidateMAG": 100, "DualH2COAutotrophicCandidateMAG": 34}


def _tmp_data(tmp_path):
    d = tmp_path / "data" / "v03"
    d.mkdir(parents=True)
    for f in DATA.iterdir():
        if f.is_file():
            shutil.copy(f, d / f.name)
    (tmp_path / "data" / "ontology").mkdir()
    shutil.copy(OFN, tmp_path / "data" / "ontology" / ont.OFN_NAME)
    return d


def test_mapping_table_covers_every_mapping_and_is_verified():
    rows = list(csv.DictReader((DATA / "ontology_mappings.tsv").open(), delimiter="\t"))
    by_key = {(r["rak_term"], r["target_iri"]) for r in rows}
    for term, rel, local, label in ont.MAPPINGS:
        assert (f"rak:{term}", ont._mapping_iri(local)) in by_key
    for r in rows:
        assert r["verified_label"] and r["retrieved_ols4"].startswith("2026-")
    assert {r["source_ontology"] for r in rows} >= {"go", "chebi", "envo", "uo", "sio", "pato", "ncbitaxon"}
    labels = {r["target_curie"]: r["verified_label"] for r in rows}
    assert labels["GO:0019253"] == "reductive pentose-phosphate cycle"
    assert labels["CHEBI:18276"] == "dihydrogen"
    assert labels["ENVO:00005801"] == "rhizosphere"


def test_ofn_contains_mappings_labels_and_definitions():
    text = OFN.read_text()
    assert "skos:relatedMatch rak:RuBisCOMAG <http://purl.obolibrary.org/obo/GO_0016984>" in text
    assert "rak:unit rak:ArealH2FluxEndpoint <http://purl.obolibrary.org/obo/UO_0000041>" in text
    annotated = set(ont.class_annotations())
    declared = set(re.findall(r"^Declaration\(Class\(rak:([A-Za-z0-9_]+)\)\)", text, flags=re.M))
    declared = {c for c in declared if not c.startswith("MAG_") and not c.startswith("ObsProbe")}
    assert not (declared - annotated), declared - annotated
    for name in ("hasKO",):
        assert f"AnnotationAssertion(obo:IAO_0000115 rak:{name}" in text
    assert "at least one gene assigned to that KO" in text
    # a MAG class is never asserted to be a GO process
    assert not re.search(r"SubClassOf\(rak:\w+ obo:GO_", text)


def test_ontology_is_current_and_measurement_partition_is_disjoint():
    assert ("\n".join(ont._tbox_lines()).split("\n")[0]) in OFN.read_text()
    for ep in ont.ALL_ENDPOINTS:
        assert f"SubClassOf(rak:{ep} rak:MeasuredQuantity)" in OFN.read_text()
    for ep in ont.NEW_ENDPOINTS:
        assert f"rak:unitUCUM rak:{ep}" in OFN.read_text()


def test_entailment_counts_unchanged():
    report = json.loads((DATA / "ontology_report.json").read_text())
    assert {k: report["entailed_class_counts"][k] for k in COUNTS} == COUNTS
    classes = load_classification(DATA)["classes_by_mag"]
    for cls, n in COUNTS.items():
        assert sum(cls in v for v in classes.values()) == n


@pytest.mark.parametrize("a,b,status,map_name", [
    ("GravimetricWaterContentEndpoint", "VolumetricWaterContentEndpoint", "compatible_via_map", "GravimetricToVolumetricMap"),
    ("VolumetricWaterContentEndpoint", "GravimetricWaterContentEndpoint", "compatible_via_map", "GravimetricToVolumetricMap"),
    ("FirstOrderRateConstantEndpoint", "MassNormalizedBottleUptakeRateEndpoint", "compatible_via_map", "FirstOrderToBottleUptakeMap"),
    ("MatricWaterPotentialEndpoint", "TotalWaterPotentialEndpoint", "compatible_via_map", "MatricToTotalPotentialMap"),
    ("MassNormalizedBottleUptakeRateEndpoint", "ArealH2FluxEndpoint", "incompatible", None),
    ("ArealH2FluxEndpoint", "MassNormalizedBottleUptakeRateEndpoint", "incompatible", None),
    ("GenomicCapacityEndpoint", "MappedReadCompositionEndpoint", "incompatible", None),
    ("GenomeEquivalentsEndpoint", "MappedReadShareEndpoint", "incompatible", None),
    ("MappedReadShareEndpoint", "MappedReadCompositionEndpoint", "compatible", None),
    ("ArealH2FluxEndpoint", "ArealH2FluxEndpoint", "compatible", None),
])
def test_observation_compatible(a, b, status, map_name):
    r = observation_compatible(a, b, DATA)
    assert r["status"] == status
    assert r["map"] == map_name
    assert r["justification"]


def test_observation_assumptions_and_undeclared_map():
    r = observation_compatible("GravimetricWaterContentEndpoint", "VolumetricWaterContentEndpoint", DATA)
    assert r["assumptions"] == ["BulkDensityAssumption"]
    r = observation_compatible("FirstOrderRateConstantEndpoint", "MassNormalizedBottleUptakeRateEndpoint", DATA)
    assert set(r["assumptions"]) == {"HeadspaceGasAmountInput", "SoilMassInput"}
    r = observation_compatible("MassNormalizedBottleUptakeRateEndpoint", "ArealH2FluxEndpoint", DATA)
    assert r["would_require"] == ["TransportGeometryAssumption"] and "Not declared" in r["justification"]
    with pytest.raises(ValueError):
        observation_compatible("ArealH2FluxEndpoint", "NoSuchEndpoint", DATA)


def test_reasoner_probe_agrees_with_exported_table(tmp_path):
    pytest.importorskip("pandas")
    try:
        ont.find_km()
    except FileNotFoundError:
        pytest.skip("pinned reasoner not installed")
    fresh = ont.build_observation_table(tmp_path)
    exported = json.loads((DATA / "observation_maps.json").read_text())
    assert {k: v["relation"] for k, v in fresh["pairs"].items()} == {k: v["relation"] for k, v in exported["pairs"].items()}


def test_observation_table_tamper_rejected(tmp_path):
    d = _tmp_data(tmp_path)
    t = json.loads((d / "observation_maps.json").read_text())
    t["pairs"]["GenomicCapacityEndpoint|MappedReadCompositionEndpoint"]["relation"] = "subsumed"
    (d / "observation_maps.json").write_text(json.dumps(t))
    with pytest.raises(ValueError, match="checksum"):
        observation_compatible("ArealH2FluxEndpoint", "ArealH2FluxEndpoint", d)


def test_runtime_digest_checks(tmp_path):
    d = _tmp_data(tmp_path)
    assert load_classification(d)["ontology_sha256"] == ont.sha256(d.parent / "ontology" / ont.OFN_NAME)
    # tampered shipped ontology
    ofn = d.parent / "ontology" / ont.OFN_NAME
    ofn.write_text(ofn.read_text() + "\n")
    with pytest.raises(ValueError, match="Ontology digest"):
        load_classification(d)
    ofn.write_bytes(OFN.read_bytes())
    # missing ontology
    ofn.unlink()
    with pytest.raises(FileNotFoundError):
        load_classification(d)
    ofn.write_bytes(OFN.read_bytes())
    # reasoner digest differs from pin (consistent re-signing of both artefacts and manifest)
    for name in ("entailed_memberships.json", "ontology_report.json"):
        x = json.loads((d / name).read_text())
        x["reasoner_binary_sha256"] = "0" * 64
        (d / name).write_text(json.dumps(x))
    m = json.loads((d / "manifest.json").read_text())
    for name in ("entailed_memberships.json", "ontology_report.json"):
        m["artifacts"][name] = ont.sha256(d / name)
    (d / "manifest.json").write_text(json.dumps(m))
    with pytest.raises(ValueError, match="Reasoner digest"):
        load_classification(d)


def test_pinned_reasoner_constant_matches_install_script():
    src = (ROOT / "scripts/install_reasoner.py").read_text()
    assert re.search(rf"SHA256 = '{ont.KM_SHA256}'", src)
    assert re.search(rf"VERSION = '{re.escape(ont.KM_VERSION)}'", src)
    assert json.loads((DATA / "ontology_report.json").read_text())["reasoner_binary_sha256"] == ont.KM_SHA256


def test_wheel_ships_ontology_file():
    assert '"data/ontology" = "eq_biodt/v03/ontology"' in (ROOT / "pyproject.toml").read_text()


def test_all_justifications_exported():
    j = json.loads((DATA / "justifications.json").read_text())
    assert j["counts"] == {"H2AutotrophicCandidateMAG": 136, "COCarrierMAG": 198}
    assert j["all_verified_subset_minimal"] is True
    assert j["ontology_sha256"] == json.loads((DATA / "entailed_memberships.json").read_text())["ontology_sha256"]
    classes = load_classification(DATA)["classes_by_mag"]
    for cls, members in j["classes"].items():
        assert set(members) == {m for m, v in classes.items() if cls in v}
    some = next(iter(j["classes"]["COCarrierMAG"].values()))
    assert any("RAK_GUILD_COXL_FORMI" in a for a in some["mag_axioms"])
    h2 = next(iter(j["classes"]["H2AutotrophicCandidateMAG"].values()))
    assert len(h2["tbox_axiom_ids"]) + len(h2["mag_axioms"]) > 10
