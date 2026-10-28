"""Independent source and counterexample regressions for the publication audit."""
from dataclasses import replace
import json
from pathlib import Path
import sys

import pandas as pd
import pytest
from fastapi.testclient import TestClient

from eq_biodt.v03.api import app, DATA
from eq_biodt.v03.evidence import source_directory
from eq_biodt.v03.process import Parameters, State, simulate
from eq_biodt.v03.ontology import load_classification

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "reproducibility"))
from independent_reference import satisfied, module_definition


def test_all_pathway_calls_against_independent_pinned_kegg_expression():
    source = source_directory()
    expression = module_definition((source / "M00165.txt").read_text())
    raw = pd.read_csv(source / "rep_ko_v2.tsv.gz", sep="\t", names=["rep", "kos"], keep_default_na=False)
    expected = raw.set_index("rep").kos.map(lambda s: satisfied(expression, set(s.split(","))))
    actual = pd.read_csv(DATA / "mag_evidence.tsv", sep="\t").set_index("rep")
    assert expected.reindex(actual.index).equals(actual.cbb_all_steps_detected)
    assert int(expected.sum()) == 232
    assert int((actual.cbb_all_steps_detected & ~actual.sbpase_detected).sum()) == 149
    transaldolase_only = actual.loc[actual.cbb_all_steps_detected & ~actual.sbpase_detected]
    assert transaldolase_only.transaldolase_detected.all()


def test_mag_without_co_evidence_cannot_use_co_for_energy():
    mags = pd.read_csv(DATA / "mag_evidence.tsv", sep="\t")
    row = mags.loc[mags.h2_cbb_all_steps & ~mags.co_form1].iloc[0]
    request = dict(mag_id=row.rep, h2=1e-9, co=1e-8, oxygen=1e-8,
                   external_co2=1e-8, capture_fraction=.25, maintenance_w_m2=1e-5,
                   synthesis_j_mol_c=500000)
    response = TestClient(app).post("/v3/resource-ceiling", json=request)
    assert response.status_code == 200
    result = response.json()
    assert result["fluxes"]["co_respiration"] == 0
    assert result["fluxes"]["h2_respiration"] > 0
    assert result["genomic_evidence"]["ontology_sha256"]
    assert TestClient(app).post("/v3/resource-ceiling", json=request | {"co_oxidation_supported": True}).status_code == 422


def test_subdaily_gas_integral_resolves_short_wetting_pulses():
    p = Parameters(k_h2_s=.001, k_co_s=.001)
    f = [{"water_mm":12/14, "evaporation_mm":1, "temp_c":25} for _ in range(14)]
    default = simulate(State(1.8), f, p)["summary"]["h2_uptake_mol_m2"]
    fine = simulate(State(1.8), f, replace(p, integration_steps_per_day=384))["summary"]["h2_uptake_mol_m2"]
    assert default > 0
    assert default == pytest.approx(fine, rel=0.001)


def test_tampered_reasoner_export_rejected(tmp_path):
    for name in ("mag_evidence.tsv", "manifest.json", "entailed_memberships.json", "ontology_report.json"):
        (tmp_path / name).write_bytes((DATA / name).read_bytes())
    artifact = json.loads((tmp_path / "entailed_memberships.json").read_text())
    artifact["classes_by_mag"][next(iter(artifact["classes_by_mag"]))].append("COCarrierMAG")
    (tmp_path / "entailed_memberships.json").write_text(json.dumps(artifact))
    with pytest.raises(ValueError, match="checksum"):
        load_classification(tmp_path)


def test_interventions_match_flat_calculation_and_preserve_invariant_evidence():
    from copy import deepcopy
    from eq_biodt.v03.interventions import intervene
    p = Parameters(k_h2_s=.001, k_co_s=.001)
    forcing = [{'water_mm':0, 'evaporation_mm':1, 'temp_c':25} for _ in range(3)]
    original = deepcopy(forcing)
    result = intervene(State(1.8), forcing, p, {'forcing':{'0':{'water_mm':12}}}, DATA)
    flat = deepcopy(forcing); flat[0]['water_mm'] = 12
    assert result['intervened'] == simulate(State(1.8), flat, p)
    assert forcing == original
    assert result['trace']['definitions_and_annotations_invariant']
    assert result['difference']['net_co2_removal_g_m2'] is None
    assert result['difference']['h2_uptake_mol_m2'] > 0
    for invalid in ({'h2_uptake_mol_m2':1}, {'parameters':{'cbb_all_steps_detected':True}},
                    {'forcing':{'3':{'water_mm':1}}}, {'parameters':{'porosity':.001}}):
        with pytest.raises(ValueError):
            intervene(State(1.8), forcing, p, invalid, DATA)
    request = {'parameters':{'k_h2_s':.001,'k_co_s':.001},'initial':{'water_mm':1.8},
               'forcing':forcing,'intervention':{'forcing':{'0':{'water_mm':12}}}}
    response = TestClient(app).post('/v3/intervene', json=request)
    assert response.status_code == 200
    assert response.json()['intervened'] == result['intervened']
