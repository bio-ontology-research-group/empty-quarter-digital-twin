"""v0.6 causal-ontological coupling: graph, translation, community factors and invariance checks."""
from dataclasses import replace
import json

import pytest
from fastapi.testclient import TestClient

from eq_biodt.v03 import community as comm
from eq_biodt.v03 import theory
from eq_biodt.v03.api import app, DATA
from eq_biodt.v03.interventions import intervene
from eq_biodt.v03.planning import compare, parameter_ensemble
from eq_biodt.v03.process import Parameters, State, simulate

P = Parameters(k_h2_s=.001, k_co_s=.001)
DRY = [{'water_mm': 0, 'evaporation_mm': 1, 'temp_c': 25} for _ in range(3)]
STILL = [{'water_mm': 0, 'evaporation_mm': 0, 'temp_c': 25} for _ in range(3)]
GOLDEN_SCENARIO_H2 = 0.00028263236248340614      # data/v03/example_scenario.json before v0.6
GOLDEN_INTERVENTION_DH2 = 0.0011023333247692896  # data/v03/example_intervention.json before v0.6


def example(name):
    return json.loads((DATA/f'example_{name}.json').read_text())


def test_reference_coupling_reproduces_previous_outputs_exactly():
    v = example('scenario')
    s = simulate(State(**v['initial']), v['forcing'], Parameters(**v['parameters']))
    assert s['summary']['h2_uptake_mol_m2'] == GOLDEN_SCENARIO_H2
    v = example('intervention')
    r = intervene(State(**v['initial']), v['forcing'], Parameters(**v['parameters']), v['intervention'], DATA)
    assert r['difference']['h2_uptake_mol_m2'] == GOLDEN_INTERVENTION_DH2
    assert r['community']['mode'] == 'reference' and r['deduced'] == []
    assert all(r['checks'].values()) and r['trace']['definitions_and_annotations_invariant']
    explicit = Parameters(**v['parameters'], gamma_h2=1., gamma_co=1.)
    assert simulate(State(**v['initial']), v['forcing'], explicit)['trajectory'] == \
        simulate(State(**v['initial']), v['forcing'], Parameters(**v['parameters']))['trajectory']
    with pytest.raises(ValueError):
        Parameters(k_h2_s=.001, k_co_s=.001, gamma_h2=-1.)


def test_gate_closes_co_for_zero_share_sample_and_contrast_is_exactly_zero():
    closed = comm.closed_gates()
    assert closed['co_closed'], 'expected samples with no COCarrierMAG weight'
    sample = closed['co_closed'][0]
    s = comm.carrier_shares(sample)
    assert s['shares']['COCarrierMAG'] == 0 and comm.gamma(s['shares'], 'gate')['co'] == 0
    r = intervene(State(1.8), DRY, P, {'forcing': {'0': {'water_mm': 12}}, 'community': {'sample': sample}, 'coupling': 'gate'}, DATA)
    assert r['gas_uptake']['co_mol_m2']['baseline'] == 0 and r['gas_uptake']['co_mol_m2']['difference'] == 0
    assert r['gas_uptake']['h2_mol_m2']['difference'] > 0
    d, = r['deduced']
    assert d['gas'] == 'co' and d['effect'] == 'zero_by_deduction' and d['share'] == 0 and sample in d['samples']
    ref = intervene(State(1.8), DRY, P, {'forcing': {'0': {'water_mm': 12}}}, DATA)
    assert ref['gas_uptake']['co_mol_m2']['difference'] != 0


def test_gamma_logic_on_synthetic_shares():
    shares = {c: 0. for c in comm.CLASSES}
    shares['H2CarrierMAG'] = .46
    assert comm.gamma(shares, 'reference') == {'h2': 1., 'co': 1.}
    assert comm.gamma(shares, 'gate') == {'h2': 1., 'co': 0.}
    assert comm.gamma(shares, 'proportional')['h2'] == pytest.approx(.46/comm.medians()['H2CarrierMAG'])
    with pytest.raises(ValueError, match='mode'):
        comm.gamma(shares, 'bogus')


def test_proportional_mode_is_monotone_in_share():
    t = comm._load()['table'].sort_values('H2CarrierMAG')
    lo, hi = t.index[len(t)//10], t.index[-len(t)//10]
    up = {}
    for sample in (lo, hi):
        r = intervene(State(1.8), DRY, P, {'forcing': {'0': {'water_mm': 12}}, 'community': {'sample': sample}, 'coupling': 'proportional'}, DATA)
        up[sample] = r['gas_uptake']['h2_mol_m2']['intervened']
    assert t.H2CarrierMAG[lo] < t.H2CarrierMAG[hi] and up[lo] < up[hi]
    u = [simulate(State(5.), DRY, replace(P, gamma_h2=g))['summary']['h2_uptake_mol_m2'] for g in (0., .5, 1., 2.)]
    assert u[0] == 0 and u == sorted(u) and len(set(u)) == 4


def test_visit_averages_replicates_and_lifted_consistency():
    v = comm.carrier_shares(community={'site': 10, 'compartment': 'Deep', 'trip': 1})
    assert v['samples'] and set(v['shares']) == set(comm.CLASSES) and v['evidence_sha256'] and v['read_mapping_sha256']
    checks = comm.lifted_consistency()
    assert checks and all(checks.values())
    with pytest.raises(ValueError):
        comm.carrier_shares(community={'site': 10})


@pytest.mark.parametrize('target,fragment', [
    ('theta', 'constitutive'), ('J_h2', 'causal outcome mechanism without physical handle'),
    ('carrier_share', 'deduced from evidence and ontology'), ('membership', 'observation update'),
    ('psi', 'version')])
def test_inadmissible_targets_raise_with_zone_message(target, fragment):
    with pytest.raises(ValueError, match=fragment):
        theory.check_admissible([target])
    with pytest.raises(ValueError, match=fragment):
        intervene(State(1.8), DRY, P, {'forcing': {'0': {target: .1}}}, DATA)
    assert theory.check_admissible(['water_mm', 'porosity', 'W', 'N']) == ['I', 'phi', 'W', 'N']


def test_state_targets_translate_and_validate():
    a = intervene(State(5.), DRY, P, {'state': {'1': {'theta': .1}}}, DATA)
    b = intervene(State(5.), DRY, P, {'state': {'1': {'water_mm': 10.}}}, DATA)
    assert a['intervened']['trajectory'][1]['water_mm'] == pytest.approx(b['intervened']['trajectory'][1]['water_mm'], rel=1e-12)
    assert a['intervened']['trajectory'][1]['state_override']['water_mm'] == pytest.approx(10.)
    assert a['baseline']['trajectory'] == b['baseline']['trajectory'] and all(a['checks'].values())
    for bad in ({'theta': .5}, {'psi_mpa': -5.}, {'psi_mpa': -5., 'version': 'x'}, {'psi_mpa': -5., 'version': 'matric', 'water_mm': 3.},
                {'J_h2': 1.}, {'temp_c': 20.}, {'water_mm': 1e6}):
        with pytest.raises(ValueError):
            intervene(State(5.), DRY, P, {'state': {'0': bad}}, DATA)
    with pytest.raises(ValueError, match='matric'):
        intervene(State(5.), DRY, P, {'state': {'0': {'psi_mpa': -5.}}}, DATA)
    with pytest.raises(ValueError):
        theory.translate_constituted('psi', -.0001, State(20.), P, 'matric')     # above air-entry potential


def test_psi_versions_agree_on_activation_and_differ_on_uptake():
    psi = -.3
    out = {}
    for version in ('matric', 'osmotic'):
        out[version] = intervene(State(20.), STILL, P, {'state': {'0': {'psi_mpa': psi, 'version': version}}}, DATA)
    for version, r in out.items():
        row = r['intervened']['trajectory'][0]
        assert row['gas']['total_water_potential_mpa'] == pytest.approx(psi, rel=1e-9)
        assert r['state_interventions'][0]['resolved']['psi_mpa']['handle'] == {'matric': 'W', 'osmotic': 'N'}[version]
    a = [out[v]['intervened']['trajectory'][0]['gas']['activation'] for v in out]
    j = [out[v]['intervened']['trajectory'][0]['gas']['h2_mol_m2_s'] for v in out]
    assert a[0] == pytest.approx(a[1], rel=1e-9) and 0 < a[0] < 1
    assert abs(j[0]-j[1]) > 1e-3*abs(j[0])
    rel = out['matric']['version_relevance']
    assert rel['A']['relevant'] is False and rel['J_h2']['relevant'] is True
    assert ['W', 'theta', 'D_h2', 'J_h2'] in rel['J_h2']['bypass_paths']


def test_composition_setting_input_to_its_baseline_value_changes_nothing():
    r = intervene(State(1.8), DRY, P, {'forcing': {'1': {'water_mm': 0, 'temp_c': 25}}, 'parameters': {'porosity': .4}}, DATA)
    assert r['intervened'] == r['baseline'] and all(v == 0 for v in r['difference'].values() if v is not None)
    r = intervene(State(1.8), DRY, P, {'state': {'0': {'water_mm': 1.8, 'nacl_mol_m2': 0.}}}, DATA)
    assert r['intervened']['trajectory'][0]['water_mm'] == r['baseline']['trajectory'][0]['water_mm']
    assert r['difference']['h2_uptake_mol_m2'] == 0


def test_invariance_checks_all_true_and_failure_raises(monkeypatch):
    v = example('intervention')
    r = intervene(State(**v['initial']), v['forcing'], Parameters(**v['parameters']), v['intervention'], DATA)
    assert r['checks'] and all(r['checks'].values()) and all(r['subsumption_checks'].values())
    assert {'memberships_digest_equal_across_arms', 'gamma_equal_across_arms', 'intervened.power_identity'} <= set(r['checks'])
    monkeypatch.setattr(comm, 'SUBSUMPTIONS', comm.SUBSUMPTIONS + (('RuBisCOMAG', 'H2AutotrophicCandidateMAG'),))
    with pytest.raises(ValueError, match='Invariance checks failed'):
        intervene(State(1.8), DRY, P, {'forcing': {'0': {'water_mm': 3}}}, DATA)


def test_graph_is_acyclic_and_well_formed():
    g = theory.graph()
    json.dumps(g, allow_nan=False)
    ids = {n['id'] for n in g['nodes']}
    assert len(ids) == len(g['nodes']) and theory.slice_is_acyclic()
    assert {n['zone'] for n in g['nodes']} == set(theory.ZONES)
    for n in g['nodes']:
        assert set(n['parents']) <= ids and n['equation'] and n['unit'] is not None
    for e in g['edges']:
        assert e['source'] in ids and e['target'] in ids and e['kind'] in {'mechanism', 'definition', 'deduction', 'lag'}
        assert (e['kind'] == 'lag') == (e['lag'] == 1)
    assert {(e['source'], e['target']) for e in g['edges'] if e['kind'] == 'lag'} == {('W', 'W'), ('N', 'N')}
    assert theory.version_relevance('psi', 'A')['relevant'] is False
    assert theory.version_relevance('psi', 'J_h2')['relevant'] is True
    assert theory.version_relevance('theta', 'J_h2')['unique_translation'] is True


def test_treatment_versions_and_community_in_compare():
    pulse = [{'water_mm': 12, 'evaporation_mm': 1, 'temp_c': 25}] + [{'water_mm': 0, 'evaporation_mm': 1, 'temp_c': 25}]*3
    split = [{'water_mm': 3, 'evaporation_mm': 1, 'temp_c': 25} for _ in range(4)]
    dry = [{'water_mm': 0, 'evaporation_mm': 1, 'temp_c': 25}]*4
    scenarios = [{'name': 'control', 'forcing': dry, 'disturbance_m2': 0.}, {'name': 'pulse', 'forcing': pulse, 'disturbance_m2': 0.},
                 {'name': 'split', 'forcing': split, 'disturbance_m2': 0.}]
    members = parameter_ensemble(P)[:4]
    r = compare(State(1.8), scenarios, members)
    tv = r['treatment_versions']
    g, = tv['groups']
    assert g['versions'] == ['pulse', 'split'] and tv['control'] == 'control' and len(g['members']) == 4
    assert tv['version_relevant'] is True and tv['effect_vs_control_range_mol_m2'][0] <= tv['effect_vs_control_range_mol_m2'][1]
    assert r['community']['mode'] == 'reference'
    closed = comm.closed_gates()['co_closed'][0]
    gated = compare(State(1.8), scenarios, members, {'sample': closed}, 'gate')
    assert gated['community']['gamma']['co'] == 0 and all(p['gamma_co'] == 0 for p in gated['parameters'])


def test_api_graph_community_and_version_relevance():
    client = TestClient(app)
    assert client.get('/v3/graph').json()['nodes']
    c = client.get('/v3/community/10Dr2').json()
    assert set(c['gamma']) == set(comm.MODES) and all(c['consistency'].values())
    r = client.post('/v3/version-relevance', json={'target': 'psi', 'outcome': 'J_h2'}).json()
    assert r['relevant'] is True
    body = {'parameters': {'k_h2_s': .001, 'k_co_s': .001}, 'initial': {'water_mm': 1.8}, 'forcing': DRY,
            'intervention': {'forcing': {'0': {'water_mm': 12}}}, 'community': {'sample': comm.closed_gates()['co_closed'][0]}, 'coupling': 'gate'}
    out = client.post('/v3/intervene', json=body)
    assert out.status_code == 200 and out.json()['deduced'][0]['gas'] == 'co'
    bad = client.post('/v3/intervene', json={**body, 'intervention': {'forcing': {'0': {'theta': .2}}}, 'community': None, 'coupling': None})
    assert bad.status_code == 422 and 'constitutive' in bad.json()['detail']
