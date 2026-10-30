"""Paired, version-specific regimes under the stratified causal-ontological model.

An intervention replaces exogenous inputs, prescribed soil parameters or the causal states W and N
(directly, or through the translation of a constituted target theta / psi). Definitions, evidence
and deduced memberships are never touched; the invariance checks below recompute them in each arm.
See docs/causal-ontological-theory.md. This is an assumed structural model, not observational
causal identification.
"""
from copy import deepcopy
from dataclasses import replace
import hashlib
import json
from pathlib import Path
from math import isfinite

from . import community as comm
from . import theory
from .contract import contract
from .ontology import load_classification
from .process import Parameters, State, simulate

FORCING_TARGETS = {'water_mm', 'nacl_mol_l', 'evaporation_mm', 'temp_c'}
PARAMETER_TARGETS = {'porosity', 'psi_entry_mpa', 'retention_b'}
STATE_TARGETS = {'water_mm', 'nacl_mol_m2', 'theta', 'psi_mpa'}
TOP_KEYS = {'forcing', 'parameters', 'state', 'community', 'coupling'}


def fingerprint(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, allow_nan=False).encode()).hexdigest()


def _number(v, what):
    if isinstance(v, bool) or not isinstance(v, (int, float)) or not isfinite(v):
        raise ValueError(f'{what} must be a finite number')


def _day(key, n):
    if not isinstance(key, str) or not key.isdecimal() or str(int(key)) != key or not 0 <= int(key) < n:
        raise ValueError('Intervention day must be a zero-based day in the baseline schedule')
    return int(key)


def _close(a, b):
    return abs(a - b) <= 1e-9 * max(1., abs(a), abs(b))


def _state_override(spec: dict, temp_c: float, params: Parameters, record: list):
    """Build the State -> State map for one day from a validated state intervention."""
    version = spec.get('version')
    keys = set(spec) - {'version'}
    for k in keys:
        if k not in STATE_TARGETS:
            theory.check_admissible([k])               # raises the zone message for non-manipulable targets
            raise ValueError(f'{k} is an exogenous input, not a state; set it under forcing or parameters')
        _number(spec[k], f'State target {k}')
    if not keys:
        raise ValueError('Empty state intervention')
    if 'psi_mpa' in keys and keys != {'psi_mpa'}:
        raise ValueError('psi_mpa is set through W or N; do not combine it with another state target')
    if {'water_mm', 'theta'} <= keys:
        raise ValueError('Specify water_mm or theta, not both: theta = W/(1000 L)')
    if version is not None and 'psi_mpa' not in keys:
        raise ValueError('version applies to psi_mpa targets only')

    def apply(state: State) -> State:
        water, salt = state.water_mm, state.salt_mol_m2
        resolved = {}
        if 'psi_mpa' in keys:
            t = theory.translate_constituted('psi', spec['psi_mpa'], state, params, version, temp_c)
            water, salt = t['water_mm'], t['salt_mol_m2']
            resolved['psi_mpa'] = t
        if 'theta' in keys:
            t = theory.translate_constituted('theta', spec['theta'], state, params)
            water = t['water_mm']
            resolved['theta'] = t
        if 'water_mm' in keys:
            water = spec['water_mm']
        if 'nacl_mol_m2' in keys:
            salt = spec['nacl_mol_m2']
        if water < 0 or salt < 0:
            raise ValueError('State intervention must be non-negative')
        record.append({'spec': spec, 'resolved': resolved, 'state_before': {'water_mm': state.water_mm, 'salt_mol_m2': state.salt_mol_m2},
                       'state_after': {'water_mm': water, 'salt_mol_m2': salt}})
        return State(water, salt)
    return apply


def _gas_totals(result):
    rows = result['trajectory']
    if any(r['gas'] is None for r in rows):
        return {'h2_mol_m2': None, 'co_mol_m2': None}
    return {'h2_mol_m2': sum(r['gas']['h2_mol_m2_s']*86400 for r in rows),
            'co_mol_m2': sum(r['gas']['co_mol_m2_s']*86400 for r in rows)}


def _arm_checks(arm: str, result: dict, p: Parameters) -> dict:
    """Constitutive identities along one trajectory (theta, W_max, balances, power)."""
    wmax = 1000*p.depth_m*p.porosity
    rows = result['trajectory']
    gas = [r['gas'] for r in rows if r['gas'] is not None]
    return {
        f'{arm}.water_within_pore_storage': all(0 <= r['water_mm'] <= wmax*(1+1e-12) for r in rows),
        f'{arm}.theta_identity': all(0 <= r['water_mm']/(1000*p.depth_m) <= p.porosity*(1+1e-12) for r in rows),
        f'{arm}.water_balance_closure': all(abs(r['water_balance_residual_mm']) < 1e-9 for r in rows),
        f'{arm}.salt_balance_closure': all(abs(r['salt_balance_residual_mol_m2']) < 1e-9 for r in rows),
        f'{arm}.power_identity': all(_close(g['chemical_power_w_m2'], 202000*g['h2_mol_m2_s'] + 243000*g['co_mol_m2_s']) for g in gas),
        f'{arm}.captured_power_identity': all(_close(g['captured_power_w_m2'], p.energy_capture_fraction*g['chemical_power_w_m2']) for g in gas),
    }


def community_factors(community: dict | None, coupling: str | None) -> dict:
    """Resolve the evidence-selecting community and coupling mode to gamma (conditioning, not intervention)."""
    if community is None:
        if coupling not in (None, 'reference'):
            raise ValueError(f'Coupling mode {coupling!r} needs a community (sample or site visit)')
        return {'mode': 'reference', 'community': None, 'gamma': {'h2': 1., 'co': 1.}, 'shares': None, 'samples': [], 'provenance': None}
    mode = coupling or 'gate'
    s = comm.carrier_shares(community=community)
    return {'mode': mode, 'community': community, 'gamma': comm.gamma(s['shares'], mode), 'shares': s['shares'],
            'samples': s['samples'],
            'provenance': {k: s[k] for k in ('evidence_sha256', 'ontology_sha256', 'read_mapping_sha256')}}


def intervene(initial: State, forcing: list[dict], p: Parameters, intervention: dict, data: Path) -> dict:
    """Replace selected inputs, prescribed soil properties or causal states in the intervened arm only."""
    if not forcing:
        raise ValueError('At least one forcing day is required')
    if not isinstance(intervention, dict) or set(intervention) - TOP_KEYS:
        raise ValueError('Interventions may target forcing, prescribed soil parameters or causal state (W, N) only; '
                         'derived endpoints, evidence and ontology are not intervention targets')
    changes = intervention.get('parameters', {})
    days = intervention.get('forcing', {})
    state_spec = intervention.get('state', {})
    if not all(isinstance(x, dict) for x in (changes, days, state_spec)):
        raise ValueError('Intervention forcing, parameters and state must be objects')
    theory.check_admissible(list(changes))
    if set(changes) - PARAMETER_TARGETS:
        raise ValueError('Unsupported parameter intervention; specify porosity and, if measured, retention parameters')
    for v in changes.values():
        _number(v, 'Physical parameter interventions')
    modified = deepcopy(forcing)
    targets, nodes = [], set(theory.check_admissible(list(changes)))
    for key, replacements in days.items():
        d = _day(key, len(forcing))
        if not isinstance(replacements, dict):
            raise ValueError('Forcing intervention must be an object')
        nodes |= set(theory.check_admissible(list(replacements)))
        if set(replacements) - FORCING_TARGETS:
            raise ValueError('Only water, specified NaCl, evaporation and temperature can be clamped')
        for v in replacements.values():
            _number(v, 'Forcing interventions')
        modified[d].update(replacements)
        targets.extend(f'forcing[{key}].{k}' for k in sorted(replacements))
    targets.extend('parameters.'+k for k in sorted(changes))
    # community selects evidence; it conditions both arms identically
    factors = community_factors(intervention.get('community'), intervention.get('coupling'))
    base = replace(p, gamma_h2=factors['gamma']['h2'], gamma_co=factors['gamma']['co']) if factors['community'] else p
    changed_p = replace(base, **changes)
    records, overrides, psi_used = [], {}, False
    for key, spec in state_spec.items():
        d = _day(key, len(forcing))
        if not isinstance(spec, dict):
            raise ValueError('State intervention must be an object')
        psi_used |= 'psi_mpa' in spec
        overrides[d] = _state_override(spec, modified[d]['temp_c'], changed_p, rec := [])
        records.append((d, rec))
        for k in sorted(set(spec) - {'version'}):
            targets.append(f'state[{key}].{k}')
            nodes.add({'water_mm': 'W', 'nacl_mol_m2': 'N'}.get(k, k))
    classified = load_classification(data)
    baseline = simulate(initial, forcing, base)
    changed = simulate(initial, modified, changed_p, overrides or None)
    delta = {}
    for key, value in baseline['summary'].items():
        other = changed['summary'][key]
        if isinstance(value, (int, float)) and isinstance(other, (int, float)):
            delta[key] = other - value
        elif value is None or other is None:
            delta[key] = None
    # uptake of each gas, with the contrast
    gb, gc = _gas_totals(baseline), _gas_totals(changed)
    uptake = {g: {'baseline': gb[g], 'intervened': gc[g],
                  'difference': None if gb[g] is None or gc[g] is None else gc[g] - gb[g]} for g in gb}
    # deduced block: gamma = 0 closes the gas in every admissible regime
    deduced = []
    for gas, cls in comm.GAS_CLASS.items():
        if factors['community'] and factors['gamma'][gas] == 0:
            d = uptake[gas+'_mol_m2']['difference']
            deduced.append({'gas': gas, 'effect': 'zero_by_deduction', 'class': cls, 'share': factors['shares'][cls],
                            'samples': factors['samples'], 'simulated_contrast_mol_m2': d,
                            'derivation': f'q_s({cls}) = 0 under the reasoner memberships, so gamma_{gas} = 0, k_{gas} = 0 and J_{gas} = 0 '
                                          'in every admissible regime; absence of entailment is not read as biological absence.'})
    # invariance checks, recomputed per arm
    arms = {}
    for arm in ('baseline', 'intervened'):
        arm_class = load_classification(data)
        arm_factors = community_factors(intervention.get('community'), intervention.get('coupling'))
        arms[arm] = {'digest': fingerprint(arm_class['classes_by_mag']), 'shares': arm_factors['shares'], 'gamma': arm_factors['gamma'],
                     'classes': arm_class['classes_by_mag']}
    checks = {'memberships_digest_equal_across_arms': arms['baseline']['digest'] == arms['intervened']['digest'],
              'community_shares_equal_across_arms': arms['baseline']['shares'] == arms['intervened']['shares'],
              'gamma_equal_across_arms': (arms['baseline']['gamma'] == arms['intervened']['gamma']
                                          and (baseline['parameters']['gamma_h2'], baseline['parameters']['gamma_co'])
                                          == (changed['parameters']['gamma_h2'], changed['parameters']['gamma_co'])
                                          == (arms['intervened']['gamma']['h2'], arms['intervened']['gamma']['co'])
                                          if factors['community'] else True),
              'initial_state_shared': baseline['initial'] == changed['initial'],
              'untargeted_forcing_shared': all(forcing[i] == modified[i] for i in range(len(forcing)) if str(i) not in days)}
    sub = comm.membership_consistency(arms['baseline']['classes'])
    if factors['community']:
        sub.update(comm.lifted_consistency(factors['shares']))
    checks['subsumption_inequalities_hold'] = all(sub.values())
    checks.update(_arm_checks('baseline', baseline, base))
    checks.update(_arm_checks('intervened', changed, changed_p))
    checks['deduced_null_effects_match_simulation'] = all(x['simulated_contrast_mol_m2'] in (0., None) for x in deduced)
    failed = sorted(k for k, v in checks.items() if not v)
    if failed:
        raise ValueError(f'Invariance checks failed: {failed}')
    relevance = ({o: theory.version_relevance('psi', o) for o in ('A', 'J_h2', 'J_co')} if psi_used else None)
    return {'baseline': baseline, 'intervened': changed, 'difference': delta,
            'intervention': deepcopy(intervention),
            'gas_uptake': uptake, 'deduced': deduced,
            'community': {k: factors[k] for k in ('mode', 'community', 'gamma', 'shares', 'samples', 'provenance')},
            'state_interventions': [{'day': d, **r} for d, rec in records for r in rec],
            'version_relevance': relevance,
            'checks': checks, 'subsumption_checks': sub,
            'trace': {'targets_replaced': targets, 'admissible_nodes': sorted(nodes),
                      'contract_sha256': fingerprint(contract()),
                      'ontology_sha256': classified['ontology_sha256'],
                      'evidence_sha256': classified['evidence_sha256'],
                      'process_code_sha256': hashlib.sha256(Path(__file__).with_name('process.py').read_bytes()).hexdigest(),
                      'graph_sha256': fingerprint(theory.graph()),
                      'initial_state_shared': checks['initial_state_shared'],
                      'definitions_and_annotations_invariant': all(checks.values()),
                      'endpoint': 'conditional integrated H2 uptake, not genomic capacity or carbon removal',
                      'scope': 'Specified regimes under assumed mechanisms: exogenous inputs, prescribed soil parameters and causal states W, N. '
                               'Constituted targets are translated; evidence, ontology and deduced memberships are not targets. '
                               'No observational identification or transfer of Boolean-calculus soundness theorems.'}}
