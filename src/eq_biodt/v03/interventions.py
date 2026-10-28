"""Paired, version-specific forcing regimes with invariant evidence and definitions.

This evaluates an assumed structural model; it is not an implementation of
observational causal identification or the Boolean natural-deduction calculus.
"""
from copy import deepcopy
from dataclasses import replace
import hashlib
import json
from pathlib import Path
from math import isfinite

from .contract import contract
from .ontology import load_classification
from .process import Parameters, State, simulate

FORCING_TARGETS = {'water_mm', 'nacl_mol_l', 'evaporation_mm', 'temp_c'}
PARAMETER_TARGETS = {'porosity', 'psi_entry_mpa', 'retention_b'}


def fingerprint(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, allow_nan=False).encode()).hexdigest()


def intervene(initial: State, forcing: list[dict], p: Parameters, intervention: dict, data: Path) -> dict:
    """Replace selected day-level inputs and specified post-disturbance soil properties."""
    if not forcing:
        raise ValueError('At least one forcing day is required')
    if set(intervention) - {'forcing', 'parameters'}:
        raise ValueError('Interventions may target forcing or prescribed soil parameters only; derived endpoints and evidence are invariant')
    changes = intervention.get('parameters', {})
    days = intervention.get('forcing', {})
    if not isinstance(changes, dict) or not isinstance(days, dict):
        raise ValueError('Intervention forcing and parameters must be objects')
    if set(changes) - PARAMETER_TARGETS:
        raise ValueError('Unsupported parameter intervention; specify porosity and, if measured, retention parameters')
    if any(isinstance(v, bool) or not isinstance(v, (int,float)) or not isfinite(v) for v in changes.values()):
        raise ValueError("Physical parameter interventions must be finite numbers")
    modified = deepcopy(forcing)
    targets = []
    for key, replacements in days.items():
        if not isinstance(key, str) or not key.isdecimal() or str(int(key)) != key or not 0 <= int(key) < len(forcing):
            raise ValueError('Forcing intervention day must be a zero-based day in the baseline schedule')
        if not isinstance(replacements, dict) or set(replacements) - FORCING_TARGETS:
            raise ValueError('Only water, specified NaCl, evaporation and temperature can be clamped')
        if any(isinstance(v, bool) or not isinstance(v, (int,float)) or not isfinite(v) for v in replacements.values()):
            raise ValueError("Forcing interventions must be finite numbers")
        modified[int(key)].update(replacements)
        targets.extend(f'forcing[{key}].{k}' for k in sorted(replacements))
    targets.extend('parameters.'+k for k in sorted(changes))
    classified = load_classification(data)
    baseline = simulate(initial, forcing, p)
    changed = simulate(initial, modified, replace(p, **changes))
    delta = {}
    for key, value in baseline['summary'].items():
        other = changed['summary'][key]
        if isinstance(value, (int, float)) and isinstance(other, (int, float)):
            delta[key] = other - value
        elif value is None or other is None:
            delta[key] = None
    return {'baseline': baseline, 'intervened': changed, 'difference': delta,
            'intervention': deepcopy(intervention),
            'trace': {'targets_replaced': targets, 'contract_sha256': fingerprint(contract()),
                      'ontology_sha256': classified['ontology_sha256'],
                      'evidence_sha256': classified['evidence_sha256'],
                      'process_code_sha256': hashlib.sha256(Path(__file__).with_name('process.py').read_bytes()).hexdigest(),
                      'initial_state_shared': True, 'definitions_and_annotations_invariant': True,
                      'endpoint': 'conditional integrated H2 uptake, not genomic capacity or carbon removal',
                      'scope': 'Specified forcing regimes under assumed mechanisms. No do-calculus identification, endogenous-state clamps or transfer of Boolean-calculus soundness theorems.'}}
