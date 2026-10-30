# Physical interventions with invariant evidence

The formal statement (three kinds of change, the stratified model, admissibility,
treatment versions, deductive null effects, identification) is in
[causal-ontological-theory.md](causal-ontological-theory.md). This note documents
the endpoints. The implemented distinction follows the separation between invariant ontology /
constitutive knowledge and replaceable physical assignments in causal-ontological
models. In a regime, definitions retain their meaning while designated physical
inputs change. The Twin evaluates a deterministic, assumed process model; it does
not identify intervention effects from observational metagenomes.

The paired endpoint is `POST /v3/intervene`, or
`eq-biodt intervene --input data/v03/example_intervention.json`.

```json
{
  "parameters": {"k_h2_s": 0.001, "k_co_s": 0.001},
  "initial": {"water_mm": 1.8},
  "forcing": [
    {"water_mm": 0, "evaporation_mm": 1, "temp_c": 25},
    {"water_mm": 0, "evaporation_mm": 1, "temp_c": 25}
  ],
  "intervention": {"forcing": {"0": {"water_mm": 12}}}
}
```

Day indices start at zero. Forcing replacements target water addition (mm), NaCl in
added water (mol/L), evaporation demand (mm/day) or temperature (°C).
`intervention.parameters` can specify post-disturbance `porosity`,
`psi_entry_mpa` and `retention_b`. A porosity replacement alone keeps the
retention parameters fixed. This is a specified compaction scenario, not a
model of mechanical dose or soil deformation.

State interventions (`state`, keyed by day) set the causal states directly:
`water_mm`, `nacl_mol_m2`, or the constituted `theta` (translated uniquely to
`water_mm`) or `psi_mpa` (ambiguous; requires `"version": "matric"` or
`"osmotic"`, naming the handle `W` or `N` that realizes it). The response lists the
resolved values and, for `psi_mpa`, the version relevance for activation and each
gas. `POST /v3/version-relevance` (`{"target": "psi", "outcome": "h2_uptake_mol_m2"}`)
answers the same question from the graph without simulating; `GET /v3/graph` returns
the graph.

A community and a coupling mode (`reference`, `gate`, `proportional`) may be given in
the request body or inside `intervention`. They select evidence and condition both
arms identically; they are not interventions. `GET /v3/community/{sample}` returns
the reasoner-entailed carrier shares and the coupling factors for each mode.
`POST /v3/compare` accepts the same `community` and `coupling` and adds a
treatment-version block: wetting schedules with equal total water are versions of
one superclass intervention, reported per parameter member with their rankings.

The API returns baseline and intervened trajectories, their summary differences, gas
uptake contrasts, deduced null effects, and a trace naming the replaced inputs,
shared initial state and digests of fixed evidence, ontology, contract, graph and
process code. The invariance checks (memberships and shares equal across arms,
subsumption inequalities, constitutive identities, water and salt balances, power
identity) are recomputed on every call; a failed check raises an error.

Invalid requests are rejected with an error that names the zone of the target:
ontology, evidence or entailed memberships (not manipulable), constituted variables
(translate through constituents), outcome mechanisms (`A`, `k`, `J`), nonfinite
values, and regimes that violate a constitutive constraint, such as stored water
above pore storage. Uptake, relative abundance, carbon removal and annotation
membership cannot be set through this interface. A saline regime may be physically
accounted for while its gas output is unavailable. If the selected community has no
entailed carrier for a gas, its uptake is zero in every regime and is reported as
`zero_by_deduction`; this reflects the annotation evidence, not biological absence.

The resource ceiling is a distinct coupling of deduction and numerical calculation.
Only verified H2-plus-CBB membership admits the synthesis calculation. CO oxidation
additionally requires CO-carrier membership. Unentailed reactions are excluded
under an evidence-admission policy; absence of entailment is not asserted biological
absence. Supplied substrate bounds and kinetics remain independent assumptions.

## What this contributes

The two-zone regime interpretation in the causal-deductive work provides a useful
design principle: carry invariant definitions into every scenario and replace only
the designated physical assignments. It explains why changing an aggregate endpoint
or a pathway label is not a physical intervention. Explicit treatment versions also
explain why “add 12 mm water” is incomplete without a schedule: single and fractionated
pulses reverse their uptake ranking across the evaluated parameter grid.

The Twin does not claim the Boolean natural-deduction soundness results or graphical
do-calculus identification for its continuous structural equations; these would
require a separate formal semantics. The executable checks above are assertions about
this implementation, and interventional quantities are model-derived rather than
identified from observational metagenomes (theory note, section 5). Independent flat
calculations with the same information agree; predictive superiority from OWL
representation is not an outcome of this study.
