# Physical interventions with invariant evidence

The implemented distinction follows the separation between invariant ontology /
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

Day indices start at zero. Replacements target water addition (mm), NaCl in
added water (mol/L), evaporation demand (mm/day), or temperature (°C).
`intervention.parameters` can specify post-disturbance `porosity`,
`psi_entry_mpa`, and `retention_b`. A porosity replacement alone keeps the
retention parameters fixed. This is a specified compaction scenario, not a
model of mechanical dose or soil deformation.

The API returns baseline and intervened trajectories, their summary differences,
and a trace naming the replaced inputs, shared initial state, and digests of
fixed evidence, ontology, contract and process code. Invalid targets, nonfinite
inputs and states exceeding pore storage are rejected. Uptake, relative abundance,
carbon removal and annotation membership cannot be clamped through this interface.
A saline regime may be physically accounted for while its gas output is unavailable.

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

The Twin does not import the full natural-deduction engine, graphical do-calculus
identification, or a Boolean soundness theorem. Its continuous structural equations
and numerical integration would require a separate formal semantics and verification
for such a claim. Current evidence establishes executable admission rules,
version-specific scenarios, and preserved definitions. Independent flat calculations
with the same information agree; predictive superiority from OWL representation is
not an outcome of this study.
