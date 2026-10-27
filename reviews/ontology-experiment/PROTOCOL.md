# Ontology experiment: genome-resolved composition of carbon-fixation queries

Protocol date: 2026-10-27. Status: **designed; evaluation not run**. A descriptive
[feasibility audit](design-audit.json) has been run on the existing v0.3 evidence.
These data have already been inspected, so this is a prospective analysis plan
for a retrospective benchmark, not a blinded biological validation or a public
preregistration. Commit this protocol before implementing the comparisons.

## Question and claim being tested

Does explicit reasoning about entity identity, pathway membership, overlapping
carrier classes and measurement types correctly connect metagenomic observations
to the digital twin's conditional carbon-fixation calculation?

The target benefit is **valid model composition and query answering**. The
decisive query is: which reference MAGs have both a screened hydrogenase and all
required CBB steps in that same MAG, and how much mapped-read weight belongs to
their union in a specified sample? This determines which MAG-specific resource
ceiling requests the current API admits. It also distinguishes an answer supported
by genome linkage from an interval supported by marginal observations alone.

The primary comparison is a generic rule-based implementation against an
independently written, equally informed flat implementation. Equal answers are
expected when both implement the same mathematics. A successful experiment would
establish that the ontology preserves and exposes the required semantics across
these queries. Predictive superiority or improved biological carbon fixation
would require additional outcomes and a different experiment.

## Existing evidence and feasibility

The frozen source is v0.3.0, model commit
`0657fe9bdf29d5ad1470a90e5f4c08b34692323f`, with exact input digests in
[design-audit.json](design-audit.json) and the
[source manifest](../../data/v03/manifest.json). The audit uses the bundled
MAG and sample evidence; independent evaluation must reconstruct its reference
answers from the underlying corrected KO, hydrogenase-classification and mapping
tables, rather than copying the bundled derived Boolean columns.

| Observed annotation stratum | MAGs |
| --- | ---: |
| Screened hydrogenase and all CBB steps | 23 |
| Screened hydrogenase, incomplete CBB screen | 1,071 |
| All CBB steps, no screened hydrogenase detected | 60 |
| Neither positive screen | 4,284 |

There are 47 MAGs with multiple screened hydrogenase categories and 995 with
RuBisCO detected but an incomplete full CBB screen. These are annotation states;
an incomplete screen permits several explanations, including genome incompleteness
and annotation failure. The stricter quality sensitivity retains the existing
completeness ≥90% and contamination <5% rule.

Of 551 samples, 550 have usable mapped-read measurements. F5Dr1 has unavailable
mapping and remains in the denominator as an unavailable query. In 543 of the
550 usable samples, hydrogenase and CBB marginal shares leave a nonzero interval
for their intersection. The median interval width is 0.02019 of mapped-read
weight. These counts establish that the proposed semantic contrasts have real
instances; they are not performance estimates.

For the first sample in lexical order, 10Dr2, hydrogenase-carrier share is
0.259474 and full-CBB share is 0.012289. Those two marginals permit an intersection
between 0 and 0.012289. The genome-linked observed intersection is 0.000635.
This example was chosen by identifier order, not by disagreement magnitude.

Reproduce the inventory from the repository root:

```bash
uv run python reproducibility/v03/ontology_design_audit.py
```

## Representation and complete experimental assumptions

The present release implements typed arithmetic and finite set membership.
The proposed experiment adds a small declarative rule evaluator that emits
derivation certificates. This evaluator is **not implemented by this protocol**.
It need not use OWL; any implementation must state its actual rule language and
semantics. The following local rule labels are protocol identifiers, not external
ontology identifiers. No new biological ontology IRIs are asserted here.

| Rule | Content and role | Existing contract |
| --- | --- | --- |
| E1 | A gene annotation, a reference MAG, a sample and an observation are distinct entity types. Exact reference identifiers link annotations to MAGs and mapping weights to that same catalogue. A MAG denotes a genome reconstruction with uncertain biological fidelity. | O1, M1 |
| E2 | All pathway-step witnesses for a MAG-level query have the same MAG subject. A conjunction of existential witnesses with a shared subject differs from a conjunction of community-level existential witnesses. | M1, M2; explicit logical formulation below |
| E3 | Gene-to-step mappings and the finite required-step list follow the frozen KEGG M00165 screen. Alternative KOs satisfy a step by disjunction; required steps combine by conjunction. The biological interpretation is detected annotation support. | M2 |
| E4 | A MAG carrying several screened hydrogenases occurs once in the carrier union. Hydrogenase and CBB carrier classes can intersect. Each carrier and its exact complement partition a fixed reference universe. | O2, O3 |
| E5 | An additive measure retains endpoint, unit, denominator, sample and source context. Mapped-read weights, gene-copy GE and absolute gas flux have distinct measurement maps. | O1–O4, M3, M4 |
| E6 | Unavailable source observations propagate to unavailable answers. An unsuccessful annotation screen returns insufficient evidence for eligibility; biological incapacity remains unasserted. Finite-table non-detection can be reported explicitly. | O5, M1, M2 |
| E7 | The eligibility conjunction controls admission to a conditional optimization with independently supplied flux and energy bounds. The simulation and resource optimization keep their existing equations and inputs. | C5, C11, C12 |

E1, E2, E4 and the compatible-measure part of E5 express identity, logical
structure and constitutive relations. E3 and the assay-specific part of E5 are
observation mappings supported by annotation conventions. E6 is an epistemic
policy. E7 is the boundary between evidence and a causal/resource model. A
successful screen permits a conditional calculation; enzyme expression, viable
biomass, turnover and measured activity require independent observations. The
ontology alone supplies none of those quantities.

Let G be the fixed MAG catalogue, K(g) its detected KOs, S the eleven required
steps, and A_s the alternative KOs for step s. Let H(g) mean that at least one
of the three screened hydrogenase categories is detected in g. Define

```text
C(g) = for every s in S, K(g) intersects A_s
E(g) = H(g) and C(g)
q(sample) = sum over g in G of w(sample,g) * indicator[E(g)]
           / sum over g in G of w(sample,g)
```

The annotation screen allows one KO to witness more than one step where the
frozen mapping explicitly permits it. K01602 remains optional in M00165.
The corresponding community-pooled expression, “for every step there exists
some MAG with a witness,” is a different query. No assertion excludes community
cross-feeding in general; a distributed pathway would require its own transport
and reaction model.

For marginal-only input, let h and c be carrier shares under the same measure.
The admitted intersection obeys

```text
max(0, h + c - 1) <= q <= min(h, c).
```

These are sharp bounds when only the two marginal measures and total are known.
With additional discrete MAG weights they are conservative outer bounds unless
an exact allocation problem is solved. Do not silently give one method hidden
MAG weights or linkage. Missing marginals yield unavailable bounds. Both rules
and flat baselines use the same tolerance, 1e-12 for shares.

The pathway definition was checked against the primary
[KEGG M00165 entry](https://www.kegg.jp/entry/M00165) on 2026-10-27; the release
also retains its [source audit](../v03/source-verification.json). It supports the
annotation screen. MAG integrity and the relation between detected genes and
physiology remain empirical assumptions requiring separate checks.

## Methods and information controls

Run all methods on the identical frozen query manifest. Keep the information
contrast separate from the representation contrast.

| Method | Available evidence | Computation and purpose |
| --- | --- | --- |
| Marginals with flat bounds | Sample H and C shares, common denominator and provenance | Exact marginal bounds; legitimate reduced-information baseline |
| Marginals with ontology rules | Identical marginal evidence | Same bounds with typed derivation; tests correct abstention from a point answer |
| Linked flat | Full per-MAG KOs, hydrogenase classifications, quality and mapping weights | Independently written relational joins, Boolean screen and set union; primary matched comparator |
| Linked ontology | Identical linked evidence | Generic typed rules derive step support, carrier intersection, union measure and admission certificate |
| Current v0.3 | Existing derived evidence tables and public API | Regression comparator; labels are precomputed, so report its preprocessing cost and provenance |

The gain from marginal to linked evidence is an **information gain**. Any
ontology-versus-flat difference at matched information is a representation,
implementation or computation effect that needs diagnosis. Both linked methods
receive the same mappings and quality filters; neither receives extra expression
or physiological data. There is no statistical model fitting or hyperparameter
selection in the primary experiment.

Run three ablations of the linked ontology implementation, keeping raw evidence
fixed: (i) remove the shared-MAG requirement, (ii) replace the carrier union with
addition of category totals, (iii) remove quantity/denominator checks. These are
deliberately incomplete implementations, not competitive baselines. They measure
which answers depend on each rule. Also permute gene-to-MAG links within source
assembly with seed 42, preserving gene labels and counts, as an information
corruption diagnostic. This corruption changes evidence, so it cannot establish
the benefit of an ontology representation. Leave unchanged cases in all reports.

## Query manifest and independent reference

Freeze four query families before running any method:

1. **Admission:** all 5,438 MAG identifiers, with two evidence conditions: full
   linked observations and marginal-only observations. Request the evidence state
   and permission to calculate a MAG-specific ceiling. A marginal-only sample
   cannot identify which named MAG passes. Report all four screen strata.
2. **Composition:** all 551 samples, requesting hydrogenase union share,
   H∩C share and its mapped-universe complement. Score point answers for linked
   input and intervals for marginal input. Keep the missing sample explicit.
3. **Measurement compatibility:** for every sample, attempt to substitute its
   observed total hydrogenase GE and its all-read carrier percentage for a
   mapped-read carrier share. Supply the actual endpoint/normalizer metadata.
   Include the valid all-read-to-mapped-read conversion when the mapped fraction
   is observed. These are constructed interface challenges using real values;
   report them separately from natural data queries. A conversion with an
   explicit measurement map must remain accepted.
4. **Causal interface:** for each linked MAG, pass the admission outcome to the
   unchanged resource-ceiling model, using the same fixed illustrative supplies:
   H₂ = O₂ = external CO₂ = 1e-9 mol m⁻² s⁻¹, CO = 0,
   capture fraction 0.3, maintenance 0 W m⁻², synthesis cost
   500,000 J mol-C⁻¹. Zero maintenance is a declared scenario assumption.
   An insufficient screen returns unavailable, not a biological zero. Accepted
   runs must have identical ceilings and balance certificates across methods.
   These uniform inputs test wiring; they do not estimate per-MAG field fluxes.

Construct reference answers independently from raw frozen tables with explicit
joins and separately transcribed KEGG alternatives. This evaluates fidelity to
the declared annotation contract. It is not an independent physiological gold
standard. Have two reviewers independently inspect all 23 co-carriers, all 60
CBB-only MAGs, all 47 multi-hydrogenase MAGs, and the first 30 identifiers in
each remaining screen stratum after SHA-256 sorting of `42|MAG_ID`; deduplicate
the review list. Include every missing-input case and every method disagreement.
Reviewers inspect source rows and definitions before seeing method names or
predictions; then record disagreements and their adjudication. If two independent
reviewers are unavailable, report the automated comparison separately and leave
the independent-review gate unmet. Agents reconstructing the same tables are not
a substitute for independent biological validation.

Because the current evidence and aggregate counts are already known, do not call
the catalogue held out. To claim generalization of a reusable rule engine later,
freeze this protocol and implementation, then have an independent curator supply
a new annotation release or schema with a separately signed reference map. That
extension requires a separate versioned protocol before outcomes are inspected.

## Outcomes, denominators and failure decisions

The primary outcome is **exact joint agreement of evidence state, admission and
query answer between linked ontology and independently reconstructed reference**,
reported beside linked-flat agreement on every query family. Boolean states must
match exactly; shares use absolute tolerance 1e-12. Compare solver fluxes with
absolute tolerance 1e-18 mol m⁻² s⁻¹ plus relative tolerance 1e-8, and independently
check balance/inequality residuals. Missing outputs, exceptions and invalid
certificates remain in the denominator as failures unless unavailability is the
reference outcome. An always-abstaining method fails all answerable queries.

Report false admissions, false refusals, unsupported point identification,
correct unavailability, point coverage, interval containment and interval width.
Reference answers for partial information include the correct bounds, so a
vacuous [0,1] interval is not counted as an exact answer. Report ablation changes
and certificate validity separately; require each claimed rule benefit to have
at least one independently checked nontrivial witness. A rule with no affected
natural query has no demonstrated natural-data contribution in this cohort.

This is a finite-catalogue logical/computational audit. Report exact counts and
errors rather than treating 551 repeated site samples or many queries per MAG
as independent experimental replicates. Do not claim a population significance
test from these dependent queries. Runtime is secondary: record hardware,
versions, preparation/compilation and query time separately, use one warm-up and
five measured repetitions, and report median time and peak memory. Make no
speed claim unless both implementations perform the same work.

Prespecified decisions:

- Any ontology-versus-flat discrepancy at equal information requires a source,
  rule or implementation diagnosis. A larger answer set alone is not a benefit.
- Any unsupported admission, incorrect missing-to-zero conversion or failed
  certificate fails the semantic-correctness claim, even if aggregate error is low.
- Equal correct ontology and flat answers support faithful semantic compilation.
  Report no incremental accuracy gain from the representation.
- Narrower answers after adding valid MAG linkage support the value of that
  measurement relation. Report the information source explicitly.
- If independent review finds wrong MAG assignments or annotation mappings,
  freeze the failed evaluation, revise the evidence in a new version and rerun
  all methods. Keep the original failure report.

## Deliverables, execution order and biological follow-up

The next implementation should produce a frozen query manifest, local rule
file with source links, independent reference answers, method outputs with
premise-level certificates, adjudication forms, per-query errors and a compact
comparison table. Tests should exercise the generic rule interpreter and the
resource interface independently. Keep code and full provenance beside results.

Execution order: (1) freeze protocol and source hashes; (2) implement independent
flat reference and generic rules; (3) freeze query manifest and rule versions;
(4) obtain blinded source review; (5) run every method and diagnostic once under
the frozen specification; (6) publish all disagreements and denominators. No
field campaign date is needed for this computational experiment.

A subsequent biological test would require linking hydrogen-supported carbon
assimilation to identified organisms under controlled conditions, for example
isotope-resolved assimilation with independently measured H₂ uptake and genome
identity. Whole-community uptake alone leaves organism-level co-carriage and
cross-feeding unresolved. Such a study needs an assay-specific protocol,
appropriate killed/background controls and replication based on pilot variance.
It is a separate proposed validation, with no claimed data or power here.

The manuscript currently reports the completed v0.3 analyses. This protocol may
be cited in Discussion as a planned test. Its proposed rule evaluator, review
outcomes and comparisons enter Methods/Results only after implementation and
evaluation, with all experimental assumptions above retained.
