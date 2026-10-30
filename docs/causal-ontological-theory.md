# Causal–ontological theory of EQ-BioDT (v0.6)

This note states how ontological deduction and causal intervention interact in
the twin. It is the reference for the code in `src/eq_biodt/v03/theory.py`,
`community.py`, `interventions.py` and `planning.py`, and for the Methods section
of the manuscript. The construction follows the stratified causal–ontological
models of Bochman and Lifschitz (2015), in which non-causal constraints restrict
the solutions of structural equations and `do()` replaces only structural
equations. Here it is specialized to a deterministic, time-discretized soil
process model.

## 1. Three kinds of change

A digital twin changes in three distinct ways, and the model keeps them apart.

| Change | What changes | Example | Interface |
|---|---|---|---|
| Observation update | evidence `E` | a new campaign, a new metagenome | evidence build, `freeze` |
| Model revision | ontology `K` or mechanism `F` | a new KEGG module rule, a new retention law | new release |
| Intervention | exogenous inputs `u` or designated mechanisms | add 12 mm water, compact the soil | `/v3/intervene`, `/v3/compare` |

An intervention never changes `E` or `K`. Choosing a community (a sample or
site visit) selects which evidence instantiates the twin; it conditions on
evidence and is not an intervention.

## 2. The stratified model

The twin is a tuple `M = (K, E, V_ex, V_const, V_caus, Φ, F)`:

* **Ontological zone `K`.** The OWL 2 TBox: KEGG step classes, the M00165
  conjunction, marker-guild hierarchy, carrier classes and endpoint classes.
* **Evidence `E`.** MAG KO annotations, marker guilds, MAG quality, mapped-read
  weights `w_ig` and genome equivalents.
* **Deduced memberships `D_K(E)`.** Reasoner entailments such as
  `g ∈ H2CarrierMAG`, `g ∈ COCarrierMAG`, `g ∈ CBBCompleteMAG`.
* **Exogenous variables `V_ex`.** Day-level forcing `I_t` (added water, mm),
  `c_t` (NaCl in added water, mol/L), `Ev_t` (evaporation demand) and `T_t`
  (temperature), plus soil parameters `φ, ψ_e, b`, layer depth `L`, drainage
  constants and kinetic constants.
* **Constituted variables `V_const` with maps `Φ`.** These are definitions, not
  mechanisms: `θ = W/(1000L)`, `W_max = 1000Lφ`, `ψ_m = ψ_e(θ/φ)^(-b)`,
  `ψ_o = -2ϕ(N/W)RT·10⁻³`, `ψ = ψ_m + ψ_o`, `D_g = D_0,g(φ-θ)^(10/3)/φ²`,
  carrier shares `q_s(C) = Σ_{g∈C} w_sg / Σ_g w_sg`, the community factor `γ_g`,
  chemical power, and the declared observation maps `H`.
* **Causal variables `V_caus` with mechanisms `F`.** Water storage `W_{t+1}`,
  salt inventory `N_{t+1}`, activation `A_t = A(ψ_t)`, the kinetic coefficient
  `k_{g,t} = k_{g,ref} γ_g A_t Q10^((T_t-T_ref)/10)` and uptake
  `J_{g,t} = C_a √(D_g k_g) tanh(L√(k_g/D_g))`, with proportional oxygen
  limitation.

The graph `G` has a directed edge from every argument of a mechanism or
constitutive map to its output. It is acyclic after unrolling over time.
`theory.graph()` exports the template for one time slice, with zone,
manipulability and the equation for every node.

### Coupling between deduction and mechanism

Deduction enters the mechanisms through `γ_g`, which is defined from the
reasoner's memberships, not from table flags:

* `reference` (default, no community selected): `γ_g = 1`. This reproduces the
  reference scenarios.
* `gate`: `γ_g = 1` if `q_s(C_g) > 0`, otherwise `0`, with
  `C_H2 = H2CarrierMAG` and `C_CO = COCarrierMAG`.
* `proportional`: `γ_g = q_s(C_g) / median_s' q_s'(C_g)`. This assumes that the
  kinetic coefficient scales with carrier share. The assumption is stated, not
  validated.

The same memberships admit reactions in the resource ceiling.

## 3. Interventions

An intervention `a = do(X_I = x_I)` is **admissible** only when both of the
following hold:

1. Every target is a manipulable exogenous input or a causal state with a
   physical handle (`W_t`, `N_t`). Targets in `K`, `E`, `D_K(E)` or `V_const`
   cannot be clamped, and neither can the outcome mechanisms `A`, `k` or `J`.
2. The regime is feasible: some state satisfies every constitutive map and every
   unsevered mechanism. For example, `do(φ = φ')` together with an initial water
   `W_0 > 1000Lφ'` is inadmissible because it violates `W ≤ W_max`.

Two kinds of target need special handling:

* **Constituted targets.** An intervention on a constituted variable is defined
  through its constituents.
  * `do(θ = θ*)` has the unique translation `do(W = 1000Lθ*)`, because the map
    is injective for fixed `L`.
  * `do(ψ = ψ*)` is an **ambiguous manipulation**: many `(W, N)` pairs give the
    same `ψ`, through either the matric or the osmotic component. The effect on
    `A` is version-independent, because `A` depends on `ψ` alone. The effect on
    `J` is version-dependent, because the path `W → θ → D_g → J` bypasses `ψ`.
    `theory.version_relevance(target, outcome)` decides this from the graph.
    A request must name a version (`matric` or `osmotic`), and the response
    reports the bypass paths.
* **Treatment versions.** "Add 12 mm over 14 days" is a superclass of the
  schedules that realize it (one pulse, four pulses, daily). Its effect is
  defined per version, or as a mixture under a stated version policy. The
  comparison reports the effect set across versions and parameter members, and
  flags version relevance when the ranking of versions changes.

## 4. Results the implementation enforces or computes

1. **Regime invariance of deduced quantities.** `D_K(E)`, `q_s(C)` and `γ_g`
   have no causal ancestors, so they are identical in every admissible regime.
   Every ontological consequence holds in every regime:
   `H2AutotrophicCandidateMAG ⊑ H2CarrierMAG` implies
   `q_s(H2Autotrophic) ≤ q_s(H2Carrier)`. The intervention trace recomputes
   these quantities in each arm and checks equality and the subsumption
   inequalities. It no longer records a constant `true`.
2. **Constitutive identities in every regime.** Interventions replace no map in
   `Φ`, so `θ = W/1000L`, `W ≤ W_max`, salt mass balance and the power identity
   hold along both trajectories. The trace checks these numerically.
3. **Deductive null effects.** If `γ_g(s) = 0`, then `J_g ≡ 0` in every
   admissible regime. Every interventional contrast on `J_g` is therefore zero,
   and the twin derives this from the closure without simulating. In the other
   direction, measured uptake of gas `g` in a community whose gate is closed
   would show that the evidence is incomplete: annotation or binning missed a
   carrier. It would not refute a definition in `K`. Absence of entailment is
   never read as biological absence.
4. **Paired regimes are unit-level counterfactuals.** The model is
   deterministic, and the two arms share the initial state, the non-targeted
   forcing, the parameters, `E` and `K`. So `Y_a − Y_a'` is a counterfactual
   contrast for the same specified unit. Composition holds: setting an input to
   the value it already takes leaves the trajectory unchanged (tested).

## 5. Identification and the role of external data

* **The twin's interventional quantities are model-derived.** They follow from
  `F` and `Φ`, not from the observational metagenomes, which do not identify
  any effect of water on uptake.
* **The compartment contrasts and the ridge forecast are associational.**
* **The Namib and Sonoran hydration microcosms are interventional data**, of
  the form `do(θ)` on bottled soil. They test the predicted interventional
  response directly.
* **The Sonoran field campaigns are observational.** Moisture and temperature
  co-vary across seasons, so `Q10` and the moisture response are not separately
  identified from the four campaigns. The boundary selection `Q10 = 1` must be
  read with that in mind.
* **Transfer to another desert is a transportability question** (Pearl and
  Bareinboim, 2014). Selection nodes point at the soil parameters `φ, ψ_e, b`
  and the activation thresholds. Calibration re-estimates exactly those
  mechanisms, and the external tests show which of them differ.

When a prediction fails, the stratification says where the fault is. `K` and
`Φ` are definitional, so a failed uptake prediction implicates a mechanism or an
assumed parameter in `F`, or an observation map. The failed Namib dry-activity
prediction implicates the activation threshold and the retention relation.

## 6. Identification with deductive closure

`src/eq_biodt/v03/identification.py` implements the interaction between the
do-calculus and deduction from the context-specific ontological d-separation
theorem. It is symbolic and deterministic: it checks graph conditions and
derives entailed values; it estimates no probability.

**Graphs.** `unrolled_graph(days, latent)` unrolls the slice of `theory.graph()`
over days (`W@t`, with `W@t -> W@t+1` and `N@t -> N@t+1`) and keeps the
within-slice edges exactly. An optional latent node, for the field setting a
season `S` with edges `S -> I@t, T@t, Ev@t`, is declared unobserved.
`d_separated` is d-separation on the moral graph of the ancestral set;
`d_connecting_path` returns a witness path.

**Closure.** `closure(observed, coupling)` returns `Cl_O`, the values entailed
by the invariant core (`K`, the release constants) and the mechanisms of the
non-admissible nodes. `check_admissible` refuses `do()` on `A`, `k`, `J`, `q`,
`gamma`, the constitutive nodes, `K`, `E` and the constants, so no admissible
regime severs their mechanisms. An entailment that uses only those mechanisms
and observed values therefore holds in every admissible regime, which is the
regime-universality the theorem needs. Manipulable nodes are never closed by
rule. The rules are (a) determination: a non-admissible node with all parents
observed or closed is closed (`determined`); (b) zero propagation: `q_g = 0`
gives `gamma_g = 0` under gate or proportional coupling, then `k_g = 0`, then
`J_g = 0` (the flux `C_a sqrt(D k) tanh(...)` vanishes at `k = 0`, and the
shared oxygen limit only scales it), then `P = P_cap = 0`; (c) `psi <= psi_d`
gives `A = 0` and the same chain. Under `reference` coupling `gamma = 1` and
`q = 0` closes nothing.

**Rules 1-3.** `rule1`, `rule2`, `rule3` test the usual conditions on the
mutilated graphs with the conditioning set `Z* = Z u Cl(observed)`. Each returns
`holds`, the plain `graph_test`, the test with `Z*`, the test with `Z*` after
removing the edges that the closure makes inert (`J_g` with `k_g = 0` is the
constant 0, so its inputs and the oxygen coupling edges `D_g, k_g -> J_other`
are pruned), whether the outcome is a deductive constant, `Z_star`,
`closed_by_deduction` and a reason. A closure node that is a collider or a
descendant of one opens paths when conditioned on; plain d-separation already
handles this.

**Queries** (`identify`; `POST /v3/identify`, `eq-biodt identify --query`;
`GET /v3/closure?q_co=0&coupling=gate`, `POST /v3/closure`):

* `deductive_null` (Q1). With `q_CO = 0`, `P(J_CO | do(I@0), q_CO = 0) =
  P(J_CO | q_CO = 0)` by Rule 3. The plain test fails: `I@0 -> W -> theta -> D
  -> J_CO` is open. The closure fixes `gamma_CO = k_CO = J_CO = 0`; a deductive
  constant is independent of every intervention, so Rule 3 holds. The same
  holds for every admissible time-0 target (`I, c, Ev, T, W, N, phi, psi_e,
  b`) and under gate and proportional coupling. Under reference coupling it does
  not hold.
* `oxygen_context` (Q2). `J_H2` depends on the CO channel (`D_CO`, `k_CO`)
  only through the shared oxygen demand. In the stratum `q_CO = 0`, `k_CO` is
  closed and `J_CO = 0` for every `D_CO`, so the probe `do(D_CO)` leaves
  `J_H2` unchanged: Rule 3 fails with `Z*` alone (the direct edge remains) and
  holds after the closure prunes the inert edges. Without the stratum it fails.
  `D_CO` is not an admissible regime; it is a formal probe of the CO channel.
  Water still changes `J_H2`. A check with `process.gas_exchange` and a binding
  oxygen limit confirms it numerically.
* `field_backdoor` (Q3). In observational field data (`I, T, Ev, W, J_H2`
  measured, `phi` a measured site property, season `S` latent), `P(J_H2@t |
  do(W@t))` is identified by back-door. For `days = 2` the minimal sets are
  `{I@0, I@1, T@1, W@0, phi}` and `{Ev@0, Ev@1, I@0, I@1, W@0, phi}`. Without
  `T` the second set still works; there is no back-door set if `T` and `Ev`
  are both unmeasured, or if `phi` is unmeasured. The effect of `T` needs the
  moisture state or its drivers (`{Ev@0, Ev@1, I@0, I@1}` or `{I@0, I@1, W@0,
  W@1, phi}`). The identification is graphical. The Sonoran campaigns have
  near-collinear moisture and temperature, so the strata of these sets are
  nearly empty; the positivity violation, not the graph, is why the `Q10`
  question stays unresolved.

Reproduce with `reproducibility/v06/identification.py`, which writes
`reviews/v06/identification.json` and `paper/tables/v06_identification.tex`.

**Limits.** The reasoning is deterministic and symbolic. The latent season is
an assumed structure, not an estimate. Q3 uses back-door search only; failing to
find a set is not a proof of non-identifiability. The value rules cover the
zero chain and the dry-threshold chain, not every entailment of the continuous
equations.

## 7. Scope

The do-calculus identification of section 6 is implemented for the three queries above. The calculus here is the deterministic, single-unit case: point-mass exogenous
inputs, with paired regimes treated as counterfactuals. Probabilistic `do()`
over uncertain parameters is represented as a sensitivity set, not as a
probability distribution. The Boolean natural-deduction results of the
causal–ontological calculus motivate the two-zone design, but the twin does not
claim them as theorems for its continuous equations. The checks above are
executable assertions about this implementation.
