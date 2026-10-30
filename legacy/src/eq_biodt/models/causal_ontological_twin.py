"""Two-Layer Causal-Ontological Engine for the Empty Quarter Biodiversity Digital Twin.

Directly implements the mathematical framework from:
- `~/Documents/papers/causal-knowledge-graphs` (Two-Layer Causal-Ontological Model
  M = (T_ont, V_caus, V_const, Phi_const, M_caus, P_U), Ontological Polytope C(T_ont),
  Ontological Consistency Score OCS, Version-Ambiguity Bounds [L(c), U(c)], and
  Agreement-Weighted Class Transport Proposition 2).
- `~/Documents/papers/onto-causal` (Phylogenetic deconfounding of microbial trait-habitat
  associations via phyloglm attenuation weights).
- `CKG_USEFULNESS.md` (strict separation of constitutive aggregation maps Phi_const from
  causal structural equations F_caus, and separation of encoded genomic capacity from
  realized thermodynamic activity/flux).
- Automated OWL 2 DL semantic compilation (`rubalkhali_twin_constraints.owl` -> linear
  equality/inequality polyhedron and summing matrix S) and Shrinkage-MinT / Nonnegative-MinT
  hierarchical reconciliation (Wickramasuriya et al. 2019; Ledoit & Wolf 2004).
"""

from __future__ import annotations

import xml.etree.ElementTree as ET
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
from scipy.optimize import minimize


@dataclass(frozen=True)
class SubsumptionRelation:
    """Axiom `child_var <= parent_var` induced by `ChildClass rdfs:subClassOf ParentClass`."""

    child_var: str
    parent_var: str
    ontology_axiom: str


@dataclass(frozen=True)
class PartitionSumConstraint:
    """Axiom `sum(children) <= parent` (or `== parent` when exhaustive)."""

    children_vars: Tuple[str, ...]
    parent_var: str
    exhaustive: bool = False
    ontology_axiom: str = ""


class OntologicalPolytope:
    """Convex polytope C(T_ont) = {y in R^p : A_ont y <= b_ont} induced by OWL 2 DL TBox axioms.

    Enforces:
    1. Completeness: every declared state variable must be present and finite.
    2. Non-negativity of abundances, genome equivalents, and fluxes.
    3. Subsumption monotonicity (`y_child <= y_parent`).
    4. Disjointness / partition bounds (`sum_{c in children} y_c <= y_parent` or `== y_parent`).
    5. Upper simplex bounds (`y_parent <= 1.0` for relative abundances).
    """

    DEFAULT_OWL_PATH = (
        Path(__file__).resolve().parents[3]
        / "data"
        / "ontology"
        / "rubalkhali_twin_constraints.owl"
    )

    def __init__(
        self,
        variable_names: List[str],
        subsumptions: Optional[List[SubsumptionRelation]] = None,
        partitions: Optional[List[PartitionSumConstraint]] = None,
        upper_bounds: Optional[Dict[str, float]] = None,
    ) -> None:
        self.variable_names = list(variable_names)
        self.var_idx = {name: i for i, name in enumerate(self.variable_names)}
        self.subsumptions = subsumptions or []
        self.partitions = partitions or []
        self.upper_bounds = upper_bounds or {}
        self.last_projection_solver_success: bool = True
        self.last_projection_used_fallback: bool = False

    @classmethod
    def from_owl_file(cls, owl_path: Optional[Path] = None) -> "OntologicalPolytope":
        """Compile an OWL 2 DL XML ontology specification directly into an OntologicalPolytope."""
        path = owl_path or cls.DEFAULT_OWL_PATH
        if not path.exists():
            raise FileNotFoundError(f"Missing OWL constraint specification: {path}")

        tree = ET.parse(path)
        root = tree.getroot()
        ns = {"owl": "http://www.w3.org/2002/07/owl#"}

        iri_to_var: Dict[str, str] = {}
        upper_bounds: Dict[str, float] = {}
        variables: List[str] = []

        # 1. Extract AnnotationAssertion bindings (hasStateVariable and upperSimplexBound)
        raw_bounds_by_iri: Dict[str, float] = {}
        for ann in root.findall("owl:AnnotationAssertion", ns):
            prop_el = ann.find("owl:AnnotationProperty", ns)
            iri_el = ann.find("owl:IRI", ns)
            lit_el = ann.find("owl:Literal", ns)
            if prop_el is None or iri_el is None or lit_el is None:
                continue
            prop_iri = prop_el.attrib.get("IRI", "")
            subj_iri = (iri_el.text or "").strip()
            lit_val = (lit_el.text or "").strip()
            if prop_iri.endswith("hasStateVariable") and subj_iri and lit_val:
                iri_to_var[subj_iri] = lit_val
                if lit_val not in variables:
                    variables.append(lit_val)
            elif prop_iri.endswith("upperSimplexBound") and subj_iri and lit_val:
                raw_bounds_by_iri[subj_iri] = float(lit_val)

        for subj_iri, bound_val in raw_bounds_by_iri.items():
            if subj_iri in iri_to_var:
                upper_bounds[iri_to_var[subj_iri]] = bound_val

        # 2. Extract SubClassOf axioms
        subsumptions: List[SubsumptionRelation] = []
        for sub_el in root.findall("owl:SubClassOf", ns):
            classes = sub_el.findall("owl:Class", ns)
            if len(classes) == 2:
                c_iri = classes[0].attrib.get("IRI", "")
                p_iri = classes[1].attrib.get("IRI", "")
                if c_iri in iri_to_var and p_iri in iri_to_var:
                    subsumptions.append(
                        SubsumptionRelation(
                            child_var=iri_to_var[c_iri],
                            parent_var=iri_to_var[p_iri],
                            ontology_axiom=f"{c_iri} SubClassOf {p_iri}",
                        )
                    )

        # 3. Extract DisjointUnion axioms (exhaustive partitions) and DisjointClasses
        partitions: List[PartitionSumConstraint] = []
        for du_el in root.findall("owl:DisjointUnion", ns):
            classes = du_el.findall("owl:Class", ns)
            if len(classes) >= 2:
                p_iri = classes[0].attrib.get("IRI", "")
                c_iris = [c.attrib.get("IRI", "") for c in classes[1:]]
                if p_iri in iri_to_var and all(c in iri_to_var for c in c_iris):
                    c_vars = tuple(iri_to_var[c] for c in c_iris)
                    partitions.append(
                        PartitionSumConstraint(
                            children_vars=c_vars,
                            parent_var=iri_to_var[p_iri],
                            exhaustive=True,
                            ontology_axiom=f"DisjointUnion({p_iri} == {', '.join(c_iris)})",
                        )
                    )

        return cls(
            variable_names=variables,
            subsumptions=subsumptions,
            partitions=partitions,
            upper_bounds=upper_bounds,
        )

    def check_violations(self, state: Dict[str, float], tol: float = 1e-6) -> List[str]:
        """Return a list of violated TBox axioms for a state dictionary.

        Fails closed when required coordinates are missing or non-finite.
        """
        if not state:
            return ["EmptyState: no state variables provided"]

        violations: List[str] = []
        for name in self.variable_names:
            if name not in state:
                violations.append(f"MissingCoordinate: {name}")
                continue
            val = float(state[name])
            if not np.isfinite(val):
                violations.append(f"NonFinite: {name}={val}")
                continue
            if val < -tol:
                violations.append(f"NonNegativity: {name}={val:.6f} < 0")
            if name in self.upper_bounds and val > self.upper_bounds[name] + tol:
                violations.append(
                    f"UpperBound: {name}={val:.6f} > {self.upper_bounds[name]:.6f}"
                )
        if any(
            v.startswith("MissingCoordinate:") or v.startswith("NonFinite:")
            for v in violations
        ):
            return violations

        for sub in self.subsumptions:
            c_val = float(state[sub.child_var])
            p_val = float(state[sub.parent_var])
            if c_val > p_val + tol:
                violations.append(
                    f"Subsumption({sub.ontology_axiom}): {sub.child_var}={c_val:.6f} > {sub.parent_var}={p_val:.6f}"
                )

        for part in self.partitions:
            c_sum = sum(float(state[c]) for c in part.children_vars)
            p_val = float(state[part.parent_var])
            if part.exhaustive:
                if abs(c_sum - p_val) > tol:
                    violations.append(
                        f"ExhaustivePartition({part.ontology_axiom}): sum({part.children_vars})={c_sum:.6f} != {part.parent_var}={p_val:.6f}"
                    )
            else:
                if c_sum > p_val + tol:
                    violations.append(
                        f"DisjointPartition({part.ontology_axiom}): sum({part.children_vars})={c_sum:.6f} > {part.parent_var}={p_val:.6f}"
                    )
        return violations

    def compute_ocs(self, states: List[Dict[str, float]], tol: float = 1e-5) -> float:
        """Compute Ontological Consistency Score (OCS in [0, 1]) over a list of states."""
        if not states:
            return 0.0
        valid = sum(1 for s in states if len(self.check_violations(s, tol=tol)) == 0)
        return float(valid) / float(len(states))

    def project(
        self,
        raw_state: Dict[str, float],
        weight_matrix_inv: Optional[np.ndarray] = None,
    ) -> Dict[str, float]:
        """Oblique (W^{-1}) or Euclidean (I) projection Pi_{C(T_ont)} onto the ontological polyhedron."""
        if not raw_state:
            raise ValueError("Cannot project an empty state dictionary.")
        for k, raw_val in raw_state.items():
            if isinstance(raw_val, (int, float, np.number)) and not np.isfinite(float(raw_val)):
                raise ValueError(f"Cannot project state containing non-finite value for '{k}': {raw_val}")
        missing = set(self.variable_names) - raw_state.keys()
        if missing:
            raise ValueError(f"Cannot project missing coordinates: {sorted(missing)}")
        y0 = np.array([float(raw_state[k]) for k in self.variable_names], dtype=float)
        if not np.all(np.isfinite(y0)):
            raise ValueError("Cannot project state vector containing NaN or infinite values.")
        bounds = [
            (0.0, self.upper_bounds.get(k, None)) for k in self.variable_names
        ]

        constraints: List[Dict[str, Any]] = []
        for sub in self.subsumptions:
            ci = self.var_idx[sub.child_var]
            pi = self.var_idx[sub.parent_var]
            constraints.append(
                {
                    "type": "ineq",
                    "fun": lambda y, ci=ci, pi=pi: y[pi] - y[ci],
                }
            )
        for part in self.partitions:
            cis = [self.var_idx[c] for c in part.children_vars]
            pi = self.var_idx[part.parent_var]
            if part.exhaustive:
                constraints.append(
                    {
                        "type": "eq",
                        "fun": lambda y, cis=cis, pi=pi: y[pi] - sum(y[i] for i in cis),
                    }
                )
            else:
                constraints.append(
                    {
                        "type": "ineq",
                        "fun": lambda y, cis=cis, pi=pi: y[pi] - sum(y[i] for i in cis),
                    }
                )

        if weight_matrix_inv is not None:
            w_inv = np.asarray(weight_matrix_inv, dtype=float)
            if (w_inv.shape != (len(y0), len(y0)) or not np.isfinite(w_inv).all()
                    or not np.allclose(w_inv, w_inv.T)):
                raise ValueError("Projection metric must be finite, symmetric and correctly sized")
            try:
                np.linalg.cholesky(w_inv)
            except np.linalg.LinAlgError as exc:
                raise ValueError("Projection metric must be positive definite") from exc
            obj_fun = lambda y: 0.5 * float((y - y0) @ w_inv @ (y - y0))
            jac_fun = lambda y: w_inv @ (y - y0)
        else:
            obj_fun = lambda y: 0.5 * float(np.sum((y - y0) ** 2))
            jac_fun = lambda y: y - y0

        res = minimize(
            fun=obj_fun,
            x0=np.maximum(y0, 0.0),
            jac=jac_fun,
            bounds=bounds,
            constraints=constraints,
            method="SLSQP",
            options={"ftol": 1e-10, "maxiter": 300},
        )
        self.last_projection_solver_success = bool(res.success)
        if not res.success:
            raise RuntimeError(f"Projection solver failed: {res.message}")
        projected = res.x
        out = dict(raw_state)
        for i, k in enumerate(self.variable_names):
            ub = self.upper_bounds.get(k, None)
            val = float(max(0.0, projected[i]))
            out[k] = min(val, ub) if ub is not None else val

        # A feasibility repair is not an exact metric projection. Fail explicitly.
        if len(self.check_violations(out, tol=1e-5)) > 0:
            raise RuntimeError("Projection solver returned an infeasible state")
        self.last_projection_used_fallback = False
        return out


class MinTReconciler:
    """Optimal hierarchical reconciliation via Shrinkage-MinT and Nonnegative-MinT.

    Follows Wickramasuriya, Athanasopoulos & Hyndman (2019, JASA) and Ledoit & Wolf (2004):
    - Given base forecasts `Y_hat` (N x m) and summing matrix `S` (m x k_leaf) such that
      coherent states satisfy `y = S b`:
      1. Estimate error covariance `W_hat` strictly from inner training-fold base residuals
         `E_train = Y_hat_train - Y_train` using analytical Ledoit-Wolf shrinkage toward the
         diagonal variance target `diag(S_emp)`.
      2. Compute unconstrained MinT reconciled predictions:
         `Y_tilde = Y_hat @ W_inv @ S @ (S^T W_inv S)^{-1} @ S^T`
      3. For Nonnegative-MinT (`nonnegative=True`), solve for each row with any negative
         leaf coordinate:
         `min_{b >= 0} 0.5 * (S b - y_hat)^T W_inv (S b - y_hat)`, returning `y_tilde = S b`.
    """

    @staticmethod
    def estimate_ledoit_wolf_covariance(
        residuals: np.ndarray, ridge_floor: float = 1e-6
    ) -> Tuple[np.ndarray, float]:
        """Compute analytical Ledoit-Wolf shrinkage covariance toward diagonal variance target."""
        E = np.asarray(residuals, dtype=float)
        n, p = E.shape
        if n <= 1:
            return np.eye(p, dtype=float), 1.0
        E_centered = E - E.mean(axis=0, keepdims=True)
        S_emp = (E_centered.T @ E_centered) / float(n)
        diag_target = np.diag(np.diag(S_emp))

        # Ledoit-Wolf optimal shrinkage intensity toward diagonal variance matrix
        E2 = E_centered**2
        var_s_ij = (E2.T @ E2) / float(n) - (S_emp**2)
        off_diag_mask = ~np.eye(p, dtype=bool)
        num = float(np.sum(var_s_ij[off_diag_mask]))
        denom = float(np.sum((S_emp[off_diag_mask]) ** 2)) * float(n)
        shrinkage_lambda = float(np.clip(num / max(1e-12, denom), 0.05, 0.95))

        W_shrunk = (1.0 - shrinkage_lambda) * S_emp + shrinkage_lambda * diag_target
        W_shrunk += ridge_floor * np.eye(p, dtype=float)
        return W_shrunk, shrinkage_lambda

    @classmethod
    def reconcile_matrix(
        cls,
        y_hat: np.ndarray,
        summing_matrix_s: np.ndarray,
        w_cov: np.ndarray,
        nonnegative: bool = True,
    ) -> np.ndarray:
        """Reconcile an (N x m) base prediction matrix under covariance `w_cov` and summing matrix `S`."""
        Y = np.asarray(y_hat, dtype=float)
        S = np.asarray(summing_matrix_s, dtype=float)
        W = np.asarray(w_cov, dtype=float)
        m, k = S.shape

        W_inv = np.linalg.pinv(W, rcond=1e-10)
        middle = np.linalg.pinv(S.T @ W_inv @ S, rcond=1e-10)
        P_mint = middle @ (S.T @ W_inv)  # k x m mapping from all nodes to bottom leaves

        B_tilde = (P_mint @ Y.T).T  # N x k leaf predictions
        if not nonnegative:
            return B_tilde @ S.T

        # Nonnegative-MinT: solve strictly convex QP in leaf space for rows with any b_i < 0
        Q = S.T @ W_inv @ S
        bounds = [(0.0, None) for _ in range(k)]
        B_nn = B_tilde.copy()
        for idx in range(Y.shape[0]):
            if np.all(B_nn[idx] >= 0.0):
                continue
            y_row = Y[idx]
            lin_vec = S.T @ (W_inv @ y_row)
            b0 = np.maximum(B_nn[idx], 0.0)
            res = minimize(
                fun=lambda b: 0.5 * float(b @ Q @ b) - float(lin_vec @ b),
                x0=b0,
                jac=lambda b: Q @ b - lin_vec,
                bounds=bounds,
                method="L-BFGS-B",
            )
            B_nn[idx] = np.maximum(res.x if res.success else b0, 0.0)
        return B_nn @ S.T


class CausalOntologicalTwinEngine:
    """Two-Layer Causal-Ontological Engine M = (T_ont, V_caus, V_const, Phi_const, M_caus, P_U).

    Separates:
    - Layer 1 (Causal Structural Equations `evaluate_causal_layer`):
      Computes leaf guild genome equivalents, phylum relative abundances, thermodynamic
      trace-gas oxidation power (uW m^-3), dark-proteome coding fraction, and resistome risk
      from abiotic/biotic intervention drivers.
    - Layer 2 (Constitutive Ontological Map `evaluate_constitutive_layer`):
      Deterministically computes super-class aggregates (`total_h2ase_ge`,
      `total_trace_gas_scavenger_ge`, `total_prokaryotic_rel_abund`,
      `shannon_16s`, `hill_q1_diversity`, `viable_dna_fraction`) via Phi_const and
      projects onto the OWL-compiled polytope C(T_ont) so that OCS == 1.0 identically.
    """

    # Phylogenetic deconfounding weights from ~/Documents/papers/onto-causal:
    PHYLOGLM_DECONFOUNDING_WEIGHTS: Dict[str, float] = {
        "group_1l_nife_aridity_selection": 0.94,      # Retained under phyloglm across 4 phyla
        "form_I_coxL_aridity_selection": 0.89,        # Retained under phyloglm
        "thermoplasmatota_subsurface_refugium": 0.91, # Retained across Nitrososphaerota/Thermoplasmatota
        "halobacteriota_evaporite_selection": 0.96,   # Strong physiological osmoregulation law
        "rbcL_hyaB_plasmid_rhizosphere_rescue": 0.88, # Horizontal plasmid transfer decouples from host tree
        "naive_gc_content_aridity": 0.18,             # Heavily confounded by Actinomycetota clade share
    }

    def __init__(self, owl_path: Optional[Path] = None) -> None:
        target_owl = owl_path or OntologicalPolytope.DEFAULT_OWL_PATH
        if target_owl.exists():
            self.polytope = OntologicalPolytope.from_owl_file(target_owl)
        else:
            self.polytope = self._build_default_polytope()

    @staticmethod
    def _build_default_polytope() -> OntologicalPolytope:
        return OntologicalPolytope.from_owl_file(OntologicalPolytope.DEFAULT_OWL_PATH)

    def evaluate_constitutive_layer(self, leaf_state: Dict[str, float]) -> Dict[str, float]:
        """Apply constitutive aggregation map Phi_const and project onto C(T_ont)."""
        if not leaf_state:
            raise ValueError("Cannot evaluate constitutive layer on an empty state dictionary.")
        for k, raw_val in leaf_state.items():
            if isinstance(raw_val, (int, float, np.number)) and not np.isfinite(float(raw_val)):
                raise ValueError(f"Non-finite value in leaf_state['{k}']: {raw_val}")

        state = dict(leaf_state)

        def _nonneg_finite(key: str, default: float) -> float:
            raw_v = state.get(key, default)
            val = float(default if raw_v is None else raw_v)
            if not np.isfinite(val):
                raise ValueError(f"Non-finite value for state coordinate '{key}': {raw_v}")
            return max(0.0, val)

        # Normalize phylum partition on the unit simplex
        act = _nonneg_finite("actinomycetota_rel_abund", 0.45)
        therm = _nonneg_finite("thermoplasmatota_rel_abund", 0.04)
        halo = _nonneg_finite("halobacteriota_rel_abund", 0.02)
        known_sum = act + therm + halo
        if known_sum > 0.96:
            scale = 0.96 / known_sum
            act *= scale
            therm *= scale
            halo *= scale
            known_sum = 0.96
        other = max(0.0, 1.0 - known_sum)

        state["actinomycetota_rel_abund"] = act
        state["thermoplasmatota_rel_abund"] = therm
        state["halobacteriota_rel_abund"] = halo
        state["other_prokaryotes_rel_abund"] = other
        state["total_prokaryotic_rel_abund"] = act + therm + halo + other

        # Constitutive sum of [NiFe]-hydrogenase & trace-gas scavenger genome equivalents
        h1l = _nonneg_finite("hydB_1l_ge", 0.61)
        h1h = _nonneg_finite("hydB_1h_ge", 0.12)
        h2a = _nonneg_finite("hydB_2a_ge", 0.03)
        cox = _nonneg_finite("coxL_formI_ge", 1.15)

        state["hydB_1l_ge"] = h1l
        state["hydB_1h_ge"] = h1h
        state["hydB_2a_ge"] = h2a
        state["total_h2ase_ge"] = h1l + h1h + h2a
        state["coxL_formI_ge"] = cox
        state["total_trace_gas_scavenger_ge"] = state["total_h2ase_ge"] + cox

        # Constitutive sum of realized thermodynamic trace-gas oxidation power (uW m^-3)
        p_h2 = _nonneg_finite("h2_oxidation_power_uw_m3", 18.0)
        p_co = _nonneg_finite("co_oxidation_power_uw_m3", 2.2)
        state["h2_oxidation_power_uw_m3"] = p_h2
        state["co_oxidation_power_uw_m3"] = p_co
        state["total_trace_gas_power_uw_m3"] = p_h2 + p_co

        # Constitutive sum of PMAxx viable vs relic DNA fractions
        viable_raw = _nonneg_finite("viable_dna_fraction", 0.56)
        viable = float(np.clip(viable_raw, 0.01, 0.99))
        state["viable_dna_fraction"] = viable
        state["relic_dna_fraction"] = 1.0 - viable
        state["total_dna_fraction"] = 1.0

        # Project onto ontological polytope C(T_ont) to guarantee OCS == 1.0
        projected = self.polytope.project(state)

        # Derived biodiversity indices (Shannon & Hill q1)
        shannon = _nonneg_finite("shannon_16s", float(projected.get("shannon_16s", 4.8)))
        projected["shannon_16s"] = shannon
        projected["hill_q1_diversity"] = float(np.exp(shannon))
        projected["ocs_score"] = self.polytope.compute_ocs([projected])
        projected["projection_solver_success"] = self.polytope.last_projection_solver_success
        projected["projection_used_fallback"] = self.polytope.last_projection_used_fallback
        return projected

    @staticmethod
    def compute_version_ambiguity_interval(
        target_parent_total: float,
        leaf_sensitivities: Dict[str, float],
        min_leaf_shares: Optional[Dict[str, float]] = None,
    ) -> Tuple[float, float]:
        """Compute exact version-ambiguity interval [L_k(c), U_k(c)] for an ontology-class intervention."""
        if target_parent_total < 0.0:
            raise ValueError("target_parent_total must be non-negative.")
        if not leaf_sensitivities or target_parent_total == 0.0:
            return (0.0, 0.0)
        keys = list(leaf_sensitivities.keys())
        betas = np.array([float(leaf_sensitivities[k]) for k in keys], dtype=float)
        min_shares = min_leaf_shares or {}
        base_alloc = np.array([float(min_shares.get(k, 0.0)) for k in keys], dtype=float)
        if np.any(base_alloc < 0.0) or float(base_alloc.sum()) > 1.0 + 1e-9:
            raise ValueError("min_leaf_shares must be non-negative and sum to at most 1.0.")
        rem_share = max(0.0, 1.0 - float(base_alloc.sum()))

        base_effect = float(np.dot(base_alloc, betas)) * target_parent_total
        low = base_effect + rem_share * float(betas.min()) * target_parent_total
        high = base_effect + rem_share * float(betas.max()) * target_parent_total
        return (min(low, high), max(low, high))

    @staticmethod
    def agreement_weighted_class_transport(
        sibling_estimates: Dict[str, Tuple[float, float]],
        prior_tau: float = 0.25,
        baseline_abundances: Optional[Dict[str, float]] = None,
    ) -> Dict[str, Any]:
        """Proposition 2 (Sibling Class Effect Pooling & Abundance-Weighted Aggregate Transport)."""
        if not sibling_estimates:
            return {
                "pooled_effect": 0.0,
                "pooled_se": 0.0,
                "dersimonian_laird_tau2": 0.0,
                "cochran_q": 0.0,
                "sign_concordance": 1.0,
                "abundance_weighted_parent_log2fc": 0.0,
                "weights": {},
            }

        names = list(sibling_estimates.keys())
        betas = np.array([sibling_estimates[k][0] for k in names], dtype=float)
        ses = np.array([max(1e-6, sibling_estimates[k][1]) for k in names], dtype=float)

        pos_frac = float(np.mean(betas >= 0.0))
        sign_concordance = 2.0 * abs(pos_frac - 0.5)

        w_fe = 1.0 / (ses**2)
        beta_fe = float(np.sum(w_fe * betas) / np.sum(w_fe))
        cochran_q = float(np.sum(w_fe * ((betas - beta_fe) ** 2)))
        df = max(1, len(names) - 1)
        c_denom = float(np.sum(w_fe) - np.sum(w_fe**2) / np.sum(w_fe))
        tau2_dl = max(0.0, (cochran_q - df) / max(1e-9, c_denom))

        effective_tau2 = max(prior_tau**2, tau2_dl)
        inv_var = 1.0 / (effective_tau2 + ses**2)
        weights = inv_var / np.sum(inv_var)
        pooled_beta = float(np.dot(weights, betas))
        pooled_se = float(np.sqrt(1.0 / np.sum(inv_var)))

        if baseline_abundances:
            a_vec = np.array(
                [max(0.0, float(baseline_abundances.get(k, 0.0))) for k in names], dtype=float
            )
        else:
            a_vec = np.ones_like(betas)
        a_sum = max(1e-9, float(a_vec.sum()))
        abund_parent_log2fc = float(np.log2(np.sum(a_vec * (2.0**betas)) / a_sum))

        return {
            "pooled_effect": round(pooled_beta, 6),
            "pooled_se": round(pooled_se, 6),
            "dersimonian_laird_tau2": round(float(tau2_dl), 6),
            "effective_tau2": round(float(effective_tau2), 6),
            "cochran_q": round(cochran_q, 4),
            "sign_concordance": round(sign_concordance, 4),
            "abundance_weighted_parent_log2fc": round(abund_parent_log2fc, 6),
            "weights": {names[i]: round(float(weights[i]), 4) for i in range(len(names))},
        }
