"""Do-calculus identification with deductive value closure on the twin's causal-ontological graph.

Implements, symbolically and deterministically, the context-specific ontological d-separation of the
causal-knowledge-graph theory (Definition "deductive value closure", Theorem "context-specific ontological
D-separation"): Pearl's three rules are checked on mutilated graphs, with the conditioning set enlarged to
Z* = Z u Cl(Z = z), the set of nodes whose values are entailed by the invariant core (K, Phi) and the mechanisms
that no admissible regime severs. No probability is estimated; see docs/causal-ontological-theory.md section 6.
"""
from __future__ import annotations

from itertools import combinations

from . import theory

# ---------------------------------------------------------------- graphs


class DAG:
    """Minimal directed graph with the mutilations used by the do-calculus."""

    def __init__(self, edges=(), nodes=(), latent=(), info=None):
        self.edges = {(a, b) for a, b in edges}
        self.nodes = set(nodes) | {x for e in self.edges for x in e}
        self.latent = set(latent)
        self.info = info or {}

    def copy(self) -> "DAG":
        return DAG(self.edges, self.nodes, self.latent, self.info)

    def parents(self, n):
        return {a for a, b in self.edges if b == n}

    def children(self, n):
        return {b for a, b in self.edges if a == n}

    def _closure(self, start, step):
        seen, stack = set(), list(start)
        while stack:
            for m in step(stack.pop()):
                if m not in seen:
                    seen.add(m)
                    stack.append(m)
        return seen

    def ancestors(self, S, inclusive=True):
        S = set(S)
        out = self._closure(S, self.parents)
        return out | S if inclusive else out

    def descendants(self, S, inclusive=True):
        S = set(S)
        out = self._closure(S, self.children)
        return out | S if inclusive else out

    def remove_in(self, S) -> "DAG":
        S = set(S)
        return DAG({e for e in self.edges if e[1] not in S}, self.nodes, self.latent, self.info)

    def remove_out(self, S) -> "DAG":
        S = set(S)
        return DAG({e for e in self.edges if e[0] not in S}, self.nodes, self.latent, self.info)

    def without_edges(self, drop) -> "DAG":
        return DAG(self.edges - set(drop), self.nodes, self.latent, self.info)

    def is_acyclic(self) -> bool:
        indeg = {n: 0 for n in self.nodes}
        for _, b in self.edges:
            indeg[b] += 1
        queue, done = [n for n, d in indeg.items() if d == 0], 0
        while queue:
            n = queue.pop()
            done += 1
            for c in self.children(n):
                indeg[c] -= 1
                if indeg[c] == 0:
                    queue.append(c)
        return done == len(self.nodes)


def _as_dag(G) -> DAG:
    return G if isinstance(G, DAG) else DAG(G)


def _names(x):
    return [x] if isinstance(x, str) else list(x)


def d_separated(G, X, Y, Z=()) -> bool:
    """X and Y are d-separated by Z: separation in the moral graph of the ancestral set of X u Y u Z."""
    G = _as_dag(G)
    X, Y, Z = set(_names(X)), set(_names(Y)), set(_names(Z))
    if X & Y:
        return False
    anc = G.ancestors(X | Y | Z)
    adj = {n: set() for n in anc}
    for a, b in G.edges:
        if a in anc and b in anc:
            adj[a].add(b)
            adj[b].add(a)
    for n in anc:                                 # marry parents
        ps = [p for p in G.parents(n) if p in anc]
        for p, q in combinations(ps, 2):
            adj[p].add(q)
            adj[q].add(p)
    seen, stack = set(), [x for x in X if x not in Z]
    seen.update(stack)
    while stack:
        n = stack.pop()
        if n in Y:
            return False
        for m in adj[n]:
            if m not in Z and m not in seen:
                seen.add(m)
                stack.append(m)
    return True


def d_connecting_path(G, X, Y, Z=()):
    """A d-connecting walk from X to Y given Z (reachability with explicit directions), or None."""
    G = _as_dag(G)
    X, Y, Z = set(_names(X)), set(_names(Y)), set(_names(Z))
    anc_z = G.ancestors(Z)
    # state (node, 'up'): arrived from a child; (node, 'down'): arrived from a parent
    start = [(x, "up") for x in X]
    prev = {s: None for s in start}
    queue = list(start)
    while queue:
        node, d = queue.pop(0)
        if node in Y and node not in X:
            path, s = [], (node, d)
            while s:
                path.append(s[0])
                s = prev[s]
            return path[::-1]
        nxt = []
        if node not in Z:
            if d == "up":
                nxt += [(p, "up") for p in G.parents(node)] + [(c, "down") for c in G.children(node)]
            else:
                nxt += [(c, "down") for c in G.children(node)]
        if d == "down" and node in anc_z:        # collider (or its descendant) activated by conditioning
            nxt += [(p, "up") for p in G.parents(node)]
        for s in nxt:
            if s not in prev:
                prev[s] = (node, d)
                queue.append(s)
    return None


# time-varying nodes of the slice; every other node is static
TIME_VARYING = ("I", "c", "Evap", "T", "C_h2", "C_co", "W", "N", "theta", "psi_m", "psi_o", "psi", "D_h2", "D_co", "A",
                "k_h2", "k_co", "J_h2", "J_co", "P", "P_cap")
FIELD_LATENT = {"S": {"children": ["I", "T", "Evap"], "observed": False,
                      "meaning": "latent season / synoptic regime (assumed structure)"}}


def _t(n, t):
    return f"{n}@{t}" if n in TIME_VARYING else n


def base_of(name: str) -> str:
    return name.split("@")[0]


def time_of(name: str):
    return int(name.split("@")[1]) if "@" in name else None


def unrolled_graph(days: int = 1, latent: dict | None = None) -> DAG:
    """Time-unrolled DAG over `days` slices. Within-slice edges are exactly theory.graph()'s non-lag edges;
    lag edges (W, N) connect slice t to t+1. `latent` maps a node name to {"children": [base ids]}; each latent
    node is static, unobserved, and points to the listed variables in every slice."""
    if days < 1:
        raise ValueError("days must be >= 1")
    g = theory.graph()
    edges, info = set(), {}
    for n in g["nodes"]:
        slots = range(days) if n["id"] in TIME_VARYING else [None]
        for t in slots:
            nm = n["id"] if t is None else f"{n['id']}@{t}"
            info[nm] = {"base": n["id"], "t": t, "zone": n["zone"]}
    for e in g["edges"]:
        s, d = e["source"], e["target"]
        for t in range(days):
            if e["kind"] == "lag":
                if t + 1 < days:
                    edges.add((f"{s}@{t}", f"{d}@{t + 1}"))
            elif d in TIME_VARYING:
                edges.add((_t(s, t), f"{d}@{t}"))
        if e["kind"] != "lag" and d not in TIME_VARYING:
            edges.add((s, d))
    lat = set()
    for ln, spec in (latent or {}).items():
        lat.add(ln)
        info[ln] = {"base": ln, "t": None, "zone": "latent"}
        for ch in spec["children"]:
            for t in range(days):
                edges.add((ln, _t(ch, t)))
    G = DAG(edges, info.keys(), lat, info)
    assert G.is_acyclic()
    return G


# ---------------------------------------------------------------- deductive value closure

GAS = {"h2": "H2", "co": "CO"}


def _non_severable(base: str) -> bool:
    """True iff check_admissible refuses do() on the node, i.e. no admissible regime severs its mechanism."""
    try:
        theory.check_admissible([base])
        return False
    except ValueError:
        return True


def closure(observed: dict, coupling: str = "gate", graph: DAG | None = None,
            dry_threshold_mpa: float = -85.0) -> dict:
    """Deductive value closure Cl_O(observed).

    Cl collects the nodes whose value is entailed by the observed values together with the invariant core:
    the ontology K, the release constants (exogenous_fixed zone) and the mechanisms of non-admissible nodes.
    check_admissible refuses do() on A, k, J, q, gamma, D_K, the constitutive nodes, K, E and the constants, so
    the mechanisms of these nodes are not severed in any admissible regime; an entailment that uses only such
    mechanisms and observed values therefore holds in every admissible regime, which is the regime-universality
    that the theorem needs to add the closure to the conditioning set. Manipulable nodes (I, c, Ev, T, phi,
    psi_e, b, W, N) are never closed by rule: their mechanism may be severed.

    Rules (all fire to a fixpoint):
      (a) Geiger-Verma-Pearl determination: a non-severable node whose parents are all observed or closed is
          closed with symbolic value 'determined'.
      (b) value-specific zero propagation: q_g = 0 => gamma_g = 0 (gate: 1[q>0]; proportional: q/median) =>
          k_g = k_ref gamma_g A Q10^.. = 0 => J_g = C sqrt(D k) tanh(L sqrt(k/D)) = 0 (also after the shared
          oxygen limit, which only scales J_g) => P = P_cap = 0 when both fluxes vanish.
      (c) psi <= psi_d => A = 0 => k_h2 = k_co = 0 => J = 0.
    Observed keys are node ids (static) or 'id@t'. Returns closed nodes, values and the derivation list.
    """
    if coupling not in ("reference", "gate", "proportional"):
        raise ValueError("coupling must be reference, gate or proportional")
    G = graph or unrolled_graph(1)
    unknown = [k for k in observed if k not in G.nodes]
    if unknown:
        raise ValueError(f"Unknown nodes in observed: {unknown}")
    val: dict[str, object] = dict(observed)
    deriv: list[dict] = []

    def close(n, v, rule, prem):
        if n in val:
            return False
        val[n] = v
        deriv.append({"node": n, "value": v, "rule": rule, "premises": sorted(prem)})
        return True

    for n in sorted(G.nodes):
        z = G.info[n]["zone"]
        if n == "K":
            close(n, "invariant", "invariant core K", [])
        elif z == "exogenous_fixed":
            close(n, "constant", "release constant (not admissible, no handle)", [])
    for t in sorted({G.info[n]["t"] for n in G.nodes if G.info[n]["t"] is not None}):
        pass
    times = sorted({G.info[n]["t"] for n in G.nodes if G.info[n]["t"] is not None}) or [None]

    def nm(b, t):
        return f"{b}@{t}" if b in TIME_VARYING else b

    changed = True
    while changed:
        changed = False
        for g in GAS:
            q, gam = f"q_{g}", f"gamma_{g}"
            if coupling == "reference":
                changed |= close(gam, 1.0, "coupling=reference: gamma = 1", [])
            elif q in val and isinstance(val[q], (int, float)):
                if val[q] == 0:
                    changed |= close(gam, 0.0, f"(b) q=0 => gamma=0 under {coupling} coupling", [q])
                elif coupling == "gate":
                    changed |= close(gam, 1.0, "(b) q>0 => gamma=1 under gate coupling", [q])
            for t in times:
                k, J, A = nm(f"k_{g}", t), nm(f"J_{g}", t), nm("A", t)
                if k in G.nodes and J in G.nodes:
                    if val.get(gam) == 0.0:
                        changed |= close(k, 0.0, "(b) gamma=0 => k = k_ref gamma A Q10^.. = 0", [gam])
                    if A in val and val[A] == 0.0:
                        changed |= close(k, 0.0, "(c) A=0 => k=0", [A])
                    if val.get(k) == 0.0:
                        changed |= close(J, 0.0, "(b) k=0 => J = C sqrt(D k) tanh(..) = 0", [k])
        for t in times:
            psi, A = nm("psi", t), nm("A", t)
            if psi in val and isinstance(val[psi], (int, float)) and val[psi] <= dry_threshold_mpa and A in G.nodes:
                changed |= close(A, 0.0, f"(c) psi <= psi_d = {dry_threshold_mpa} MPa => A = 0", [psi])
            P, Pc = nm("P", t), nm("P_cap", t)
            if all(val.get(nm(f"J_{g}", t)) == 0.0 for g in GAS) and P in G.nodes:
                changed |= close(P, 0.0, "(b) J_h2 = J_co = 0 => P = 0", [nm("J_h2", t), nm("J_co", t)])
            if val.get(P) == 0.0 and Pc in G.nodes:
                changed |= close(Pc, 0.0, "(b) P = 0 => P_cap = eta P = 0", [P])
        for n in sorted(G.nodes):                      # (a) determination by closed parents
            if n in val or not G.parents(n):
                continue
            b = G.info[n]["base"]
            if _non_severable(b) and all(p in val for p in G.parents(n)):
                changed |= close(n, "determined", "(a) all parents observed or closed", G.parents(n))
    closed = {n: v for n, v in val.items() if n not in observed}
    return {"observed": dict(observed), "coupling": coupling, "closed": closed, "derivations": deriv,
            "non_trivial_closed": {n: v for n, v in closed.items() if v not in ("constant", "invariant")}}


def csi_edges(G: DAG, cl: dict) -> set:
    """Edges that are inert in the stratum: with k_g = 0 the flux J_g is the constant 0 whatever its inputs, and
    gas g adds nothing to the shared oxygen demand, so D_g and k_g do not reach the other flux."""
    val = {**cl["observed"], **cl["closed"]}
    drop = set()
    for g, o in (("h2", "co"), ("co", "h2")):
        for n in G.nodes:
            if G.info[n]["base"] == f"k_{g}" and val.get(n) == 0.0:
                t = G.info[n]["t"]
                Jg, Jo = f"J_{g}@{t}", f"J_{o}@{t}"
                drop |= {e for e in G.edges if e[1] == Jg}
                drop |= {(f"D_{g}@{t}", Jo), (n, Jo)} & G.edges
    return drop


# ---------------------------------------------------------------- the three rules

def _prep(G, Y, X, Z, W, observed, coupling, use_csi):
    G = _as_dag(G)
    Y, X, Z, W = (set(_names(v)) for v in (Y, X, Z, W))
    observed = dict(observed or {})
    cl = closure(observed, coupling, G if isinstance(G, DAG) and G.info else None) if observed else \
        {"observed": {}, "closed": {}, "derivations": [], "coupling": coupling, "non_trivial_closed": {}}
    closed_nodes = set(cl["closed"]) - set(Y) - set(Z) - set(X)
    Zs = W | set(observed) | closed_nodes
    drop = csi_edges(G, cl) if (use_csi and cl["closed"]) else set()
    return G, Y, X, Z, W, cl, Zs, drop


def _result(rule, plain, star, star_csi, det, Zs, cl, extra, reason_true, reason_false, path):
    holds = plain or star or star_csi or det
    why = ("plain d-separation" if plain else "d-separation after adding Cl to the conditioning set" if star
           else "d-separation in the stratum after pruning edges made inert by the closure" if star_csi
           else "outcome is a deductive constant" if det else "")
    return {"rule": rule, "holds": holds, "graph_test": plain, "graph_test_Z_star": star,
            "graph_test_Z_star_csi": star_csi, "deterministic_outcome": det, "Z_star": sorted(Zs),
            "closed_by_deduction": sorted(cl["closed"]), "open_path_plain": path,
            "reason": (f"holds by {why}: {reason_true}" if holds else reason_false), **extra}


def _all_rules(rule, G, Y, X, Z, W, observed, coupling, use_csi):
    G, Y, X, Z, W, cl, Zs, drop = _prep(G, Y, X, Z, W, observed, coupling, use_csi)
    Wp = W | set(observed or {})
    if rule == 1:
        H = G.remove_in(X)
        test = lambda H, S: d_separated(H, Y, Z, S)
        stmt = "P(y|do(x),z,w) = P(y|do(x),w)"
        mut = "G with edges into X removed"
    elif rule == 2:
        H = G.remove_in(X).remove_out(Z)
        test = lambda H, S: d_separated(H, Y, Z, S)
        stmt = "P(y|do(x),do(z),w) = P(y|do(x),z,w)"
        mut = "G with edges into X and out of Z removed"
    else:
        def zw(S):
            return Z - G.remove_in(X).ancestors(S - Z)
        H = None
        stmt = "P(y|do(x),do(z),w) = P(y|do(x),w)"
        mut = "G with edges into X and into Z(W) removed"
    out = {}
    for label, cond, dr in (("plain", X | Wp, set()), ("star", X | Zs, set()), ("star_csi", X | Zs, drop)):
        Gm = G.without_edges(dr)
        if rule == 3:
            Zw = Z - Gm.remove_in(X).ancestors(cond - Z)
            Hm = Gm.remove_in(X | Zw)
        else:
            Hm = H if not dr else (Gm.remove_in(X) if rule == 1 else Gm.remove_in(X).remove_out(Z))
        out[label] = d_separated(Hm, Y, Z, cond - Z) if True else None
    # a Y that is entailed to be a constant is independent of everything: the context-specific step
    vals = {**cl["observed"], **cl["closed"]}
    det = bool(Y) and all(isinstance(vals.get(y), (int, float)) and not isinstance(vals.get(y), bool) for y in Y)
    path = None
    if not out["plain"]:
        Hp = (G.remove_in(X) if rule == 1 else G.remove_in(X).remove_out(Z) if rule == 2 else
              G.remove_in(X | (Z - G.remove_in(X).ancestors(Wp | X - Z))))
        path = d_connecting_path(Hp, Y, Z, (X | Wp) - Z)
    extra = {"statement": stmt, "mutilation": mut, "csi_pruned_edges": sorted(f"{a}->{b}" for a, b in drop)}
    return _result(rule, out["plain"], out["star"], out["star_csi"], det, Zs, cl, extra,
                   f"{stmt} with Z* = {sorted(Zs)}", f"{stmt} not licensed: d-connected even given Z*", path)


def rule1(G, Y, X, Z, W=(), observed=None, coupling="gate", use_csi=True):
    """Insertion/deletion of observations: Y _||_ Z | X, W in G_{bar X}, with Z* = X u W u Cl."""
    return _all_rules(1, G, Y, X, Z, W, observed, coupling, use_csi)


def rule2(G, Y, X, Z, W=(), observed=None, coupling="gate", use_csi=True):
    """Action/observation exchange: Y _||_ Z | X, W in G_{bar X, underline Z}."""
    return _all_rules(2, G, Y, X, Z, W, observed, coupling, use_csi)


def rule3(G, Y, X, Z, W=(), observed=None, coupling="gate", use_csi=True):
    """Insertion/deletion of actions: Y _||_ Z | X, W in G_{bar X, bar Z(W)}, Z(W) = Z minus anc(W) in G_{bar X}."""
    return _all_rules(3, G, Y, X, Z, W, observed, coupling, use_csi)


# ---------------------------------------------------------------- back-door search

def backdoor_sets(G, x, y, pool, known=(), max_size=None) -> dict:
    """Minimal back-door adjustment sets for P(y | do(x)) among `pool`: no member descends from x and the set
    d-separates x and y in G with the edges out of x removed. Latent nodes cannot be in a set."""
    G = _as_dag(G)
    desc = G.descendants([x])
    cands = [p for p in pool if p not in desc and p not in G.latent and p not in (x, y)]
    Gb = G.remove_out([x])
    found: list[tuple] = []
    for r in range(0, (max_size or len(cands)) + 1):
        for S in combinations(cands, r):
            if any(set(m) <= set(S) for m in found):
                continue
            if d_separated(Gb, x, y, set(S) | set(known)):
                found.append(S)
    return {"x": x, "y": y, "pool": list(pool), "known_constants": sorted(known),
            "minimal_sets": [sorted(s) for s in found], "identified_by_backdoor": bool(found)}


# ---------------------------------------------------------------- the three concrete queries

def _measure_oxygen_invariance():
    """Executable check for Q2: with gamma_co = 0, J_h2 does not change with the CO diffusivity even when the
    shared oxygen limit binds; with gamma_co = 1 it does."""
    from .process import Parameters, gas_exchange
    out = {}
    for gco in (0.0, 1.0):
        vals = []
        for d in (1e-5, 2.1e-5, 8e-5):
            p = Parameters(k_h2_s=.001, k_co_s=.05, d_co_m2_s=d, gamma_co=gco, oxygen_supply_mol_m2_s=3e-10)
            r = gas_exchange(.15, 25., p)
            vals.append({"d_co_m2_s": d, "h2": r["h2_mol_m2_s"], "co": r["co_mol_m2_s"],
                         "oxygen_used": r["oxygen_used_mol_m2_s"]})
        out[f"gamma_co={gco:g}"] = {"runs": vals, "h2_invariant": len({v["h2"] for v in vals}) == 1}
    return out


def deductive_null(coupling: str = "gate", days: int = 1) -> dict:
    """Q1: P(J_CO | do(I@0), q_CO = 0) = P(J_CO | q_CO = 0) as Rule 3."""
    G = unrolled_graph(days)
    obs = {"q_co": 0}
    r = rule3(G, "J_co@0", [], ["I@0"], observed=obs, coupling=coupling)
    cl = closure(obs, coupling, G)
    targets = ["I@0", "c@0", "Evap@0", "T@0", "W@0", "N@0", "phi", "psi_e", "b"]
    allr = {t: rule3(G, "J_co@0", [], [t], observed=obs, coupling=coupling) for t in targets}
    ref = rule3(G, "J_co@0", [], ["I@0"], observed=obs, coupling="reference")
    return {"query": "deductive_null", "coupling": coupling, "days": days,
            "claim": "P(J_CO@0 | do(I@0), q_CO=0) = P(J_CO@0 | q_CO=0)", "rule": "Rule 3", "result": r,
            "closure": cl,
            "plain_graph_test_fails_because": r["open_path_plain"],
            "context_specific_step": ("Cl contains J_CO@0 = 0, a constant; a constant is independent of every "
                                      "intervention, so Y subset Cl gives the equality by determinism even though "
                                      "I@0 -> W -> theta -> D (D_CO or, through the shared oxygen demand, D_H2) -> J_CO stays open: D is not in Cl."),
            "all_admissible_time0_targets": {t: v["holds"] for t, v in allr.items()},
            "negative_control_reference_coupling": {"holds": ref["holds"], "closed_by_deduction": ref["closed_by_deduction"],
                                                    "note": "gamma=1 under reference coupling, so q_CO=0 closes nothing"},
            "verdict": "holds with closure, fails without" if r["holds"] and not r["graph_test"] else
                       "holds without closure" if r["holds"] else "does not hold"}


def oxygen_context(coupling: str = "gate", days: int = 1) -> dict:
    """Q2: in the stratum q_CO = 0, J_H2 does not depend on the CO channel (D_CO, k_CO)."""
    G = unrolled_graph(days)
    obs = {"q_co": 0}
    probe = "D_co@0"
    r = rule3(G, "J_h2@0", [], [probe], observed=obs, coupling=coupling)
    unstrat = rule3(G, "J_h2@0", [], [probe], observed={}, coupling=coupling)
    cl = closure(obs, coupling, G)
    other_gas = rule3(G, "J_h2@0", [], ["I@0"], observed=obs, coupling=coupling)
    return {"query": "oxygen_context", "coupling": coupling, "days": days,
            "formulation": ("Rule 3 with the CO-only probe Z = D_CO@0 (a sensitivity probe on the CO channel; "
                            "check_admissible refuses it as a regime) for Y = J_H2@0. J_H2 has D_CO and k_CO as parents "
                            "only through the shared oxygen demand .5(J_H2 + J_CO). k_CO is closed (= 0), so any path "
                            "through it is blocked, and J_CO = 0 for every D_CO, so the direct edge D_CO -> J_H2 is "
                            "inert in the stratum (context-specific independence, not graph d-separation)."),
            "claim": "P(J_H2@0 | do(D_CO@0), q_CO=0) = P(J_H2@0 | q_CO=0)",
            "stratum": r, "without_stratum": unstrat,
            "csi_pruned_edges": r["csi_pruned_edges"], "closure": cl,
            "control_water_still_affects_h2": {"holds": other_gas["holds"],
                                                "note": "do(I@0) still changes J_H2: only the CO channel is cut"},
            "executable_check": _measure_oxygen_invariance(),
            "verdict": "independent in the stratum; dependent without it" if r["holds"] and not unstrat["holds"]
            else "unexpected"}


def field_backdoor(days: int = 2) -> dict:
    """Q3: observational field data with a latent season S. Identification of P(J_H2@t | do(W@t)) by back-door."""
    if days < 2:
        raise ValueError("field_backdoor needs days >= 2 (lag edges)")
    G = unrolled_graph(days, FIELD_LATENT)
    t = days - 1
    W, J = f"W@{t}", f"J_h2@{t}"
    fixed = sorted(n for n in G.nodes if G.info[n]["zone"] == "exogenous_fixed")
    measured = [f"{b}@{s}" for s in range(days) for b in ("I", "T", "Evap", "W", "J_h2")]
    site = ["phi"]

    def go(drop=(), x=W, y=J, extra=()):
        pool = [m for m in measured + site + list(extra) if m not in drop]
        return backdoor_sets(G, x, y, pool, known=fixed)
    full = go()
    no_T = go(drop=[f"T@{s}" for s in range(days)])
    no_TEv = go(drop=[f"T@{s}" for s in range(days)] + [f"Evap@{s}" for s in range(days)])
    no_phi = go(drop=["phi"])
    paths = d_connecting_path(G.remove_out([W]), W, J, set(fixed))
    eff_T = backdoor_sets(G, f"T@{t}", J, measured + site, known=fixed)
    return {"query": "field_backdoor", "days": days, "latent": FIELD_LATENT,
            "estimand": f"P({J} | do({W}))", "assumptions": [
                "season S is a latent common cause of I, T, Ev in every slice (assumed structure, not estimated)",
                "I, T, Ev, W, J_H2 measured each day; phi measured as a site property; release constants known",
                "lag edges W@t -> W@t+1 and N@t -> N@t+1 as in theory.graph(); N and psi are not measured"],
            "open_backdoor_path_without_adjustment": paths,
            "all_measured": full, "T_unmeasured": no_T, "T_and_Ev_unmeasured": no_TEv, "phi_unmeasured": no_phi,
            "effect_of_T": eff_T,
            "formula": (f"P({J} | do({W}=w)) = sum_z P({J} | {W}=w, z) P(z) over z in each minimal set"),
            "reading": ("Back-door paths leave W@t through its parents I@t, Ev@t, W@(t-1) and W_max (phi). "
                        "S opens W <- Ev <- S -> T -> k -> J, W <- I <- S -> T -> k -> J and, through the previous slice, "
                        "W(t) <- Ev(t) <- S -> I(t-1) -> N(t-1) -> N(t) -> psi_o -> psi -> A -> k -> J. Adjusting for "
                        "T@t alone leaves the last route open; I@(t-1), I@t and W@(t-1) are also needed (W <- I -> N and "
                        "W(t-1) -> N(t-1) -> N(t)). The search finds two minimal sets: one contains T@t, the other "
                        "replaces T@t by Ev@(t-1), Ev@t. Temperature is therefore a sufficient blocker, not a necessary one: "
                        "without T the effect is still identified through the evaporation record, and it is lost by "
                        "back-door only if both T and Ev are unmeasured (or phi is)."),
            "positivity_caveat": ("Identification here is graphical. The Sonoran campaigns have near-collinear moisture "
                                  "and temperature across seasons, so the strata of the adjustment set are almost empty "
                                  "at the contrasts of interest (a positivity violation); that, not the graph, is why "
                                  "the Q10 question was not resolved."),
            "limits": "back-door only: no front-door or full ID-algorithm; absence of a set is not a proof of non-identifiability",
            "verdict": "identified by back-door given phi, I@(t-1), I@t, W@(t-1) and either T@t or Ev@(t-1), Ev@t; no back-door set otherwise"}


QUERIES = {"deductive_null": deductive_null, "oxygen_context": oxygen_context, "field_backdoor": field_backdoor}


def identify(query: str, **kw) -> dict:
    if query not in QUERIES:
        raise ValueError(f"Unknown query {query!r}; choose from {sorted(QUERIES)}")
    return QUERIES[query](**kw)
