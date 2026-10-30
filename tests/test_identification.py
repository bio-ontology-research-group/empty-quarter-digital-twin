"""Do-calculus identification with deductive value closure (v0.6)."""
import random
from dataclasses import replace

import pytest
from fastapi.testclient import TestClient

from eq_biodt.v03 import identification as ident
from eq_biodt.v03 import theory
from eq_biodt.v03.api import app
from eq_biodt.v03.process import Parameters, State, simulate

client = TestClient(app)


def test_d_separation_textbook_cases():
    chain = ident.DAG([("a", "b"), ("b", "c")])
    assert not ident.d_separated(chain, "a", "c") and ident.d_separated(chain, "a", "c", "b")
    fork = ident.DAG([("b", "a"), ("b", "c")])
    assert not ident.d_separated(fork, "a", "c") and ident.d_separated(fork, "a", "c", "b")
    col = ident.DAG([("a", "b"), ("c", "b"), ("b", "d")])
    assert ident.d_separated(col, "a", "c")
    assert not ident.d_separated(col, "a", "c", "b")
    assert not ident.d_separated(col, "a", "c", "d")          # descendant of the collider
    assert ident.d_connecting_path(col, "a", "c", "d") is not None and ident.d_connecting_path(col, "a", "c") is None


def test_moral_graph_agrees_with_path_search_on_random_dags():
    rnd = random.Random(1)
    for _ in range(150):
        order = list(range(8))
        edges = [(f"n{i}", f"n{j}") for i in order for j in order if i < j and rnd.random() < .3]
        G = ident.DAG(edges, [f"n{i}" for i in order])
        x, y = rnd.sample(G.nodes and sorted(G.nodes), 2)
        Z = {n for n in G.nodes if n not in (x, y) and rnd.random() < .3}
        assert ident.d_separated(G, x, y, Z) == (ident.d_connecting_path(G, x, y, Z) is None)


def test_unrolled_graph_keeps_slice_edges_and_lags():
    G = ident.unrolled_graph(3, ident.FIELD_LATENT)
    assert G.is_acyclic() and "S" in G.latent
    for e in theory.graph()["edges"]:
        if e["kind"] != "lag":
            s, d = e["source"], e["target"]
            for t in range(3):
                a = f"{s}@{t}" if s in ident.TIME_VARYING else s
                b = f"{d}@{t}" if d in ident.TIME_VARYING else d
                assert (a, b) in G.edges
    assert ("W@0", "W@1") in G.edges and ("N@1", "N@2") in G.edges
    assert {("S", "I@1"), ("S", "T@1"), ("S", "Evap@1")} <= G.edges
    assert ("S", "I@0") not in ident.unrolled_graph(1).edges


def test_closure_zero_propagation():
    for coupling in ("gate", "proportional"):
        cl = ident.closure({"q_co": 0}, coupling)
        assert cl["closed"]["gamma_co"] == 0 and cl["closed"]["k_co@0"] == 0 and cl["closed"]["J_co@0"] == 0
        assert "J_h2@0" not in cl["closed"] and "gamma_h2" not in cl["closed"]
        rules = {d["node"]: d["premises"] for d in cl["derivations"]}
        assert rules["k_co@0"] == ["gamma_co"] and rules["J_co@0"] == ["k_co@0"]
    assert "gamma_co" in ident.closure({"q_co": 0}, "reference")["closed"]
    assert ident.closure({"q_co": 0}, "reference")["closed"]["gamma_co"] == 1.0
    assert "k_co@0" not in ident.closure({"q_co": 0}, "reference")["closed"]


def test_closure_determination_and_activation_rule():
    cl = ident.closure({"W@0": 5.0})
    assert cl["closed"]["theta@0"] == "determined"             # theta = W/(1000 L), L a release constant
    assert "psi_m@0" not in cl["closed"]                        # phi, psi_e, b unobserved and severable
    assert "W@0" not in cl["closed"]                            # observed, not closed
    dry = ident.closure({"psi@0": -100.0})
    assert dry["closed"]["A@0"] == 0 and dry["closed"]["k_h2@0"] == 0 and dry["closed"]["J_co@0"] == 0
    assert dry["closed"]["P@0"] == 0 and dry["closed"]["P_cap@0"] == 0
    assert ident.closure({"psi@0": -10.0})["closed"]["A@0"] == "determined"
    with pytest.raises(ValueError):
        ident.closure({"nonsense": 1})


def test_closure_never_closes_severable_nodes():
    cl = ident.closure({"q_co": 0, "q_h2": 0})
    for n in cl["closed"]:
        base = ident.base_of(n)
        if base not in ("K",):
            try:
                theory.check_admissible([base])
            except ValueError:
                continue
            raise AssertionError(f"{n} is admissible but was closed")


def test_q1_deductive_null_rule3():
    r = ident.deductive_null("gate")
    res = r["result"]
    assert res["graph_test"] is False and res["graph_test_Z_star"] is False
    assert res["deterministic_outcome"] and res["holds"]
    assert {"gamma_co", "k_co@0", "J_co@0"} <= set(res["closed_by_deduction"])
    assert "theta@0" in res["open_path_plain"] and "I@0" in res["open_path_plain"]
    assert all(r["all_admissible_time0_targets"].values())
    assert r["negative_control_reference_coupling"]["holds"] is False
    assert ident.deductive_null("proportional")["result"]["holds"]


def test_q2_oxygen_context_independence():
    r = ident.oxygen_context()
    assert r["stratum"]["holds"] and r["stratum"]["graph_test_Z_star"] is False and r["stratum"]["graph_test_Z_star_csi"]
    assert "D_co@0->J_h2@0" in r["csi_pruned_edges"]
    assert r["without_stratum"]["holds"] is False
    assert r["control_water_still_affects_h2"]["holds"] is False
    chk = r["executable_check"]
    assert chk["gamma_co=0"]["h2_invariant"] and not chk["gamma_co=1"]["h2_invariant"]


def test_q3_field_backdoor():
    r = ident.field_backdoor()
    full = r["all_measured"]["minimal_sets"]
    assert ["I@0", "I@1", "T@1", "W@0", "phi"] in full
    assert r["T_unmeasured"]["identified_by_backdoor"]          # the evaporation record also blocks the season
    assert r["T_unmeasured"]["minimal_sets"] == [["Evap@0", "Evap@1", "I@0", "I@1", "W@0", "phi"]]
    assert not r["T_and_Ev_unmeasured"]["identified_by_backdoor"]
    assert not r["phi_unmeasured"]["identified_by_backdoor"]
    assert all(any(x.startswith(("W@", "Evap@")) for x in s) for s in r["effect_of_T"]["minimal_sets"])
    assert not ident.backdoor_sets(ident.unrolled_graph(2, ident.FIELD_LATENT), "W@1", "J_h2@1", [],
                                   known=r["all_measured"]["known_constants"])["identified_by_backdoor"]
    # no set contains a descendant of the treatment
    G = ident.unrolled_graph(2, ident.FIELD_LATENT)
    desc = G.descendants(["W@1"])
    assert not any(set(s) & desc for s in full)


def test_rules_on_simple_graph_and_closure_caveat():
    G = ident.DAG([("z", "x"), ("z", "y"), ("x", "y")])
    assert ident.rule2(G, "y", [], ["x"], ["z"])["holds"] and not ident.rule2(G, "y", [], ["x"])["holds"]
    assert not ident.rule3(G, "y", [], ["x"])["holds"]
    assert ident.rule1(ident.DAG([("x", "y"), ("y", "w")]), "y", ["x"], ["w"])["holds"] is False
    assert ident.rule1(ident.DAG([("x", "y"), ("x", "w")]), "y", ["x"], ["w"])["holds"]


def test_closure_conditioning_can_open_a_collider_path():
    # Conditioning on a closure node that is a collider is just conditioning: d-separation opens the path.
    G = ident.unrolled_graph(1)
    r = ident.rule1(G, "J_h2@0", [], ["D_co@0"], observed={"q_co": 0})
    assert r["Z_star"] and "k_co@0" in r["closed_by_deduction"]


def test_symbolic_null_matches_executable_model():
    base = Parameters(k_h2_s=.001, k_co_s=.001, gamma_co=0.)      # gamma_co = 0 is what q_CO = 0 closes to
    forcing = lambda w: [{"water_mm": w, "evaporation_mm": 1, "temp_c": 25} for _ in range(3)]
    runs = [simulate(State(water_mm=10), forcing(w), base) for w in (0., 20.)]
    co = [[r["gas"]["co_mol_m2_s"] for r in run["trajectory"]] for run in runs]
    assert co[0] == co[1] == [0., 0., 0.]
    h2 = [run["summary"]["h2_uptake_mol_m2"] for run in runs]
    assert h2[0] != h2[1]                                         # the H2 pathway still responds to do(I)
    live = simulate(State(water_mm=10), forcing(20.), replace(base, gamma_co=1.))
    assert live["trajectory"][0]["gas"]["co_mol_m2_s"] > 0
    assert ident.closure({"q_co": 0})["closed"]["J_co@0"] == 0


def test_api_and_cli():
    r = client.post("/v3/identify", json={"query": "deductive_null"})
    assert r.status_code == 200 and r.json()["result"]["holds"]
    assert client.post("/v3/identify", json={"query": "oxygen_context", "coupling": "proportional"}).json()["stratum"]["holds"]
    f = client.post("/v3/identify", json={"query": "field_backdoor", "days": 2}).json()
    assert f["all_measured"]["identified_by_backdoor"]
    assert client.post("/v3/identify", json={"query": "nope"}).status_code in (400, 422)
    c = client.get("/v3/closure", params={"q_co": 0, "coupling": "gate"}).json()
    assert c["closed"]["J_co@0"] == 0
    c2 = client.post("/v3/closure", json={"observed": {"q_co": 0}, "coupling": "gate", "days": 2}).json()
    assert c2["closed"]["J_co@1"] == 0
    assert client.get("/v3/graph").status_code == 200
    from eq_biodt.v03 import cli
    import sys
    sys.argv = ["eq-biodt", "identify", "--query", "deductive_null"]
    assert cli.main() == 0
