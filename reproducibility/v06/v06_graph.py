"""The v0.6 view of the causal graph.

v0.7 adds the exogenous boundary-condition nodes C_h2 and C_co (atmospheric H2/CO as forcing) and renames C_a to 'Default
atmospheric gas concentrations'. The v0.6 outputs describe the v0.6 graph, so the v0.6 scripts read it through this function
and stay byte-identical (the v0.7 graph is documented in the v0.7 composition notes)."""
from copy import deepcopy

V07_NODES = {"C_h2", "C_co"}
V06_LABELS = {"C_a": "Atmospheric gas concentrations"}


def graph(theory):
    g = deepcopy(theory.graph())
    g["nodes"] = [{**n, "label": V06_LABELS.get(n["id"], n["label"])} for n in g["nodes"] if n["id"] not in V07_NODES]
    g["edges"] = [e for e in g["edges"] if e["source"] not in V07_NODES and e["target"] not in V07_NODES]
    return g
