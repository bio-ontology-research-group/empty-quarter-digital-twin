"""Do-calculus identification with deductive closure: reviews/v06/identification.json and the LaTeX table."""
from pathlib import Path
import json
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
from eq_biodt.v03 import identification as ident  # noqa: E402


def tf(b):
    return "holds" if b else "fails"


def main():
    q1, q2, q3 = (ident.identify(q) for q in ("deductive_null", "oxygen_context", "field_backdoor"))
    out = {"note": "Symbolic, deterministic derivations; no probabilistic estimation. Latent season S is an assumed structure.",
           "deductive_null": q1, "oxygen_context": q2, "field_backdoor": q3}
    (ROOT / "reviews/v06/identification.json").write_text(json.dumps(out, indent=2, allow_nan=False) + "\n")
    r1, r2 = q1["result"], q2["stratum"]
    def sym(x):
        b, _, t = x.partition("@")
        b = {"Evap": "\\mathrm{Ev}", "phi": "\\phi", "J_h2": "J_{\\mathrm{H}_2}"}.get(b, b)
        return b + (f"_{t}" if t else "")
    sets = lambda r: " or ".join("$\\{" + ",".join(sym(x) for x in s) + "\\}$" for s in r["minimal_sets"]) or "none"
    rows = [
        ("Q1 deductive null, $P(J_{\\mathrm{CO}}\\mid do(I_0),q_{\\mathrm{CO}}{=}0)$ (Rule 3)", tf(r1["graph_test"]),
         tf(r1["holds"]), "none needed", "identified; equals $P(J_{\\mathrm{CO}}\\mid q_{\\mathrm{CO}}{=}0)=\\delta_0$"),
        ("Q2 oxygen coupling, $P(J_{\\mathrm{H}_2}\\mid do(D_{\\mathrm{CO}}),q_{\\mathrm{CO}}{=}0)$ (Rule 3)",
         tf(q2["without_stratum"]["graph_test"]), tf(r2["holds"]), "none needed", "independent in the stratum only"),
        ("Q3 field, $P(J_{\\mathrm{H}_2}\\mid do(W_t))$, latent season", "fails (open back-door)",
         "not applicable", sets(q3["all_measured"]) + "; without $T$: " + sets(q3["T_unmeasured"]),
         "identified by back-door; none if $T$ and Ev both unmeasured"),
    ]
    tex = "\n".join(" & ".join(r) + " \\\\" for r in rows) + "\n"
    (ROOT / "paper/tables/v06_identification.tex").write_text(tex)
    print("identification:", q1["verdict"], "|", q2["verdict"], "|", q3["verdict"])


if __name__ == "__main__":
    main()
