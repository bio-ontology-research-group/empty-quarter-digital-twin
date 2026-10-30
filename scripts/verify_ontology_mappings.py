"""Verify every external IRI used by the v03 ontology against the live OLS4 API and write
data/v03/ontology_mappings.tsv (rak term, axiom type, target IRI, verified label, ontology, retrieval date).

Run: uv run python scripts/verify_ontology_mappings.py   (needs network; fails on any label mismatch)
"""
import csv
import datetime as dt
import json
import sys
import urllib.parse
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from eq_biodt.v03.ontology import MAPPINGS, VOCABULARY, _mapping_iri  # noqa: E402

OLS = "https://www.ebi.ac.uk/ols4/api"


def lookup(iri: str) -> dict:
    for kind in ("terms", "properties", "individuals"):
        url = f"{OLS}/{kind}?iri={urllib.parse.quote(iri, safe='')}"
        try:
            data = json.load(urllib.request.urlopen(url, timeout=60))
        except Exception:
            continue
        found = data.get("_embedded", {}).get(kind, [])
        if found:
            term = next((t for t in found if t.get("is_defining_ontology")), found[0])
            return {"label": term["label"], "ontology": term["ontology_name"], "obo_id": term.get("obo_id") or "",
                    "obsolete": bool(term.get("is_obsolete"))}
    raise SystemExit(f"IRI not found in OLS4: {iri}")


def main():
    today = dt.date.today().isoformat()
    cache: dict[str, dict] = {}
    rows = []

    def verified(iri, expected):
        if iri not in cache:
            cache[iri] = lookup(iri)
        got = cache[iri]
        if got["label"] != expected or got["obsolete"]:
            raise SystemExit(f"Label/obsolescence mismatch for {iri}: OLS={got} expected={expected!r}")
        return got

    for term, rel, local, label in MAPPINGS:
        iri = _mapping_iri(local)
        got = verified(iri, label)
        axiom = {"SubClassOf": "SubClassOf", "(existing) child": "existing SubClassOf(external rak:parent) usage"}.get(rel, f"AnnotationAssertion({rel})")
        rows.append([f"rak:{term}", axiom, iri, got["obo_id"], got["label"], got["ontology"], today])
    for iri, label in VOCABULARY:
        got = verified(iri, label)
        rows.append(["(vocabulary)", "term used by the ontology/export", iri, got["obo_id"], got["label"], got["ontology"], today])
    out = ROOT / "data/v03/ontology_mappings.tsv"
    with out.open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["rak_term", "axiom_type", "target_iri", "target_curie", "verified_label", "source_ontology", "retrieved_ols4"])
        w.writerows(rows)
    print(f"{len(rows)} rows, {len(cache)} distinct IRIs verified -> {out}")


if __name__ == "__main__":
    main()
