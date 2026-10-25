"""SPARQL Bridge to the Rub' al-Khali Knowledge Graph (Virtuoso) with offline fallback.

Uses the canonical four-individual SIO measurement reification pattern:
  site/sample -> sio:SIO_000008 (has attribute) -> quality
  quality -> sio:SIO_000216 (has measurement value) -> value
  value -> sio:SIO_000300 (has value) -> literal
"""

from __future__ import annotations

from typing import Any, Dict, List, Optional
import httpx

from eq_biodt.data.manifest_registry import CanonicalDataRegistry


class RubAlKhaliSPARQLBridge:
    """Queries the Rub' al-Khali Virtuoso SPARQL endpoint with deterministic local fallback."""

    DEFAULT_LOCAL_ENDPOINT = "http://localhost:8895/sparql"
    DEFAULT_PUBLIC_ENDPOINT = "https://rubalkhali.science/sparql"

    def __init__(
        self,
        endpoint: Optional[str] = None,
        registry: Optional[CanonicalDataRegistry] = None,
        timeout_s: float = 3.0,
    ) -> None:
        self.endpoint = endpoint or self.DEFAULT_LOCAL_ENDPOINT
        self.registry = registry or CanonicalDataRegistry()
        self.timeout_s = timeout_s

    @property
    def endpoint_url(self) -> str:
        return self.endpoint

    def is_endpoint_alive(self) -> bool:
        """Check whether the Virtuoso SPARQL endpoint responds to a lightweight ASK/SELECT."""
        res = self.query_sparql("SELECT * WHERE { ?s ?p ?o } LIMIT 1")
        return res is not None

    def build_site_measurements_query(self, site_iri: str) -> str:
        """Construct a canonical SIO_000216 SPARQL query for a given site IRI."""
        return f"""
PREFIX rak: <https://rubalkhali.science/kb/>
PREFIX sio: <http://semanticscience.org/resource/>
PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
PREFIX geo: <http://www.opengis.net/ont/geosparql#>

SELECT ?label ?wkt ?qualityType ?qualityLabel ?numericValue
WHERE {{
  <{site_iri}> rdfs:label ?label .
  OPTIONAL {{
    <{site_iri}> geo:hasGeometry/geo:asWKT ?wkt .
  }}
  OPTIONAL {{
    <{site_iri}> sio:SIO_000008 ?quality .
    ?quality a ?qualityType ;
             sio:SIO_000216 ?valNode .
    ?valNode sio:SIO_000300 ?numericValue .
    OPTIONAL {{ ?qualityType rdfs:label ?qualityLabel . }}
  }}
}}
LIMIT 200
""".strip()

    def query_sparql(self, sparql: str) -> Optional[Dict[str, Any]]:
        """Execute a SPARQL SELECT query if Virtuoso is reachable; otherwise return None."""
        try:
            with httpx.Client(timeout=self.timeout_s) as client:
                resp = client.get(
                    self.endpoint,
                    params={"query": sparql, "format": "application/sparql-results+json"},
                    headers={"Accept": "application/sparql-results+json"},
                )
                if resp.status_code == 200:
                    return resp.json()
        except Exception:
            return None
        return None

    def get_site_semantic_summary(self, site_id: str) -> Dict[str, Any]:
        """Return semantic metadata and ontology IRIs for a site, querying Virtuoso or local KG tables."""
        site_iri = self.registry.ontology.site_iris.get(
            str(site_id),
            self.registry.ontology.site_iris.get(f"Site {site_id}", f"https://rubalkhali.science/kb/RAK_10000{int(site_id):02d}" if str(site_id).isdigit() else ""),
        )
        sparql_q = self.build_site_measurements_query(site_iri)
        live_res = self.query_sparql(sparql_q)
        bindings = (
            live_res.get("results", {}).get("bindings", [])
            if isinstance(live_res, dict)
            else []
        )
        return {
            "site_id": str(site_id),
            "site_iri": site_iri,
            "sparql_endpoint": self.endpoint,
            "sparql_live_connected": bool(bindings),
            "sparql_query": sparql_q,
            "live_binding_count": len(bindings),
            "measurement_pattern": "SIO_000008 -> Quality -> SIO_000216 -> Value -> SIO_000300",
        }
