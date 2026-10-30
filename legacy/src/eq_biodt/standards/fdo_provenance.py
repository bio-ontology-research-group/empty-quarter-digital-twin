"""BioDT FAIR Digital Object (FDO) & Workflow Run RO-Crate (wfrun/0.4) provenance generator.

Conforms to `BioDT/biodt-DTOs` and WorkflowHub RO-Crate v1.3 / Process Run Crate v0.4:
- Wraps every Digital Twin counterfactual simulation or multi-objective optimization
  in a verifiable `CreateAction` entity with SHA-256 digests of input datasets,
  ontology constraint specification, software version, and Ontological Consistency Score (`OCS`).
"""

from __future__ import annotations

import datetime as dt
import hashlib
import json
import uuid
from pathlib import Path
from typing import Any, Dict, Optional, Tuple


class FDOProvenanceRecorder:
    """Builds FAIR Digital Object (FDO) and Workflow Run RO-Crate records for twin runs."""

    WFRUN_PROFILE = "https://w3id.org/ro/wfrun/process/0.4"
    RO_CRATE_CONTEXT = "https://w3id.org/ro/crate/1.1/context"
    SOFTWARE_VERSION = "eq-biodt/0.2.0"
    _REPO_ROOT = Path(__file__).resolve().parents[3]
    _DIGEST_CACHE: Dict[str, str] = {}

    @classmethod
    def _file_sha256(cls, rel_path: str) -> str:
        if rel_path in cls._DIGEST_CACHE:
            return cls._DIGEST_CACHE[rel_path]
        target = cls._REPO_ROOT / rel_path
        if target.exists() and target.is_file():
            digest = f"sha256:{hashlib.sha256(target.read_bytes()).hexdigest()}"
            cls._DIGEST_CACHE[rel_path] = digest
            return digest
        return f"unresolved:{rel_path}"

    @classmethod
    def _default_input_digests(cls) -> Dict[str, str]:
        candidates = (
            "data/processed/eq_compartment_state_cube.parquet",
            "data/processed/eq_sites_cube.parquet",
            "data/processed/eq_era5_h3_corridor_grid.parquet",
        )
        return {rel: cls._file_sha256(rel) for rel in candidates}

    @classmethod
    def build_run_crate(
        cls,
        intervention_name: str,
        site_iri: str,
        parameters: Dict[str, Any],
        outputs: Dict[str, Any],
        ocs_score: float = 1.0,
        version_ambiguity_interval: Optional[Tuple[float, float]] = None,
        input_dataset_digests: Optional[Dict[str, str]] = None,
        ontology_constraint_digest: Optional[str] = None,
    ) -> Dict[str, Any]:
        """Generate a JSON-LD Workflow Run RO-Crate + FDO kernel metadata dictionary."""
        now_iso = dt.datetime.now(dt.timezone.utc).isoformat()
        input_digests = input_dataset_digests or cls._default_input_digests()
        onto_digest = ontology_constraint_digest or cls._file_sha256(
            "data/ontology/rubalkhali_twin_constraints.owl"
        )

        payload_bytes = json.dumps(
            {
                "software_version": cls.SOFTWARE_VERSION,
                "intervention": intervention_name,
                "site_iri": site_iri,
                "parameters": parameters,
                "outputs": outputs,
                "input_dataset_digests": input_digests,
                "ontology_constraint_digest": onto_digest,
            },
            sort_keys=True,
            default=str,
        ).encode("utf-8")
        digest = hashlib.sha256(payload_bytes).hexdigest()
        run_id = f"urn:uuid:{uuid.uuid5(uuid.NAMESPACE_OID, digest)}"
        fdo_pid = f"https://rubalkhali.science/fdo/eq-biodt/{digest[:16]}"

        amb_low, amb_high = version_ambiguity_interval or (0.0, 0.0)

        return {
            "@context": cls.RO_CRATE_CONTEXT,
            "fdo_kernel_metadata": {
                "fdo_pid": fdo_pid,
                "fdo_profile": "https://github.com/BioDT/biodt-DTOs#FDOProfile",
                "software_version": cls.SOFTWARE_VERSION,
                "sha256_checksum": digest,
                "ontology_constraint_digest": onto_digest,
                "input_dataset_digests": input_digests,
                "created_at": now_iso,
                "ontological_consistency_score_ocs": round(float(ocs_score), 6),
                "version_ambiguity_bounds": [round(float(amb_low), 6), round(float(amb_high), 6)],
            },
            "@graph": [
                {
                    "@id": "ro-crate-metadata.json",
                    "@type": "CreativeWork",
                    "conformsTo": [
                        {"@id": "https://w3id.org/ro/crate/1.1"},
                        {"@id": cls.WFRUN_PROFILE},
                    ],
                    "about": {"@id": "./"},
                },
                {
                    "@id": "./",
                    "@type": "Dataset",
                    "name": f"EQ-BioDT Simulation Run: {intervention_name}",
                    "version": cls.SOFTWARE_VERSION,
                    "datePublished": now_iso,
                    "mentions": [{"@id": run_id}],
                },
                {
                    "@id": run_id,
                    "@type": "CreateAction",
                    "name": f"Counterfactual {intervention_name} at {site_iri}",
                    "startTime": now_iso,
                    "endTime": now_iso,
                    "actionStatus": "http://schema.org/CompletedActionStatus",
                    "instrument": {
                        "@id": "src/eq_biodt/models/causal_ontological_twin.py",
                        "version": cls.SOFTWARE_VERSION,
                    },
                    "location": {"@id": site_iri},
                    "object": [
                        {
                            "@type": "PropertyValue",
                            "name": k,
                            "value": v,
                        }
                        for k, v in parameters.items()
                    ],
                    "result": {
                        "@id": fdo_pid,
                        "@type": "Dataset",
                        "sha256": digest,
                        "ontologicalConsistencyScore": round(float(ocs_score), 6),
                        "summary": outputs,
                    },
                },
            ],
        }

    @classmethod
    def save_run_crate(cls, crate_dict: Dict[str, Any], output_path: Path) -> Path:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        output_path.write_text(json.dumps(crate_dict, indent=2, default=str), encoding="utf-8")
        return output_path
