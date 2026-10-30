# Reproducibility (`reproducibility/`)

The current workflow is `reproducibility/v05/reproduce.py`. It verifies the bundled source digests, then runs `reproducibility/v03/reproduce.py`, the external-transfer, primary-connection, audit and figure steps, and writes the manuscript tables and figures.

```bash
uv sync
uv run python scripts/install_reasoner.py
uv run python reproducibility/v05/reproduce.py
uv run pytest -q
```

## Contents

- `v03/` - v0.3/v0.4 pipeline stages (evidence rebuild, audit, website export, ontology design audit) invoked by the v0.5 workflow.
- `v05/` - v0.5 workflow: `reproduce.py`, `external_transfer.py`, `primary_connection.py`, `audit.py`, `figures.py`, `fetch_pmc_cloud.py`.
- `data/`, `outputs/`, `independent_reference.py` - bundled reference tables and generated outputs retained from earlier versions.

## Legacy material

The earlier v0.2 exploratory model, its `reproduce_all.py` pipeline, the 205-assertion numerical verifier and the old Figures 1-3 are archived under `legacy/` for provenance only. They do not describe the model in the current manuscript and are not run in CI.
