# Audit publication record

This audit is maintained across three existing repositories:

| Repository | Target branch | Change |
|---|---|---|
| [Digital Twin](https://github.com/bio-ontology-research-group/empty-quarter-digital-twin) | `main` | Scientific audit, machine-readable evidence, diagnostics, coefficient inventory, API console and documentation |
| [Website](https://github.com/bio-ontology-research-group/empty-quarter/commit/4f232a5d2ca4b4f8298bd500d3d81b8b8a15532a) | `analysis/16s-statistical-review` | `4f232a5d2ca4b4f8298bd500d3d81b8b8a15532a`: evidence guide, corrected displays and the existing twin dependencies needed for an independent build |
| [Overleaf manuscript](https://www.overleaf.com/project/6addceec3f96e54f5b1ddd02) | `main` | `2e3c565dcac858505a8dfa8a0c925e274b7b2a1b`: explicit ontological-causal methods, implemented carbon equations, assumptions note and bibliography corrections |

The working source under `Public/software/empty-quarter-digital-twin` has no Git metadata. The canonical model checkout is `Documents/papers/empty-quarter-digital-twin`, with a separate Overleaf repository in its `paper/` directory. Publication was prepared from existing remote history in isolated checkouts under `~/tmp/eq-twin-publication-20261026`. The website's unrelated pending changes were preserved. No force push or website deployment is part of this update.

## Verification

- Original audit: 9 pytest tests and 205 numerical fixtures passed; all 23 registered source hashes matched.
- Additional contract audit: 4 passes and 7 failures, retained as findings. Python model source hashes are unchanged.
- Publication checkout: production frontend build and 7 scenario-response tests passed.
- Publication checkout: 10 browser checks passed, including a real API scenario, mobile overflow, keyboard interaction and failed-request handling.
- Overleaf checkout: both main and supplementary PDFs compiled successfully.

The [audit report](report.md) is the scientific interpretation. [Model contract results](model-contracts.json) and [browser checks](browser-checks.json) are machine-readable evidence. The [coefficient inventory](../../reproducibility/outputs/model_coefficient_inventory.json) covers 1,500 source locations in 28 Python files; it does not certify their empirical support.

## Reproduce checks

From the model repository:

```bash
uv sync
uv run pytest -q
uv run python scripts/audit_model_contracts.py
uv run python scripts/inventory_model_assumptions.py
```

The contract script records failed claims without modifying coefficients. Browser checks additionally require Chromium at `/usr/bin/chromium`, a running website and the API console:

```bash
TWIN_WEBSITE_URL=http://127.0.0.1:5188/twin \
TWIN_CONSOLE_URL=http://127.0.0.1:8098/ \
uv run --with playwright python reviews/validation-2026-10-26/browser_check.py
```

Build the manuscript with `latexmk -pdf main.tex` and `latexmk -pdf supplementary.tex` in the separate Overleaf checkout. In the website repository, run `npm ci` and `npm run build` from `website/frontend`, and `node --test tests/twin/test_scenario_impacts.mjs` from the repository root.

Before-edit backups, one-time editing scripts, live-server logs, dependency directories and LaTeX build auxiliaries are kept out of the published audit. The before/after patch and retained logs preserve the review evidence. Source publication does not resolve the seven model findings or establish causal effects.
