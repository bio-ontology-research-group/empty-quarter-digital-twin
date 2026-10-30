# Legacy v0.2 exploratory model

This directory holds the earlier v0.2 exploratory P0-P6 portfolio model: cloud seeding, revegetation, brine, megafauna and thermal-refugium interventions, Pareto OCS scores, and carbon-storage estimates, together with its API, CLI, validation suite, figure/table scripts and 205-assertion manuscript-number checks.

It is NOT the model described in the manuscript. The current twin lives in `src/eq_biodt/v03/` and is reproduced with `reproducibility/v06/reproduce.py`.

This code is kept for provenance only. It is not packaged (it is outside `src/`), not tested in CI, and not maintained. Layout mirrors the original repository paths (`legacy/src/eq_biodt/...`, `legacy/tests/test_eq_biodt.py`, `legacy/scripts/`, `legacy/reproducibility/`, `legacy/paper/`). To run it, restore the files to their original locations (see git history) and note that it also reads the bundled `reproducibility/data/` and `data/processed/` inputs.
