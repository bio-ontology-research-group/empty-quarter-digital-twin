# Twin question pages and full uncertainty analysis — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Release EQ-BioDT 0.8.0 with a full uncertainty analysis for 14 questions, and rebuild rubalkhali.science/model as a question catalogue with a plain and a scientific page per question.

**Architecture:** A question is a list of series of levels; each level is a replayable chain spec (action, upstream models, twin request, downstream models) evaluated by the existing `composition.run_chain`. The UQ module samples literature constants and twin parameters jointly (Saltelli design on a scrambled Sobol' sequence), evaluates chunks as a Slurm array on unimatrix, gathers and analyses on ws, and writes one JSON per question; the website renders those JSONs as paired plain/science pages.

**Tech Stack:** Python 3.12, numpy, scipy (`scipy.stats.qmc`), existing `eq_biodt.v03` (process, composition, external_models, community); Slurm on unimatrix; React 19 + Vite + react-bootstrap + react-router-dom 7; Playwright (Python) for browser checks; pdflatex/bibtex.

**Spec:** `docs/superpowers/specs/2026-10-31-twin-questions-site-design.md` (twin repo, branch `v0.8-questions`).

**Repos and places:**
- Twin: ws `/mnt/data1/agents/projects/eq-biodt-v07/twin`, branch `v0.8-questions` (from `main` = v0.7.0).
- Website: ws `/mnt/data1/agents/projects/eq-biodt-v07/web`, new branch `twin-v08-questions` from `analysis/16s-statistical-review` (`811ef3b8`). Paths below prefixed `web:` are relative to `web/website/frontend/`.
- UQ cluster: unimatrix `/home/hohndor/eq-biodt-uq/` (clone + uv env), scratch `$SLURM_TMPDIR`.

## Global Constraints

- All compute on ws (`ssh ws-ts`) except UQ sampling on unimatrix (`ssh unimatrix01`, Slurm partition `debug`); nothing heavy on the laptop.
- Never write to `/tmp`; on ws `TMPDIR=$HOME/tmp`; on unimatrix job I/O in `$SLURM_TMPDIR`, results small JSON under `/home/hohndor/eq-biodt-uq/results/`.
- v0.6 and v0.7 outputs (`reviews/v06`, `reviews/v07`, their paper tables/figures) must stay byte-identical; `uv run python reproducibility/v06/verify_manuscript.py` must pass after every paper change.
- Every number shown on the website comes from exported JSON; no hand-typed numbers.
- Parameter-uncertainty intervals only; every place that shows an interval says it excludes model-structure error (plain wording fixed in the spec).
- Distributions: literature constants triangular(low, value, high); twin parameters per `data/external/twin_parameter_ranges.json`.
- Saltelli design, N = 1024 base samples, ≤ 7 levels per series; seeds derived from task keys; results independent of worker count, node and task order.
- Refusal rule: refusal fraction < 5 % → intervals + Sobol' on admitted samples (flagged); ≥ 5 % → interval on admitted samples + Sobol' indices of the refusal indicator.
- Closed deductive gate → exact 0, interval [0, 0], provenance `deduced`, not sampled.
- Plain-page jargon list (must not appear in plain texts): Greek letters, `do(`, `Sobol`, `ontology`, `MAG`, `EC_`, `ECe`, `Maas`, `CBB`, `γ`, model ids (any `snake_case` token), unit symbols other than °C, %, mm, g, m², kg, dS/m.
- Commit messages end with `Claude-Session: https://claude.ai/code/session_01QpuQV6aqPiEnJ67xiyvSrj`; never `git add -A`; push only after checks pass.
- Deploy only after a local `vite preview` + Playwright check of every page type at 1280 and 390 px with 0 page errors; the bundle/question version gate must accept `0.6.x`–`0.8.x`.

## Review Focus

1. A question whose level is refused for most samples (salty water at 8 dS/m, NaCl > 0.1 M): the page must show the refusal plainly, not an empty or zero bar — test in Task 6 (`test_refused_level_rendering_fields`) and Task 10 (Playwright on `salty-water` top level).
2. A closed gate inside a sampled chain (carbon bound for 14PRr2, CO for 51PRr3): must be exact 0 in every draw with zero variance and no NaN indices — test in Task 3 (`test_closed_gate_is_exact_zero`).
3. Chunked cluster evaluation vs single-machine evaluation of the same task: byte-identical gathered JSON — test in Task 4 (`test_chunked_equals_single`).
4. Old URLs (`/model?view=research#identification`, `/model#composition`, `/model?view=whatif`) land on the matching new page — test in Task 9 (Playwright `test_redirects`).
5. Plain text when an interval spans zero ("may go up or down") and when a value is an upper bound ("at most") — test in Task 6 (`test_plain_text_rules`).

---

### Task 1: Twin parameter registry

**Files:**
- Create: `data/external/twin_parameter_ranges.json`
- Create: `src/eq_biodt/v03/uq_inputs.py`
- Test: `tests/test_uq_inputs.py`

**Interfaces:**
- Produces: `load_parameter_ranges() -> dict[str, dict]` (keys: `id, dist ∈ {"triangular","uniform","loguniform"}, low, value, high, unit, meaning, status, source{authors,year,title,doi_or_url,page_or_table}, note`); `Input` dataclass `(name: str, kind: Literal["constant","parameter"], dist: str, low: float, value: float, high: float, status: str)`; `inputs_for(constant_ids: list[str], exclude_params: set[str]) -> list[Input]` (constants from the composition registry, then parameters, sorted by name).

- [ ] **Step 1: Write the failing tests**

```python
def test_registry_has_required_parameters():
    R = load_parameter_ranges()
    assert {"k_h2_s","k_co_s","psi_entry_mpa","retention_b","wet_activation_mpa","dry_threshold_mpa",
            "activation_shape","q10","d_h2_m2_s","d_co_m2_s","drainage_fraction_day",
            "field_capacity_fraction","porosity"} <= set(R)
    for e in R.values():
        assert e["status"] in {"verified_in_source","modeller_assumption","source_not_checked"}
        assert e["low"] <= e["value"] <= e["high"]
        assert e["value"] == getattr(Parameters(k_h2_s=.001, k_co_s=.001), e["id"])

def test_existing_ensemble_bounds_reused():
    R = load_parameter_ranges()
    assert (R["k_h2_s"]["dist"], R["k_h2_s"]["low"], R["k_h2_s"]["high"]) == ("loguniform", 1e-4, 1e-2)
    assert (R["retention_b"]["low"], R["retention_b"]["high"]) == (2.0, 6.0)
    assert (R["dry_threshold_mpa"]["low"], R["dry_threshold_mpa"]["high"]) == (-100.0, -70.0)
    assert (R["activation_shape"]["low"], R["activation_shape"]["high"]) == (0.5, 2.0)
    assert (R["porosity"]["low"], R["porosity"]["high"]) == (0.35, 0.5)
    assert (R["psi_entry_mpa"]["dist"], R["psi_entry_mpa"]["low"], R["psi_entry_mpa"]["high"]) == ("loguniform", -0.01, -0.0001)
    assert (R["q10"]["low"], R["q10"]["high"]) == (1.0, 2.0)

def test_inputs_for_excludes_intervened_parameter():
    names = [i.name for i in inputs_for(["fao29_ece_per_ecw"], exclude_params={"porosity"})]
    assert "porosity" not in names and "fao29_ece_per_ecw" in names and names == sorted(names)
```

- [ ] **Step 2: Run to verify failure** — `uv run pytest tests/test_uq_inputs.py -q` → FAIL (module missing).
- [ ] **Step 3: Write the registry and loader.** Bounds above come from `src/eq_biodt/v03/planning.py:parameter_ensemble` (rates ×0.1–×10, b 3–5 widened to the transfer grid 2–6, dry −70…−100, shape 0.5–2) and `reproducibility/v05/external_transfer.py:CANDIDATES` (porosity 0.35–0.5, ψ_entry −1e-4…−1e-2, Q10 1–2); cite them in `source` with `page_or_table` = file:function and status `modeller_assumption`. Remaining parameters (diffusivities ±20 %, drainage 0.3–0.7, field capacity 0.4–0.6, wet activation −0.3…−0.03 log) get status `modeller_assumption` with a `note`. Log-uniform for negative values samples the magnitude.
- [ ] **Step 4: Run tests** → PASS.
- [ ] **Step 5: Commit** `feat(v0.8): twin parameter range registry and UQ input space`.

### Task 2: Question specifications (14 questions as replayable level specs)

**Files:**
- Create: `src/eq_biodt/v03/questions.py`
- Test: `tests/test_questions.py`

**Interfaces:**
- Consumes: `composition.run_chain`, `composition.twin_request`, `external_models.MODELS`, `reproducibility/v06/analyse.reference_schedules` (import pattern as in `reproducibility/v07/composition.py`).
- Produces:
  - `LevelSpec(label: str, x: float | str, action: dict, upstream: list[str], twin: dict, downstream: list[str], intervened_params: frozenset[str], sample: str | None)` where `twin` holds `twin_request` kwargs (forcing, params overrides, community, coupling, intervention).
  - `Series(id: str, label_plain: str, label_science: str, levels: list[LevelSpec])`.
  - `Question(slug, theme, title_plain, title_science, level_kind ∈ {"numeric","categorical","sample"}, x_label_plain, x_label_science, series: list[Series], headline: list[str], baseline_label: str)`.
  - `QUESTIONS: dict[str, Question]` with slugs: `rain-pattern, drip-or-flood, cloud-seeding, salty-water, alfalfa-bacterium, hotter, mulch-shade, solar-park, compaction, biochar, missing-co-genes, adding-microbes, air-hydrogen, carbon-fixation`; themes `water, salt, heat, soil, microbes`.
  - `evaluate(spec: LevelSpec, constants: dict[str,float] | None = None, params: dict[str,float] | None = None) -> dict[str, float | None | str]` — runs `run_chain(..., constants=constants, envelope=False)` with `Parameters` replaced by `params`; returns the headline outputs (keys as in the spec table, e.g. `h2_uptake_contrast_umol_m2`, `relative_yield_alfalfa`), `"__refused__": reason` when refused, and exact `0.0` with key suffix provenance via `__deduced__: [outputs]` for closed gates.

- [ ] **Step 1: Write the failing tests**

```python
def test_fourteen_questions_with_themes():
    assert len(QUESTIONS) == 14
    assert {q.theme for q in QUESTIONS.values()} == {"water","salt","heat","soil","microbes"}
    assert all(1 <= len(s.levels) <= 7 for q in QUESTIONS.values() for s in q.series)

def test_v07_levels_reproduce_stored_central_values():
    # each series sourced from a v0.7 experiment reproduces reviews/v07 central values exactly
    stored = json.loads((ROOT/"reviews/v07/composition/brackish_irrigation.json").read_text())
    spec = next(l for l in QUESTIONS["salty-water"].series[1].levels if l.x == 8.0)
    out = evaluate(spec)
    assert out["relative_yield_alfalfa"] == pytest.approx(stored_central_ry_alfalfa(stored, 8.0), rel=1e-12)

def test_v06_points_placed_on_axis():
    hot = QUESTIONS["hotter"]; assert any(l.x == 10.0 for l in hot.series[0].levels)          # 35 C fortnight = +10 C
    comp = QUESTIONS["compaction"]; assert any(l.twin.get("params",{}).get("porosity") == 0.3 for s in comp.series for l in s.levels)
```

- [ ] **Step 2: Run to verify failure.**
- [ ] **Step 3: Implement.** Build v0.7 series by mirroring the level construction of the corresponding `expNN_*` functions in `reproducibility/v07/composition.py` (same actions, models, levels); do not modify that file. v0.6 series (rain patterns, salty rain once with 0.05 M placed at EC via `nacl_mol_l_per_dS_m`, warmer soil only 0…+10 °C, porosity 0.3 point, CO-gate sample picker over the closed-gate samples listed in `reviews/v06/paired_regimes_summary.json["closed_gates"]` plus two open-gate samples) are chains with no external models.
- [ ] **Step 4: Run tests** → PASS; also `uv run pytest -q` (all old tests) → PASS.
- [ ] **Step 5: Commit** `feat(v0.8): 14 question specs as replayable chains`.

### Task 3: UQ core (design, evaluation, analysis)

**Files:**
- Create: `src/eq_biodt/v03/uq.py`
- Test: `tests/test_uq.py`

**Interfaces:**
- Consumes: `uq_inputs.Input`, `uq_inputs.inputs_for`, `questions.evaluate`.
- Produces:
  - `task_seed(key: str) -> int` (sha256 of key, first 8 bytes).
  - `design(inputs: list[Input], n: int, seed: int) -> np.ndarray` — Saltelli rows `[A, B, AB_1..AB_d, AB_G1..AB_G3]` (groups = statuses), scrambled `qmc.Sobol(2d, seed)` mapped through each input's inverse CDF; shape `(n*(d+2+3), d)`.
  - `evaluate_rows(spec, inputs, rows: np.ndarray) -> list[dict]` (one `evaluate` per row).
  - `analyse(inputs, n, outputs: list[dict], headline: list[str], boot: int = 500, seed: int = 0) -> dict` with per output: `median, q05, q25, q75, q95, refusal_fraction, mode ∈ {"variance","refusal_indicator","deduced"}, S1, ST (each {name: [est, lo, hi]}), group_S1, group_ST (by status), convergence {half_n_shift_frac: float, flagged: bool}`; contrasts computed per row (paired).
- Estimators: first order Saltelli 2010, total Jansen 1999; bootstrap percentile CIs with `seed`.

- [ ] **Step 1: Write the failing tests**

```python
def test_ishigami_indices():
    # analytic: S1 = [0.3139, 0.4424, 0.0], ST = [0.5576, 0.4424, 0.2437]
    res = analyse_function(ishigami, uniform_inputs(3, -np.pi, np.pi), n=4096, seed=1)
    assert np.allclose([res["S1"][k][0] for k in "x1 x2 x3".split()], [0.3139, 0.4424, 0.0], atol=0.03)
    assert np.allclose([res["ST"][k][0] for k in "x1 x2 x3".split()], [0.5576, 0.4424, 0.2437], atol=0.03)

def test_design_deterministic():
    a = design(INPUTS, 64, task_seed("q|s|0")); b = design(INPUTS, 64, task_seed("q|s|0"))
    assert np.array_equal(a, b)

def test_closed_gate_is_exact_zero():
    spec = level_of("carbon-fixation", sample="14PRr2")
    res = run_level(spec, n=32)
    out = res["outputs"]["c_fixed_upper_bound_g_m2_yr"]
    assert out["mode"] == "deduced" and out["q05"] == out["q95"] == 0.0

def test_refusal_mode_switch():
    res = analyse(INPUTS, 32, [{"y": None, "__refused__": "x"}]*10 + [{"y": 1.0}]*(32*(len(INPUTS)+5)-10), ["y"])
    assert res["y"]["refusal_fraction"] < 0.05 or res["y"]["mode"] == "refusal_indicator"
```

- [ ] **Step 2: Run to verify failure.**
- [ ] **Step 3: Implement** (`analyse_function`, `uniform_inputs`, `run_level`, `level_of` are thin test helpers in `uq.py`).
- [ ] **Step 4: Run tests** → PASS.
- [ ] **Step 5: Commit** `feat(v0.8): Saltelli/Jansen UQ with group indices, refusal handling, convergence`.

### Task 4: Task runner, Slurm array, gather, local mode

**Files:**
- Create: `reproducibility/v08/tasks.py`, `reproducibility/v08/uq_task.py`, `reproducibility/v08/gather.py`, `reproducibility/v08/slurm/uq_array.sbatch`, `reproducibility/v08/slurm/setup_unimatrix.sh`
- Test: `tests/test_uq_tasks.py`

**Interfaces:**
- Produces: `tasks.enumerate_tasks(n: int = 1024, chunk_rows: int = 2000) -> list[dict]` (`{key: "slug|series|level|chunk", slug, series, level, rows: [start, stop]}`), written to `reviews/v08/uq/tasks.json`; `uq_task.py --task-index I --tasks reviews/v08/uq/tasks.json --out DIR` writes `DIR/<key>.json` (raw outputs per row, input digests); `gather.py --raw DIR` writes `reviews/v08/uq/<slug>.json` via `uq.analyse`; `reproduce.py --local -j N` runs every task with multiprocessing then gather.
- sbatch: `#SBATCH --partition=debug --cpus-per-task=1 --mem=2G --array=0-<T-1>%250 --time=02:00:00 --output=/home/hohndor/eq-biodt-uq/logs/%A_%a.out`; runs `uv run --frozen python reproducibility/v08/uq_task.py --task-index $SLURM_ARRAY_TASK_ID ...` with `TMPDIR=$SLURM_TMPDIR`.

- [ ] **Step 1: Write the failing test**

```python
def test_chunked_equals_single(tmp_path):
    tasks = [t for t in enumerate_tasks(n=16, chunk_rows=40) if t["slug"] == "biochar" and t["level"] == 1]
    for i, _ in enumerate(tasks): run_task(tasks, i, tmp_path/"chunks")
    run_task(merge_chunks(tasks), 0, tmp_path/"single")
    assert gather_bytes(tmp_path/"chunks") == gather_bytes(tmp_path/"single")
```

- [ ] **Step 2–4:** fail, implement, pass.
- [ ] **Step 5: Commit** `feat(v0.8): UQ task array, Slurm script, gather and local mode`.

### Task 5: Run the UQ on unimatrix (pilot, then full)

**Files:** results only: `reviews/v08/uq/*.json`, `reviews/v08/uq/run_record.json` (job ids, nodes, wall time, cost per run, pilot comparison).

- [ ] **Step 1: Set up** — `ssh unimatrix01`; install `uv` to `~/.local/bin` (official installer), clone the twin repo at the branch commit into `/home/hohndor/eq-biodt-uq/twin` (git bundle from ws if GitHub access is missing), `uv sync --frozen`, install km reasoner if tests need it (UQ itself does not).
- [ ] **Step 2: Pilot** — tasks of `biochar` level 1 and `salty-water` top level only; submit array; when done `rsync` raw results to ws; gather; compare with `reproduce.py --local` for the same tasks → identical bytes. Record cost per run.
- [ ] **Step 3: Full run** — submit all tasks (throttle `%250`); monitor with `squeue`; resubmit failed indices only.
- [ ] **Step 4: Gather on ws**; convergence report: list outputs with `flagged: true`; if any, rerun those levels at N = 2048 and note it in `run_record.json`.
- [ ] **Step 5: Commit** results `data(v0.8): UQ results (N=1024, unimatrix job <id>)`.

### Task 6: Question files with plain and scientific texts

**Files:**
- Create: `reproducibility/v08/questions.py`, `reproducibility/v08/texts.py` (templates)
- Test: `tests/test_question_files.py`
- Output: `reviews/v08/questions/<slug>.json`, `reviews/v08/questions/index.json`

**Interfaces:**
- Produces JSON schema (version `"question-1"`): `{slug, theme, title_plain, title_science, level_kind, x_label_plain, x_label_science, baseline_label, series:[{id,label_plain,label_science, levels:[{label, x, scene:{water_mm, theta, salt_mol_l, temp_c, porosity, active_days, activity_by_day[14]}, outcomes:{<output>:{value, q05,q25,q75,q95, unit, unit_plain, provenance, mode, refusal_fraction}}, trajectories:{baseline:{W,N,psi,A,J_h2,J_co}, intervened:{…}}}]}], answer_plain, answer_science, chain:[{id, name_plain, name_science, kind}], limits_plain, limits_science, methods:{models:[{id,name,equation_tex,constants}], data:[…], gates:[…], admissibility}, uncertainty:{S1,ST,group_S1,group_ST,inputs:[{name,dist,low,value,high,status,source}], convergence, oat:{…}}, reproduce:{files, commands, digests, paper_section}}`.
- Scene and trajectories come from the central run of each level (`evaluate` with no overrides, `keep_twin_result=True`).

- [ ] **Step 1: Write the failing tests**

```python
def test_schema_and_count(): assert len(load_index()) == 14 and all(validate(q) for q in load_all())
def test_plain_text_rules():
    for q in load_all():
        txt = " ".join([q["answer_plain"], q["limits_plain"], q["title_plain"]])
        assert not JARGON.search(txt), q["slug"]
        assert "the models themselves could also be wrong" in q["limits_plain"]
    assert "at most" in load("carbon-fixation")["answer_plain"]
def test_plain_matches_science():
    for q in load_all():
        assert plain_number(q) == round_plain(science_central(q))   # same core number
def test_refused_level_rendering_fields():
    lvl = top_level(load("salty-water"), series="daily-irrigation")
    assert lvl["outcomes"]["h2_uptake_contrast_umol_m2"]["refusal_fraction"] > 0
```

- [ ] **Step 2–4:** fail, implement, pass. Template rules: interval spanning 0 → "may go up or down"; bound → "at most"; refusal ≥ 5 % → "the model cannot tell for such <condition>".
- [ ] **Step 5: Commit** `feat(v0.8): question files with plain and scientific texts`.

### Task 7: Paper, verifier, reproduce, release 0.8.0

**Files:** `paper/main.tex`, `paper/supplementary.tex`, `paper/figures/v08_sensitivity.pdf` (+ `reproducibility/v08/figures.py`), `paper/tables/v08_uq.tex`, `reproducibility/v06/verify_manuscript.py` (load v0.8 outputs), `reproducibility/v08/reproduce.py`, `pyproject.toml`/`CITATION.cff`/`ro-crate-metadata.json`/`src/eq_biodt/v03/__init__.py` → 0.8.0, `README.md`, `CHANGELOG.md`, `docs/causal-ontological-theory.md` (short UQ note).

- [ ] **Step 1:** Methods paragraph (inputs, distributions, Saltelli/Jansen, N, bootstrap, convergence, refusal rule, scope = parameter uncertainty only); Results: where v0.7 text reports one-at-a-time ranges for headline numbers, report joint 90 % intervals and the top two Sobol' contributors; supplement: parameter registry table, convergence summary, sensitivity figure (ST per question, grouped by status).
- [ ] **Step 2:** `uv run python reproducibility/v06/verify_manuscript.py` → `manuscript numbers verified`; `uv run pytest -q` → all pass; PDFs build with 0 undefined references.
- [ ] **Step 3:** `reproducibility/v08/reproduce.py` (questions + figures + digests; UQ gather from stored raw results; `--local` recomputes) run twice → identical sha256 of `reviews/v08/**/*.json`, table, figure.
- [ ] **Step 4: Commit** `release: EQ-BioDT 0.8.0 (question pages data, full uncertainty analysis)`; tag `v0.8.0` only after Task 13 passes.

### Task 8: Web export of question files

**Files:** `scripts/web/export_web_bundle.py`, `tests/test_web_bundle.py`.

- [ ] **Step 1: Failing test** — `test_question_files_exported`: after export, `<out>/questions/index.json` lists 14 slugs and every file's sha256 equals `reviews/v08/questions/<slug>.json`; manifest records them.
- [ ] **Step 2–4:** fail, implement (copy files verbatim into `<out>/questions/`, add to manifest), pass.
- [ ] **Step 5: Commit** `feat(web): export question files`.

### Task 9: Website routing, catalogue, redirects

**Files:** `web:src/App.jsx`, create `web:src/pages/twin/Catalogue.jsx`, `web:src/pages/twin/QuestionPage.jsx` (shell), `web:src/pages/twin/HowPage.jsx` (shell), `web:src/pages/twin/TryPage.jsx`, `web:src/components/twin08/ModeSwitch.jsx`, `web:src/components/twin08/data.js`, `web:src/components/twin08/redirects.js`; modify `web:src/pages/Twin.jsx` (becomes Catalogue route or is replaced); Test: `web:e2e/test_site.py` (Playwright, Python).

**Interfaces:**
- Routes: `/model` (catalogue), `/model/q/:slug`, `/model/q/:slug/science`, `/model/how/:topic`, `/model/how/:topic/science`, `/model/try`, `/model/guide` (existing `TwinGuide`), `/twin` → `/model`.
- `data.js`: `loadIndex(): Promise<Index>`, `loadQuestion(slug): Promise<Question>` from `/twin-data/eq_biodt_v06/questions/`; version gate accepts `/^0\.[678]\./`.
- `ModeSwitch({mode, onChange})`; mode persisted in localStorage key `eqbiodt.twin.mode` (try/catch), URL is authoritative.
- `redirects.js`: `redirectFor(location) -> string | null` mapping `?view=research` → `/model/how/mechanisms/science`, `?view=whatif` → `/model`, section ids: `footprint→/model/how/samples`, `genomes→/model/how/deduction`, `coupling→/model/how/community`, `graph`/`identification`/`versions`→`/model/how/mechanisms`, `run→/model/try`, `transfer→/model/how/other-deserts`, `ontology→/model/how/definitions`, `composition→/model` (science mode), `research→/model/how/mechanisms/science`.

- [ ] **Step 1: Failing e2e test** (`test_catalogue_lists_14_questions_in_5_themes`, `test_mode_switch_keeps_question`, `test_redirects` for the mapping above) against `vite preview`.
- [ ] **Step 2–4:** fail, implement, pass (`npx eslint src/pages src/components/twin08` 0 errors; `npm run build`).
- [ ] **Step 5: Commit** `feat(web): question catalogue, routes, mode switch, redirects`.

### Task 10: Question page core (scene, slider, bars, chain, limits)

**Files:** create `web:src/components/twin08/{QuestionCore.jsx, Scene.jsx, LevelSlider.jsx, OutcomeBars.jsx, ChainStrip.jsx}`; reuse/extend `web:src/components/twin06/SoilScene.jsx` (props `water_mm, theta, salt_mol_l, temp_c, porosity, activity`, reduced motion); test in `web:e2e/test_site.py`.

**Interfaces:** `QuestionCore({question, mode})` renders sections in this order with ids `#question #answer #picture #how #limits`; `LevelSlider` snaps to level indices (numeric, categorical, sample picker); `OutcomeBars` shows `value` (plain: rounded, `unit_plain`; science: value + `[q05, q95]`), source label from `provenance`, refusal text when `refusal_fraction ≥ 0.05`.

- [ ] **Step 1: Failing e2e tests**: `test_all_questions_render_both_modes` (28 pages, 1280 and 390 px, 0 page errors, no horizontal scroll), `test_same_core_number` (plain and science `#answer` contain the same rounded headline from the JSON), `test_salty_top_level_shows_refusal`.
- [ ] **Step 2–4:** fail, implement, pass.
- [ ] **Step 5: Commit** `feat(web): question page core with slider and soil scene`.

### Task 11: Science extras (methods, uncertainty, trajectories, reproduce)

**Files:** create `web:src/components/twin08/{Methods.jsx, Uncertainty.jsx, Trajectories.jsx, Reproduce.jsx}`; charts via existing `charts.jsx`/`MiniChart` style (inline SVG).

**Interfaces:** sections after the core with ids `#methods #uncertainty #trajectories #reproduce`, collapsible except `#uncertainty`. `Uncertainty` shows: curve across levels with 50 % and 90 % bands, ST/S1 bars with CI whiskers (top 8 + "other"), group-by-status bar, convergence note, refusal fraction, OAT view under a disclosure, inputs table (filter by status).

- [ ] **Step 1: Failing e2e test** `test_science_sections_present` (all four ids on every science page; none on plain pages) and `test_uncertainty_scope_sentence` (text "excludes … model structure" present).
- [ ] **Step 2–4:** fail, implement, pass.
- [ ] **Step 5: Commit** `feat(web): science sections — methods, uncertainty, trajectories, reproduce`.

### Task 12: "How the twin works" pages

**Files:** create `web:src/pages/twin/how/{samples,deduction,community,mechanisms,other-deserts,definitions}.jsx` reusing `Footprint`, `Genomes`, `Coupling`, `GraphView`+`Identification`+`Versions`, `Transfer`, `OntologyView`; each topic exports `{Plain, Science}` with the same core order (question → answer → picture → how we know → limits), plain intros written for the public.

- [ ] **Step 1: Failing e2e test** `test_how_pages_render_both_modes` (12 pages, 0 errors, no horizontal scroll) and jargon check on plain how-pages.
- [ ] **Step 2–4:** fail, implement, pass.
- [ ] **Step 5: Commit** `feat(web): how-the-twin-works pages (plain/science)`.

### Task 13: Review, release, deploy

- [ ] **Step 1:** Independent jargon/contradiction review (fresh subagent) of all plain pages against the science pages and the JSONs; fix findings.
- [ ] **Step 2:** Twin: push `v0.8-questions`, fast-forward `main`, tag `v0.8.0`, push (SSH). Re-export bundle from the tagged commit (`tree_dirty: false`).
- [ ] **Step 3:** Website: full e2e suite on `vite preview` (ws) → pass; commit; push `twin-v08-questions` and to `analysis/16s-statistical-review`.
- [ ] **Step 4:** Deploy release `twin08-YYYYMMDD` with the overlay procedure (`build*.sh`/`accept*.sh` adapted from `~/tmp/twin06e-deploy/`, `deploy_vm.py stage`), add `/model/q/rain-pattern`, `/model/q/rain-pattern/science`, `/model/how/mechanisms`, `/model/try` to the acceptance URL list; the staged port is VM-local, so run the full Playwright suite on ws against `vite preview` of the exact `dist/` that went into the overlay tarball (same sha256) before `activate`; after activation run it against https://rubalkhali.science and roll back immediately on any failure; record the release and rollback command in `website/README.md`.
