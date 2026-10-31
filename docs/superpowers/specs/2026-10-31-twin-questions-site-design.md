# EQ-BioDT question pages: plain and scientific, with full uncertainty analysis

Status: design agreed with R. Hoehndorf on 2026-10-31 (layout A, picture style A, full uncertainty analysis);
approved 2026-10-31 with one change (UQ compute on unimatrix). Implementation runs on ws (`/mnt/data1/agents/projects/eq-biodt-v07/{twin,web}`),
never on the laptop.

## Intent

The public twin page (rubalkhali.science/model) has grown into one long page with a "What if…?" view and a "Research
view". It is cluttered and things are hard to find; several questions the twin answers have no picture, and the two
views are not parallel. Robert wants **one page per question the twin answers**, in **two versions**: a plain,
public-facing page and a scientific page for experts. The two share the same core (same question, same answer, same
picture, same order); the scientific page adds **methods** (how the answer was determined), **exact numbers**,
**uncertainty quantification** and careful hedging. Uncertainty becomes a **full analysis** (joint sampling of all
uncertain inputs and global sensitivity), replacing the one-at-a-time ranges as the primary uncertainty statement.

Success: a member of the public finds their question from the start page in one click and understands the answer from
the picture and two sentences; an expert switches to the scientific page of the same question and finds the method,
the joint uncertainty interval, which inputs drive it, and how to reproduce it.

## Goals

1. Twin repo (release 0.8.0): an uncertainty-quantification (UQ) module and analyses for every question; per-question
   data files for the website; paper update.
2. Website: a question catalogue start page, 14 question pages × 2 versions, "How the twin works" pages × 2 versions,
   the live twin and the guide kept; old URLs keep working.
3. Everything reproducible, verified, built and tested on ws; deployed with a pre-deploy browser check.

## Non-goals and anti-goals

- No new process science: the twin's mechanisms, the external models and their central values stay as in 0.7.0
  (except where UQ needs ranges that do not exist yet; those are added as labelled ranges, see below).
- No claims beyond the outputs: plant yields and carbon remain the external models' numbers.
- Anti-goal: a plain page that hedges in every sentence, or a science page that says something different from the
  plain page. The core must agree; the science page refines, it does not contradict.
- Anti-goal: hand-typed numbers in the site. Every number comes from the exported data.
- Anti-goal: breaking the live site again. Every deploy is preceded by a local build + headless-browser check of all
  pages on ws (the twin07 incident: a version gate rejected the new bundle).

## Part 1 — Twin repo: uncertainty quantification (0.8.0)

### Uncertain inputs

- **Literature constants**: every entry of `data/external/composition_constants.json` used by a chain, with its
  low/central/high. Distribution: triangular(low, central, high). Constants with `status: modeller_assumption` or
  `source_not_checked` keep the same distribution but are flagged in all outputs.
- **Twin parameters** (new registry `data/external/twin_parameter_ranges.json`, same schema and status labels as the
  constants registry): at least `k_h2_s`, `k_co_s` (log-uniform over one order of magnitude around the reference
  0.001 s⁻¹, status `modeller_assumption`), `psi_entry_mpa`, `retention_b`, `wet_activation_mpa`,
  `dry_threshold_mpa`, `activation_shape`, `q10`, `d_h2_m2_s`, `d_co_m2_s`, `drainage_fraction_day`,
  `field_capacity_fraction`, `evaporation` scaling and porosity where porosity is not the intervened variable.
  Ranges come from the paper's existing sensitivity grid / scenario ensemble where it exists (reuse those bounds and
  cite the paper section), otherwise from literature with source, otherwise labelled modeller assumption. Parameters
  that an experiment intervenes on are not sampled in that experiment.
- Admissibility and validity still apply to every sample: a sample whose chain is refused (validity domain, twin
  abstention above 0.1 M NaCl) is counted and reported as a refusal fraction, never silently dropped.

### Analyses (per question, per level, per headline output; also for paired contrasts vs baseline)

- **Joint propagation**: scrambled Sobol' quasi-random sampling (`scipy.stats.qmc.Sobol`, fixed seed), N ≥ 1024
  samples; report median, 5–95 % and 25–75 % intervals, refusal fraction. Contrasts (intervened − baseline) are
  computed per sample on the same draw (paired), not as differences of marginals.
- **Global sensitivity**: first-order and total Sobol' indices (`scipy.stats.sobol_indices` or an equivalent
  Saltelli/Jansen estimator) with bootstrap confidence intervals; report the top contributors per output.
- **Convergence check**: indices and intervals at N/2 vs N; flag outputs whose 90 % interval endpoints move by more
  than 5 % of their width.
- **One-at-a-time ranges** from 0.7.0 are kept as a secondary, clearly labelled "sensitivity to each constant" view.
- Deductive nulls stay exact: where a gate is closed the output is 0 in every sample and is reported as such
  (interval [0, 0], provenance deduced), not as a sampled quantity.

### Questions covered (14) and their data

The question set merges the v0.6 public what-ifs and the v0.7 chained experiments (one question per action):

| Theme | Question (plain) | Levels (slider) | Sources |
|---|---|---|---|
| Water & rain | One heavy rain, four showers or daily drizzle? | none, single, four, daily (categorical) | v0.6 paired regimes / treatment versions |
| | Drip or flood irrigation? | every 1, 3, 10 days | v0.7 drip_vs_flood |
| | What if cloud seeding adds rain? | +0…+20 % | v0.7 cloud_seeding |
| Salt | What if the water is salty? | EC 0.5…8 dS/m (+ the v0.6 0.05 M salty-rain case as a point) | v0.6 salty + v0.7 brackish_irrigation |
| | Can a desert bacterium help alfalfa in salty water? | EC 0.5…12 dS/m | v0.7 inoculated_alfalfa_salinity |
| Heat & shade | What if it gets hotter? | +0…+10 °C (v0.6 hot fortnight 35 °C and v0.7 warming as points on one axis) | v0.6 hot + v0.7 warming |
| | Mulch or shade nets? | cover 0…1 | v0.7 mulch_shade |
| | What if we build a solar park? | cover 0…1 | v0.7 solar_park |
| Soil | What if vehicles compact the sand? | passes 0…50 (v0.6 porosity 0.3 as a point) | v0.6 compact + v0.7 compaction |
| | What if we add biochar? | 0…5 % w/w | v0.7 biochar |
| Microbes & air | What if the soil lacks the CO genes? | sample picker (closed vs open gate) | v0.6 closed CO gates |
| | What if we add microbes? | share 0…5 %, gate/proportional | v0.7 inoculation |
| | What if air hydrogen rises? | 1…10× | v0.7 h2_leakage |
| | How much carbon can the hydrogen eaters fix? | sample picker | v0.7 carbon_fixation_bound |

Where a merged question puts v0.6 and v0.7 levels on one axis, both must be computed with the same baseline and
twin version (rerun the v0.6 cases through the 0.8 pipeline; do not mix stored numbers from different baselines).

### Outputs

- `reproducibility/v08/uq.py` (analyses), `reproducibility/v08/questions.py` (one JSON per question:
  `reviews/v08/questions/<slug>.json` with: plain and scientific texts generated from templates filled with numbers,
  levels, per-level scene state for the picture (water mm and θ, salt mol/L, temperature, porosity, active days,
  activity per day), headline outputs with central value, joint intervals, OAT ranges, Sobol indices, refusal
  fraction, chain/methods description, constants and parameters used with status, file paths and paper section),
  `reproducibility/v08/reproduce.py` with digests, run twice → byte-identical.
- Paper: Methods paragraph on UQ; Results numbers switched to joint intervals where the text reports ranges
  (one-at-a-time ranges remain in the supplement); a supplementary section with the parameter registry, the
  convergence check and a sensitivity figure. Verifier extended; all tests pass; PDFs build without warnings.
- Web export: `scripts/web/export_web_bundle.py` writes the question files next to the bundle (or into it).


### Decisions settled after the builder check (supervisor's calls, marked as such)

- **Headline outputs per question** (*supervisor's call*). UQ (joint intervals + Sobol' indices) is computed for these
  outputs at every level; all other outputs get joint intervals only.

  | Question | Headline outputs | Level kind |
  |---|---|---|
  | Rain pattern | 14-day H₂ uptake; contrast vs dry | categorical |
  | Drip or flood | 14-day H₂ uptake contrast vs rain-only | numeric (interval days) |
  | Cloud seeding | H₂ uptake contrast | numeric |
  | Salty water | H₂ uptake contrast; crop relative yield (4 crops) | numeric, two series |
  | Alfalfa + SA187 | inoculated relative yield (fresh), per root-zone model A/B/C | numeric, three series |
  | Hotter | H₂ uptake contrast; active days | numeric, two series |
  | Mulch / shade | H₂ uptake contrast | numeric, two series |
  | Solar park | H₂ uptake contrast | numeric |
  | Compaction | H₂ uptake contrast; porosity | numeric |
  | Biochar | H₂ uptake contrast | numeric |
  | Missing CO genes | CO uptake (exact 0 where the gate is closed) | sample picker |
  | Adding microbes | CO uptake; contrast | numeric (share) × coupling mode |
  | Air hydrogen | H₂ uptake contrast | numeric |
  | Carbon | upper bound on C fixed per year | sample picker |

  A **sample** (metagenome) and the **coupling mode** (gate/proportional) are levels or display toggles, never sampled
  inputs. A closed gate gives an exact 0 with provenance `deduced` in every draw (interval [0, 0]). The carbon output is
  an upper bound; its interval is the parameter uncertainty of the bound, labelled "upper bound" on both pages.
- **What the intervals mean** (*supervisor's call*): they are **parameter-uncertainty intervals** over the constants
  and twin parameters only. They exclude the external models' structural or predictive error and the twin's
  structural error. Both pages say so (plain: "the range shows how much the answer moves when uncertain numbers in the
  models are varied; the models themselves could also be wrong"); the paper's Methods says the same.
- **Merged axes** (*supervisor's call*): a question can have several **series** on one axis when the actions differ;
  every series is computed by the 0.8 pipeline from the same baseline. Salty water: series "salty rain once"
  (single 12 mm pulse; the v0.6 0.05 M case sits on it, converted to EC with the registry factor) and "daily
  irrigation with salty water" (v0.7). Hotter: series "warmer soil only" (v0.6 style, temperature +ΔT, the
  35 °C fortnight is +10 °C) and "warmer climate" (v0.7, with evaporation scaling). Compaction: one series over
  vehicle passes, with the resulting porosity shown on a secondary scale; the v0.6 porosity-0.3 case is a marked
  point at the porosity it implies, recomputed. Rain pattern: categorical levels. The question JSON has one schema:
  `series[] → levels[]`.
- **Compute** (*Robert's decision, 2026-10-31: use unimatrix where possible*): the UQ runs execute on the
  **unimatrix** Slurm cluster (partition `debug`; nodes node002–node008, about 280 schedulable CPUs; node005/node008
  are CPU-only in Slurm). Layout: a clone of the twin repo at the release commit and a `uv` environment under
  `/home/hohndor/eq-biodt-uq/` (NFS, shared by all nodes; small files only); per-task scratch in
  `$SLURM_TMPDIR` (node-local `/scratch`); one **Slurm array job** over (question, series, level) tasks, each task
  writing one small JSON of raw sample outputs; a gather step on ws merges them (after `rsync` back) into the
  question files. Each task's Sobol' draw and seed are fixed by its task key, so results do not depend on node,
  task order or worker count. Unimatrix only runs the sampling; tests, reproduction checks, LaTeX, website builds
  and agents stay on ws. A small pilot array (one question, two levels) runs first to measure cost per run on the
  cluster and to check that the gathered result equals the same computation on ws byte for byte. At most 7 levels
  per series; Saltelli design with N = 1024 base samples (raise only if the convergence check fails). `reproduce.py`
  also supports a single-machine mode (multiprocessing) that yields identical files, so the analysis does not
  depend on the cluster.
- **Refused samples** (*supervisor's call*): report the refusal fraction per level. If it is < 5 %, compute intervals
  and Sobol' indices on the admitted samples and flag the result. If it is ≥ 5 %, report the joint interval on the
  admitted samples plus Sobol' indices of the refusal indicator (which inputs push the chain out of its valid
  domain), not variance indices of the output.
- **Assumption-driven uncertainty** (*supervisor's call*): besides the per-input indices, report **group Sobol'
  indices by status** (`verified_in_source`, `modeller_assumption`, `source_not_checked`), so the science page states
  how much of each interval comes from assumed ranges. No separate "sourced-only" interval.
- **Slider** (*supervisor's call*): the slider snaps to the computed levels; nothing is interpolated between levels.
- **Old URLs** (*supervisor's call*): redirects are client-side in the React app (hashes never reach the server); the
  builder maps every existing section id in `Twin.jsx` and the `?view=` parameter; `/twin` keeps its existing route.
- **Plain-language texts** (*supervisor's call*): the builder writes them as templates; the supervisor reviews them,
  and an independent reviewer checks for jargon against a word list (Greek letters, model ids, "do()", "Sobol",
  "ontology", "MAG", unit symbols beyond °C, %, mm, g, m²), and for contradictions between plain and science
  versions. Where an interval spans zero, the plain text says "may go up or down"; where results are refused, it
  says "the model cannot tell for such salty water".

## Part 2 — Website

### Information architecture

- `/model` — start page: one-line intro, a **Plain | Scientific** switch (remembered per viewer), and the question
  catalogue grouped in the five themes above. Each tile: the question, a small still of its scene, the headline
  number. A sixth group, **How the twin works**, links the background pages. Search/filter box over question titles.
- `/model/q/<slug>` (plain) and `/model/q/<slug>/science` (scientific); the switch keeps the question and the open
  section. Prev/next question links; breadcrumb back to the catalogue.
- `/model/how/<topic>` and `/model/how/<topic>/science` for background: where the samples come from (footprint);
  which genomes can do what (ontology deduction, candidates); how community composition sets rates (coupling); what can
  be changed and what follows (causal–ontological graph, admissibility, identification with deductive closure);
  does it work in other deserts (transfer tests); definitions (ontology mappings).
- `/model/try` — the live twin (current Pyodide paired-regime runner), unchanged in function.
- `/model/guide` — the existing guide and glossary, kept; info buttons on all pages link into it.
- Old links keep working: `/model?view=research`, `/model?view=whatif`, `/model#<section id>` and `/twin` redirect
  to the matching new page (map every existing section id).

### Question page — shared core (same order on both versions)

1. The question (plain wording on both; the science page adds a one-line technical subtitle).
2. The short answer: the same headline number(s). Plain: rounded, everyday units, one or two sentences. Science:
   exact value with joint 90 % interval and units.
3. The picture: **slider + animated soil scene** (style A): a dose slider (or sample picker / categorical steps),
   the sand cross-section showing water, salt, heat, compaction and active microbes for the selected level, and
   outcome bars next to it labelled with their source ("the twin", "another model: …", "follows by deduction").
   Science: the same scene, plus the interval under each bar.
4. How we worked it out: the chain as a strip of steps (plain: everyday names; science: model ids).
5. What we can't say: plain — one short box; science — the full, precise limits.

### Science page — additional sections (after the core)

- **Methods**: the models in the chain with their equations, data used (e.g. the SA187 field trial and its
  digitisation), the ontology classes that gate the result, why the intervention is admissible.
- **Uncertainty**: joint intervals per level (curve with 50 % and 90 % bands), Sobol' sensitivity bar chart (first
  order and total, with CIs), convergence note, refusal fraction, the one-at-a-time view as secondary, and the table
  of every sampled input with its distribution, source and status.
- **Trajectories**: the twin's 14-day W, salt, ψ, activity and uptake for the selected level (baseline vs intervened).
- **Reproduce**: result files, commands, code paths, digests, paper section.

### Visual and interaction quality

- One idea per screen; generous spacing; no more than one primary visual above the fold; secondary material in
  collapsible sections on the science page.
- Reuse and extend the existing `SoilScene`, `MiniChart`, `Info`, `glossary.js` and styles; light/dark themes,
  keyboard access, `prefers-reduced-motion` (scene animation off), phone width 390 px without horizontal page scroll.
- All texts on plain pages pass a jargon check (no Greek letters, model ids or unit symbols beyond common ones).

## Constraints

- Compute: UQ sampling on unimatrix (Slurm array); everything else (tests, reproductions, LaTeX, npm, browsers,
  agents) on ws; nothing heavy on the laptop.
- No `/tmp`; `TMPDIR=$HOME/tmp` on ws.
- Commit messages end with the Claude-Session line; push only after checks pass.
- Deploy via the existing image overlay + `deploy_vm.py stage/activate`, with a new release name; staged release is
  checked with a headless browser on all page types before activation; rollback command recorded in the website
  README.

## Left to the builder

Component structure, routing library details, exact Sobol' sample size above the minimum (subject to the convergence
check), plot styling within the existing design language, the search box implementation, and how the question JSONs
are loaded (one bundle or lazy per question).
