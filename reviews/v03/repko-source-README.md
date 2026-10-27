# repko_clean — corrected per-representative KO table

> ## STATUS: DONE (2026-08-12). `results/rep_ko_v2.tsv` built and consumed by §2.2.
> **Result:** 5,438 rows / 0 unannotated (was 168) / 4,580 identical to original / 690 reps shed at least one
> spurious KO (8,134 total). Downstream: §2.2 recomputed on v2 and shipped — Overleaf `3e0ff2f`,
> GitHub main `e01203d` (see `notes/CORRECTIONS.md` 2026-08-12 entry for the full v1→v2 diff table).
> **Runs (all COMPLETED):** step 1a Bakta `50376640` (124/124); step 1b eggNOG `50378049` (123/124, one caught
> by the `## Rate:` trailer guard) + resubmit `50379295`; step 2 `repko_2_submap.py` (1,010/1,010 corrected);
> step 3 `repko_3_build.py`.
> **Untouched by policy:** PI's `04_functional/refs/rep_ko.tsv`; `reproducibility/claims.tsv` — see
> `CLAUDE.md` "Working records and precedence".
> **KEY LESSON for future annotation jobs:** eggNOG-mapper needs `--dbmem` on WekaFS — without it, emapper
> random-reads the 41 GB `eggnog.db` and concurrent tasks deadlock in D-state on `wekafs_lock_helper`.

Fixes three defects in the per-representative KO table used by all of Part 2 (§2.1–§2.4).
Source of the problem: `scripts/eco_probes/gather_rep_ko.py` (hohndor, 2026-06-11) →
`04_functional/refs/rep_ko.tsv`. **That original file is never modified**; we build
`results/rep_ko_v2.tsv` alongside it. M = /ibex/scratch/projects/c2014/rob/metagenomics.

## The defects (all VERIFIED 2026-08-11, full 5,438 reconciliation)
| reps | state | cause |
|---|---|---|
| 4,304 | annotated, contig-identical to their raw bin | ✅ correct, untouched |
| 966 | annotated `_sub` reps, inherited the PARENT's annotation | ⚠️ superset — 915 retain ≥95% of parent, **51 retain <95%** (worst 25%). Measured on one rep: 422 of 1,521 KOs (28%) spurious |
| 44 | unannotated — `48PRr1__SemiBin_N`, no `single_`/`self_` token | ❌ resolver splits on the token and fails; **the annotation exists** |
| 80 | unannotated — raw bin was LQ pre-refinement, never on the HQ+MQ list | ❌ refinement promoted them to HQ/MQ afterwards; nothing re-annotated them |
| 44 | unannotated `_sub` whose PARENT was also never annotated (0/44 on the list) | ❌ parent-fallback had nothing to fall back to |
| **5,438** | | 168 unannotated total; 44+44 are disjoint sets (0 overlap) |

Root cause of the missing ones: annotation ran on the **raw, pre-refinement** HQ+MQ bins
(`03_annotation/mag_list.txt`, 15,469 rows) while the catalogue is built from **post-refinement**
bins — refinement both created new bins (`_sub`) and promoted LQ bins.

## Pipeline
1. **`repko_1_annotate.sh`** — Bakta + eggNOG for the **124** never-annotated raw bins
   (`results/annotate_targets.tsv`: 80 promoted_LQ + 44 parent_of_sub; all single_sample, all
   FASTAs present). Writes into the per-sample tree at the canonical paths
   (`<sample>/wgs/annotation/<mode>/{bakta,eggnog}/<bin>/`) so everything stays uniformly organised.
   **Identical stack to the originals** (`03_bakta_array.sh` / `04_eggnog_array.sh`):
   bakta db `/ibex/reference/KSL/bakta/6.0/db` `--skip-plot --force`;
   emapper data_dir `/ibex/reference/KSL/eggnog-mapper/2.1.12` `--itype proteins --override`.
   Both tools run in ONE task and **fail loudly** — the original eggNOG array `exit 0`-ed when Bakta
   output was absent, which is precisely how bins went silently unannotated.
2. **`repko_2_submap.py`** — correct KO sets for all **1,010 `_sub`** reps by keeping only genes on
   contigs the sub bin retained. Bakta renames contigs to `contig_N` but preserves input order, so
   `contig_N → original name` is recovered from the raw bin and **verified position-by-position by
   sequence length**; any bin failing that check is REFUSED, never guessed. No tool is re-run —
   same gene calls, same versions. Run again after step 1 to pick up the 44 new parents.
3. **`repko_3_build.py`** — rebuild all 5,438 rows: subset KOs for `_sub`, fixed resolver for the
   48PRr1 naming, newly annotated bins picked up automatically → `results/rep_ko_v2.tsv`, plus a
   diff against the original (identical / gained / spurious-removed / still-unannotated).

## Consumers of rep_ko.tsv (blast radius)
Ours: `06_atmospheric/lifeonair_clean/scripts/lifeonair_{1,2}_*`, `06_metabolism/uv_survival_clean` (§2.3).
PI's: `scripts/eco_probes/{aggregate_e6.py,aggregate_e6_comp.py,prep_e6_refs.sh}` — left on the
original table; his call, registered in `notes/CORRECTIONS.md` rather than changed silently.

## After the rebuild
Recompute §2.2 once (`lifeonair_2_type_guild.sh` + `lifeonair_3_abundance.py`), let §2.3 run on
v2, and register the old → new numbers in `notes/CORRECTIONS.md`.
