"""Tests for main-text and supplementary v0.8 figures and tables."""
from __future__ import annotations

from pathlib import Path
import re
import sys

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from reproducibility.v08 import figures, paper_figures  # noqa: E402

FIG_DIR = ROOT / "paper" / "figures"
TAB_DIR = ROOT / "paper" / "tables"

FORBIDDEN_IDENTIFIERS = (
    "k_h2_s",
    "k_co_s",
    "retention_b",
    "psi_entry_mpa",
    "dry_threshold_mpa",
    "wet_activation_mpa",
    "activation_shape",
    "drainage_fraction_day",
    "field_capacity_fraction",
    "d_h2_m2_s",
    "d_co_m2_s",
    "cloud_seeding_rain_increase_frac",
    "solar_panel_evap_reduction_frac",
    "solar_panel_soil_temp_reduction_c",
    "evap_scaling_per_c",
    "sa187_mock_ratio_fresh",
    "sa187_mock_ratio_dry",
    "sa187_benefit_b",
    "modeller_assumption",
    "verified_in_source",
    "source_not_checked",
    "C_fitted_fresh",
    "A_fao_fixed_leaching",
    "B_twin_physics",
    "h2_uptake_change_pct",
    "co_uptake_change_pct",
    "rain-pattern",
    "drip-or-flood",
    "cloud-seeding",
    "salty-water",
    "alfalfa-bacterium",
    "mulch-shade",
    "solar-park",
    "missing-co-genes",
    "adding-microbes",
    "air-hydrogen",
    "carbon-fixation",
)

# Matches any snake_case identifier outside LaTeX math ($...$)
MATH_SPAN_RE = re.compile(r"\$[^$]*\$")
SNAKE_ID_RE = re.compile(r"\b[A-Za-z0-9]+_[A-Za-z0-9_]+\b")


def _strip_math_and_comments(tex_or_label: str) -> str:
    lines = [
        re.sub(r"(?<!\\)%.*$", "", ln)
        for ln in tex_or_label.splitlines()
    ]
    no_comments = "\n".join(lines)
    return MATH_SPAN_RE.sub(" ", no_comments)


def test_paper_figures_script_runs_and_creates_outputs() -> None:
    paper_figures.main()
    figures.main()

    expected_files = [
        FIG_DIR / "v08_composition_main.pdf",
        FIG_DIR / "v08_composition_main.png",
        FIG_DIR / "v08_answers.pdf",
        FIG_DIR / "v08_answers.png",
        FIG_DIR / "v08_sensitivity.pdf",
        FIG_DIR / "v08_sensitivity.png",
        TAB_DIR / "v08_answers.tex",
        TAB_DIR / "v08_uq.tex",
        TAB_DIR / "v08_parameters.tex",
    ]
    for p in expected_files:
        assert p.exists() and p.stat().st_size > 0, f"Missing or empty output: {p}"


def test_answers_table_has_14_data_rows() -> None:
    tex_path = TAB_DIR / "v08_answers.tex"
    content = tex_path.read_text(encoding="utf-8")
    assert "\\begin{tabularx}{\\textwidth}" in content
    assert "\\toprule" in content and "\\midrule" in content and "\\bottomrule" in content

    body = content.split("\\midrule", 1)[1].split("\\bottomrule", 1)[0]
    data_rows = [ln.strip() for ln in body.splitlines() if "&" in ln]
    assert len(data_rows) == 14, f"Expected 14 data rows, got {len(data_rows)}"
    for row in data_rows:
        cols = [c.strip() for c in row.rstrip("\\").split("&")]
        assert len(cols) == 6, f"Expected 6 columns in row: {row!r}"


def test_no_underscore_identifiers_in_tables_or_figure_labels() -> None:
    # 1. Check generated LaTeX tables
    for tab_name in ("v08_answers.tex", "v08_uq.tex", "v08_parameters.tex"):
        raw = (TAB_DIR / tab_name).read_text(encoding="utf-8")
        assert "\\texttt" not in raw, f"{tab_name} contains \\texttt identifier"
        assert "\\_" not in raw, f"{tab_name} contains escaped underscore \\_"
        for ident in FORBIDDEN_IDENTIFIERS:
            assert ident not in raw, f"{tab_name} contains raw identifier {ident!r}"
        stripped = _strip_math_and_comments(raw)
        assert "_" not in stripped, f"{tab_name} contains underscore outside math mode: {stripped!r}"
        assert not SNAKE_ID_RE.search(stripped), f"{tab_name} contains snake_case token"

    # 2. Check all readable label dictionaries and figure label specs
    all_labels: list[str] = [
        *paper_figures.READABLE_NAMES.values(),
        *paper_figures.QUESTION_TITLES_SHORT.values(),
        *paper_figures.HEADLINE_LABELS_READABLE.values(),
        *paper_figures.HEATMAP_COLS,
    ]
    for grp_name, items in paper_figures.FIG_B_GROUPS:
        all_labels.append(grp_name)
        for row_label, _, _, _ in items:
            all_labels.append(row_label)
    for _, q_plain, intv_label, endpoint_label, _ in paper_figures.TABLE_SPECS:
        all_labels.extend([q_plain, intv_label, endpoint_label])

    for lbl in all_labels:
        for ident in FORBIDDEN_IDENTIFIERS:
            assert ident not in lbl, f"Label {lbl!r} contains raw identifier {ident!r}"
        stripped = _strip_math_and_comments(lbl)
        assert "_" not in stripped, f"Label {lbl!r} contains underscore outside math mode"
        assert not SNAKE_ID_RE.search(stripped), f"Label {lbl!r} contains snake_case token"
