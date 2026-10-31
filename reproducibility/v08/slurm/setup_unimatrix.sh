#!/usr/bin/env bash
set -euo pipefail

export PATH="$HOME/.local/bin:$PATH"
export TMPDIR="$HOME/tmp"
mkdir -p "$TMPDIR" /home/hohndor/eq-biodt-uq/logs /home/hohndor/eq-biodt-uq/raw

if ! command -v uv >/dev/null 2>&1; then
  curl -LsSf https://astral.sh/uv/install.sh | sh
fi

cd /home/hohndor/eq-biodt-uq/twin
uv sync --frozen
uv run --frozen python reproducibility/v08/tasks.py --out reviews/v08/uq/tasks.json
