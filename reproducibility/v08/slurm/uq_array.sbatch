#!/usr/bin/env bash
#SBATCH --job-name=eq-biodt-v08-uq
#SBATCH --partition=debug
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --array=0-1218%250
#SBATCH --time=02:00:00
#SBATCH --output=/home/hohndor/eq-biodt-uq/logs/%A_%a.out

set -euo pipefail

export PATH="$HOME/.local/bin:$PATH"
export TMPDIR="${SLURM_TMPDIR:-$HOME/tmp}"
mkdir -p "$TMPDIR" /home/hohndor/eq-biodt-uq/logs /home/hohndor/eq-biodt-uq/raw

cd /home/hohndor/eq-biodt-uq/twin

TASKS_FILE="${UQ_TASKS_FILE:-reviews/v08/uq/tasks.json}"
OUT_DIR="${UQ_RAW_DIR:-/home/hohndor/eq-biodt-uq/raw}"

uv run --frozen python reproducibility/v08/uq_task.py \
  --task-index "${SLURM_ARRAY_TASK_ID}" \
  --tasks "${TASKS_FILE}" \
  --out "${OUT_DIR}"
