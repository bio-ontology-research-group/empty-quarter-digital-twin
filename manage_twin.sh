#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")"

CMD="${1:-help}"

case "$CMD" in
  verify)
    uv run eq-biodt verify-data
    ;;
  build)
    uv run eq-biodt build-cubes --output-dir data
    ;;
  test)
    uv run pytest -v
    ;;
  serve)
    PORT="${2:-8095}"
    uv run eq-biodt serve --host 0.0.0.0 --port "$PORT"
    ;;
  *)
    echo "Usage: ./manage_twin.sh {verify|build|test|serve [port]}"
    exit 1
    ;;
esac
