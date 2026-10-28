#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
case "${1:-help}" in
  verify) uv run pytest -q tests/test_repairs.py ;;
  build) uv run python scripts/install_reasoner.py
         uv run python reproducibility/v03/reproduce.py --paper paper ;;
  test) uv run pytest -q ;;
  serve) uv run eq-biodt serve --host 127.0.0.1 --port "${2:-8095}" ;;
  *) echo 'Usage: ./manage_twin.sh {verify|build|test|serve [port]}'; exit 1 ;;
esac
