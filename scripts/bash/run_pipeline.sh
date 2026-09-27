#!/bin/bash

# ============================================================
# Automated High-Throughput In-Silico Drug Screening Pipeline
# ============================================================

set -e

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"

cd "$PROJECT_ROOT"

echo "============================================================"
echo "AUTOMATED DRUG SCREENING PIPELINE"
echo "============================================================"

echo
echo "Project root:"
echo "$PROJECT_ROOT"

echo
echo "============================================================"
echo "STEP 1 — PREPARE LIGANDS"
echo "============================================================"

python scripts/python/prepare_ligands.py

echo
echo "============================================================"
echo "STEP 2 — RUN DOCKING"
echo "============================================================"

python scripts/python/run_docking.py

echo
echo "============================================================"
echo "STEP 3 — COLLECT RESULTS"
echo "============================================================"

python scripts/python/collect_results.py

echo
echo "============================================================"
echo "STEP 4 — ANALYZE RESULTS"
echo "============================================================"

python scripts/python/analyze_results.py

echo
echo "============================================================"
echo "PIPELINE COMPLETE"
echo "============================================================"

echo
echo "Final screening summary:"
echo "results/ranked/screening_summary.csv"

echo
echo "Docking results:"
echo "results/ranked/docking_results.csv"

echo
echo "Figure:"
echo "results/figures/best_docking_scores.png"

echo
echo "============================================================"
