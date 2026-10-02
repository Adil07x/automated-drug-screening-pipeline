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
echo "STEP 1 — REGISTER NEW COMPOUNDS"
echo "============================================================"

python scripts/python/register_ligands.py

echo
echo "============================================================"
echo "STEP 2 — CALCULATE COMPOUND PROPERTIES"
echo "============================================================"

python scripts/python/calculate_compound_properties.py

echo
echo "============================================================"
echo "STEP 3 — RUN ADMET-AI"
echo "============================================================"

conda run --no-capture-output -n admet-ai python scripts/python/run_admet_ai.py

echo
echo "============================================================"
echo "STEP 4 — SCREEN ADMET RESULTS"
echo "============================================================"

python scripts/python/screen_admet.py

echo
echo "============================================================"
echo "STEP 5 — PREPARE ELIGIBLE LIGANDS"
echo "============================================================"

python scripts/python/prepare_ligands.py

echo
echo "============================================================"
echo "STEP 6 — RUN DOCKING"
echo "============================================================"

python scripts/python/run_docking.py

echo
echo "============================================================"
echo "STEP 7 — COLLECT DOCKING RESULTS"
echo "============================================================"

python scripts/python/collect_results.py

echo
echo "============================================================"
echo "STEP 8 — ANALYZE RESULTS"
echo "============================================================"

python scripts/python/analyze_results.py

echo
echo "============================================================"
echo "STEP 9 — INTERACTION ANALYSIS"
echo "============================================================"

python scripts/python/interaction_analysis.py

echo
echo "============================================================"
echo "PIPELINE COMPLETE"
echo "============================================================"

echo
echo "Final screening summary:"
echo "results/ranked/screening_summary.csv"

echo
echo "Docking results:"
echo "results/ranked/current_docking_results.csv"

echo
echo "Interaction results:"
echo "results/interactions/current_interactions.csv"

echo
echo "Figure:"
echo "results/figures/best_docking_scores.png"

echo
echo "============================================================"
