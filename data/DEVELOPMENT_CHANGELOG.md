# Development Change Log

## Current milestone — Incremental end-to-end pipeline

The pipeline was expanded from a docking/results workflow into a reusable incremental compound-screening workflow.

### Added

- `scripts/python/screen_admet.py`
- incremental ADMET screening
- docking eligibility generation
- canonical-SMILES structure hashing for duplicate detection
- incremental compound registration
- incremental compound-property calculation
- SCScore integration
- ADMET-AI database integration
- ADMET-AI engine wrapper
- incremental ADMET prediction
- MDAnalysis/ProLIF interaction workflow
- SQLite interaction storage
- CSV interaction reporting
- complete nine-stage `run_pipeline.sh` orchestration

### Modified

- `scripts/bash/run_pipeline.sh`
- `scripts/python/calculate_compound_properties.py`
- `scripts/python/run_admet_ai.py`
- `scripts/python/analyze_results.py`
- related database/result-generation code used by the incremental workflow

### Environment architecture

Two Conda environments are currently used:

- `drug-screening-pipeline` — core pipeline, RDKit, docking, MDAnalysis, ProLIF, etc.
- `admet-ai` — ADMET-AI and its model dependencies.

The main pipeline invokes ADMET-AI explicitly through:

```bash
conda run --no-capture-output -n admet-ai python scripts/python/run_admet_ai.py
```

### Current screening rules

The validated screening layer currently uses:

- AMES >= 0.50
- DILI >= 0.50
- ClinTox >= 0.50
- hERG >= 0.50
- Bioavailability_Ma < 0.50
- Caco2_Wang < -5.50
- Solubility_AqSolDB < -4.00
- QED < 0.30

Decision:

- 0 toxicity flags → PASS
- 1 toxicity flag → REVIEW
- 2 or more toxicity flags → FAIL

Docking eligibility:

- PASS → ELIGIBLE
- REVIEW → REVIEW
- FAIL → EXCLUDED

SCScore, SA Score, Lipinski, Veber, PAINS, and Brenk are calculated/stored as compound-level information, but they are not all current decision criteria.

### Validation

The current database contains:

- 10 compounds
- 10 property records
- 1,040 ADMET-AI records
- 10 screening records
- 10 docking eligibility records
- 20 docking pose records
- 13 interaction records

The nine-stage pipeline completed successfully for the current test set.

### Known operational limitation

SDF acquisition from the external source is currently performed manually because automated retrieval is restricted/unreliable in the current environment.

The downstream pipeline is designed to remain independent of that acquisition mechanism.
