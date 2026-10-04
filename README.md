# Automated High-Throughput In-Silico Drug Screening Pipeline

A reusable, incremental, database-driven computational drug-screening pipeline built around a local SDF compound library, SQLite as the source of truth, RDKit, SCScore, ADMET-AI, AutoDock Vina, MDAnalysis and ProLIF.

The project is designed for **repeated screening of growing compound libraries**, not a single one-off docking run. It runs on Linux/WSL and is controlled through Python and Bash scripts.

---

## Table of contents

- [1. Scientific objective](#1-scientific-objective)
- [2. Design principles](#2-design-principles)
- [3. Inputs and outputs](#3-inputs-and-outputs)
- [4. Pipeline workflow](#4-pipeline-workflow)
- [5. Current verified run](#5-current-verified-run)
- [6. Compound ingestion and duplicate detection](#6-compound-ingestion-and-duplicate-detection)
- [7. Compound properties](#7-compound-properties)
- [8. ADMET-AI](#8-admet-ai)
- [9. Screening rules and docking eligibility](#9-screening-rules-and-docking-eligibility)
- [10. Protein, docking configuration and ligand preparation](#10-protein-docking-configuration-and-ligand-preparation)
- [11. Docking, result collection and analysis](#11-docking-result-collection-and-analysis)
- [12. Interaction analysis](#12-interaction-analysis)
- [13. Reproducibility and experiment identity](#13-reproducibility-and-experiment-identity)
- [14. Database](#14-database)
- [14A. Database schema and SQL audit reference](#14a-database-schema-and-sql-audit-reference)
- [15. Scripts and orchestration](#15-scripts-and-orchestration)
- [16. Environments and technology stack](#16-environments-and-technology-stack)
- [17. Installation and running](#17-installation-and-running)
- [18. Project structure](#18-project-structure)
- [19. Verification and idempotency](#19-verification-and-idempotency)
- [20. Known limitations](#20-known-limitations)
- [21. Scientific interpretation](#21-scientific-interpretation)
- [22. Future development](#22-future-development)
- [23. Troubleshooting](#23-troubleshooting)
- [24. Reproducibility checklist](#24-reproducibility-checklist)
- [25. Development history and status](#25-development-history-and-status)
- [25A. Git and repository hygiene](#25a-git-and-repository-hygiene)
- [26. Contributing](#26-contributing)
- [27. References](#27-references)
- [27A. Final operational checklist](#27a-final-operational-checklist)
- [28. License, author, contact and disclaimer](#28-license-author-contact-and-disclaimer)


---

## 1. Scientific objective

> **Which compounds in a screening library are computationally prioritized as potential binders of a selected protein target?**

The pipeline accepts local compound structures in SDF format, registers their chemical identities, calculates reusable molecular properties, predicts ADMET endpoints, applies project-defined screening rules, prepares eligible ligands, docks them against a configured target, collects and analyzes the results, and finally performs protein–ligand interaction analysis.

It does not claim experimental binding, biological activity, therapeutic efficacy or clinical safety. It is a computational **prioritization** workflow.

---

## 2. Design principles

### 2.1 SQLite is the source of truth

```text
data/database/screening_database.sqlite
```

The database stores compounds, targets, properties, ADMET predictions, screening decisions, docking configurations, experiments, docking results and interaction results, each linked to the inputs and settings that produced it. CSV files and figures are generated outputs for inspection and exchange; an exported CSV is never the authoritative record.

### 2.2 The compound library is persistent and reusable

New SDF files are added to an existing library without removing previously registered compounds. A compound that fails a screening rule stays in the master library and can still be selected manually (see [section 9](#manual-compound-selection)).

### 2.3 Incremental processing

Only new or changed records are processed; existing results are reused.

```text
Existing database: 10 compounds
        |
Add 100 new SDF structures
        |
        v
Registration
   +--> existing structures -> skip
   +--> new structures      -> register
                                |
                                v
                  properties -> ADMET-AI -> screening
                                |
                                v
                  (eligible only) docking -> interactions
```

### 2.4 Compound-level vs target-level work

| Reusable across targets (compound-level) | Target-specific |
|---|---|
| SDF ingestion, structure identity | Ligand preparation for a docking setup |
| Molecular properties, Lipinski, PAINS, Brenk, SA Score, SCScore | Docking |
| ADMET predictions (once per engine/version) | Docking-result analysis |
| ADMET screening | Protein–ligand interaction analysis |

```text
100,000 compounds
   +--> properties calculated once
   +--> ADMET predicted once per engine/version
   +--> Target A docking
   +--> Target B docking
   +--> Target C docking
```

In the data model, **a ligand is a reusable chemical entity** and **an experiment is one ligand × target × docking configuration**.

### 2.5 Calculation is separate from screening

A calculated property (for example SCScore) is not automatically an exclusion rule. The screening layer decides which stored values currently affect progression. This allows rules to change without rebuilding the compound database.

---

## 3. Inputs and outputs

### Inputs

1. **Protein target** — a three-dimensional, docking-ready structure (PDB-compatible). A FASTA sequence alone is not sufficient; automatic structure prediction is not part of the pipeline. The docking site is defined by a configurable docking box.
2. **Compound library** — local SDF files in `data/raw/compounds/`. Each compound receives an internal ligand ID; canonical SMILES, structure hash and external identifiers such as PubChem CID are retained where available.
3. **Screening and docking configuration** — property/ADMET screening criteria, docking engine and parameters, preparation settings and docking-box parameters.

### Outputs

- Molecular descriptors, Lipinski calculations, structural alerts, SA Score and SCScore.
- ADMET-AI endpoint predictions (104 endpoints per compound).
- Screening decisions (`PASS` / `REVIEW` / `FAIL`) and docking-eligibility records.
- Prepared 3D ligand structures and PDBQT files; receptor files.
- Docking poses and predicted scores.
- ProLIF protein–ligand interaction records.
- Experiment and configuration identifiers.
- Ranked reports, CSV exports and figures; SQLite records as the structured source of truth.

---

## 4. Pipeline workflow

```text
Local SDF compound library
        |
        v
1. Compound registration
        |
        v
2. Compound properties
   - molecular weight, LogP, HBD / HBA, Lipinski violations
   - TPSA, rotatable bonds, Veber
   - PAINS, Brenk, SA Score, SCScore
        |
        v
3. ADMET-AI (104 endpoints, ADMET-AI 2.0.1)
   absorption, distribution, metabolism, excretion, toxicity
        |
        v
4. ADMET screening
   toxicity / absorption / physicochemical flags, QED
   PASS / REVIEW / FAIL  ->  docking eligibility
        |
        v
5. Ligand preparation
   SMILES validation, hydrogens, 3D embedding,
   MMFF optimization, SDF + PDBQT (Meeko)
        |
        v
6. AutoDock Vina docking (multiple poses)
        |
        v
7. Docking-result collection
        |
        v
8. Docking analysis
        |
        v
9. MDAnalysis + ProLIF interaction analysis
        |
        v
SQLite + CSV + figures
```

Individual stages remain independently executable so failed or incomplete stages can be resumed without repeating earlier calculations.

---

## 5. Current verified run

The current validated test database snapshot (documented from the latest verified run; update these counts when the database changes):

| Component | Count |
|---|---:|
| Registered compounds | 10 |
| Compound-property records | 10 |
| ADMET-AI endpoint records | 1,040 (104 × 10) |
| ADMET screening records | 10 |
| Docking-eligibility records | 10 |
| Compounds eligible and docked | 2 |
| Docking pose records | 20 (10 per experiment) |
| Interaction records | 13 |

The library contains Azadirachtin-D, -F, -H, -I, Caffeine, Oleic acid, Quercetin, Quinic acid and Riboflavin, plus Aspirin, added later as an incremental-processing test.

### Docking results

| Compound | Screening decision | Best Vina score (kcal/mol) | Poses | Interaction events |
|---|---|---:|---:|---:|
| Oleic acid | PASS | see current docking CSV/DB | 10 | 7 |
| Quinic acid | PASS | see current docking CSV/DB | 10 | 6 |

Oleic acid contacts were hydrophobic and van der Waals; Quinic acid contacts were van der Waals. These are computational predictions only. Exact docking-score values should be read from the current SQLite/CSV snapshot rather than treated as immutable README values.

### Stage status

```text
1. Compound registration       PASS
2. Compound properties         PASS
3. ADMET-AI                    PASS
4. ADMET screening             PASS
5. Ligand preparation          PASS
6. Docking                     PASS
7. Result collection           PASS
8. Docking analysis            PASS
9. Interaction analysis        PASS
```

---

## 6. Compound ingestion and duplicate detection

### Input location and workflow

```text
data/raw/compounds/
```

SDF files are currently downloaded **manually**, because automated retrieval from the external compound source is affected by site/network restrictions. `scripts/python/download_compounds.py` (and PubChem REST API acquisition) exists and can be re-enabled without changing downstream architecture.

```text
Download SDF manually -> place *.sdf in data/raw/compounds/ -> run pipeline
```

### Chemical duplicate detection

`register_ligands.py`:

```text
SDF -> RDKit molecule -> canonical SMILES -> SHA-256 structure hash
         |
         +-- hash exists in ligands.structure_hash -> skip
         +-- new hash                              -> register
```

Identity is derived from the molecular structure, **not** the filename, compound name or source identifier. A copied Caffeine SDF with a different filename or identifier therefore maps to the existing compound rather than creating a second ligand. The filename is used only as the initial compound name.

### SDF requirements

- Each SDF must contain a readable, chemically interpretable structure; malformed molecules are skipped.
- The current registration code processes the **first molecule** in an SDF supplier. Multi-record SDFs must be handled deliberately.
- Keep source SDFs after registration; downstream ADMET and interaction code can use the registered source path.
- Large libraries (tens or hundreds of thousands of structures) should be treated as a persistent library, added incrementally.

### Planned ingestion improvements

Explicit multi-record SDF ingestion, stronger normalization, salt/solvate handling, stereochemistry-aware identity rules, changed-file detection, validation reports, duplicate audit reports and ingestion manifests.

---

## 7. Compound properties

Scripts: `calculate_compound_properties.py`, `calculate_lipinski.py`, `scscore_runner.py`.

Stored fields: molecular weight, LogP, H-bond donors/acceptors, molar refractivity, Lipinski violations, TPSA, rotatable bonds, Veber pass, PAINS alert count/details, Brenk alert count/details, SA Score, SCScore, and QED (QED is also used by the screening layer).

Property calculation is incremental: compounds that already have property records are not recalculated.

### SCScore

SCScore (Synthetic Complexity Score) is calculated locally using the NumPy-compatible model under:

```text
external/scscore/models/full_reaxys_model_1024bool/model.ckpt-10654.as_numpy.json.gz
```

> **SCScore is calculated and stored, but is not currently a PASS/REVIEW/FAIL screening threshold.** Do not describe current screening decisions as SCScore-based.

---

## 8. ADMET-AI

Current engine: **ADMET-AI 2.0.1**, run on CPU in a dedicated Conda environment (`admet-ai`), separate from the main `drug-screening-pipeline` environment because its dependency stack differs. Activating the main environment does **not** make `admet_ai` importable; the pipeline invokes it explicitly:

```bash
conda run --no-capture-output -n admet-ai \
    python scripts/python/run_admet_ai.py
```

Individual endpoint predictions are stored (not collapsed into one score), covering absorption, distribution, metabolism, excretion, toxicity and drug-likeness/physicochemical properties.

### Replacing or adding an ADMET model

The `admet_results` table stores `engine`, `engine_version`, `endpoint`, `value`, `unit`, `prediction_type`, `uncertainty`, `status`, `raw_result` and `calculated_at`, so engines are not hard-coded into the schema. To add a future model:

1. Create or update the engine wrapper (`scripts/python/admet/`).
2. Record engine name and version.
3. Define its endpoint set.
4. Store predictions in `admet_results`.
5. Update the screening layer only if the rules need endpoints the new model provides.
6. Recalculate only compounds lacking complete results for that engine/version.

Do not overwrite old model results because a newer model exists; this preserves provenance and allows model comparison.

---

## 9. Screening rules and docking eligibility

Script: `scripts/python/screen_admet.py`. The rules were reconstructed from the project's earlier A29L screening methodology and validated against the existing database. They are project-defined computational rules, not universal biological or regulatory cutoffs.

| Category | Endpoint | Flag if |
|---|---|---:|
| Toxicity | AMES | ≥ 0.50 |
| Toxicity | DILI | ≥ 0.50 |
| Toxicity | ClinTox | ≥ 0.50 |
| Toxicity | hERG | ≥ 0.50 |
| Absorption | Bioavailability_Ma | < 0.50 |
| Absorption | Caco2_Wang | < -5.50 |
| Physicochemical | Solubility_AqSolDB | < -4.00 |
| Physicochemical | QED | < 0.30 |

Decision logic (based on the toxicity flag count; absorption and physicochemical flags are recorded):

```text
toxicity_flags >= 2  -> FAIL
toxicity_flags == 1  -> REVIEW
toxicity_flags == 0  -> PASS
```

Docking eligibility mapping:

```text
PASS   -> ELIGIBLE
REVIEW -> REVIEW
FAIL   -> EXCLUDED
```

Screening is incremental and idempotent. These categories organize compounds for downstream computation and are not clinical, regulatory or experimental safety classifications; neither `PASS` nor `FAIL` demonstrates safety or toxicity.

### Manual compound selection

A user can deliberately include compounds that would not be selected automatically: ones passing the filters, requiring review, failing a rule, or a hand-picked list. The scientific reason for overriding an automated decision should be documented. Lipinski-style rules are heuristics and do not determine biological activity.

---

## 10. Protein, docking configuration and ligand preparation

### Target

The current demonstration target is **PDB 1IEP, chain A**, prepared as a docking-ready PDBQT receptor (`data/processed/protein/1IEP_chainA.pdbqt`) and registered with `register_target.py`.

### Docking configuration

The docking site is an explicit configuration, not hard-coded into the docking script. `config/docking_config.json`:

```json
{
  "target": "1IEP_chainA",
  "receptor": "data/processed/protein/1IEP_chainA.pdbqt",
  "docking_box": {
    "center_x": 15.190, "center_y": 53.903, "center_z": 16.917,
    "size_x": 19.0, "size_y": 27.0, "size_z": 24.0
  },
  "docking_parameters": { "exhaustiveness": 8, "num_modes": 10 }
}
```

The configuration is registered in `docking_configs` (box, exhaustiveness, number of modes, Vina version, SHA-256 `config_hash`). Changing an important parameter produces a distinguishable configuration and experiment rather than silently replacing a result.

### Ligand preparation

Script: `prepare_ligands.py`. Only ligands whose docking eligibility permits it are prepared. Steps: SMILES validation → hydrogen addition → 3D conformer generation → MMFF geometry optimization → SDF output → PDBQT via Meeko. Prepared structures are tied to their internal ligand IDs.

---

## 11. Docking, result collection and analysis

### Docking — `run_docking.py`

AutoDock Vina is the current docking engine. For each selected ligand the pipeline records configuration, target, ligand, Vina version, docking log, poses and predicted scores. Multiple poses are retained (10 per experiment currently). Vina scores are scoring-function estimates, not measured affinities; within the same setup, a more negative score is generally more favorable.

### Collection — `collect_results.py`

Reads Vina logs, extracts pose number, affinity and RMSD bounds, stores results in SQLite, and writes:

```text
results/ranked/current_docking_results.csv
results/ranked/current_screening_results.csv
```

The second file is exported from the `current_screening_results` view and combines screening decisions, docking scores and interaction summaries.

### Analysis — `analyze_results.py`

Selects the best pose per compound, ranks compounds by affinity, and writes:

```text
results/ranked/screening_summary.csv
results/figures/best_docking_scores.png
```

---

## 12. Interaction analysis

Script: `interaction_analysis.py`, using MDAnalysis, ProLIF and RDKit. The stage:

1. identifies completed docking experiments;
2. selects pose 1 (the top-ranked pose);
3. reads the docking pose;
4. reconstructs the docked ligand from the original SDF plus docking coordinates;
5. runs ProLIF;
6. extracts protein–ligand interaction events;
7. writes records to SQLite;
8. exports `results/interactions/current_interactions.csv`.

Recorded information includes protein residue, residue number, chain (where available), interaction type, distance (where available), ligand/protein atom information, and analysis engine and version.

Interaction detection describes contacts in a computational pose; it does not establish that the complex is experimentally stable or functional.

**Warnings.** A successful run emits MDAnalysis warnings about a deprecated topology import path and bond guessing when the AtomGroup lacks explicit bond information. These are not failures but should be reviewed on dependency upgrades.

---

## 13. Reproducibility and experiment identity

SHA-256 hashes identify ligand structures, receptor structures and docking configurations (engine, box centre and size, exhaustiveness, number of poses, engine version).

```text
Target + Ligand + Docking configuration = Experiment
```

(`experiments` enforces `UNIQUE(target_id, ligand_id, config_id)`.) If inputs and configuration are unchanged, previous results can be identified for reference; a change to protein, ligand or variant, docking site or box, engine or parameters creates a distinguishable experiment.

Results depend on exact input structures, software versions, models and settings, so results from different configurations should not be compared without accounting for those differences. When you change SDF files, thresholds, ADMET engines or versions, docking configurations, target structures or dependency versions, record the change in Git and update this README if the workflow changes. Do not silently change historical screening criteria.

The database should preserve: compound identity, source, structure, properties, ADMET engine/version, screening decision, docking configuration, docking results and interaction results.

```text
Input Structure -> Software / Model -> Configuration -> Experiment -> Computational Result
```

---

## 14. Database

Tables:

```text
targets               ligands                compound_properties
admet_results         admet_screening        docking_eligibility
docking_configs       experiments            docking_results
interaction_results
```

View: `current_screening_results` (integrated reporting across screening, experiments, docking scores and configuration, ADMET decisions and interactions).

Relationships: `targets → docking_configs → experiments → docking_results`, with each experiment also referencing a ligand. The database preserves processing history and provenance rather than only the final winners.

---


## 14A. Database schema and SQL audit reference

The current database is:

```text
data/database/screening_database.sqlite
```

The major tables are:

```text
targets
ligands
compound_properties
admet_results
admet_screening
docking_eligibility
docking_configs
experiments
docking_results
interaction_results
```

The database also contains the reporting view:

```text
current_screening_results
```

### Core schema

The following reflects the current major schema. `docking_results` and `interaction_results` contain additional fields and should be inspected with SQLite before writing ad-hoc queries.

```sql
CREATE TABLE targets (
   target_id INTEGER PRIMARY KEY AUTOINCREMENT,
   name TEXT NOT NULL,
   pdb_id TEXT,
   chain TEXT,
   receptor_file TEXT,
   receptor_hash TEXT,
   created_at TEXT DEFAULT CURRENT_TIMESTAMP
);

CREATE TABLE ligands (
   ligand_id INTEGER PRIMARY KEY AUTOINCREMENT,
   name TEXT NOT NULL,
   source TEXT,
   source_id TEXT,
   smiles TEXT,
   structure_hash TEXT NOT NULL,
   source_file TEXT,
   created_at TEXT DEFAULT CURRENT_TIMESTAMP
);

CREATE TABLE compound_properties (
   property_id INTEGER PRIMARY KEY AUTOINCREMENT,
   ligand_id INTEGER NOT NULL,
   molecular_weight REAL,
   logp REAL,
   hbd INTEGER,
   hba INTEGER,
   lipinski_violations INTEGER,
   admet_data TEXT,
   calculated_at TEXT DEFAULT CURRENT_TIMESTAMP,
   tpsa REAL,
   rotatable_bonds INTEGER,
   pains_alert_count INTEGER,
   pains_alerts TEXT,
   brenk_alert_count INTEGER,
   brenk_alerts TEXT,
   sa_score REAL,
   sc_score REAL,
   veber_pass INTEGER,
   FOREIGN KEY (ligand_id) REFERENCES ligands(ligand_id)
);

CREATE TABLE admet_results (
   admet_result_id INTEGER PRIMARY KEY AUTOINCREMENT,
   ligand_id INTEGER NOT NULL,
   engine TEXT NOT NULL,
   engine_version TEXT,
   endpoint TEXT NOT NULL,
   value REAL,
   unit TEXT,
   prediction_type TEXT,
   uncertainty REAL,
   status TEXT,
   raw_result TEXT,
   calculated_at TEXT DEFAULT CURRENT_TIMESTAMP,
   FOREIGN KEY (ligand_id) REFERENCES ligands(ligand_id)
);

CREATE TABLE admet_screening (
   screening_id INTEGER PRIMARY KEY AUTOINCREMENT,
   ligand_id INTEGER NOT NULL,
   engine TEXT NOT NULL,
   decision TEXT NOT NULL,
   toxicity_flags INTEGER DEFAULT 0,
   absorption_flags INTEGER DEFAULT 0,
   physicochemical_flags INTEGER DEFAULT 0,
   qed_value REAL,
   rationale TEXT,
   screened_at TEXT DEFAULT CURRENT_TIMESTAMP,
   FOREIGN KEY (ligand_id) REFERENCES ligands(ligand_id)
);

CREATE TABLE docking_eligibility (
   eligibility_id INTEGER PRIMARY KEY AUTOINCREMENT,
   ligand_id INTEGER NOT NULL,
   screening_id INTEGER NOT NULL,
   eligibility TEXT NOT NULL,
   reason TEXT,
   created_at TEXT DEFAULT CURRENT_TIMESTAMP,
   FOREIGN KEY (ligand_id) REFERENCES ligands(ligand_id)
);
```

`docking_results`, `interaction_results`, `docking_configs`, `experiments`, and `current_screening_results` should be inspected directly from the live database when exact columns are required.

There are currently no SQLite triggers in the verified database design; `current_screening_results` is a view.

### Open SQLite

```bash
sqlite3 data/database/screening_database.sqlite
```

### Table counts

```sql
SELECT 'ligands' AS table_name, COUNT(*) AS count FROM ligands
UNION ALL
SELECT 'compound_properties', COUNT(*) FROM compound_properties
UNION ALL
SELECT 'admet_results', COUNT(*) FROM admet_results
UNION ALL
SELECT 'admet_screening', COUNT(*) FROM admet_screening
UNION ALL
SELECT 'docking_eligibility', COUNT(*) FROM docking_eligibility
UNION ALL
SELECT 'docking_results', COUNT(*) FROM docking_results
UNION ALL
SELECT 'interaction_results', COUNT(*) FROM interaction_results;
```

### Compound audit

```sql
SELECT
    ligand_id,
    name,
    source,
    source_id,
    smiles,
    structure_hash,
    source_file
FROM ligands
ORDER BY ligand_id;
```

### Property completeness

```sql
SELECT
    l.ligand_id,
    l.name,
    CASE
        WHEN cp.ligand_id IS NULL THEN 'MISSING'
        ELSE 'COMPLETE'
    END AS property_status
FROM ligands l
LEFT JOIN compound_properties cp
    ON cp.ligand_id = l.ligand_id
ORDER BY l.ligand_id;
```

### ADMET endpoint completeness

```sql
SELECT
    l.name,
    COUNT(ar.admet_result_id) AS endpoint_count
FROM ligands l
JOIN admet_results ar
    ON l.ligand_id = ar.ligand_id
WHERE ar.engine = 'ADMET-AI'
GROUP BY l.ligand_id, l.name
ORDER BY l.ligand_id;
```

A complete current ADMET-AI record contains 104 endpoints.

### Screening audit

```sql
SELECT
    l.ligand_id,
    l.name,
    s.engine,
    s.decision,
    s.toxicity_flags,
    s.absorption_flags,
    s.physicochemical_flags,
    s.qed_value,
    s.rationale
FROM admet_screening s
JOIN ligands l
    ON l.ligand_id = s.ligand_id
ORDER BY l.ligand_id;
```

### Docking eligibility audit

```sql
SELECT
    l.ligand_id,
    l.name,
    de.eligibility,
    de.reason
FROM docking_eligibility de
JOIN ligands l
    ON l.ligand_id = de.ligand_id
ORDER BY l.ligand_id;
```

### Inspect the reporting view

```sql
SELECT *
FROM current_screening_results;
```

### Inspect docking schema before writing SQL

```sql
.schema docking_results
```

The verified docking-affinity field is:

```text
docking_affinity_kcal_mol
```

A previous query that assumed a column named `affinity` produced a `no such column: affinity` error. That was a query/schema mismatch, not evidence of database corruption.

Likewise, inspect:

```sql
.schema interaction_results
```

before assuming an interaction primary-key or column name.

Exit SQLite with:

```sql
.quit
```

---

## 15. Scripts and orchestration

| Area | Script | Role |
|---|---|---|
| Database | `create_database.py` | Create SQLite schema |
| Target | `register_target.py` | Register target/receptor |
| Acquisition | `download_compounds.py` | Automated compound download (currently not the dependable path) |
| Ingestion | `register_ligands.py` | SDF registration and duplicate detection |
| Properties | `calculate_compound_properties.py` | Descriptors, alerts, SA Score, SCScore |
| Properties | `calculate_lipinski.py` | Lipinski calculations |
| Properties | `scscore_runner.py` | SCScore integration |
| ADMET | `run_admet_ai.py` | Incremental ADMET-AI prediction |
| ADMET | `admet/admet_ai_engine.py` | Engine wrapper (replaceable) |
| Screening | `screen_admet.py` | Flags, decisions, docking eligibility |
| Docking | `prepare_ligands.py`, `run_docking.py` | Preparation and Vina docking |
| Results | `collect_results.py`, `analyze_results.py`, `interaction_analysis.py` | Collection, analysis, interactions |


### Python file reference

| File | Purpose |
|---|---|
| `create_database.py` | Creates/initializes the SQLite schema |
| `register_target.py` | Registers a target/receptor |
| `download_compounds.py` | Optional PubChem-based acquisition; not the dependable current acquisition path |
| `register_ligands.py` | Reads SDF input, derives structure identity and registers new ligands |
| `calculate_compound_properties.py` | Calculates incremental molecular properties, alerts, SA Score and SCScore |
| `calculate_lipinski.py` | Lipinski-related calculations |
| `scscore_runner.py` | SCScore integration |
| `admet/admet_ai_engine.py` | ADMET-AI engine wrapper |
| `run_admet_ai.py` | Incremental ADMET-AI execution |
| `screen_admet.py` | Current ADMET flagging, PASS/REVIEW/FAIL logic and docking eligibility |
| `prepare_ligands.py` | 3D ligand preparation and PDBQT generation |
| `run_docking.py` | AutoDock Vina execution |
| `collect_results.py` | Parses docking output and stores pose/score records |
| `analyze_results.py` | Generates ranked docking/reporting outputs |
| `interaction_analysis.py` | MDAnalysis/ProLIF interaction analysis |

All Python scripts live under `scripts/python/`.

### Bash runner — `scripts/bash/run_pipeline.sh`

```text
1. register_ligands.py
2. calculate_compound_properties.py
3. run_admet_ai.py          (via the admet-ai environment)
4. screen_admet.py
5. prepare_ligands.py
6. run_docking.py
7. collect_results.py
8. analyze_results.py
9. interaction_analysis.py
```

The script uses `set -e`, so it stops when any stage fails.

---

## 16. Environments and technology stack

### Environments

| Environment | Purpose | Validated components |
|---|---|---|
| `drug-screening-pipeline` (main) | Registration, properties, SCScore, ligand/protein preparation, docking, interaction analysis | Python 3.12, RDKit 2025.09.5, Open Babel 3.2.1, Meeko 0.8.0, AutoDock Vina, MDAnalysis 2.10.0, ProLIF, pandas, NumPy, SQLite |
| `admet-ai` | ADMET prediction | Python 3.12, ADMET-AI 2.0.1, PyTorch (CPU; CUDA not available), NumPy, pandas, RDKit |

Interaction analysis now runs from the main environment. An earlier separate `prolif-env` (ProLIF 2.2.2, MDAnalysis 2.10.0) should only be retained if a future dependency conflict requires it. Do not merge everything into one environment for convenience; ADMET-AI is kept separate because its dependency stack differs.

### Technology stack

| Technology | Role |
|---|---|
| Python, Bash | Pipeline implementation and orchestration |
| SQLite | Persistent experiment and result database |
| RDKit | Structures, descriptors, Lipinski, QED, 3D chemistry |
| SCScore | Synthetic-complexity estimation |
| ADMET-AI (PyTorch) | ADMET prediction |
| AutoDock Vina | Molecular docking |
| Meeko | Ligand preparation, PDBQT generation |
| Open Babel | Structure conversion |
| MDAnalysis, ProLIF | Structure handling and interaction analysis |
| PubChem REST API | Compound acquisition / external information |
| Pandas, NumPy, Matplotlib | Tabular processing, numerics, figures |
| Conda, Git/GitHub | Environment management, version control |
| SHA-256 | Structure, receptor and configuration identity |

**Docking engine scope:** only AutoDock Vina is part of the verified pipeline. GNINA, smina, QVina and CNN-based scoring are possible future engines, not current components.

---

## 17. Installation and running

The project targets Linux; the development setup is Ubuntu on WSL2 (Windows).

```bash
git clone git@github.com:Adil07x/automated-drug-screening-pipeline.git
cd automated-drug-screening-pipeline
```

Environment definitions are under `environment/environment.yml`.

### Verify the environments

```bash
conda activate drug-screening-pipeline
python --version
vina --version
obabel -V
python -c "from rdkit import Chem; print('RDKit OK')"
python -c "import prolif, MDAnalysis; print('ProLIF:', prolif.__version__, 'MDAnalysis:', MDAnalysis.__version__)"

conda run -n admet-ai python -c "import admet_ai; print('admet_ai: OK')"
conda run -n admet-ai python -c "import torch; print(torch.__version__, torch.cuda.is_available())"
```

### Run the complete pipeline

```bash
conda activate drug-screening-pipeline
cd ~/Bioinformatics/automated-drug-screening-pipeline
# place new *.sdf files in data/raw/compounds/
bash scripts/bash/run_pipeline.sh
```

Existing and new compounds are processed incrementally.

### Run individual stages (recommended for debugging)

```bash
python scripts/python/register_ligands.py
python scripts/python/calculate_compound_properties.py
conda run --no-capture-output -n admet-ai python scripts/python/run_admet_ai.py
python scripts/python/screen_admet.py
python scripts/python/prepare_ligands.py
python scripts/python/run_docking.py
python scripts/python/collect_results.py
python scripts/python/analyze_results.py
python scripts/python/interaction_analysis.py
```

Avoid repeatedly rerunning expensive docking while debugging.

---

## 18. Project structure

```text
automated-drug-screening-pipeline/
├── config/
│   └── docking_config.json
├── data/
│   ├── database/screening_database.sqlite
│   ├── ligands/{3d,pdbqt}/
│   ├── metadata/compound_library.csv
│   ├── processed/{ligands,protein}/
│   └── raw/{compounds,protein}/
├── environment/environment.yml
├── external/scscore/
├── results/
│   ├── docking/
│   ├── figures/
│   ├── interactions/
│   └── ranked/
├── scripts/
│   ├── bash/run_pipeline.sh
│   └── python/
│       ├── admet/admet_ai_engine.py
│       ├── analyze_results.py
│       ├── calculate_compound_properties.py
│       ├── calculate_lipinski.py
│       ├── collect_results.py
│       ├── create_database.py
│       ├── download_compounds.py
│       ├── interaction_analysis.py
│       ├── prepare_ligands.py
│       ├── register_ligands.py
│       ├── register_target.py
│       ├── run_admet_ai.py
│       ├── run_docking.py
│       ├── scscore_runner.py
│       ├── screen_admet.py
│       └── test_experiment.py
├── .gitignore
├── LICENSE
└── README.md
```

| Directory | Role |
|---|---|
| `config/` | Reproducible computational settings (docking box, Vina parameters) |
| `data/raw/` | Source compound and protein structures |
| `data/processed/` | Prepared protein and ligand structures |
| `data/ligands/` | Generated 3D ligands and docking-ready PDBQT files |
| `data/database/` | SQLite source-of-truth database |
| `external/scscore/` | Local SCScore source and model resources (large payloads excluded via `.gitignore`) |
| `results/docking/` | Vina outputs and logs |
| `results/interactions/` | ProLIF interaction results |
| `results/ranked/`, `results/figures/` | Ranked CSV reports and figures |

---

## 19. Verification and idempotency

Tested against the 10-compound library. Rerunning after registration and property calculation produced:

```text
New ligands registered: 0
Compounds requiring property calculation: 0
Compounds requiring ADMET prediction: 0
No new compounds require ADMET screening.
```

Verified behavior:

- Existing compounds are recognized, not re-registered; duplicate Caffeine structures with different identifiers map to the existing compound.
- Adding Aspirin created a new ligand ID and added downstream calculations without recalculating existing compounds.
- Existing properties and complete ADMET results are retained; screening rows are not recreated unnecessarily.
- Only eligible compounds proceed to preparation and docking (2 docked; 2 interaction analyses completed).

A rerun with no new input performs little or no expensive computation.

---

## 20. Known limitations

**Execution scope is global, not per experiment.** The database is experiment-aware (target → configuration → experiment → results), but the execution scripts are not yet fully isolated by experiment/target:

- `run_docking.py` operates on the shared ligand PDBQT folder rather than only the ligands selected for a given target run.
- `collect_results.py` reads all completed experiments and looks for their logs in the shared `results/docking/` directory.
- `config/docking_config.json` is currently specific to `1IEP_chainA`.
- Older documentation/code paths referenced `results/ranked/docking_results.csv`, while the current collection output is `current_docking_results.csv`. Verify the current `analyze_results.py` input path after future changes; do not treat the historical filename as a required current output.
- `test_experiment.py` demonstrates experiment records and is not connected to the production pipeline.

Until this is addressed, a run against a different target (for example a historical A29L validation) risks mixing ligands and results from other runs. The planned fix is to scope docking, collection and analysis to the ligands and experiments selected for the target/configuration at hand, while still drawing compounds from the shared library (no separate "campaign" schema is required).

**Other limitations**

- First-molecule-only SDF registration; no changed-file detection yet.
- Automated compound download is unreliable because of external restrictions.
- SCScore is not part of the screening thresholds.
- Docking uses a single engine (Vina) and interaction analysis uses only the top-ranked pose.
- Throughput depends on CPU resources, molecular complexity and library size; the architecture is intended to scale beyond the small demonstration set.

---

## 21. Scientific interpretation

The pipeline is a computational prioritization framework:

```text
Compound Library -> Molecular / ADMET Screening -> Docking
   -> Interaction Analysis -> Computational Prioritization -> Experimental Validation
```

- **Docking scores** are estimates from a scoring function, not experimental binding affinities or predictions of activity.
- **ADMET categories** are workflow decision rules; they are not clinical, regulatory or experimental safety classifications.
- **Lipinski and physicochemical filters** are heuristics. Failing compounds can remain in the library and be selected manually when scientifically justified.
- **Receptor structure** — docking requires an appropriate 3D structure and a defined binding region.
- **Interaction analysis** describes contacts in the analyzed pose only.

Final determination of binding, activity, pharmacokinetics, toxicity and therapeutic potential requires experimental and/or validated external evidence.

---

## 22. Future development

Future work aims to increase automation, scalability and analysis depth without changing the experiment-tracking architecture.

- **Ingestion:** multi-record SDFs, stronger normalization, salt/solvate and stereochemistry policies, changed-file detection, manifests, validation and duplicate-audit reports, batch ingestion of very large libraries.
- **Experiment isolation:** target/experiment-scoped docking, collection and analysis; target-specific configuration files.
- **Orchestration:** parallel processing, resumable runs, tracking of failed compounds and errors, checkpointing, resource-aware and batch execution, GPU methods where supported.
- **Screening:** versioned, configurable rules layer (Lipinski, Veber, SCScore, SA Score, PAINS, Brenk, QED, solubility, absorption, toxicity, other ADMET endpoints); every rule documented and versioned.
- **ADMET:** multiple engines and model versions side by side, old predictions preserved, model comparison, automatic detection of compounds needing prediction per engine/version.
- **Docking:** GNINA, smina, QVina variants, alternative Vina scoring configurations, machine-learning rescoring (each treated as a separate method), stronger pose-quality checks.
- **Pose and interaction analysis:** configurable top-N poses, clustering, interaction fingerprints and frequency, residue-level comparison, comparative binding-site analysis, explicit bond information handling, dependency compatibility updates.
- **Ranking:** documented multi-signal ranking combining ADMET status, properties, SCScore, docking score, interaction count/type, contacted residues and pose consistency.
- **Reporting and tools:** automated reports (compound, ADMET, ranking, interaction, configuration, experiment history, QC), Excel reports, dashboards, command-line queries, filtered exports, experiment comparison, compound audit trail, screening-rule and model provenance.
- **Experimental validation integration:** store externally generated binding measurements, assay results, literature evidence and validation status, clearly separated from computational predictions.
- **Method provenance:** any new method records software, version, model, input structure, configuration, parameters, execution status, output and experiment identity.

---


### Important verified error history

- **ADMET-AI import error:** `admet_ai` is installed in the dedicated `admet-ai` environment, not the main environment.
- **MDAnalysis import error:** MDAnalysis was initially absent from the main environment; the current validated installation is MDAnalysis 2.10.0 in `drug-screening-pipeline`.
- **ProLIF environment confusion:** the current successful interaction-analysis workflow runs in the main environment; the historical `prolif-env` is not required by the current runner.
- **MDAnalysis warnings:** the current successful interaction run emits a deprecated topology-import warning and a bond-guessing warning. These are warnings, not observed pipeline failures.
- **Docking SQL column error:** a prior `affinity` query failed because the current field is `docking_affinity_kcal_mol`. Inspect `.schema docking_results` before constructing audit SQL.

## 23. Troubleshooting

**`ModuleNotFoundError: No module named 'rdkit'`** — you are probably in the `base` Conda environment. Activate the project environment (`conda activate drug-screening-pipeline`) and verify with `python -c "import rdkit; print(rdkit.__version__)"`. A failed import means registration never ran, so an unchanged ligand count does not confirm duplicate detection.

**`ModuleNotFoundError: No module named 'admet_ai'`** — the main environment does not contain ADMET-AI. Verify with `conda run -n admet-ai python -c "import admet_ai; print('OK')"`; the pipeline already launches it in the `admet-ai` environment.

**`ModuleNotFoundError: No module named 'MDAnalysis'`** — install it in `drug-screening-pipeline` and verify with `conda run -n drug-screening-pipeline python -c "import MDAnalysis; print(MDAnalysis.__version__)"`.

**ProLIF interaction analysis fails** — check `import prolif` in the main environment and run `python scripts/python/interaction_analysis.py` directly.

**No new ADMET predictions** — normally expected when every registered compound already has complete results for the current engine/version.

**No eligible ligands** — inspect `admet_screening` and `docking_eligibility`. This can be the correct outcome of the configured rules, not a software failure.

**`analyze_results.py` reports missing input** — inspect the current script and the files under `results/ranked/`; older versions referenced `docking_results.csv`, while the current collection output is `current_docking_results.csv`.

---

## 24. Reproducibility checklist

**Inputs**
- [ ] Target protein structure and chain identified
- [ ] Compound identity and structure/SMILES stored; external identifier recorded when available
- [ ] Docking site and configuration defined

**Methods**
- [ ] Properties, SCScore and ADMET prediction completed
- [ ] Screening decision and docking eligibility recorded
- [ ] Ligand and protein preparation succeeded
- [ ] Docking completed; poses and scores available
- [ ] Interaction analysis completed for selected poses

**Experiment tracking**
- [ ] Software and model versions recorded
- [ ] Docking configuration identified and its identity preserved
- [ ] Experiment identity stored in SQLite and linked to the correct target and ligand
- [ ] Changed parameters create distinguishable experiments

**Results**
- [ ] SQLite holds the persistent result; CSV reports come from the current run
- [ ] Historical results are not mixed with the current run
- [ ] Scores are treated as predictions, not measurements
- [ ] Manually selected compounds have a documented rationale

**Repository hygiene**
- [ ] No sensitive or private information; no large local-only payloads, temporary files or caches
- [ ] Configuration files and required scripts included
- [ ] Documentation matches the implemented workflow; claims limited to verified components

---

## 25. Development history and status

The project evolved from a target-specific A29L natural-compound screening workflow (in silico screening against Monkeypox A29L) into a reusable incremental pipeline.

- **Ingestion:** reusable local SDF ingestion, canonical-structure duplicate detection, incremental registration; tested with duplicate Caffeine and an Aspirin addition.
- **Properties:** expanded descriptors, PAINS/Brenk, SA Score, SCScore; incremental calculation.
- **ADMET:** ADMET-AI integration, engine wrapper, engine/version tracking, 104-endpoint storage, incremental prediction, dedicated environment.
- **Screening:** `screen_admet.py`, reconstructed and validated thresholds, flags, PASS/REVIEW/FAIL decisions, docking eligibility, incremental and idempotent behavior.
- **Docking:** eligibility integrated with ligand preparation, docking, result collection and analysis.
- **Interactions:** MDAnalysis and ProLIF integration, docked-ligand reconstruction from SDF plus pose, SQLite storage, CSV reporting.
- **Orchestration:** the runner expanded from the later docking stages to the full nine-stage workflow.

**Current stage:** a functional end-to-end research and portfolio prototype in a hardening phase, not a production pharmaceutical platform. Next priorities, in order: large-library SDF ingestion, change detection, configurable screening rules, ADMET engine/version replacement, provenance and reproducibility, validation of docking and interaction outputs, reporting.

---


## 25A. Git and repository hygiene

Do not commit transient or machine-specific artifacts merely because they appear in `git status`.

Examples that should normally remain uncommitted unless intentionally required include:

```text
*.before_*
temporary backups
local caches
large generated docking payloads
local environment files
```

Generated scientific results may be committed when they are intentionally part of the public demonstration, but the repository should distinguish source/configuration from large reproducible outputs.

Before committing:

```bash
git status --short --branch
git diff --stat
git diff -- README.md
```

Then stage only the intended files:

```bash
git add README.md
git add DEVELOPMENT_CHANGELOG.md
```

Use a focused commit message and push only after the final diff has been reviewed.

## 26. Contributing

Contributions that improve reproducibility, reliability, documentation or usability are welcome, for example: automation, validated methods, database queries and reporting, ADMET and docking workflows, interaction analysis, reproducibility checks, error handling and logging, tests, and efficiency.

Contributions should describe the change clearly, preserve reproducibility, document important parameters, identify software/model versions, avoid presenting unvalidated methods as experimentally established, avoid committing large generated files, and keep current implementation distinct from future features.

A new prediction or analysis method should document its name, version, model, input requirements, configuration, output format, how results are stored, and its effect on experiment identity.

**Reporting issues:** include the OS/WSL environment, Conda environment, Python and relevant software versions, command executed, error message, configuration, input type, and expected vs observed behaviour.

---

## 27. References

The project builds on established open-source software, databases and methods: RDKit, PubChem, Open Babel, SCScore, ADMET-AI, AutoDock Vina, Meeko, ProLIF, MDAnalysis, SQLite, Python, Conda, Git and GitHub. Consult each project's documentation and publications, and cite the original methods, software, models and data sources in derived work rather than treating this repository as their origin.

---


## 27A. Final operational checklist

Before treating a run as complete:

- [ ] Input SDF files are readable and intentionally retained.
- [ ] New/duplicate compound registration was verified.
- [ ] Compound properties are complete for new compounds.
- [ ] ADMET-AI results contain the expected engine/version and endpoint count.
- [ ] Screening flags and decisions are recorded.
- [ ] Docking eligibility matches the screening state.
- [ ] Eligible ligands were prepared successfully.
- [ ] Docking completed without unhandled errors.
- [ ] Docking results were collected into SQLite.
- [ ] Ranked outputs were generated.
- [ ] Interaction analysis completed for the intended poses.
- [ ] SQLite counts and important relationships were audited.
- [ ] Current configuration and software versions are known.
- [ ] No stale or temporary files are being mistaken for authoritative results.
- [ ] README claims match the implemented code and latest verified run.

## 28. License, author, contact and disclaimer

### License

Released under the **MIT License** (see `LICENSE`). Integrated third-party components (RDKit, AutoDock Vina, Meeko, Open Babel, ProLIF, MDAnalysis, ADMET-AI, SCScore, PubChem resources) remain under their own licenses and terms; consult upstream projects for licensing, citation and usage requirements. If this project or methodology is used in research, cite the relevant upstream software, models, databases and methods in addition to this repository.

### Author

**Sk. Adil Siraj** — developed as a computational biology and drug-discovery portfolio project demonstrating molecular informatics, bioinformatics, cheminformatics, ADMET prediction, structure-based virtual screening, molecular docking, protein–ligand interaction analysis, Python-based scientific software development, SQLite experiment tracking and reproducible computational workflows.

### Contact

- Repository: https://github.com/Adil07x/automated-drug-screening-pipeline
- GitHub: https://github.com/Adil07x
- Portfolio: https://adil07x.github.io/

For questions, collaboration or discussion, use the contact information on the GitHub profile and portfolio.

### Disclaimer

This project is for computational research, education and demonstration of reproducible scientific software. Its outputs must not be interpreted as experimental confirmation of binding, proof of biological activity, clinical evidence, drug-safety assessment, therapeutic efficacy, or regulatory approval or recommendation. Docking scores, ADMET predictions, physicochemical properties, SCScore values and interaction predictions depend on the underlying models, structures, parameters and software. Experimental validation and appropriate scientific review are required before drawing biological or therapeutic conclusions.
