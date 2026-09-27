# Automated High-Throughput In-Silico Drug Screening Pipeline

## Project Overview

This project is a reproducible, database-driven computational drug screening pipeline designed to prioritize small molecules for further investigation against a protein target.

The pipeline integrates compound acquisition, molecular property calculation, synthetic-complexity assessment, ADMET prediction, rule-based screening, ligand preparation, molecular docking, protein–ligand interaction analysis, experiment tracking, and result reporting into a structured workflow.

The system is designed to support repeatable screening experiments rather than one-off docking runs. Molecular structures, screening results, docking configurations, experiment metadata, computational results, and interaction data are persisted in a SQLite database and linked to the corresponding input structures and computational settings.

The current implementation runs on Linux/WSL and is primarily controlled through Python and Bash scripts.

## Scientific Objective

The central scientific question is:

> **Which compounds in a screening library are computationally prioritized as potential binders of a selected protein target?**

The pipeline does not claim experimental binding, biological activity, therapeutic efficacy, or clinical safety. Instead, it provides a computational prioritization workflow that combines multiple filtering and structure-based analysis stages before compounds are selected for further investigation.

## Inputs

The pipeline is designed around three primary inputs:

1. **Protein target**
   - A suitable three-dimensional protein structure for docking.
   - Protein structures can be represented using PDB-compatible structural data.
   - A FASTA sequence alone is not sufficient for molecular docking; a suitable three-dimensional structure is required before docking.
   - The docking site is defined through a configurable docking-box configuration.

2. **Compound library**
   - Small-molecule structures represented using formats such as SDF and SMILES.
   - Compounds are registered in the project database with an internal compound ID.
   - Canonical or isomeric SMILES and external identifiers such as PubChem CID can be retained when available.
   - New compounds can be added to the existing master library without permanently removing previously registered compounds.

3. **Screening and docking configuration**
   - Molecular-property and filtering criteria.
   - ADMET screening configuration.
   - Docking engine and parameters.
   - Protein/ligand preparation settings.
   - Binding-site and docking-box parameters.

## Outputs

A screening run can produce:

- Molecular descriptors and physicochemical properties.
- Lipinski-rule calculations.
- SCScore synthetic-complexity estimates.
- ADMET prediction results.
- Computational screening decisions such as `PASS`, `REVIEW`, and `EXCLUDED`.
- Docking eligibility records.
- Prepared ligand structures and PDBQT files.
- Protein preparation records and receptor files.
- Docking poses and predicted docking scores.
- Protein–ligand interaction records generated with ProLIF.
- Ranked docking and screening reports.
- Persistent experiment records and configuration identifiers.
- CSV reports for analysis and data exchange.
- SQLite records providing the structured source of truth for the screening workflow.

The resulting data can therefore be used to identify compounds for subsequent computational analysis or experimental investigation.
## End-to-End Workflow

The pipeline follows a staged screening architecture:

```text
Compound Library
       │
       ▼
Compound Acquisition / Registration
       │
       ▼
Molecular Structure Validation
       │
       ▼
RDKit Molecular Properties
       │
       ├── Lipinski descriptors
       ├── Physicochemical properties
       ├── QED
       └── Structural filters
       │
       ▼
SCScore
       │
       ▼
ADMET-AI Prediction
       │
       ├── Absorption
       ├── Distribution
       ├── Metabolism
       ├── Excretion
       └── Toxicity
       │
       ▼
Screening Decision
(PASS / REVIEW / EXCLUDED)
       │
       ▼
Docking Eligibility
       │
       ▼
3D Ligand Preparation
       │
       ├── Hydrogen addition
       ├── 3D coordinate generation
       ├── Geometry optimization
       └── PDBQT preparation
       │
       ▼
Protein Preparation
       │
       ▼
AutoDock Vina
       │
       ├── Multiple poses
       └── Predicted docking scores
       │
       ▼
ProLIF Interaction Analysis
       │
       ▼
SQLite Experiment Database
       │
       ▼
Ranked Results / CSV Reports
## Implemented Capabilities

### Compound Management

The pipeline supports:

- Registration of compounds in a persistent SQLite database.
- Internal compound identifiers.
- Storage of molecular structures and SMILES.
- SHA-256 structure hashing.
- Detection of previously registered molecular structures.
- Addition of new compounds without permanently removing compounds from the master library.
- Compound acquisition using the PubChem REST API where applicable.

### Molecular Property Calculation

RDKit is used for molecular structure processing and property calculation, including:

- Molecular weight.
- Hydrogen-bond donors.
- Hydrogen-bond acceptors.
- LogP.
- Molar refractivity.
- Lipinski rule-of-five properties.
- Lipinski rule violations.
- QED.
- Additional molecular descriptors used by the screening workflow.
### Structure-Based Screening

Compounds selected for downstream analysis can proceed through:

1. Ligand preparation.
2. Protein preparation.
3. Molecular docking.
4. Docking-pose collection.
5. Docking-score ranking.
6. Protein–ligand interaction analysis.

### Persistent Results

The computational results are stored in a structured SQLite database rather than relying only on terminal output.

This keeps compounds, targets, configurations, experiments, docking results, ADMET predictions, and interaction results associated with their corresponding computational records.
## Molecular Complexity, ADMET and Screening

### SCScore

The pipeline incorporates **SCScore (Synthetic Complexity Score)** as part of compound-property evaluation.

SCScore provides a computational estimate of molecular synthetic complexity and is calculated locally using the SCScore model stored under:

```text
external/scscore/
```

The current implementation uses the NumPy-compatible SCScore model:

external/scscore/models/full_reaxys_model_1024bool/model.ckpt-10654.as_numpy.json.gz

SCScore is integrated into the compound-property calculation workflow and is executed locally using the available CPU environment.

### ADMET Prediction

The current ADMET prediction engine is **ADMET-AI**.

The pipeline stores individual ADMET endpoint predictions rather than reducing the entire prediction to a single score.

The current implementation evaluates a broad range of endpoints covering:

- Absorption
- Distribution
- Metabolism
- Excretion
- Toxicity
- Drug-likeness and physicochemical properties

The current demonstration evaluates **104 ADMET endpoints per compound**.

ADMET-AI runs in a dedicated Conda environment using a CPU-based PyTorch installation.

### Screening Categories

The pipeline combines selected ADMET and molecular-property signals into project-defined screening categories:

- `PASS`
- `REVIEW`
- `EXCLUDED`

These categories represent computational workflow decisions and are not clinical, regulatory, or experimental safety classifications.

The screening logic considers selected toxicity, absorption, physicochemical, and drug-likeness indicators.

The thresholds used by the project are computational screening rules and should not be interpreted as universal biological or regulatory cutoffs.

### Docking Eligibility

Screening results are translated into a separate docking-eligibility record.

This creates a clear separation between:

```text
ADMET / molecular-property assessment
                |
                v
       Screening decision
                |
                v
        Docking eligibility

The master compound library is not permanently modified when a compound fails a screening criterion.

Manual Compound Selection

The workflow also supports deliberate manual inclusion of compounds that would otherwise not be selected automatically.

A user can select:

Compounds passing the configured filters.
Compounds requiring review.
Compounds failing a particular screening rule.
A specific compound of interest.
A manually selected list of compounds.

Manual inclusion can be documented separately together with the scientific reason for overriding an automated screening decision.

This is important because computational filters such as Lipinski's rules are heuristics and should not be treated as absolute rules for biological activity or drug development.
## Reproducibility and Experiment History

A major design goal of the project is to make computational screening experiments reproducible and distinguishable from one another.

### SQLite as the Source of Truth

The primary structured data store is:

    data/database/screening_database.sqlite

The database stores information about:

- Compounds and ligands.
- Protein targets.
- Molecular properties.
- ADMET predictions.
- Screening decisions.
- Docking configurations.
- Experiments.
- Docking results.
- Protein–ligand interaction results.

This allows computational results to remain associated with the exact compound, target, configuration, and experiment that produced them.

### Configuration Identity

Important computational configurations are identified using SHA-256 hashes.

The project uses hashing for:

- Ligand molecular structures.
- Protein/receptor structures.
- Docking configurations.

The configuration identity includes important computational parameters such as:

- Docking engine.
- Docking-box centre.
- Docking-box dimensions.
- Exhaustiveness.
- Number of output poses.
- Engine/version information.

A configuration hash allows the pipeline to distinguish an existing configuration from a materially changed configuration.

### Experiment Identity

An experiment links the relevant:

    Target
       +
    Ligand
       +
    Docking Configuration
       =
    Experiment

This provides a persistent record of which target, compound, and computational configuration were used together.

### Repeatability

The workflow is designed to distinguish repeated calculations from calculations performed with changed inputs or parameters.

If the relevant inputs and configuration remain unchanged, previously recorded results can be identified for reference and comparison.

A changed input or important computational parameter creates a distinguishable experimental configuration.

Examples include changes to:

- Protein structure.
- Ligand structure.
- Ligand variant.
- Docking site.
- Docking-box dimensions.
- Docking engine.
- Docking parameters.

## Database Architecture

The SQLite database is organized into separate tables for the major stages of the screening workflow.

The current database includes tables for:

- Compounds and ligands.
- Protein targets.
- Molecular properties.
- ADMET predictions.
- ADMET screening decisions.
- Docking configurations.
- Docking experiments.
- Docking results.
- Docking eligibility.
- Protein–ligand interaction results.

This structure separates raw computational results from screening decisions and experimental records while maintaining relationships between them.

### Current Reporting Layer

A database view named `current_screening_results` provides an integrated reporting layer.

The view combines information from:

- Screening results.
- Docking experiments.
- Docking scores.
- Docking configuration.
- ADMET decisions.
- Protein–ligand interaction analysis.

The resulting report can therefore connect a compound's screening decision with its docking score and interaction summary.

### Data Exchange

The pipeline also exports selected results as CSV files for analysis and reporting.

Current reporting outputs include:

    results/ranked/current_docking_results.csv
    results/ranked/current_screening_results.csv
    results/interactions/current_interactions.csv

SQLite remains the structured source of truth, while CSV files provide convenient data exchange and human-readable reporting.

## Protein and Ligand Preparation

### Protein Structure

The current demonstration uses the protein structure:

    PDB ID: 1IEP
    Chain: A

The workflow uses a prepared three-dimensional receptor structure for docking.

A FASTA sequence alone is not treated as a docking-ready receptor. A suitable three-dimensional protein structure is required before structure-based docking can be performed.

### Docking Site Configuration

The docking site is represented by an explicit configuration rather than being permanently hard-coded into the docking script.

The current demonstration configuration is stored in:

    config/docking_config.json

The current configured docking box is:

    Center:
    X = 15.190
    Y = 53.903
    Z = 16.917

    Size:
    X = 19 Å
    Y = 27 Å
    Z = 24 Å

The current docking parameters include:

    Exhaustiveness = 8
    Number of poses = 10

The configured docking box is reused for the current screening demonstration.

### Ligand Preparation

Eligible ligands are prepared using RDKit and Meeko.

The preparation workflow includes:

1. Molecular structure validation.
2. Hydrogen addition.
3. Three-dimensional coordinate generation.
4. Geometry optimization.
5. SDF generation.
6. PDBQT preparation.

Prepared structures are associated with their corresponding internal ligand identifiers.

## Molecular Docking and Interaction Analysis

### AutoDock Vina

The current structure-based docking engine is AutoDock Vina.

For each selected ligand, the pipeline records:

- Docking configuration.
- Target.
- Ligand.
- Vina version.
- Docking log.
- Docked poses.
- Predicted docking scores.

Multiple poses can be retained for each docking experiment.

Docking scores are computational estimates generated by the Vina scoring function. They are not experimental binding affinities.

Within the same docking setup, a more negative predicted score is generally interpreted as more favorable according to the scoring function.

### Protein–Ligand Interaction Analysis

The pipeline uses ProLIF for post-docking protein–ligand interaction analysis.

The current interaction workflow uses:

- ProLIF 2.2.2.
- MDAnalysis 2.10.0.
- RDKit.

The current implementation analyzes the top-ranked docking pose and records detected interaction events in the SQLite database.

Recorded interaction information includes:

- Protein residue.
- Residue number.
- Protein chain where available.
- Interaction type.
- Distance where available.
- Ligand and protein atom information.
- Analysis engine and version.

This provides an additional structural interpretation layer beyond the docking score alone.
## Current Demonstration

The current demonstration uses a small compound set to demonstrate the complete computational workflow.

### Screening Dataset

The demonstration contains:

- 9 compounds.
- 104 ADMET endpoints per compound.
- 936 stored ADMET prediction records.

After the configured screening stage:

- 2 compounds were classified as currently eligible for docking.
- 2 docking experiments were completed.
- 10 poses were generated for each completed docking experiment.
- The top-ranked pose from each experiment was analyzed using ProLIF.
- 17 protein–ligand interaction events were recorded.

### Current Docking Results

| Compound | Screening Decision | Best Vina Score (kcal/mol) | Docking Poses | Unique Interaction Residues |
|---|---|---:|---:|---:|
| Oleic acid | PASS | -7.143 | 10 | 8 |
| Quinic acid | PASS | -5.967 | 10 | 7 |

### Interaction Summary

Oleic acid:

- 10 interaction events.
- 8 unique protein residues.
- Detected interaction categories include hydrophobic and van der Waals contacts.

Quinic acid:

- 7 interaction events.
- 7 unique protein residues.
- Detected interaction category: van der Waals contacts.

The interaction results are stored in the SQLite database and exported through:

    results/interactions/current_interactions.csv

### Current Reporting Outputs

The current workflow generates:

    results/ranked/current_docking_results.csv
    results/ranked/current_screening_results.csv
    results/interactions/current_interactions.csv

The integrated screening report combines screening, docking, and interaction information for the current completed experiments.

### Interpretation

These results demonstrate that the computational workflow can progress from compound screening through docking and interaction analysis.

They should not be interpreted as experimental confirmation of binding or biological activity.

The docking scores represent computational predictions, while the interaction analysis describes contacts detected in the selected computational docking poses.

## Software Environments

The project uses separate Conda environments for components with different software dependencies.

### Core Pipeline Environment

    Environment: drug-screening-pipeline

The core environment contains the principal cheminformatics and docking tools used by the pipeline.

Verified runtime:

- Python 3.12
- RDKit 2025.09.5
- NumPy 2.5.3
- Pandas 3.0.6
- AutoDock Vina f458505-mod
- Open Babel 3.2.1
- Meeko 0.8.0
- SCScore local NumPy-compatible implementation

### ADMET Environment

    Environment: admet-ai

The ADMET environment contains ADMET-AI and its machine-learning dependencies.

Verified runtime:

- PyTorch 2.14.0+cpu
- CUDA: not available in the current runtime

ADMET-AI is executed in a separate environment because its dependency stack is distinct from the core docking environment.

### Interaction Analysis Environment

    Environment: prolif-env

The interaction-analysis environment contains the molecular interaction analysis stack.

Verified runtime:

- ProLIF 2.2.2
- MDAnalysis 2.10.0
- RDKit 2026.03.1
- NumPy 2.5.3
- Pandas 3.0.6
- SciPy 1.18.1
- Gemmi 0.7.5

ProLIF is maintained in a separate environment to avoid dependency conflicts with the core docking workflow.

### Environment Separation

The environments are intentionally separated:

- `drug-screening-pipeline` → compound processing, SCScore, ligand preparation, protein preparation, and docking
- `admet-ai` → ADMET-AI prediction
- `prolif-env` → ProLIF interaction analysis

This separation keeps the main pipeline reproducible while allowing specialized tools to use their own compatible dependency stacks.

## Technology Stack

### Core Technologies

| Technology | Role |
|---|---|
| Python | Pipeline implementation and scientific processing |
| Bash | Command-line workflow and environment operations |
| SQLite | Persistent experiment and result database |
| RDKit | Molecular structures, descriptors, Lipinski rules, QED and 3D chemistry |
| SCScore | Synthetic-complexity estimation |
| ADMET-AI | ADMET prediction |
| AutoDock Vina | Molecular docking |
| Meeko | Ligand preparation and PDBQT generation |
| ProLIF | Protein–ligand interaction analysis |
| MDAnalysis | Molecular structure handling for interaction analysis |

### Supporting Technologies

| Technology | Role |
|---|---|
| PyTorch | Machine-learning runtime used by ADMET-AI |
| Open Babel | Molecular structure conversion and cheminformatics support |
| PubChem REST API | Compound acquisition and external compound information |
| Pandas | Tabular data processing and CSV reporting |
| NumPy | Numerical computation and SCScore support |
| Matplotlib | Scientific visualization |
| Conda | Environment and dependency management |
| Git / GitHub | Version control and project distribution |
| SHA-256 | Molecular, receptor and configuration identity |

### Docking Engine Scope

The current reproducible docking workflow uses AutoDock Vina.

GNINA, smina, QVina and CNN-based docking scoring are not currently part of the verified reproducible pipeline.

They may be considered as future alternative docking and scoring engines rather than being represented as currently implemented components.
## Installation

The project is intended to run in a Linux environment. The development setup used for the current implementation is Ubuntu under WSL2 on Windows.

### Clone the Repository

    git clone git@github.com:Adil07x/automated-drug-screening-pipeline.git
    cd automated-drug-screening-pipeline

### Conda Environments

The project uses separate Conda environments for the core pipeline, ADMET prediction, and interaction analysis.

Activate the core environment:

    conda activate drug-screening-pipeline

Verify the principal tools:

    python --version
    vina --version
    obabel -V

Verify RDKit:

    python -c "from rdkit import Chem; print('RDKit OK')"

### Verify the ADMET Environment

    conda activate admet-ai

Verify PyTorch:

    python -c "import torch; print('PyTorch:', torch.__version__); print('CUDA available:', torch.cuda.is_available())"

The current validated ADMET environment uses CPU execution.

### Verify the Interaction Environment

    conda activate prolif-env

Verify ProLIF and MDAnalysis:

    python -c "import prolif, MDAnalysis; print('ProLIF:', prolif.__version__); print('MDAnalysis:', MDAnalysis.__version__)"

## Running the Pipeline

The pipeline is composed of individual reproducible stages.

Typical processing follows:

    Compound registration
            |
            v
    Molecular property calculation
            |
            v
    SCScore
            |
            v
    ADMET-AI
            |
            v
    Screening / docking eligibility
            |
            v
    Ligand preparation
            |
            v
    Protein preparation
            |
            v
    Docking
            |
            v
    Interaction analysis
            |
            v
    Result collection

Individual stages are implemented through the Python scripts under:

    scripts/python/

The SQLite database should be treated as the persistent record of the computational experiment.

Generated scientific artifacts should remain associated with the corresponding compound, target, configuration, and experiment.

## Reproducible Docking Configuration

The current demonstration uses:

    config/docking_config.json

Important parameters include:

- Receptor structure.
- Ligand preparation.
- Docking-box centre.
- Docking-box dimensions.
- Exhaustiveness.
- Number of poses.
- Docking engine and version.

Changing an important docking parameter should result in a distinguishable configuration and experiment rather than silently replacing an existing result.

## Project Structure

The repository is organized around reproducible computational screening, persistent experiment records, and generated scientific results.

    automated-drug-screening-pipeline/
    |
    ├── config/
    │   └── docking_config.json
    |
    ├── data/
    │   ├── database/
    │   │   └── screening_database.sqlite
    │   ├── ligands/
    │   │   ├── 3d/
    │   │   └── pdbqt/
    │   ├── metadata/
    │   │   └── compound_library.csv
    │   ├── processed/
    │   │   ├── ligands/
    │   │   └── protein/
    │   └── raw/
    │       ├── compounds/
    │       └── protein/
    |
    ├── external/
    │   └── scscore/
    |
    ├── results/
    │   ├── docking/
    │   ├── interactions/
    │   └── ranked/
    |
    ├── scripts/
    │   └── python/
    │       ├── calculate_compound_properties.py
    │       ├── calculate_lipinski.py
    │       ├── collect_results.py
    │       ├── create_database.py
    │       ├── interaction_analysis.py
    │       ├── prepare_ligands.py
    │       ├── register_ligands.py
    │       ├── run_admet_ai.py
    │       ├── run_docking.py
    │       ├── scscore_runner.py
    │       └── ...
    |
    ├── .gitignore
    ├── environment/environment.yml
    ├── LICENSE
    └── README.md

### Directory Roles

`config/`

Contains reproducible configuration files defining computational settings such as the docking box and AutoDock Vina parameters.

`data/`

Contains structured project data, molecular structures, prepared inputs, and the SQLite experiment database.

`data/raw/`

Contains source protein and compound structures before downstream processing.

`data/processed/`

Contains prepared protein and ligand structures used by computational stages.

`data/ligands/`

Contains generated three-dimensional ligand structures and docking-ready PDBQT files.

`data/database/`

Contains the SQLite source-of-truth database for compounds, molecular properties, ADMET results, screening decisions, docking experiments, docking results, and protein–ligand interaction results.

`external/scscore/`

Contains the locally integrated SCScore source and supporting model resources required for SCScore calculation.

Large SCScore model and dataset payloads are excluded from version control where appropriate through `.gitignore`.

`results/`

Contains generated computational results and reports.

`results/docking/`

Contains docking output files and logs associated with computational experiments.

`results/interactions/`

Contains protein–ligand interaction analysis results generated using ProLIF.

`results/ranked/`

Contains current ranked docking and integrated screening CSV reports.

`scripts/python/`

Contains the Python components implementing the computational workflow.

### Data Flow

    Raw Structures
          |
          v
    Compound / Target Registration
          |
          v
    Molecular Properties
          |
          v
    SCScore
          |
          v
    ADMET-AI
          |
          v
    Screening and Docking Eligibility
          |
          v
    Ligand / Protein Preparation
          |
          v
    AutoDock Vina
          |
          v
    ProLIF Interaction Analysis
          |
          v
    SQLite Experiment Database
          |
          v
    Ranked Results / CSV Reports

## Scientific Interpretation and Limitations

The results produced by this pipeline are computational screening results and should be interpreted as prioritization evidence rather than experimental validation.

### Docking Scores

AutoDock Vina docking scores are computational estimates generated by a scoring function.

Within the same target, ligand preparation procedure, docking configuration, and scoring setup, a more negative predicted score is generally interpreted as more favorable according to that scoring function.

Docking scores should not be treated as experimentally measured binding affinities or as direct predictions of biological activity.

### ADMET Screening

The ADMET screening categories used by this project are computational decision rules defined for this workflow.

The categories:

- `PASS`
- `REVIEW`
- `EXCLUDED`

are intended to organize compounds for downstream computational screening.

They are not clinical, regulatory, or experimental safety classifications.

A compound classified as `PASS` has not been demonstrated to be safe, effective, bioavailable, or therapeutically useful.

Similarly, a compound classified as `EXCLUDED` has not been experimentally demonstrated to be toxic or biologically inactive.

### Lipinski and Physicochemical Filters

Lipinski-style rules and physicochemical thresholds are used as practical screening heuristics.

These rules help prioritize compounds for computational analysis but do not determine whether a compound is biologically active.

Compounds that fail a screening rule can therefore be retained in the master compound library and can be manually selected for analysis when scientifically justified.

### Protein Structure Requirements

Molecular docking requires an appropriate three-dimensional receptor structure.

A FASTA sequence represents the protein sequence and is not, by itself, a docking-ready three-dimensional receptor.

The current workflow therefore expects a suitable prepared 3D protein structure and a defined binding region.

Automatic protein structure prediction is not currently part of the implemented pipeline.

### Interaction Analysis

ProLIF interaction analysis provides a computational description of interactions detected in the analyzed docking pose.

The current demonstration analyzes the top-ranked docking pose for each completed experiment.

Interaction detection does not establish that the predicted complex is experimentally stable or biologically functional.

### Current Demonstration Scope

The current demonstration contains a relatively small compound set compared with the intended high-throughput use of the architecture.

The demonstrated workflow currently includes:

- 9 registered compounds.
- 104 ADMET endpoints per compound.
- 936 stored ADMET endpoint results.
- 2 compounds currently eligible for docking.
- 2 completed docking experiments.
- 10 docking poses generated per experiment.
- Top-pose interaction analysis using ProLIF.
- 17 recorded interaction events.

The architecture is designed to scale to larger compound libraries, but practical throughput depends on available CPU/GPU resources, molecular complexity, software settings, and the number of compounds processed.

### Reproducibility

Computational results depend on the exact input structures, software versions, models, configuration parameters, and docking settings used during an experiment.

The project therefore stores experiment and configuration information in SQLite so that computational runs can be distinguished when important inputs or parameters change.

Results from different configurations should not be directly compared without considering the changes in computational setup.

### Overall Interpretation

The pipeline should be viewed as a computational prioritization framework:

    Compound Library
          |
          v
    Molecular / ADMET Screening
          |
          v
    Docking
          |
          v
    Interaction Analysis
          |
          v
    Computational Prioritization
          |
          v
    Experimental Validation

The final determination of binding, biological activity, pharmacokinetics, toxicity, and therapeutic potential requires appropriate experimental and/or validated external evidence.
## Future Development

The current implementation provides a reproducible computational screening workflow. Future development will focus on increasing automation, scalability, analysis depth, and usability without changing the core experiment-tracking architecture.

### Automated End-to-End Execution

The pipeline can be extended toward a single orchestrated workflow that automatically moves compounds through:

    Compound Ingestion
          |
          v
    Property Calculation
          |
          v
    SCScore
          |
          v
    ADMET Screening
          |
          v
    Docking Eligibility
          |
          v
    Ligand Preparation
          |
          v
    Protein Preparation
          |
          v
    Molecular Docking
          |
          v
    Interaction Analysis
          |
          v
    Ranking and Reporting

The existing individual stages are intended to remain independently executable so that failed or incomplete stages can be resumed without unnecessarily repeating previous calculations.

### Larger Compound Libraries

The architecture can be applied to substantially larger compound collections.

Future improvements may include:

- Batch compound ingestion.
- Incremental processing of newly added compounds.
- Automatic detection of previously processed compounds.
- Parallel processing where appropriate.
- Resumable screening runs.
- Automated tracking of failed compounds and processing errors.

### Expanded ADMET Engine Support

ADMET-AI is the current implemented ADMET engine.

The project can later support additional prediction engines through the existing engine-oriented architecture.

Different ADMET engines can be compared while preserving the identity of the computational method used for each result.

### Advanced Docking Engines

AutoDock Vina is the currently used docking engine.

Future versions may evaluate alternative docking or scoring approaches such as:

- GNINA.
- smina.
- QVina variants.
- Alternative Vina scoring configurations.
- Machine-learning-based rescoring methods.

These will be treated as separate computational methods and will not replace the current validated Vina workflow unless independently verified.

### Expanded Pose Analysis

The current interaction-analysis workflow focuses on the top-ranked docking pose.

Future development can extend this to:

- Multiple selected poses.
- Configurable top-N pose analysis.
- Pose clustering.
- Interaction-frequency analysis.
- Residue-level comparison between compounds.
- Interaction fingerprints.
- Comparative binding-site analysis.

### Improved Ranking

Future ranking can combine multiple computational signals rather than relying on docking score alone.

Potential ranking features include:

- ADMET screening status.
- Physicochemical properties.
- SCScore.
- Docking score.
- Interaction count.
- Interaction type.
- Contacted binding-site residues.
- Pose consistency.

Any combined ranking methodology should remain explicitly documented so that the resulting prioritization is reproducible and interpretable.

### Reporting and Visualization

Future versions can provide automated reports containing:

- Compound summaries.
- ADMET summaries.
- Docking rankings.
- Interaction fingerprints.
- Binding-site residue summaries.
- Configuration information.
- Experiment history.
- Quality-control information.

Additional visualization may include docking-score distributions, interaction maps, compound-property plots, and comparative screening summaries.

### Query and Export Tools

The SQLite database provides a foundation for future analytical queries and automated exports.

Potential additions include:

- Command-line result queries.
- Filtered CSV exports.
- Excel reports.
- Summary dashboards.
- Experiment comparison reports.
- Compound history reports.
- Target-specific screening reports.

### Scalability and Compute Optimization

Future work can investigate:

- Parallel docking execution.
- GPU-enabled computational methods where supported.
- Batch scheduling.
- Larger compound libraries.
- Resource-aware execution.
- Checkpointing and resumable workflows.

The actual computational capacity will depend on the available hardware and the requirements of the selected software and models.

### Experimental Validation Integration

The long-term objective is to make computational prioritization easier to connect with experimental evidence.

Future database extensions could store externally generated validation information such as:

- Experimental binding measurements.
- Assay results.
- Biological activity measurements.
- Literature evidence.
- Validation status.

Such information would remain clearly distinguished from computational predictions.

### Reproducible Method Expansion

Any new computational method should be recorded with its relevant:

- Software name.
- Software version.
- Model version.
- Input structure.
- Configuration.
- Parameters.
- Execution status.
- Output.
- Experiment identity.

This preserves the central principle of the project: computational results should remain traceable to the exact method and inputs that produced them.
## License

This project is released under the MIT License.

The MIT License permits reuse, modification, distribution, and use of the software, subject to the conditions specified in the license.

See the `LICENSE` file in the repository for the complete license text.

### Scientific Software and External Components

This project integrates or uses several open-source scientific software packages and external resources.

These components remain subject to their respective licenses and terms of use.

Examples include:

- RDKit
- AutoDock Vina
- Meeko
- Open Babel
- ProLIF
- MDAnalysis
- ADMET-AI
- SCScore
- PubChem resources

Users should consult the respective upstream projects for their licensing, citation, and usage requirements.

### Data and Scientific Results

Molecular structures, computational predictions, docking results, and interaction-analysis results generated by this project should be interpreted according to the scientific limitations described in this README.

Computational results are provided for research, educational, and portfolio purposes and should not be interpreted as experimental or clinical evidence.

### Citation

If this project or its methodology is used in research or derivative work, users should cite the relevant upstream scientific software, models, databases, and methods in addition to referencing this repository.
## Project Status

The core computational workflow is currently implemented and demonstrated on a small compound set.

### Currently Implemented

- Compound registration and persistent compound identity.
- Molecular property calculation using RDKit.
- Lipinski-style property evaluation.
- SCScore calculation using the integrated local SCScore implementation.
- ADMET prediction using ADMET-AI.
- Project-defined ADMET screening categories.
- Docking eligibility tracking.
- Manual compound selection for computational screening.
- 3D ligand preparation using RDKit.
- Docking-ready ligand preparation using Meeko.
- Protein preparation for docking.
- Configurable AutoDock Vina docking.
- Persistent docking experiment records.
- Docking pose and score storage.
- ProLIF-based interaction analysis.
- Interaction-result storage in SQLite.
- Integrated screening-result reporting.
- Current docking and screening CSV exports.
- Reproducible configuration tracking using configuration identity and experiment records.

### Current Demonstration

The current validated demonstration contains:

- 9 registered compounds.
- 936 ADMET endpoint results.
- 2 compounds currently eligible for docking.
- 2 completed docking experiments.
- 20 total docking poses.
- 2 top poses analyzed using ProLIF.
- 17 recorded interaction events.
- Persistent results stored in SQLite.

### Development Stage

The project is currently a functional computational research and portfolio prototype rather than a production pharmaceutical screening platform.

The architecture is intentionally designed so that additional compounds, targets, computational methods, and analysis stages can be incorporated without replacing the existing experiment-tracking system.

Further development will focus on automation, larger-scale screening, expanded interaction analysis, additional computational methods, reporting, and validation-oriented workflows.
## Author

**Sk. Adil Siraj**

This project was developed as a computational biology and drug-discovery portfolio project to demonstrate the practical integration of:

- Molecular informatics
- Bioinformatics
- Cheminformatics
- ADMET prediction
- Structure-based virtual screening
- Molecular docking
- Protein–ligand interaction analysis
- Python-based scientific software development
- SQLite-based experiment tracking
- Reproducible computational workflows

The project reflects an implementation-focused approach to applying computational methods to biological and pharmaceutical research problems.

### Repository

GitHub:

    https://github.com/Adil07x/automated-drug-screening-pipeline

The repository is named `automated-drug-screening-pipeline`.

### Portfolio

    https://adil07x.github.io/

The portfolio provides additional information about the author's computational biology, data analysis, programming, and research-oriented projects.
## Disclaimer

This project is intended for computational research, education, and demonstration of reproducible scientific software development.

The computational predictions generated by this pipeline should not be interpreted as:

- Experimental confirmation of molecular binding.
- Proof of biological activity.
- Clinical evidence.
- Drug-safety assessment.
- Therapeutic efficacy.
- Regulatory approval or recommendation.

Docking scores, ADMET predictions, physicochemical properties, SCScore values, and protein–ligand interaction predictions are computational outputs whose interpretation depends on the underlying models, structures, parameters, and software implementations.

Experimental validation and appropriate scientific review are required before drawing biological or therapeutic conclusions from computational predictions.
## Contact

For questions, collaboration, or discussion related to the project, please use the contact information available through the author's GitHub and professional portfolio.

### GitHub

    https://github.com/Adil07x

### Portfolio

    https://adil07x.github.io/

The repository contains the source code, configuration, documentation, and reproducible computational workflow for the project.
## Reproducibility Checklist

Before considering a screening run complete, the following items should be identifiable and traceable.

### Input Verification

- [ ] Target protein structure is identified.
- [ ] Protein chain used for docking is identified.
- [ ] Compound identity is stored in the database.
- [ ] Molecular structure or SMILES is available.
- [ ] External compound identifier is recorded when available.
- [ ] Docking site and configuration are defined.

### Computational Method Verification

- [ ] Molecular-property calculation completed.
- [ ] SCScore calculation completed where applicable.
- [ ] ADMET prediction completed where applicable.
- [ ] Screening decision recorded.
- [ ] Docking eligibility recorded.
- [ ] Ligand preparation completed successfully.
- [ ] Protein preparation completed successfully.
- [ ] Docking completed successfully.
- [ ] Docking poses and scores are available.
- [ ] Interaction analysis completed for selected poses.

### Experiment Tracking

- [ ] Software and model versions are recorded where applicable.
- [ ] Docking configuration is identified.
- [ ] Configuration identity is preserved.
- [ ] Experiment identity is stored in SQLite.
- [ ] Results are associated with the correct target and ligand.
- [ ] Changed computational parameters create distinguishable experiment records.

### Result Verification

- [ ] SQLite contains the persistent result.
- [ ] Current CSV reports were generated from the current results.
- [ ] Historical results are not accidentally mixed with the current run.
- [ ] Computational scores are interpreted as predictions rather than experimental measurements.
- [ ] Any manually selected compound has an appropriate selection rationale.

### Publication and Repository Hygiene

- [ ] Sensitive or private information is excluded.
- [ ] Large local-only model/data payloads are excluded where appropriate.
- [ ] Temporary files and caches are excluded.
- [ ] Reproducible configuration files are included.
- [ ] Required scripts are included.
- [ ] Documentation describes the actual implemented workflow.
- [ ] Claims about software or methods are limited to components that have been verified in the project.
## Contributing

Contributions that improve the scientific reproducibility, computational reliability, documentation, or usability of the project are welcome.

### Areas for Contribution

Potential areas include:

- Improving pipeline automation.
- Adding validated computational methods.
- Improving database queries and reporting.
- Extending ADMET analysis.
- Improving docking workflows.
- Expanding interaction analysis.
- Adding reproducibility checks.
- Improving error handling and logging.
- Adding tests for individual pipeline components.
- Improving documentation.
- Improving computational efficiency.

### Contribution Principles

Contributions should:

- Clearly describe the change being introduced.
- Preserve reproducibility where possible.
- Document important computational parameters.
- Identify software or model versions when relevant.
- Avoid presenting unvalidated computational methods as experimentally established.
- Avoid committing unnecessary large generated files or local development artifacts.
- Preserve the distinction between current implementation and future experimental features.

### Scientific Method Changes

When introducing a new prediction or analysis method, the implementation should document:

- Method or software name.
- Version.
- Model, when applicable.
- Input requirements.
- Configuration parameters.
- Output format.
- How results are stored.
- How the method affects experiment identity and reproducibility.

### Reporting Issues

When reporting a problem, include enough information to reproduce it where possible, such as:

- Operating system or WSL environment.
- Conda environment.
- Python version.
- Relevant software versions.
- Command executed.
- Error message.
- Relevant configuration.
- Input type.
- Expected behaviour.
- Observed behaviour.
## References and Resources

The project builds on established open-source scientific software, databases, and computational methods.

### Molecular and Chemical Informatics

- RDKit
- PubChem
- Open Babel

### Synthetic Complexity

- SCScore

### ADMET Prediction

- ADMET-AI

### Molecular Docking

- AutoDock Vina
- Meeko

### Protein–Ligand Interaction Analysis

- ProLIF
- MDAnalysis

### Data and Reproducibility

- SQLite
- Python
- Conda
- Git
- GitHub

Users reproducing or extending this project should consult the official documentation and scientific publications associated with each software package, model, and database.

Where appropriate, derived research outputs should cite the original methods, software, models, and public data sources rather than treating this repository as the original source of those underlying technologies.

### Reproducibility Principle

The project aims to make computational results traceable to:

    Input Structure
          |
          v
    Software / Model
          |
          v
    Configuration
          |
          v
    Experiment
          |
          v
    Computational Result

This allows future analyses to distinguish between results generated by different inputs, software versions, models, configurations, and experimental records.