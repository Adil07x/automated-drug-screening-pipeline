from pathlib import Path
import sqlite3

PROJECT_ROOT = Path(__file__).resolve().parents[2]
DB_FILE = PROJECT_ROOT / "data/database/screening_database.sqlite"

DB_FILE.parent.mkdir(parents=True, exist_ok=True)

connection = sqlite3.connect(DB_FILE)
cursor = connection.cursor()

cursor.executescript("""
CREATE TABLE IF NOT EXISTS targets (
    target_id INTEGER PRIMARY KEY AUTOINCREMENT,
    name TEXT NOT NULL,
    pdb_id TEXT,
    chain TEXT,
    receptor_file TEXT,
    receptor_hash TEXT,
    created_at TEXT DEFAULT CURRENT_TIMESTAMP
);

CREATE TABLE IF NOT EXISTS ligands (
    ligand_id INTEGER PRIMARY KEY AUTOINCREMENT,
    name TEXT NOT NULL,
    source TEXT,
    source_id TEXT,
    smiles TEXT,
    structure_hash TEXT NOT NULL,
    source_file TEXT,
    created_at TEXT DEFAULT CURRENT_TIMESTAMP
);

CREATE TABLE IF NOT EXISTS compound_properties (
    property_id INTEGER PRIMARY KEY AUTOINCREMENT,
    ligand_id INTEGER NOT NULL,
    molecular_weight REAL,
    logp REAL,
    hbd INTEGER,
    hba INTEGER,
    lipinski_violations INTEGER,
    admet_data TEXT,
    calculated_at TEXT DEFAULT CURRENT_TIMESTAMP,
    FOREIGN KEY (ligand_id) REFERENCES ligands(ligand_id)
);

CREATE TABLE IF NOT EXISTS docking_configs (
    config_id INTEGER PRIMARY KEY AUTOINCREMENT,
    target_id INTEGER NOT NULL,
    center_x REAL NOT NULL,
    center_y REAL NOT NULL,
    center_z REAL NOT NULL,
    size_x REAL NOT NULL,
    size_y REAL NOT NULL,
    size_z REAL NOT NULL,
    exhaustiveness INTEGER,
    num_modes INTEGER,
    vina_version TEXT,
    config_hash TEXT NOT NULL,
    created_at TEXT DEFAULT CURRENT_TIMESTAMP,
    FOREIGN KEY (target_id) REFERENCES targets(target_id)
);

CREATE TABLE IF NOT EXISTS experiments (
    experiment_id INTEGER PRIMARY KEY AUTOINCREMENT,
    target_id INTEGER NOT NULL,
    ligand_id INTEGER NOT NULL,
    config_id INTEGER NOT NULL,
    status TEXT NOT NULL,
    started_at TEXT,
    completed_at TEXT,
    FOREIGN KEY (target_id) REFERENCES targets(target_id),
    FOREIGN KEY (ligand_id) REFERENCES ligands(ligand_id),
    FOREIGN KEY (config_id) REFERENCES docking_configs(config_id),
    UNIQUE(target_id, ligand_id, config_id)
);

CREATE TABLE IF NOT EXISTS docking_results (
    result_id INTEGER PRIMARY KEY AUTOINCREMENT,
    experiment_id INTEGER NOT NULL,
    pose_number INTEGER NOT NULL,
    affinity_kcal_mol REAL,
    rmsd_lb REAL,
    rmsd_ub REAL,
    pose_file TEXT,
    FOREIGN KEY (experiment_id) REFERENCES experiments(experiment_id)
);
""")

connection.commit()
connection.close()

print(f"Database created: {DB_FILE}")
