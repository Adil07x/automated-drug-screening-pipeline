from pathlib import Path
import sqlite3
import hashlib
import json

PROJECT_ROOT = Path(__file__).resolve().parents[2]
DB_FILE = PROJECT_ROOT / "data/database/screening_database.sqlite"

# Same target and ligand, but the X dimension of the box is changed.
config = {
    "target": "1IEP_chainA",
    "center_x": 15.190,
    "center_y": 53.903,
    "center_z": 16.917,
    "size_x": 20.0,
    "size_y": 27.0,
    "size_z": 24.0,
    "exhaustiveness": 8,
    "num_modes": 10,
    "vina_version": "AutoDock Vina f458505-mod"
}

config_text = json.dumps(config, sort_keys=True)

config_hash = hashlib.sha256(
    config_text.encode("utf-8")
).hexdigest()

connection = sqlite3.connect(DB_FILE)
cursor = connection.cursor()

# Find the target
cursor.execute("""
SELECT target_id
FROM targets
WHERE name = ?
""", ("1IEP_chainA",))

target_id = cursor.fetchone()[0]

# Check whether this exact configuration already exists
cursor.execute("""
SELECT config_id
FROM docking_configs
WHERE config_hash = ?
""", (config_hash,))

existing = cursor.fetchone()

if existing:
    config_id = existing[0]
    print(f"Configuration already exists: ID={config_id}")
else:
    cursor.execute("""
    INSERT INTO docking_configs (
        target_id,
        center_x,
        center_y,
        center_z,
        size_x,
        size_y,
        size_z,
        exhaustiveness,
        num_modes,
        vina_version,
        config_hash
    )
    VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
    """, (
        target_id,
        config["center_x"],
        config["center_y"],
        config["center_z"],
        config["size_x"],
        config["size_y"],
        config["size_z"],
        config["exhaustiveness"],
        config["num_modes"],
        config["vina_version"],
        config_hash
    ))

    config_id = cursor.lastrowid
    connection.commit()

    print(f"New configuration created: ID={config_id}")

print(f"Configuration hash: {config_hash}")

# Find the Quercetin ligand
cursor.execute("""
SELECT ligand_id
FROM ligands
WHERE name = ?
""", ("Quercetin",))

ligand_id = cursor.fetchone()[0]

# Check whether this target + ligand + configuration
# combination already exists.
cursor.execute("""
SELECT experiment_id
FROM experiments
WHERE target_id = ?
  AND ligand_id = ?
  AND config_id = ?
""", (target_id, ligand_id, config_id))

existing_experiment = cursor.fetchone()

if existing_experiment:
    print(f"Experiment already exists: ID={existing_experiment[0]}")
else:
    cursor.execute("""
    INSERT INTO experiments (
        target_id,
        ligand_id,
        config_id,
        status
    )
    VALUES (?, ?, ?, ?)
    """, (
        target_id,
        ligand_id,
        config_id,
        "pending"
    ))

    experiment_id = cursor.lastrowid
    connection.commit()

    print(f"New experiment created: ID={experiment_id}")

print("\nAll experiments:")

for row in cursor.execute("""
    SELECT
        e.experiment_id,
        t.name AS target,
        l.name AS ligand,
        e.config_id,
        dc.size_x,
        dc.size_y,
        dc.size_z,
        e.status
    FROM experiments e
    JOIN targets t
        ON e.target_id = t.target_id
    JOIN ligands l
        ON e.ligand_id = l.ligand_id
    JOIN docking_configs dc
        ON e.config_id = dc.config_id
    ORDER BY e.experiment_id
"""):
    print(row)

connection.close()
