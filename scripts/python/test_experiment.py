from pathlib import Path
import sqlite3

PROJECT_ROOT = Path(__file__).resolve().parents[2]
DB_FILE = PROJECT_ROOT / "data/database/screening_database.sqlite"

connection = sqlite3.connect(DB_FILE)
cursor = connection.cursor()

target_id = 1
ligand_id = 1
config_id = 1

# Check whether this exact experiment already exists
cursor.execute("""
SELECT experiment_id, status
FROM experiments
WHERE target_id = ?
  AND ligand_id = ?
  AND config_id = ?
""", (target_id, ligand_id, config_id))

existing = cursor.fetchone()

if existing:
    print(f"Experiment already exists: ID={existing[0]}, status={existing[1]}")
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
        "completed"
    ))

    experiment_id = cursor.lastrowid
    connection.commit()

    print(f"New experiment created: ID={experiment_id}")

print("\nExperiments:")

for row in cursor.execute("""
    SELECT
        e.experiment_id,
        t.name AS target,
        l.name AS ligand,
        e.config_id,
        e.status
    FROM experiments e
    JOIN targets t ON e.target_id = t.target_id
    JOIN ligands l ON e.ligand_id = l.ligand_id
    ORDER BY e.experiment_id
"""):
    print(row)

connection.close()
