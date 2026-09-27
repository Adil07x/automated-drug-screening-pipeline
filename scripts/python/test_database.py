from pathlib import Path
import sqlite3
import hashlib

PROJECT_ROOT = Path(__file__).resolve().parents[2]
DB_FILE = PROJECT_ROOT / "data/database/screening_database.sqlite"
LIGAND_FILE = PROJECT_ROOT / "data/raw/compounds/Quercetin.sdf"

connection = sqlite3.connect(DB_FILE)
cursor = connection.cursor()

# --------------------------------------------------
# 1. Add the target
# --------------------------------------------------

receptor_file = PROJECT_ROOT / "data/processed/protein/1IEP_chainA.pdbqt"

receptor_hash = hashlib.sha256(
    receptor_file.read_bytes()
).hexdigest()

cursor.execute("""
INSERT INTO targets (
    name,
    pdb_id,
    chain,
    receptor_file,
    receptor_hash
)
VALUES (?, ?, ?, ?, ?)
""", (
    "1IEP_chainA",
    "1IEP",
    "A",
    str(receptor_file.relative_to(PROJECT_ROOT)),
    receptor_hash
))

target_id = cursor.lastrowid

# --------------------------------------------------
# 2. Add the ligand
# --------------------------------------------------

ligand_hash = hashlib.sha256(
    LIGAND_FILE.read_bytes()
).hexdigest()

cursor.execute("""
INSERT INTO ligands (
    name,
    source,
    source_id,
    structure_hash,
    source_file
)
VALUES (?, ?, ?, ?, ?)
""", (
    "Quercetin",
    "local",
    None,
    ligand_hash,
    str(LIGAND_FILE.relative_to(PROJECT_ROOT))
))

ligand_id = cursor.lastrowid

connection.commit()

print(f"Target inserted: ID={target_id}")
print(f"Ligand inserted: ID={ligand_id}")

# --------------------------------------------------
# 3. Display the records
# --------------------------------------------------

print("\nTargets:")
for row in cursor.execute("""
    SELECT target_id, name, pdb_id, chain
    FROM targets
"""):
    print(row)

print("\nLigands:")
for row in cursor.execute("""
    SELECT ligand_id, name, source, source_file
    FROM ligands
"""):
    print(row)

connection.close()
