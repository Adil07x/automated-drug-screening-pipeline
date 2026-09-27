from pathlib import Path
import sqlite3
import hashlib

PROJECT_ROOT = Path(__file__).resolve().parents[2]
DB_FILE = PROJECT_ROOT / "data/database/screening_database.sqlite"

TARGET_NAME = "1IEP_chainA"
PDB_ID = "1IEP"
CHAIN = "A"

RECEPTOR_FILE = (
    PROJECT_ROOT / "data/processed/protein/1IEP_chainA.pdbqt"
)

if not RECEPTOR_FILE.exists():
    raise FileNotFoundError(
        f"Receptor file not found: {RECEPTOR_FILE}"
    )

# Calculate a fingerprint of the actual receptor file.
receptor_hash = hashlib.sha256(
    RECEPTOR_FILE.read_bytes()
).hexdigest()

connection = sqlite3.connect(DB_FILE)
cursor = connection.cursor()

# Check whether this exact receptor is already registered.
cursor.execute("""
SELECT target_id, receptor_hash
FROM targets
WHERE name = ?
""", (TARGET_NAME,))

existing = cursor.fetchone()

if existing:
    target_id, existing_hash = existing

    if existing_hash == receptor_hash:
        print(
            f"Target already registered: ID={target_id}"
        )
    else:
        print(
            f"Target name exists, but receptor file has changed."
        )
        print(
            "A new target version will be registered."
        )

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
            TARGET_NAME,
            PDB_ID,
            CHAIN,
            str(RECEPTOR_FILE.relative_to(PROJECT_ROOT)),
            receptor_hash
        ))

        target_id = cursor.lastrowid
        connection.commit()

        print(f"New target version registered: ID={target_id}")

else:
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
        TARGET_NAME,
        PDB_ID,
        CHAIN,
        str(RECEPTOR_FILE.relative_to(PROJECT_ROOT)),
        receptor_hash
    ))

    target_id = cursor.lastrowid
    connection.commit()

    print(f"New target registered: ID={target_id}")

print(f"Receptor hash: {receptor_hash}")

connection.close()
