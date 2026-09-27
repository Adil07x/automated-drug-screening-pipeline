from pathlib import Path
import sqlite3
import hashlib

from rdkit import Chem

PROJECT_ROOT = Path(__file__).resolve().parents[2]
DB_FILE = PROJECT_ROOT / "data/database/screening_database.sqlite"
LIGAND_DIR = PROJECT_ROOT / "data/raw/compounds"

if not LIGAND_DIR.exists():
    raise FileNotFoundError(
        f"Ligand directory not found: {LIGAND_DIR}"
    )

connection = sqlite3.connect(DB_FILE)
cursor = connection.cursor()

sdf_files = sorted(LIGAND_DIR.glob("*.sdf"))

if not sdf_files:
    raise RuntimeError("No SDF files found.")

print(f"Found {len(sdf_files)} SDF files.\n")

registered = 0
skipped = 0

for sdf_file in sdf_files:

    supplier = Chem.SDMolSupplier(
        str(sdf_file),
        removeHs=False
    )

    mol = supplier[0]

    if mol is None:
        print(f"SKIPPED: {sdf_file.name} — invalid molecule")
        skipped += 1
        continue

    # Generate canonical SMILES from the molecular structure.
    smiles = Chem.MolToSmiles(
        mol,
        canonical=True
    )

    # Use the canonical structure as the identity.
    structure_hash = hashlib.sha256(
        smiles.encode("utf-8")
    ).hexdigest()

    compound_name = sdf_file.stem

    # Check whether this exact structure already exists.
    cursor.execute("""
    SELECT ligand_id, name
    FROM ligands
    WHERE structure_hash = ?
    """, (structure_hash,))

    existing = cursor.fetchone()

    if existing:
        print(
            f"EXISTS: {compound_name} "
            f"→ ligand ID {existing[0]}"
        )
        skipped += 1
        continue

    cursor.execute("""
    INSERT INTO ligands (
        name,
        source,
        source_id,
        smiles,
        structure_hash,
        source_file
    )
    VALUES (?, ?, ?, ?, ?, ?)
    """, (
        compound_name,
        "local_sdf",
        None,
        smiles,
        structure_hash,
        str(sdf_file.relative_to(PROJECT_ROOT))
    ))

    ligand_id = cursor.lastrowid

    connection.commit()

    print(
        f"REGISTERED: {compound_name} "
        f"→ ligand ID {ligand_id}"
    )

    registered += 1

print("\n--------------------------------")
print(f"New ligands registered: {registered}")
print(f"Already existing/skipped: {skipped}")
print("--------------------------------")

print("\nDatabase ligands:")

for row in cursor.execute("""
    SELECT ligand_id, name
    FROM ligands
    ORDER BY ligand_id
"""):
    print(row)

connection.close()
