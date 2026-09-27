from pathlib import Path
import sqlite3

from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen, Lipinski

PROJECT_ROOT = Path(__file__).resolve().parents[2]
DB_FILE = PROJECT_ROOT / "data/database/screening_database.sqlite"

connection = sqlite3.connect(DB_FILE)
cursor = connection.cursor()

cursor.execute("""
SELECT ligand_id, name, source_file
FROM ligands
ORDER BY ligand_id
""")

ligands = cursor.fetchall()

if not ligands:
    raise RuntimeError("No ligands registered in the database.")

for ligand_id, name, source_file in ligands:

    sdf_file = PROJECT_ROOT / source_file

    supplier = Chem.SDMolSupplier(
        str(sdf_file),
        removeHs=False
    )

    mol = supplier[0]

    if mol is None:
        print(f"SKIPPED: {name} — invalid molecule")
        continue

    molecular_weight = Descriptors.MolWt(mol)
    logp = Crippen.MolLogP(mol)
    hbd = Lipinski.NumHDonors(mol)
    hba = Lipinski.NumHAcceptors(mol)

    # Lipinski Rule of Five violations:
    # MW > 500
    # LogP > 5
    # HBD > 5
    # HBA > 10
    violations = 0

    if molecular_weight > 500:
        violations += 1

    if logp > 5:
        violations += 1

    if hbd > 5:
        violations += 1

    if hba > 10:
        violations += 1

    # Remove an existing property record for this ligand.
    cursor.execute("""
    DELETE FROM compound_properties
    WHERE ligand_id = ?
    """, (ligand_id,))

    cursor.execute("""
    INSERT INTO compound_properties (
        ligand_id,
        molecular_weight,
        logp,
        hbd,
        hba,
        lipinski_violations
    )
    VALUES (?, ?, ?, ?, ?, ?)
    """, (
        ligand_id,
        molecular_weight,
        logp,
        hbd,
        hba,
        violations
    ))

    connection.commit()

    print(
        f"{name}: "
        f"MW={molecular_weight:.2f}, "
        f"LogP={logp:.2f}, "
        f"HBD={hbd}, "
        f"HBA={hba}, "
        f"Violations={violations}"
    )

connection.close()

print("\nLipinski calculations stored in database.")
