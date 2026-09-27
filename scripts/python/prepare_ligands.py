#!/usr/bin/env python3

from pathlib import Path
import sqlite3

from rdkit import Chem
from rdkit.Chem import AllChem
from meeko import MoleculePreparation


# --------------------------------------------------
# Configuration
# --------------------------------------------------

PROJECT_ROOT = Path(__file__).resolve().parents[2]

DATABASE_FILE = (
    PROJECT_ROOT /
    "data/database/screening_database.sqlite"
)

OUTPUT_3D_DIR = (
    PROJECT_ROOT /
    "data/ligands/3d"
)

OUTPUT_PDBQT_DIR = (
    PROJECT_ROOT /
    "data/ligands/pdbqt"
)


# --------------------------------------------------
# Database
# --------------------------------------------------

def get_eligible_ligands():
    """
    Retrieve ligands that passed the automated
    ADMET screening eligibility layer.
    """

    connection = sqlite3.connect(DATABASE_FILE)

    query = """
        SELECT
            l.ligand_id,
            l.name,
            l.smiles
        FROM ligands l
        JOIN docking_eligibility de
            ON l.ligand_id = de.ligand_id
        WHERE de.eligibility = 'ELIGIBLE'
        ORDER BY l.ligand_id
    """

    try:
        cursor = connection.execute(query)
        rows = cursor.fetchall()
    finally:
        connection.close()

    return rows


# --------------------------------------------------
# Utility
# --------------------------------------------------

def safe_filename(name):
    """
    Convert a compound name into a filesystem-safe name.
    """

    safe = "".join(
        character if character.isalnum() or character in "._-"
        else "_"
        for character in name
    )

    return safe.strip("_")


# --------------------------------------------------
# Ligand preparation
# --------------------------------------------------

def prepare_ligand(ligand_id, name, smiles):
    """
    Convert a SMILES string into a 3D molecule and
    prepare it as PDBQT using Meeko.
    """

    print(f"[Ligand {ligand_id}] {name}")

    # ----------------------------------------------
    # Step 1: Validate SMILES
    # ----------------------------------------------

    mol = Chem.MolFromSmiles(smiles)

    if mol is None:
        raise ValueError(
            "RDKit could not parse the SMILES"
        )

    print(
        f"  SMILES validation: OK "
        f"({mol.GetNumHeavyAtoms()} heavy atoms)"
    )

    # ----------------------------------------------
    # Step 2: Add hydrogens
    # ----------------------------------------------

    mol = Chem.AddHs(mol)

    # ----------------------------------------------
    # Step 3: Generate 3D coordinates
    # ----------------------------------------------

    embedding_status = AllChem.EmbedMolecule(
        mol,
        randomSeed=42
    )

    if embedding_status != 0:
        raise ValueError(
            "RDKit failed to generate 3D coordinates"
        )

    print("  3D embedding: OK")

    # ----------------------------------------------
    # Step 4: Optimize geometry
    # ----------------------------------------------

    optimization_status = AllChem.MMFFOptimizeMolecule(
        mol
    )

    if optimization_status == -1:
        print(
            "  MMFF optimization: unavailable"
        )
    else:
        print(
            "  MMFF optimization: OK"
        )

    # ----------------------------------------------
    # Step 5: Save 3D SDF
    # ----------------------------------------------

    filename = safe_filename(name)

    sdf_file = (
        OUTPUT_3D_DIR /
        f"{ligand_id}_{filename}.sdf"
    )

    writer = Chem.SDWriter(str(sdf_file))
    writer.write(mol)
    writer.close()

    print(
        f"  3D SDF written: {sdf_file.name}"
    )

    # ----------------------------------------------
    # Step 6: Prepare PDBQT with Meeko
    # ----------------------------------------------

    preparation = MoleculePreparation()

    preparation.prepare(mol)

    pdbqt_file = (
        OUTPUT_PDBQT_DIR /
        f"{ligand_id}_{filename}.pdbqt"
    )

    preparation.write_pdbqt_file(
        str(pdbqt_file)
    )

    print(
        f"  PDBQT written: {pdbqt_file.name}"
    )

    return sdf_file, pdbqt_file


# --------------------------------------------------
# Main
# --------------------------------------------------

def main():

    OUTPUT_3D_DIR.mkdir(
        parents=True,
        exist_ok=True
    )

    OUTPUT_PDBQT_DIR.mkdir(
        parents=True,
        exist_ok=True
    )

    print("=" * 60)
    print("Automated Ligand Preparation")
    print("=" * 60)

    print(
        f"Database: {DATABASE_FILE}"
    )

    print()

    ligands = get_eligible_ligands()

    if not ligands:
        raise RuntimeError(
            "No ELIGIBLE ligands found in docking_eligibility"
        )

    print(
        f"Eligible ligands found: {len(ligands)}"
    )

    print()

    successful = 0
    failed = 0

    for ligand_id, name, smiles in ligands:

        try:

            prepare_ligand(
                ligand_id,
                name,
                smiles
            )

            successful += 1

        except Exception as error:

            print(
                f"  ERROR: {error}"
            )

            failed += 1

        print()

    print("=" * 60)
    print("Ligand preparation summary")
    print("=" * 60)

    print(
        f"Successful: {successful}"
    )

    print(
        f"Failed:     {failed}"
    )

    print(
        f"3D output:  {OUTPUT_3D_DIR}"
    )

    print(
        f"PDBQT output: {OUTPUT_PDBQT_DIR}"
    )


if __name__ == "__main__":
    main()