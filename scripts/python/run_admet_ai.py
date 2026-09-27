import json
import sqlite3
from pathlib import Path

from rdkit import Chem

from admet.admet_ai_engine import ADMETAIEngine


# ------------------------------------------------------------
# Project paths
# ------------------------------------------------------------

PROJECT_ROOT = Path(__file__).resolve().parents[2]

DB_PATH = PROJECT_ROOT / "data/database/screening_database.sqlite"


# ------------------------------------------------------------
# Database helpers
# ------------------------------------------------------------

def get_connection():
    """Create a connection to the screening database."""
    connection = sqlite3.connect(DB_PATH)
    connection.execute("PRAGMA foreign_keys = ON")
    return connection


def get_ligands(connection):
    """
    Return all ligands currently stored in the database.
    """
    query = """
        SELECT ligand_id, name, source_file
        FROM ligands
        ORDER BY ligand_id
    """

    return connection.execute(query).fetchall()


# ------------------------------------------------------------
# Molecular structure helper
# ------------------------------------------------------------

def smiles_from_sdf(sdf_path):
    """
    Read the first molecule from an SDF file and convert it
    to an isomeric SMILES string.
    """

    if not sdf_path.exists():
        raise FileNotFoundError(f"SDF file not found: {sdf_path}")

    supplier = Chem.SDMolSupplier(
        str(sdf_path),
        removeHs=False,
    )

    molecule = supplier[0]

    if molecule is None:
        raise ValueError(f"Could not read molecule from: {sdf_path}")

    smiles = Chem.MolToSmiles(
        molecule,
        isomericSmiles=True,
    )

    return smiles


# ------------------------------------------------------------
# Database insertion
# ------------------------------------------------------------

def replace_admet_results(
    connection,
    ligand_id,
    engine_name,
    engine_version,
    records,
):
    """
    Replace existing results from a particular engine for
    one ligand.

    This makes the pipeline safely re-runnable.
    """

    connection.execute(
        """
        DELETE FROM admet_results
        WHERE ligand_id = ?
          AND engine = ?
        """,
        (ligand_id, engine_name),
    )

    insert_query = """
        INSERT INTO admet_results (
            ligand_id,
            engine,
            engine_version,
            endpoint,
            value,
            unit,
            prediction_type,
            uncertainty,
            status,
            raw_result
        )
        VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
    """

    for record in records:
        connection.execute(
            insert_query,
            (
                ligand_id,
                engine_name,
                engine_version,
                record["endpoint"],
                record["value"],
                record["unit"],
                record["prediction_type"],
                record["uncertainty"],
                record["status"],
                record["raw_result"],
            ),
        )

    connection.commit()


# ------------------------------------------------------------
# Main pipeline
# ------------------------------------------------------------

def main():

    print("=" * 70)
    print("ADMET-AI DATABASE PIPELINE")
    print("=" * 70)

    print(f"Database: {DB_PATH}")
    print()

    # Load ADMET-AI once.
    # This is important because loading the model for every
    # compound would be unnecessarily slow.
    engine = ADMETAIEngine()

    connection = get_connection()

    try:
        ligands = get_ligands(connection)

        print()
        print(f"Compounds found in database: {len(ligands)}")
        print()

        if not ligands:
            print("No ligands found.")
            return

        successful = 0
        failed = 0

        for index, (ligand_id, name, source_file) in enumerate(
            ligands,
            start=1,
        ):

            print("-" * 70)
            print(f"[{index}/{len(ligands)}] {name}")
            print(f"Ligand ID: {ligand_id}")
            print(f"Source: {source_file}")

            try:

                # Convert database source path to an absolute path.
                sdf_path = PROJECT_ROOT / source_file

                # Read molecule and generate SMILES.
                smiles = smiles_from_sdf(sdf_path)

                print(f"SMILES: {smiles}")

                # Run ADMET-AI.
                predictions = engine.predict(smiles)

                print(
                    f"ADMET-AI endpoints returned: "
                    f"{len(predictions)}"
                )

                # Convert predictions to database records.
                records = engine.to_database_records(
                    predictions
                )

                # Replace previous ADMET-AI results.
                replace_admet_results(
                    connection=connection,
                    ligand_id=ligand_id,
                    engine_name="ADMET-AI",
                    engine_version=engine.engine_version,
                    records=records,
                )

                print(
                    f"Stored in database: "
                    f"{len(records)} endpoint records"
                )

                print("Status: SUCCESS")

                successful += 1

            except Exception as error:

                print(f"Status: FAILED")
                print(f"Error: {error}")

                failed += 1

        print()
        print("=" * 70)
        print("ADMET-AI PIPELINE COMPLETE")
        print("=" * 70)

        print(f"Successful: {successful}")
        print(f"Failed:     {failed}")

    finally:
        connection.close()


if __name__ == "__main__":
    main()
