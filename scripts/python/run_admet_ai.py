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


def get_ligands(connection, engine_name, engine_version):
    """
    Return ligands that require ADMET-AI prediction.

    A ligand requires processing when:
    1. It has no ADMET-AI results.
    2. Its stored ADMET-AI version differs from the current version.
    3. Its current-version results are incomplete.

    The current complete endpoint count is determined from the
    existing database state rather than being hard-coded.
    """

    complete_endpoint_count = connection.execute(
        """
        SELECT MAX(endpoint_count)
        FROM (
            SELECT
                ligand_id,
                COUNT(DISTINCT endpoint) AS endpoint_count
            FROM admet_results
            WHERE engine = ?
              AND engine_version = ?
            GROUP BY ligand_id
        )
        """,
        (engine_name, engine_version),
    ).fetchone()[0]

    if complete_endpoint_count is None:
        complete_endpoint_count = 0

    query = """
        SELECT
            l.ligand_id,
            l.name,
            l.source_file
        FROM ligands l
        LEFT JOIN (
            SELECT
                ligand_id,
                engine_version,
                COUNT(DISTINCT endpoint) AS endpoint_count
            FROM admet_results
            WHERE engine = ?
            GROUP BY ligand_id, engine_version
        ) a
            ON l.ligand_id = a.ligand_id
        WHERE
            a.ligand_id IS NULL
            OR a.engine_version != ?
            OR a.endpoint_count < ?
        ORDER BY l.ligand_id
    """

    ligands = connection.execute(
        query,
        (
            engine_name,
            engine_version,
            complete_endpoint_count,
        ),
    ).fetchall()

    return ligands, complete_endpoint_count


def get_admet_status(connection, engine_name, engine_version):
    """
    Return the current ADMET-AI status for every ligand.

    This is used for reporting and validation.
    """

    query = """
        SELECT
            l.ligand_id,
            l.name,
            COALESCE(a.engine_version, 'NONE') AS engine_version,
            COALESCE(a.endpoint_count, 0) AS endpoint_count
        FROM ligands l
        LEFT JOIN (
            SELECT
                ligand_id,
                engine_version,
                COUNT(DISTINCT endpoint) AS endpoint_count
            FROM admet_results
            WHERE engine = ?
            GROUP BY ligand_id, engine_version
        ) a
            ON l.ligand_id = a.ligand_id
        ORDER BY l.ligand_id
    """

    return connection.execute(
        query,
        (engine_name,),
    ).fetchall()


# ------------------------------------------------------------
# Molecular structure helper
# ------------------------------------------------------------

def smiles_from_sdf(sdf_path):
    """
    Read the first molecule from an SDF file and convert it
    to an isomeric SMILES string.
    """

    if not sdf_path.exists():
        raise FileNotFoundError(
            f"SDF file not found: {sdf_path}"
        )

    supplier = Chem.SDMolSupplier(
        str(sdf_path),
        removeHs=False,
    )

    molecule = supplier[0]

    if molecule is None:
        raise ValueError(
            f"Could not read molecule from: {sdf_path}"
        )

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
        (
            ligand_id,
            engine_name,
        ),
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
    print("ADMET-AI INCREMENTAL DATABASE PIPELINE")
    print("=" * 70)

    print(f"Database: {DB_PATH}")
    print()

    # Load ADMET-AI once.
    engine = ADMETAIEngine()

    engine_name = "ADMET-AI"
    engine_version = engine.engine_version

    connection = get_connection()

    try:

        # Determine which compounds actually require processing.
        ligands, complete_endpoint_count = get_ligands(
            connection,
            engine_name,
            engine_version,
        )

        print()
        print(f"Engine: {engine_name}")
        print(f"Engine version: {engine_version}")
        print(
            "Expected complete endpoint count: "
            f"{complete_endpoint_count}"
        )
        print()

        print(
            "Compounds requiring ADMET prediction: "
            f"{len(ligands)}"
        )
        print()

        if not ligands:
            print(
                "All compounds already have complete "
                "ADMET-AI results."
            )
            return

        successful = 0
        failed = 0

        for index, (ligand_id, name, source_file) in enumerate(
            ligands,
            start=1,
        ):

            print("-" * 70)
            print(
                f"[{index}/{len(ligands)}] "
                f"{name}"
            )
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
                    "ADMET-AI endpoints returned: "
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
                    engine_name=engine_name,
                    engine_version=engine_version,
                    records=records,
                )

                print(
                    "Stored in database: "
                    f"{len(records)} endpoint records"
                )

                print("Status: SUCCESS")

                successful += 1

            except Exception as error:

                print("Status: FAILED")
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
