#!/usr/bin/env python3

"""
Protein-ligand interaction analysis using ProLIF.

Workflow:
    SQLite
      ↓
    completed docking experiments
      ↓
    best docking pose
      ↓
    Vina PDBQT MODEL 1
      ↓
    reconstruct chemically correct ligand SDF
      ↓
    ProLIF interaction analysis
      ↓
    results/interactions/current_interactions.csv

The original ligand SDF and receptor PDB are never modified.
"""

from pathlib import Path
import csv
import sqlite3
import tempfile

import MDAnalysis as mda
import prolif as plf
from rdkit import Chem


# ---------------------------------------------------------------------
# Project paths
# ---------------------------------------------------------------------

PROJECT_ROOT = Path(__file__).resolve().parents[2]

DATABASE = PROJECT_ROOT / "data/database/screening_database.sqlite"

PROTEIN_FILE = (
    PROJECT_ROOT / "data/processed/protein/1IEP_chainA.pdb"
)

LIGAND_SDF_DIR = PROJECT_ROOT / "data/ligands/3d"

OUTPUT_DIR = PROJECT_ROOT / "results/interactions"

OUTPUT_CSV = OUTPUT_DIR / "current_interactions.csv"


# ---------------------------------------------------------------------
# Database
# ---------------------------------------------------------------------

def get_completed_experiments():
    """
    Return the best docking pose for every completed experiment.
    """

    connection = sqlite3.connect(DATABASE)
    connection.row_factory = sqlite3.Row

    query = """
        SELECT
            e.experiment_id,
            e.ligand_id,
            l.name AS ligand_name,
            dr.pose_number,
            dr.affinity_kcal_mol,
            dr.pose_file
        FROM experiments e
        JOIN ligands l
            ON e.ligand_id = l.ligand_id
        JOIN docking_results dr
            ON e.experiment_id = dr.experiment_id
        WHERE e.status = 'completed'
          AND dr.pose_number = 1
        ORDER BY e.experiment_id
    """

    rows = connection.execute(query).fetchall()
    connection.close()

    return rows


# ---------------------------------------------------------------------
# File helpers
# ---------------------------------------------------------------------

def ligand_filename(ligand_id, ligand_name):
    """
    Build the expected original ligand SDF filename.
    """

    safe_name = ligand_name.replace(" ", "_")

    return LIGAND_SDF_DIR / f"{ligand_id}_{safe_name}.sdf"


def extract_model_one(pdbqt_file):
    """
    Extract MODEL 1 from a multi-model Vina PDBQT file.
    """

    with open(pdbqt_file, "r") as handle:
        lines = handle.readlines()

    model_lines = []
    inside_model = False

    for line in lines:

        if line.startswith("MODEL 1"):
            inside_model = True
            continue

        if inside_model and line.startswith("ENDMDL"):
            break

        if inside_model:
            model_lines.append(line)

    if not model_lines:
        raise RuntimeError(
            f"MODEL 1 was not found in {pdbqt_file}"
        )

    return model_lines, lines


def get_vina_smiles(all_lines):
    """
    Extract the SMILES stored by Vina/Meeko in the PDBQT.
    """

    for line in all_lines:
        if line.startswith("REMARK SMILES "):
            return line.strip()[14:]

    raise RuntimeError(
        "REMARK SMILES was not found in the PDBQT."
    )


def get_vina_heavy_atom_coordinates(model_lines):
    """
    Extract coordinates for heavy atoms from MODEL 1.

    PDBQT hydrogen atom names may appear as H while the AutoDock
    atom type is HD, so hydrogen identification is based on the
    atom-name field rather than the element field.
    """

    atoms = []

    for line in model_lines:

        if not line.startswith("ATOM"):
            continue

        atom_name = line[12:16].strip()

        # Skip explicit hydrogen atoms.
        if atom_name.upper() == "H":
            continue

        x = float(line[30:38])
        y = float(line[38:46])
        z = float(line[46:54])

        atoms.append(
            {
                "atom_name": atom_name,
                "x": x,
                "y": y,
                "z": z,
            }
        )

    return atoms


# ---------------------------------------------------------------------
# Ligand reconstruction
# ---------------------------------------------------------------------

def reconstruct_docked_ligand(
    original_sdf,
    pdbqt_file,
    output_sdf,
):
    """
    Preserve the original ligand chemistry from SDF while transferring
    the docked heavy-atom coordinates from Vina MODEL 1.

    This avoids relying on PDBQT for bond orders and explicit-H chemistry.
    """

    # -------------------------------------------------------------
    # Read original ligand SDF
    # -------------------------------------------------------------

    supplier = Chem.SDMolSupplier(
        str(original_sdf),
        removeHs=False,
    )

    template = supplier[0]

    if template is None:
        raise RuntimeError(
            f"Could not read ligand SDF: {original_sdf}"
        )

    # -------------------------------------------------------------
    # Read Vina pose
    # -------------------------------------------------------------

    model_lines, all_lines = extract_model_one(pdbqt_file)

    vina_atoms = get_vina_heavy_atom_coordinates(
        model_lines
    )

    vina_smiles = get_vina_smiles(all_lines)

    # -------------------------------------------------------------
    # Validate atom counts
    # -------------------------------------------------------------

    template_heavy = Chem.RemoveHs(template)

    template_heavy_count = template.GetNumHeavyAtoms()
    vina_heavy_count = len(vina_atoms)

    if vina_heavy_count != template_heavy_count:
        raise RuntimeError(
            f"Heavy atom count mismatch for {original_sdf.name}: "
            f"Vina={vina_heavy_count}, "
            f"SDF={template_heavy_count}"
        )

    # -------------------------------------------------------------
    # Parse Vina SMILES
    # -------------------------------------------------------------

    vina_mol = Chem.MolFromSmiles(vina_smiles)

    if vina_mol is None:
        raise RuntimeError(
            f"Could not parse Vina SMILES: {vina_smiles}"
        )

    # -------------------------------------------------------------
    # Map Vina atom order to original SDF heavy-atom order
    # -------------------------------------------------------------

    matches = vina_mol.GetSubstructMatches(
        template_heavy
    )

    if len(matches) != 1:
        raise RuntimeError(
            f"Expected exactly one atom mapping for "
            f"{original_sdf.name}; found {len(matches)}"
        )

    mapping = matches[0]

    # -------------------------------------------------------------
    # Copy original chemically correct molecule
    # -------------------------------------------------------------

    docked = Chem.Mol(template)

    conformer = docked.GetConformer()

    # mapping[vina_atom_index] gives the corresponding
    # original-SDF heavy-atom index.
    for vina_index, template_index in enumerate(mapping):

        coordinates = vina_atoms[vina_index]

        conformer.SetAtomPosition(
            template_index,
            (
                coordinates["x"],
                coordinates["y"],
                coordinates["z"],
            ),
        )

    # -------------------------------------------------------------
    # Write reconstructed SDF
    # -------------------------------------------------------------

    writer = Chem.SDWriter(str(output_sdf))
    writer.write(docked)
    writer.close()

    return docked


# ---------------------------------------------------------------------
# ProLIF
# ---------------------------------------------------------------------

def calculate_interactions(docked_sdf):
    """
    Calculate ProLIF interaction metadata for one docked ligand.
    """

    # -------------------------------------------------------------
    # Protein
    # -------------------------------------------------------------

    protein_universe = mda.Universe(
        str(PROTEIN_FILE)
    )

    # The receptor PDB has no explicit hydrogens.
    protein = plf.Molecule.from_mda(
        protein_universe,
        "protein",
        inferrer=None,
    )

    # -------------------------------------------------------------
    # Ligand
    # -------------------------------------------------------------

    supplier = Chem.SDMolSupplier(
        str(docked_sdf),
        removeHs=False,
    )

    ligand_mol = supplier[0]

    if ligand_mol is None:
        raise RuntimeError(
            f"Could not read reconstructed ligand: {docked_sdf}"
        )

    ligand = plf.Molecule.from_rdkit(
        ligand_mol
    )

    # -------------------------------------------------------------
    # ProLIF fingerprint
    # -------------------------------------------------------------

    fingerprint = plf.Fingerprint()

    interaction_fingerprint = fingerprint.generate(
        ligand,
        protein,
        metadata=True,
    )

    return interaction_fingerprint

# ---------------------------------------------------------------------
# SQLite output
# ---------------------------------------------------------------------

def write_interactions_to_database(output_rows):
    """
    Replace current ProLIF interaction results in SQLite.

    The CSV remains the export file, while SQLite stores the same
    structured interaction events for downstream analysis.
    """

    connection = sqlite3.connect(DATABASE)

    try:
        connection.execute("PRAGMA foreign_keys = ON")

        # Remove previous results generated by this analysis.
        connection.execute(
            "DELETE FROM interaction_results"
        )

        insert_query = """
            INSERT INTO interaction_results (
                experiment_id,
                ligand_id,
                pose_number,
                protein_residue,
                protein_residue_number,
                protein_chain,
                interaction_type,
                distance_angstrom,
                ligand_atom_indices,
                protein_atom_indices,
                ligand_parent_indices,
                protein_parent_indices,
                engine,
                engine_version
            )
            VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        """

        rows = []

        for row in output_rows:
            rows.append(
                (
                    row["experiment_id"],
                    row["ligand_id"],
                    1,
                    row["protein_residue"],
                    row["protein_residue_number"],
                    row["protein_chain"],
                    row["interaction_type"],
                    row["distance_angstrom"],
                    row["ligand_atom_indices"],
                    row["protein_atom_indices"],
                    row["ligand_parent_indices"],
                    row["protein_parent_indices"],
                    "ProLIF",
                    plf.__version__,
                )
            )

        connection.executemany(
            insert_query,
            rows,
        )

        connection.commit()

        return len(rows)

    except Exception:
        connection.rollback()
        raise

    finally:
        connection.close()
# ---------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------

def main():

    print("=" * 70)
    print("ProLIF Protein-Ligand Interaction Analysis")
    print("=" * 70)

    print(f"Database: {DATABASE}")
    print(f"Protein:  {PROTEIN_FILE}")
    print()

    # -------------------------------------------------------------
    # Find completed experiments
    # -------------------------------------------------------------

    experiments = get_completed_experiments()

    print(
        f"Completed experiments with pose 1: "
        f"{len(experiments)}"
    )
    print()

    if not experiments:
        print("No completed docking experiments found.")
        return

    OUTPUT_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    # CSV rows
    output_rows = []

    successful = 0
    failed = 0

    # -------------------------------------------------------------
    # Process each experiment
    # -------------------------------------------------------------

    for experiment in experiments:

        experiment_id = experiment["experiment_id"]
        ligand_id = experiment["ligand_id"]
        ligand_name = experiment["ligand_name"]
        affinity = experiment["affinity_kcal_mol"]

        print("-" * 70)
        print(
            f"Experiment {experiment_id}: "
            f"{ligand_name} "
            f"(ligand_id={ligand_id})"
        )

        try:

            # -----------------------------------------------------
            # Locate files
            # -----------------------------------------------------

            pose_file = PROJECT_ROOT / experiment["pose_file"]

            original_sdf = ligand_filename(
                ligand_id,
                ligand_name,
            )

            if not pose_file.exists():
                raise FileNotFoundError(
                    f"Docking pose not found: {pose_file}"
                )

            if not original_sdf.exists():
                raise FileNotFoundError(
                    f"Original ligand SDF not found: "
                    f"{original_sdf}"
                )

            print(f"  Docking pose: {pose_file.name}")
            print(f"  Original SDF: {original_sdf.name}")
            print(
                f"  Best affinity: "
                f"{affinity:.3f} kcal/mol"
            )

            # -----------------------------------------------------
            # Temporary reconstructed SDF
            # -----------------------------------------------------

            with tempfile.TemporaryDirectory() as temp_dir:

                reconstructed_sdf = (
                    Path(temp_dir)
                    / f"{ligand_id}_docked_model1.sdf"
                )

                docked_molecule = reconstruct_docked_ligand(
                    original_sdf,
                    pose_file,
                    reconstructed_sdf,
                )

                print(
                    f"  Reconstructed ligand: "
                    f"{docked_molecule.GetNumAtoms()} atoms"
                )

                # -------------------------------------------------
                # ProLIF
                # -------------------------------------------------

                ifp = calculate_interactions(
                    reconstructed_sdf
                )

            interaction_count = 0

            # -----------------------------------------------------
            # Convert ProLIF metadata into flat CSV rows
            # -----------------------------------------------------

            for residue_pair, interactions in ifp.items():

                ligand_residue = residue_pair[0]
                protein_residue = residue_pair[1]

                for interaction_name, metadata_entries in (
                    interactions.items()
                ):

                    for metadata in metadata_entries:

                        distance = metadata.get(
                            "distance"
                        )

                        ligand_indices = metadata[
                            "indices"
                        ].get("ligand", ())

                        protein_indices = metadata[
                            "indices"
                        ].get("protein", ())

                        parent_ligand_indices = metadata[
                            "parent_indices"
                        ].get("ligand", ())

                        parent_protein_indices = metadata[
                            "parent_indices"
                        ].get("protein", ())

                        output_rows.append(
                            {
                                "experiment_id": experiment_id,
                                "ligand_id": ligand_id,
                                "ligand_name": ligand_name,
                                "docking_affinity_kcal_mol": affinity,
                                "ligand_residue": str(
                                    ligand_residue
                                ),
                                "protein_residue": (
                                    protein_residue.name
                                ),
                                "protein_residue_number": (
                                    protein_residue.number
                                ),
                                "protein_chain": (
                                    protein_residue.chain
                                ),
                                "interaction_type": (
                                    interaction_name
                                ),
                                "distance_angstrom": distance,
                                "ligand_atom_indices": ",".join(
                                    map(
                                        str,
                                        ligand_indices,
                                    )
                                ),
                                "protein_atom_indices": ",".join(
                                    map(
                                        str,
                                        protein_indices,
                                    )
                                ),
                                "ligand_parent_indices": ",".join(
                                    map(
                                        str,
                                        parent_ligand_indices,
                                    )
                                ),
                                "protein_parent_indices": ",".join(
                                    map(
                                        str,
                                        parent_protein_indices,
                                    )
                                ),
                            }
                        )

                        interaction_count += 1

            print(
                f"  Residue pairs: {len(ifp)}"
            )
            print(
                f"  Interaction events: "
                f"{interaction_count}"
            )

            successful += 1

        except Exception as error:

            failed += 1

            print(
                f"  ERROR: {error}"
            )

    # -------------------------------------------------------------
    # Write CSV
    # -------------------------------------------------------------

    fieldnames = [
        "experiment_id",
        "ligand_id",
        "ligand_name",
        "docking_affinity_kcal_mol",
        "ligand_residue",
        "protein_residue",
        "protein_residue_number",
        "protein_chain",
        "interaction_type",
        "distance_angstrom",
        "ligand_atom_indices",
        "protein_atom_indices",
        "ligand_parent_indices",
        "protein_parent_indices",
    ]

    with open(
        OUTPUT_CSV,
        "w",
        newline="",
    ) as handle:

        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
        )

        writer.writeheader()
        writer.writerows(output_rows)

    # -------------------------------------------------------------
    # Write SQLite
    # -------------------------------------------------------------

    database_rows = write_interactions_to_database(
        output_rows
    )

    print(
        f"SQLite interaction rows written: "
        f"{database_rows}"
    )
    # -------------------------------------------------------------
    # Summary
    # -------------------------------------------------------------

    print()
    print("=" * 70)
    print("Interaction analysis complete")
    print("=" * 70)
    print(f"Successful experiments: {successful}")
    print(f"Failed experiments:     {failed}")
    print(f"Interaction rows:       {len(output_rows)}")
    print(f"CSV written:            {OUTPUT_CSV}")
    print()


if __name__ == "__main__":
    main()
