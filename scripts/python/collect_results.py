from pathlib import Path
import re
import csv
import sqlite3


# ============================================================
# Configuration
# ============================================================

PROJECT_ROOT = Path(__file__).resolve().parents[2]

DB_FILE = PROJECT_ROOT / "data/database/screening_database.sqlite"

DOCKING_DIR = PROJECT_ROOT / "results/docking"
OUTPUT_DIR = PROJECT_ROOT / "results/ranked"

# IMPORTANT:
# Keep the historical docking_results.csv untouched.
OUTPUT_FILE = OUTPUT_DIR / "current_docking_results.csv"

# Integrated current-run report from SQLite VIEW.
SCREENING_OUTPUT_FILE = OUTPUT_DIR / "current_screening_results.csv"


# ============================================================
# Extract docking results from one Vina log
# ============================================================

def parse_vina_log(log_file):

    results = []

    with open(log_file, "r") as file:
        lines = file.readlines()

    for line in lines:

        match = re.match(
            r"\s*(\d+)\s+"
            r"(-?\d+(?:\.\d+)?)\s+"
            r"(-?\d+(?:\.\d+)?)\s+"
            r"(-?\d+(?:\.\d+)?)",
            line
        )

        if match:

            mode = int(match.group(1))
            affinity = float(match.group(2))
            rmsd_lb = float(match.group(3))
            rmsd_ub = float(match.group(4))

            results.append({
                "mode": mode,
                "affinity_kcal_mol": affinity,
                "rmsd_lb": rmsd_lb,
                "rmsd_ub": rmsd_ub
            })

    return results


# ============================================================
# Get current completed docking experiments
# ============================================================

def get_current_experiments():

    connection = sqlite3.connect(DB_FILE)

    cursor = connection.cursor()

    cursor.execute("""
        SELECT
            e.experiment_id,
            l.ligand_id,
            l.name
        FROM experiments e
        JOIN ligands l
            ON e.ligand_id = l.ligand_id
        WHERE e.status = 'completed'
        ORDER BY e.experiment_id
    """)

    experiments = cursor.fetchall()

    connection.close()

    return experiments


# ============================================================
# Export integrated screening results from SQLite VIEW
# ============================================================

def export_screening_results():

    connection = sqlite3.connect(DB_FILE)

    cursor = connection.cursor()

    cursor.execute("""
        SELECT
            experiment_id,
            ligand_id,
            ligand_name,
            target_id,
            target_name,
            pdb_id,
            target_chain,
            admet_engine,
            admet_decision,
            toxicity_flags,
            absorption_flags,
            physicochemical_flags,
            qed_value,
            vina_version,
            pose_number,
            docking_affinity_kcal_mol,
            interaction_event_count,
            unique_interaction_residue_count,
            interaction_types,
            contacted_residues
        FROM current_screening_results
        ORDER BY ligand_id
    """)

    rows = cursor.fetchall()

    columns = [
        description[0]
        for description in cursor.description
    ]

    connection.close()

    if not rows:
        print("WARNING: current_screening_results view returned no rows.")
        return 0

    with open(
        SCREENING_OUTPUT_FILE,
        "w",
        newline=""
    ) as file:

        writer = csv.writer(file)

        writer.writerow(columns)
        writer.writerows(rows)

    print()
    print("Integrated screening report:")
    print(f"  Rows exported: {len(rows)}")
    print(f"  CSV written:   {SCREENING_OUTPUT_FILE}")

    return len(rows)


# ============================================================
# Main
# ============================================================

def main():

    print("=" * 60)
    print("COLLECTING CURRENT DOCKING RESULTS")
    print("=" * 60)

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    experiments = get_current_experiments()

    if not experiments:
        print("ERROR: No completed docking experiments found.")
        return

    print()
    print(f"Completed experiments found: {len(experiments)}")

    all_results = []

    processed = 0

    for experiment_id, ligand_id, ligand_name in experiments:

        # Docking filenames are based on ligand ID + name.
        safe_name = ligand_name.replace(" ", "_")
        compound_name = f"{ligand_id}_{safe_name}"

        log_file = DOCKING_DIR / f"{compound_name}_docking.log"

        print()
        print(
            f"Experiment {experiment_id}: "
            f"{ligand_name} (ligand_id={ligand_id})"
        )

        if not log_file.exists():

            print(
                f"  WARNING: Log file not found: "
                f"{log_file.name}"
            )

            continue

        results = parse_vina_log(log_file)

        if not results:

            print("  WARNING: No docking results found.")

            continue

        for result in results:

            all_results.append({
                "experiment_id": experiment_id,
                "ligand_id": ligand_id,
                "compound": ligand_name,
                **result
            })

        processed += 1

        print(f"  Log: {log_file.name}")
        print(f"  Poses found: {len(results)}")
        print(
            f"  Best affinity: "
            f"{results[0]['affinity_kcal_mol']:.3f} kcal/mol"
        )

    if not all_results:

        print()
        print("ERROR: No docking results collected.")
        return

    # --------------------------------------------------------
    # Sort by affinity
    # --------------------------------------------------------

    all_results.sort(
        key=lambda x: x["affinity_kcal_mol"]
    )

    # --------------------------------------------------------
    # Write current-run docking CSV
    # --------------------------------------------------------

    fieldnames = [
        "experiment_id",
        "ligand_id",
        "compound",
        "mode",
        "affinity_kcal_mol",
        "rmsd_lb",
        "rmsd_ub"
    ]

    with open(
        OUTPUT_FILE,
        "w",
        newline=""
    ) as file:

        writer = csv.DictWriter(
            file,
            fieldnames=fieldnames
        )

        writer.writeheader()
        writer.writerows(all_results)

    # --------------------------------------------------------
    # Export integrated SQLite report
    # --------------------------------------------------------

    screening_rows = export_screening_results()

    # --------------------------------------------------------
    # Summary
    # --------------------------------------------------------

    print()
    print("=" * 60)
    print("RESULT COLLECTION SUMMARY")
    print("=" * 60)

    print(f"Experiments found:  {len(experiments)}")
    print(f"Logs processed:     {processed}")
    print(f"Total poses:        {len(all_results)}")
    print(f"Docking CSV:        {OUTPUT_FILE}")
    print(f"Screening rows:     {screening_rows}")
    print(f"Screening CSV:      {SCREENING_OUTPUT_FILE}")

    print()
    print("Best docking results:")

    best_by_compound = {}

    for result in all_results:

        compound = result["compound"]

        if compound not in best_by_compound:

            best_by_compound[compound] = result

    for result in sorted(
        best_by_compound.values(),
        key=lambda x: x["affinity_kcal_mol"]
    ):

        print(
            f"  {result['compound']}: "
            f"{result['affinity_kcal_mol']:.3f} kcal/mol"
        )


if __name__ == "__main__":
    main()
