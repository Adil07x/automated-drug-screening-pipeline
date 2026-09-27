from pathlib import Path
import subprocess
import sys
import json


# ============================================================
# Configuration
# ============================================================

PROJECT_ROOT = Path(__file__).resolve().parents[2]

CONFIG_FILE = PROJECT_ROOT / "config/docking_config.json"

LIGAND_DIR = PROJECT_ROOT / "data/ligands/pdbqt"
OUTPUT_DIR = PROJECT_ROOT / "results/docking"

VINA = "vina"

# Reference ligand used for validation only.
REFERENCE_LIGAND = "STI_reference.pdbqt"


# ============================================================
# Load docking configuration
# ============================================================

if not CONFIG_FILE.exists():
    print("ERROR: Configuration file not found:")
    print(CONFIG_FILE)
    sys.exit(1)

try:
    with CONFIG_FILE.open("r", encoding="utf-8") as file:
        CONFIG = json.load(file)

except json.JSONDecodeError as error:
    print("ERROR: Invalid JSON configuration file.")
    print(error)
    sys.exit(1)


# ============================================================
# Target and receptor
# ============================================================

TARGET = CONFIG["target"]

RECEPTOR = PROJECT_ROOT / CONFIG["receptor"]


# ============================================================
# Docking box
# ============================================================

BOX = CONFIG["docking_box"]

CENTER_X = BOX["center_x"]
CENTER_Y = BOX["center_y"]
CENTER_Z = BOX["center_z"]

SIZE_X = BOX["size_x"]
SIZE_Y = BOX["size_y"]
SIZE_Z = BOX["size_z"]


# ============================================================
# Vina parameters
# ============================================================

DOCKING_PARAMETERS = CONFIG["docking_parameters"]

EXHAUSTIVENESS = DOCKING_PARAMETERS["exhaustiveness"]
NUM_MODES = DOCKING_PARAMETERS["num_modes"]


# ============================================================
# Dock one ligand
# ============================================================

def dock_ligand(ligand_file):

    ligand_name = ligand_file.stem

    output_file = OUTPUT_DIR / f"{ligand_name}_docked.pdbqt"
    log_file = OUTPUT_DIR / f"{ligand_name}_docking.log"

    command = [
        VINA,

        "--receptor", str(RECEPTOR),
        "--ligand", str(ligand_file),

        "--center_x", str(CENTER_X),
        "--center_y", str(CENTER_Y),
        "--center_z", str(CENTER_Z),

        "--size_x", str(SIZE_X),
        "--size_y", str(SIZE_Y),
        "--size_z", str(SIZE_Z),

        "--exhaustiveness", str(EXHAUSTIVENESS),
        "--num_modes", str(NUM_MODES),

        "--out", str(output_file),
    ]

    print()
    print("=" * 60)
    print(f"Ligand: {ligand_name}")
    print("=" * 60)

    print("Running:")
    print(" ".join(command))

    try:

        result = subprocess.run(
            command,
            capture_output=True,
            text=True,
            check=False
        )

        # Save both stdout and stderr to the log
        log_file.write_text(
            result.stdout +
            "\n\n===== STDERR =====\n\n" +
            result.stderr
        )

        if result.returncode != 0:

            print(f"FAILED: {ligand_name}")
            print(f"See log: {log_file}")

            return False

        print(f"SUCCESS: {ligand_name}")
        print(f"Output: {output_file}")
        print(f"Log:    {log_file}")

        return True

    except Exception as error:

        print(f"ERROR: {ligand_name}")
        print(error)

        return False


# ============================================================
# Main pipeline
# ============================================================

def main():

    print("=" * 60)
    print("AUTOMATED DOCKING")
    print("=" * 60)

    print(f"Target: {TARGET}")

    # --------------------------------------------------------
    # Show docking configuration
    # --------------------------------------------------------

    print()
    print("Docking configuration:")
    print(
        f"  Center: "
        f"{CENTER_X}, {CENTER_Y}, {CENTER_Z}"
    )

    print(
        f"  Size:   "
        f"{SIZE_X}, {SIZE_Y}, {SIZE_Z}"
    )

    print(
        f"  Exhaustiveness: "
        f"{EXHAUSTIVENESS}"
    )

    print(
        f"  Number of modes: "
        f"{NUM_MODES}"
    )

    # --------------------------------------------------------
    # Check receptor
    # --------------------------------------------------------

    if not RECEPTOR.exists():

        print()
        print("ERROR: Receptor not found:")
        print(RECEPTOR)

        sys.exit(1)

    print()
    print(f"Receptor: {RECEPTOR}")

    # --------------------------------------------------------
    # Create output directory
    # --------------------------------------------------------

    OUTPUT_DIR.mkdir(
        parents=True,
        exist_ok=True
    )

    # --------------------------------------------------------
    # Find ligand PDBQT files
    # --------------------------------------------------------

    ligand_files = sorted(
        LIGAND_DIR.glob("*.pdbqt")
    )

    # Exclude the reference ligand
    ligand_files = [
        ligand for ligand in ligand_files
        if ligand.name != REFERENCE_LIGAND
    ]

    print()
    print(
        f"Screening ligands found: "
        f"{len(ligand_files)}"
    )

    if not ligand_files:

        print("No screening ligands found.")
        sys.exit(1)

    # --------------------------------------------------------
    # Dock all ligands
    # --------------------------------------------------------

    successful = 0
    failed = 0

    for ligand_file in ligand_files:

        if dock_ligand(ligand_file):
            successful += 1
        else:
            failed += 1

    # --------------------------------------------------------
    # Summary
    # --------------------------------------------------------

    print()
    print("=" * 60)
    print("DOCKING SUMMARY")
    print("=" * 60)

    print(f"Total ligands: {len(ligand_files)}")
    print(f"Successful:    {successful}")
    print(f"Failed:        {failed}")

    print("=" * 60)

    if failed > 0:
        sys.exit(1)


# ============================================================
# Script entry point
# ============================================================

if __name__ == "__main__":
    main()