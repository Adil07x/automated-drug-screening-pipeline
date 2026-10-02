from pathlib import Path
import sqlite3
import sys


PROJECT_ROOT = Path(__file__).resolve().parents[2]
DB_FILE = PROJECT_ROOT / "data/database/screening_database.sqlite"

ENGINE = "ADMET-AI"


# ============================================================
# VALIDATED SCREENING RULES
# ============================================================

TOXICITY_THRESHOLDS = {
    "AMES": 0.50,
    "DILI": 0.50,
    "ClinTox": 0.50,
    "hERG": 0.50,
}

ABSORPTION_THRESHOLDS = {
    "Bioavailability_Ma": 0.50,
    "Caco2_Wang": -5.50,
}

PHYSICOCHEMICAL_THRESHOLDS = {
    "Solubility_AqSolDB": -4.00,
    "QED": 0.30,
}


# ============================================================
# HELPERS
# ============================================================

def flag_high(value, threshold):
    """Flag when a value is >= threshold."""
    return int(value >= threshold)


def flag_low(value, threshold):
    """Flag when a value is < threshold."""
    return int(value < threshold)


def get_decision(toxicity_flags):
    if toxicity_flags >= 2:
        return "FAIL"

    if toxicity_flags == 1:
        return "REVIEW"

    return "PASS"


def get_docking_eligibility(decision):
    if decision == "PASS":
        return "ELIGIBLE"

    if decision == "REVIEW":
        return "REVIEW"

    return "EXCLUDED"


# ============================================================
# MAIN
# ============================================================

def main():

    if not DB_FILE.exists():
        print(f"ERROR: Database not found: {DB_FILE}")
        sys.exit(1)

    connection = sqlite3.connect(DB_FILE)
    connection.row_factory = sqlite3.Row

    try:
        cursor = connection.cursor()

        # ----------------------------------------------------
        # Determine the ADMET-AI version currently in database
        # ----------------------------------------------------

        cursor.execute(
            """
            SELECT DISTINCT engine_version
            FROM admet_results
            WHERE engine = ?
              AND engine_version IS NOT NULL
            ORDER BY engine_version
            """,
            (ENGINE,),
        )

        versions = [row["engine_version"] for row in cursor.fetchall()]

        if not versions:
            print("ERROR: No ADMET-AI results found.")
            sys.exit(1)

        if len(versions) > 1:
            print(
                "ERROR: Multiple ADMET-AI versions found in database:"
            )
            for version in versions:
                print(f"  - {version}")

            print(
                "\nPlease resolve the version situation before screening."
            )
            sys.exit(1)

        engine_version = versions[0]

        print("=" * 70)
        print("ADMET SCREENING")
        print("=" * 70)
        print(f"Engine:         {ENGINE}")
        print(f"Engine version: {engine_version}")
        print()

        # ----------------------------------------------------
        # Required endpoints
        # ----------------------------------------------------

        required_endpoints = set(
            list(TOXICITY_THRESHOLDS.keys())
            + list(ABSORPTION_THRESHOLDS.keys())
            + list(PHYSICOCHEMICAL_THRESHOLDS.keys())
        )

        placeholders = ",".join("?" for _ in required_endpoints)

        # ----------------------------------------------------
        # Find compounds with ADMET results but no screening
        # ----------------------------------------------------

        cursor.execute(
            """
            SELECT
                l.ligand_id,
                l.name
            FROM ligands l
            JOIN admet_results ar
                ON ar.ligand_id = l.ligand_id
            LEFT JOIN admet_screening s
                ON s.ligand_id = l.ligand_id
               AND s.engine = ?
            WHERE ar.engine = ?
              AND s.screening_id IS NULL
            GROUP BY l.ligand_id, l.name
            ORDER BY l.ligand_id
            """,
            (ENGINE, ENGINE),
        )

        candidates = cursor.fetchall()

        if not candidates:
            print("No new compounds require ADMET screening.")
            print()
            print("Existing screening records were left unchanged.")
            return

        print(f"Compounds requiring screening: {len(candidates)}")
        print()

        screened = 0
        skipped = 0

        for compound in candidates:

            ligand_id = compound["ligand_id"]
            name = compound["name"]

            # ------------------------------------------------
            # Fetch the required endpoint values
            # ------------------------------------------------

            cursor.execute(
                f"""
                SELECT endpoint, value
                FROM admet_results
                WHERE ligand_id = ?
                  AND engine = ?
                  AND endpoint IN ({placeholders})
                """,
                (
                    ligand_id,
                    ENGINE,
                    *required_endpoints,
                ),
            )

            endpoint_values = {
                row["endpoint"]: row["value"]
                for row in cursor.fetchall()
            }

            missing = sorted(
                required_endpoints - set(endpoint_values.keys())
            )

            null_values = sorted(
                endpoint
                for endpoint in required_endpoints
                if endpoint in endpoint_values
                and endpoint_values[endpoint] is None
            )

            if missing or null_values:

                print(f"SKIPPED: {name} (ligand {ligand_id})")

                if missing:
                    print(
                        "  Missing endpoints: "
                        + ", ".join(missing)
                    )

                if null_values:
                    print(
                        "  NULL endpoints: "
                        + ", ".join(null_values)
                    )

                skipped += 1
                continue

            # ------------------------------------------------
            # Calculate toxicity flags
            # ------------------------------------------------

            toxicity_flags = sum(
                flag_high(
                    endpoint_values[endpoint],
                    threshold,
                )
                for endpoint, threshold
                in TOXICITY_THRESHOLDS.items()
            )

            # ------------------------------------------------
            # Calculate absorption flags
            # ------------------------------------------------

            absorption_flags = (
                flag_low(
                    endpoint_values["Bioavailability_Ma"],
                    ABSORPTION_THRESHOLDS["Bioavailability_Ma"],
                )
                +
                flag_low(
                    endpoint_values["Caco2_Wang"],
                    ABSORPTION_THRESHOLDS["Caco2_Wang"],
                )
            )

            # ------------------------------------------------
            # Calculate physicochemical flags
            # ------------------------------------------------

            physicochemical_flags = (
                flag_low(
                    endpoint_values["Solubility_AqSolDB"],
                    PHYSICOCHEMICAL_THRESHOLDS[
                        "Solubility_AqSolDB"
                    ],
                )
                +
                flag_low(
                    endpoint_values["QED"],
                    PHYSICOCHEMICAL_THRESHOLDS["QED"],
                )
            )

            # ------------------------------------------------
            # Final decision
            # ------------------------------------------------

            decision = get_decision(toxicity_flags)

            docking_eligibility = get_docking_eligibility(
                decision
            )

            rationale = (
                f"Toxicity flags: {toxicity_flags}/4; "
                f"Absorption flags: {absorption_flags}/2; "
                f"Physicochemical flags: "
                f"{physicochemical_flags}/2; "
                f"Decision based on validated toxicity rule: "
                f"{decision}"
            )

            # ------------------------------------------------
            # Insert screening result
            # ------------------------------------------------

            cursor.execute(
                """
                INSERT INTO admet_screening (
                    ligand_id,
                    engine,
                    decision,
                    toxicity_flags,
                    absorption_flags,
                    physicochemical_flags,
                    qed_value,
                    rationale
                )
                VALUES (?, ?, ?, ?, ?, ?, ?, ?)
                """,
                (
                    ligand_id,
                    ENGINE,
                    decision,
                    toxicity_flags,
                    absorption_flags,
                    physicochemical_flags,
                    endpoint_values["QED"],
                    rationale,
                ),
            )

            screening_id = cursor.lastrowid

            # ------------------------------------------------
            # Insert docking eligibility
            # ------------------------------------------------

            cursor.execute(
                """
                INSERT INTO docking_eligibility (
                    ligand_id,
                    screening_id,
                    eligibility,
                    reason
                )
                VALUES (?, ?, ?, ?)
                """,
                (
                    ligand_id,
                    screening_id,
                    docking_eligibility,
                    (
                        f"ADMET decision = {decision}; "
                        f"docking eligibility = "
                        f"{docking_eligibility}"
                    ),
                ),
            )

            connection.commit()

            print(f"SCREENED: {name}")
            print(f"  Ligand ID:             {ligand_id}")
            print(f"  Toxicity flags:        {toxicity_flags}")
            print(f"  Absorption flags:      {absorption_flags}")
            print(
                f"  Physicochemical flags: "
                f"{physicochemical_flags}"
            )
            print(f"  Decision:              {decision}")
            print(
                f"  Docking eligibility:   "
                f"{docking_eligibility}"
            )
            print()

            screened += 1

        # ----------------------------------------------------
        # Summary
        # ----------------------------------------------------

        print("=" * 70)
        print("SCREENING SUMMARY")
        print("=" * 70)
        print(f"New compounds screened: {screened}")
        print(f"Skipped/incomplete:     {skipped}")
        print("=" * 70)

    finally:
        connection.close()


if __name__ == "__main__":
    main()
