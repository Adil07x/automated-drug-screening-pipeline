#!/usr/bin/env python3

import time
from pathlib import Path

import pandas as pd
import requests


# --------------------------------------------------
# Configuration
# --------------------------------------------------

PROJECT_ROOT = Path(__file__).resolve().parents[2]

INPUT_FILE = PROJECT_ROOT / "data/metadata/compound_library.csv"
OUTPUT_DIR = PROJECT_ROOT / "data/raw/compounds"

PUBCHEM_URL = (
    "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/"
    "cid/{cid}/SDF"
)

MAX_RETRIES = 5
INITIAL_WAIT = 5


# --------------------------------------------------
# Main
# --------------------------------------------------

def main():

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(INPUT_FILE)

    required_columns = {
        "compound_id",
        "name",
        "pubchem_cid",
        "role",
    }

    missing = required_columns - set(df.columns)

    if missing:
        raise ValueError(
            f"Missing required columns: {sorted(missing)}"
        )

    print(f"Compounds to download: {len(df)}")
    print(f"Output directory: {OUTPUT_DIR}")
    print()

    session = requests.Session()

    session.headers.update({
        "User-Agent": (
            "Automated-Drug-Screening-Pipeline/1.0 "
            "(research workflow)"
        )
    })

    successful = 0
    failed = 0

    for _, row in df.iterrows():

        compound_id = row["compound_id"]
        name = row["name"]
        cid = int(row["pubchem_cid"])

        output_file = OUTPUT_DIR / f"{compound_id}_{cid}.sdf"

        print(
            f"[{compound_id}] {name} "
            f"(PubChem CID {cid})"
        )

        if output_file.exists():
            print("  Already exists - skipping")
            successful += 1
            continue

        url = PUBCHEM_URL.format(cid=cid)

        downloaded = False

        for attempt in range(1, MAX_RETRIES + 1):

            try:

                response = session.get(
                    url,
                    timeout=30
                )

                if response.status_code == 200:

                    output_file.write_text(
                        response.text,
                        encoding="utf-8"
                    )

                    print(
                        f"  Downloaded -> {output_file.name}"
                    )

                    successful += 1
                    downloaded = True
                    break

                elif response.status_code == 503:

                    wait_time = INITIAL_WAIT * (2 ** (attempt - 1))

                    print(
                        f"  PubChem server busy "
                        f"(attempt {attempt}/{MAX_RETRIES})"
                    )

                    if attempt < MAX_RETRIES:
                        print(
                            f"  Waiting {wait_time} seconds..."
                        )
                        time.sleep(wait_time)

                else:

                    print(
                        f"  HTTP error: "
                        f"{response.status_code}"
                    )
                    break

            except requests.RequestException as error:

                wait_time = INITIAL_WAIT * (2 ** (attempt - 1))

                print(
                    f"  Connection error "
                    f"(attempt {attempt}/{MAX_RETRIES}): "
                    f"{error}"
                )

                if attempt < MAX_RETRIES:
                    print(
                        f"  Waiting {wait_time} seconds..."
                    )
                    time.sleep(wait_time)

        if not downloaded:
            print("  FAILED")

            failed += 1

        # Keep request rate conservative
        time.sleep(2)

    print()
    print("Download summary")
    print("----------------")
    print(f"Successful: {successful}")
    print(f"Failed:    {failed}")


if __name__ == "__main__":
    main()