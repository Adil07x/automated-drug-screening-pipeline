from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt


# ============================================================
# Configuration
# ============================================================

PROJECT_ROOT = Path(__file__).resolve().parents[2]

INPUT_FILE = PROJECT_ROOT / "results/ranked/docking_results.csv"

OUTPUT_DIR = PROJECT_ROOT / "results/ranked"
FIGURE_DIR = PROJECT_ROOT / "results/figures"

SUMMARY_FILE = OUTPUT_DIR / "screening_summary.csv"
FIGURE_FILE = FIGURE_DIR / "best_docking_scores.png"


# ============================================================
# Main analysis
# ============================================================

def main():

    print("=" * 60)
    print("DOCKING RESULT ANALYSIS")
    print("=" * 60)

    # --------------------------------------------------------
    # Check input
    # --------------------------------------------------------

    if not INPUT_FILE.exists():
        print("ERROR: Docking results CSV not found:")
        print(INPUT_FILE)
        return

    # --------------------------------------------------------
    # Load results
    # --------------------------------------------------------

    df = pd.read_csv(INPUT_FILE)

    print(f"Input rows: {len(df)}")
    print(f"Compounds:  {df['compound'].nunique()}")

    # --------------------------------------------------------
    # Select best pose for each compound
    # --------------------------------------------------------

    best_results = (
        df.loc[
            df.groupby("compound")["affinity_kcal_mol"].idxmin()
        ]
        .copy()
    )

    # Rank compounds by docking affinity
    best_results = best_results.sort_values(
        "affinity_kcal_mol"
    ).reset_index(drop=True)

    best_results.insert(
        0,
        "rank",
        range(1, len(best_results) + 1)
    )

    # --------------------------------------------------------
    # Keep summary columns
    # --------------------------------------------------------

    summary = best_results[
        [
            "rank",
            "compound",
            "mode",
            "affinity_kcal_mol",
            "rmsd_lb",
            "rmsd_ub"
        ]
    ]

    # --------------------------------------------------------
    # Save summary
    # --------------------------------------------------------

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)

    summary.to_csv(
        SUMMARY_FILE,
        index=False
    )

    # --------------------------------------------------------
    # Print summary
    # --------------------------------------------------------

    print()
    print("Screening summary:")
    print(summary.to_string(index=False))

    print()
    print(f"Summary written to:")
    print(SUMMARY_FILE)

    # --------------------------------------------------------
    # Create figure
    # --------------------------------------------------------

    plt.figure(figsize=(8, 5))

    plt.bar(
        summary["compound"],
        summary["affinity_kcal_mol"]
    )

    plt.xlabel("Compound")
    plt.ylabel("Best docking affinity (kcal/mol)")
    plt.title("Best Docking Scores")

    plt.xticks(
        rotation=45,
        ha="right"
    )

    plt.tight_layout()

    plt.savefig(
        FIGURE_FILE,
        dpi=300
    )

    plt.close()

    print()
    print(f"Figure written to:")
    print(FIGURE_FILE)

    print()
    print("=" * 60)
    print("ANALYSIS COMPLETE")
    print("=" * 60)


if __name__ == "__main__":
    main()
