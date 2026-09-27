from admet.admet_ai_engine import ADMETAIEngine


SMILES = "O=C(O)c1ccccc1"


def main():
    engine = ADMETAIEngine()

    predictions = engine.predict(SMILES)

    print()
    print("=" * 60)
    print("ADMET-AI ENGINE TEST")
    print("=" * 60)

    print(f"Number of endpoints: {len(predictions)}")

    print()
    print("First 10 endpoints:")

    for endpoint, value in list(predictions.items())[:10]:
        print(f"{endpoint}: {value}")

    print()
    print("=" * 60)
    print("TEST PASSED")
    print("=" * 60)


if __name__ == "__main__":
    main()
