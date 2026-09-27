import json
from importlib.metadata import PackageNotFoundError, version

from admet_ai import ADMETModel


class ADMETAIEngine:
    """
    Wrapper around ADMET-AI.

    The engine:
    1. Loads the ADMET-AI model once.
    2. Accepts a SMILES string.
    3. Returns the complete ADMET-AI prediction dictionary.
    """

    def __init__(self):
        print("Loading ADMET-AI model...")
        self.model = ADMETModel()
        self.engine_version = self._get_version()
        print(f"ADMET-AI version: {self.engine_version}")
        print("ADMET-AI model loaded.")

    @staticmethod
    def _get_version():
        """Return the installed ADMET-AI package version."""
        try:
            return version("admet-ai")
        except PackageNotFoundError:
            return "unknown"

    @staticmethod
    def _is_numeric(value):
        """Check whether a prediction value can be stored as a number."""
        return isinstance(value, (int, float)) and not isinstance(value, bool)

    @staticmethod
    def classify_prediction_type(endpoint, value):
        """
        Classify an ADMET-AI output for database/export purposes.

        This is a storage classification, not a biological interpretation.
        """

        # Basic molecular descriptors returned by ADMET-AI
        descriptor_endpoints = {
            "molecular_weight",
            "logP",
            "hydrogen_bond_acceptors",
            "hydrogen_bond_donors",
            "Lipinski",
            "QED",
            "stereo_centers",
            "tpsa",
            "PAINS_alert",
            "BRENK_alert",
            "NIH_alert",
        }

        if endpoint in descriptor_endpoints:
            return "descriptor"

        # ADMET-AI DrugBank percentile outputs
        if endpoint.endswith("_drugbank_approved_percentile"):
            return "percentile"

        # Known probability/classification-style endpoints
        probability_endpoints = {
            "AMES",
            "BBB_Martins",
            "Bioavailability_Ma",
            "CYP1A2_Veith",
            "CYP2C19_Veith",
            "CYP2C9_Substrate_CarbonMangels",
            "CYP2C9_Veith",
            "CYP2D6_Substrate_CarbonMangels",
            "CYP2D6_Veith",
            "CYP3A4_Substrate_CarbonMangels",
            "CYP3A4_Veith",
            "Carcinogens_Lagunin",
            "ClinTox",
            "DILI",
            "HIA_Hou",
            "NR-AR-LBD",
            "NR-AR",
            "NR-AhR",
            "NR-Aromatase",
            "NR-ER-LBD",
            "NR-ER",
            "NR-PPAR-gamma",
            "PAMPA_NCATS",
            "Pgp_Broccatelli",
            "SR-ARE",
            "SR-ATAD5",
            "SR-HSE",
            "SR-MMP",
            "SR-p53",
            "Skin_Reaction",
            "hERG",
        }

        if endpoint in probability_endpoints:
            return "probability"

        # Everything else is treated as a continuous model output.
        return "continuous"

    def predict(self, smiles):
        """
        Run ADMET-AI prediction for one molecule.

        Parameters
        ----------
        smiles : str
            Molecular SMILES.

        Returns
        -------
        dict
            Complete ADMET-AI prediction dictionary.
        """

        if not smiles or not isinstance(smiles, str):
            raise ValueError("A valid SMILES string is required.")

        predictions = self.model.predict(smiles=smiles)

        if not isinstance(predictions, dict):
            raise TypeError(
                f"Unexpected ADMET-AI output type: {type(predictions)}"
            )

        return predictions

    def to_database_records(self, predictions):
        """
        Convert the ADMET-AI prediction dictionary into rows suitable
        for insertion into the admet_results table.
        """

        raw_result = json.dumps(
            predictions,
            ensure_ascii=False,
            sort_keys=True,
        )

        records = []

        for endpoint, value in predictions.items():

            numeric_value = None

            if self._is_numeric(value):
                numeric_value = float(value)

            prediction_type = self.classify_prediction_type(
                endpoint,
                value,
            )

            records.append(
                {
                    "endpoint": endpoint,
                    "value": numeric_value,
                    "unit": None,
                    "prediction_type": prediction_type,
                    "uncertainty": None,
                    "status": "success",
                    "raw_result": raw_result,
                }
            )

        return records
