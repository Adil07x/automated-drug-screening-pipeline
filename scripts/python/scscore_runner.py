import math
import sys

from pathlib import Path

import numpy as np
from rdkit import Chem


# Project root
PROJECT_ROOT = Path(__file__).resolve().parents[2]

# SCScore source
SCSCORE_ROOT = PROJECT_ROOT / "external" / "scscore"

sys.path.insert(0, str(SCSCORE_ROOT))

from scscore.standalone_model_numpy import SCScorer


class CompatibleSCScorer(SCScorer):
    """
    Compatibility wrapper for the original SCScore NumPy implementation.

    The original implementation can return a one-element NumPy array from
    the final layer. Python's math.exp() requires a scalar, so we convert
    that value explicitly.
    """

    def apply(self, x):
        if not self._restored:
            raise ValueError("Must restore model weights!")

        for i in range(0, len(self.vars), 2):
            last_layer = (i == len(self.vars) - 2)

            W = self.vars[i]
            b = self.vars[i + 1]

            x = np.matmul(x, W) + b

            if not last_layer:
                x = x * (x > 0)

        x = np.asarray(x).item()

        return 1 + (self.score_scale - 1) * (
            1 / (1 + math.exp(-x))
        )


def calculate_scscore(smiles):
    """
    Calculate SCScore for a SMILES string.
    """

    model_path = (
        SCSCORE_ROOT
        / "models"
        / "full_reaxys_model_1024bool"
        / "model.ckpt-10654.as_numpy.json.gz"
    )

    model = CompatibleSCScorer()
    model.restore(str(model_path))

    canonical_smiles, score = model.get_score_from_smi(smiles)

    return canonical_smiles, float(score)


if __name__ == "__main__":

    ligand_file = PROJECT_ROOT / "data" / "raw" / "compounds" / "Quercetin.sdf"

    mol = Chem.SDMolSupplier(
        str(ligand_file),
        removeHs=False
    )[0]

    if mol is None:
        raise ValueError("Could not read Quercetin.sdf")

    smiles = Chem.MolToSmiles(
        mol,
        isomericSmiles=True
    )

    canonical_smiles, score = calculate_scscore(smiles)

    print("SCScore compatibility test: OK")
    print("Quercetin SMILES:", canonical_smiles)
    print("Quercetin SCScore:", round(score, 3))
