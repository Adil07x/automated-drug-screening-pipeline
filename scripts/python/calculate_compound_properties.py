import json
import math
import sqlite3
import sys
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import Crippen
from rdkit.Chem import Descriptors
from rdkit.Chem import FilterCatalog
from rdkit.Chem import Lipinski
from rdkit.Chem import rdMolDescriptors


# ============================================================
# Project paths
# ============================================================

PROJECT_ROOT = Path(__file__).resolve().parents[2]

DB_PATH = PROJECT_ROOT / "data" / "database" / "screening_database.sqlite"

COMPOUND_DIR = PROJECT_ROOT / "data" / "raw" / "compounds"

SCSCORE_ROOT = PROJECT_ROOT / "external" / "scscore"

SA_SCORE_PATH = (
    Path(sys.prefix)
    / "share"
    / "RDKit"
    / "Contrib"
    / "SA_Score"
)

sys.path.insert(0, str(SA_SCORE_PATH))

import sascorer


# SCScore compatibility wrapper
sys.path.insert(0, str(SCSCORE_ROOT))

from scscore.standalone_model_numpy import SCScorer


class CompatibleSCScorer(SCScorer):
    """
    Compatibility wrapper for the original SCScore NumPy model.

    The original implementation may return a one-element NumPy array
    from the final layer. We convert that value to a scalar before
    passing it to math.exp().
    """

    def apply(self, x):

        if not self._restored:
            raise ValueError("Must restore model weights!")

        for i in range(0, len(self.vars), 2):

            last_layer = (i == len(self.vars) - 2)

            W = self.vars[i]
            b = self.vars[i + 1]

            x = self._matmul(x, W) + b

            if not last_layer:
                x = x * (x > 0)

        x = float(x.item())

        return 1 + (self.score_scale - 1) * (
            1 / (1 + math.exp(-x))
        )

    @staticmethod
    def _matmul(x, W):
        import numpy as np
        return np.matmul(x, W)


# ============================================================
# Initialize PAINS and Brenk catalogs
# ============================================================

pains_params = FilterCatalog.FilterCatalogParams()

pains_params.AddCatalog(
    FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS
)

pains_catalog = FilterCatalog.FilterCatalog(pains_params)


brenk_params = FilterCatalog.FilterCatalogParams()

brenk_params.AddCatalog(
    FilterCatalog.FilterCatalogParams.FilterCatalogs.BRENK
)

brenk_catalog = FilterCatalog.FilterCatalog(brenk_params)


# ============================================================
# SCScore model
# ============================================================

SC_MODEL_PATH = (
    SCSCORE_ROOT
    / "models"
    / "full_reaxys_model_1024bool"
    / "model.ckpt-10654.as_numpy.json.gz"
)

scscore_model = CompatibleSCScorer()

scscore_model.restore(str(SC_MODEL_PATH))


# ============================================================
# Helper functions
# ============================================================

def get_alerts(catalog, mol):

    matches = catalog.GetMatches(mol)

    alerts = []

    for match in matches:
        description = match.GetDescription()

        if description:
            alerts.append(description)
        else:
            alerts.append(str(match))

    return alerts


def calculate_lipinski_violations(mol):

    violations = 0

    if Descriptors.MolWt(mol) > 500:
        violations += 1

    if Crippen.MolLogP(mol) > 5:
        violations += 1

    if Lipinski.NumHDonors(mol) > 5:
        violations += 1

    if Lipinski.NumHAcceptors(mol) > 10:
        violations += 1

    return violations


def calculate_properties(mol):

    molecular_weight = Descriptors.MolWt(mol)

    logp = Crippen.MolLogP(mol)

    hbd = Lipinski.NumHDonors(mol)

    hba = Lipinski.NumHAcceptors(mol)

    lipinski_violations = calculate_lipinski_violations(mol)

    tpsa = rdMolDescriptors.CalcTPSA(mol)

    rotatable_bonds = rdMolDescriptors.CalcNumRotatableBonds(mol)

    veber_pass = (
        tpsa <= 140
        and rotatable_bonds <= 10
    )

    pains_alerts = get_alerts(
        pains_catalog,
        mol
    )

    brenk_alerts = get_alerts(
        brenk_catalog,
        mol
    )

    sa_score = sascorer.calculateScore(mol)

    smiles = Chem.MolToSmiles(
        mol,
        isomericSmiles=True
    )

    _, sc_score = scscore_model.get_score_from_smi(
        smiles
    )

    return {
        "molecular_weight": molecular_weight,
        "logp": logp,
        "hbd": hbd,
        "hba": hba,
        "lipinski_violations": lipinski_violations,
        "tpsa": tpsa,
        "rotatable_bonds": rotatable_bonds,
        "veber_pass": int(veber_pass),
        "pains_alert_count": len(pains_alerts),
        "pains_alerts": json.dumps(pains_alerts),
        "brenk_alert_count": len(brenk_alerts),
        "brenk_alerts": json.dumps(brenk_alerts),
        "sa_score": sa_score,
        "sc_score": float(sc_score),
    }


# ============================================================
# Database update
# ============================================================

def main():

    conn = sqlite3.connect(DB_PATH)

    cursor = conn.cursor()

    compounds = cursor.execute(
        """
        SELECT ligand_id, name, source_file
        FROM ligands
        ORDER BY ligand_id
        """
    ).fetchall()

    print(f"Compounds found in database: {len(compounds)}")
    print()

    for ligand_id, name, source_file in compounds:

        sdf_path = PROJECT_ROOT / source_file

        if not sdf_path.exists():
            print(f"SKIP: {name}")
            print(f"Missing file: {sdf_path}")
            print()
            continue

        supplier = Chem.SDMolSupplier(
            str(sdf_path),
            removeHs=False
        )

        mol = supplier[0]

        if mol is None:
            print(f"SKIP: {name}")
            print("Could not read molecule.")
            print()
            continue

        properties = calculate_properties(mol)

        cursor.execute(
            """
            UPDATE compound_properties

            SET
                molecular_weight = ?,
                logp = ?,
                hbd = ?,
                hba = ?,
                lipinski_violations = ?,
                tpsa = ?,
                rotatable_bonds = ?,
                veber_pass = ?,
                pains_alert_count = ?,
                pains_alerts = ?,
                brenk_alert_count = ?,
                brenk_alerts = ?,
                sa_score = ?,
                sc_score = ?,
                calculated_at = CURRENT_TIMESTAMP

            WHERE ligand_id = ?
            """,
            (
                properties["molecular_weight"],
                properties["logp"],
                properties["hbd"],
                properties["hba"],
                properties["lipinski_violations"],
                properties["tpsa"],
                properties["rotatable_bonds"],
                properties["veber_pass"],
                properties["pains_alert_count"],
                properties["pains_alerts"],
                properties["brenk_alert_count"],
                properties["brenk_alerts"],
                properties["sa_score"],
                properties["sc_score"],
                ligand_id,
            )
        )

        print(
            f"{name}: "
            f"MW={properties['molecular_weight']:.2f}, "
            f"LogP={properties['logp']:.2f}, "
            f"Lipinski={properties['lipinski_violations']}, "
            f"TPSA={properties['tpsa']:.2f}, "
            f"RotB={properties['rotatable_bonds']}, "
            f"Veber={'PASS' if properties['veber_pass'] else 'FAIL'}, "
            f"PAINS={properties['pains_alert_count']}, "
            f"Brenk={properties['brenk_alert_count']}, "
            f"SA={properties['sa_score']:.3f}, "
            f"SC={properties['sc_score']:.3f}"
        )

    conn.commit()
    conn.close()

    print()
    print("Compound property calculation completed.")


if __name__ == "__main__":
    main()
