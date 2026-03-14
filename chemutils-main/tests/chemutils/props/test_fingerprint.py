import pytest
from numpy.testing import assert_equal
from rdkit import Chem

from chemutils.molecule.smiles import oemol_from_smiles
from chemutils.props.fingerprint import get_morgan_fingerprint

from .drug_smiles import DRUG_SMILES


@pytest.mark.parametrize(
    "smiles",
    [*DRUG_SMILES.values()],
)
def test_agrees_rdkit(smiles):
    oemol = oemol_from_smiles(smiles)
    rdmol = Chem.MolFromSmiles(smiles)

    assert_equal(get_morgan_fingerprint(oemol), get_morgan_fingerprint(rdmol))


@pytest.mark.parametrize(
    ["smilesA", "smilesB", "chirality", "match"],
    [
        ("C[C@H](Cl)F", "C[C@@H](Cl)F", False, True),
        ("C[C@H](Cl)F", "C[C@@H](Cl)F", True, False),
        ("CC(Cl)F", "C[C@@H](Cl)F", False, True),
        ("CC(Cl)F", "C[C@H](Cl)F", True, False),
    ],
)
def test_morgan_fingerprint_chirality(smilesA, smilesB, chirality, match):
    molA = oemol_from_smiles(smilesA)
    molB = oemol_from_smiles(smilesB)

    fpA = get_morgan_fingerprint(molA, useChirality=chirality)
    fpB = get_morgan_fingerprint(molB, useChirality=chirality)

    assert_equal(all(fpA == fpB), match)
