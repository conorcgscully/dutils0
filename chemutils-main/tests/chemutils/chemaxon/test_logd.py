import pytest

from chemutils.chemaxon import get_logd


@pytest.mark.parametrize(
    ["smiles", "expected"],
    [
        # Acetamide
        ("O=C(N)C", -1.03),
        # Methanol
        ("CO", -0.52),
        # Formic Acid
        ("O=CO", -3.26),
        # Diethyl Ether
        ("CCOCC", 0.84),
        # Dichlorobenzene
        ("Clc1ccc(Cl)cc1", 3.18),
        # Paracetemol
        ("CC(=O)Nc1ccc(O)cc1", 0.91),
    ],
)
def test_get_logd(smiles, expected):
    assert get_logd(smiles, pH=7.4) == pytest.approx(expected, abs=0.01)
