import pytest

from chemutils.chemaxon import get_logp


@pytest.mark.parametrize(
    ["smiles", "expected"],
    [
        # Acetamide
        ("O=C(N)C", -1.03),
        # Methanol
        ("CO", -0.52),
        # Formic Acid
        ("O=CO", -0.27),
        # Diethyl Ether
        ("CCOCC", 0.84),
        # Dichlorobenzene
        ("Clc1ccc(Cl)cc1", 3.18),
        # Paracetemol
        ("CC(=O)Nc1ccc(O)cc1", 0.91),
    ],
)
def test_get_logp(smiles, expected):
    assert get_logp(smiles) == pytest.approx(expected, abs=0.01)
