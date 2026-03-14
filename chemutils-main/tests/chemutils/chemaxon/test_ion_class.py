import pytest

from chemutils.chemaxon import get_ion_class, get_major_microspecies_ion_class
from chemutils.molecule import oemol_from_smiles


@pytest.mark.parametrize(
    ["smiles", "expected_as_drawn", "expected_dominant_microspecies"],
    [
        # nitro group
        ("O=[N+]([O-])C1=CC=CC=C1", "uncharged", "uncharged"),
        # nitrone
        ("[O-][N+]1=CC=CC=C1", "uncharged", "uncharged"),
        # lysine
        ("N[C@@H](CCCCN)C(O)=O", "uncharged", "base"),
        # glycine
        ("NCC(O)=O", "uncharged", "zwitterion"),
        # glycine drawn charged
        ("[NH3+]CC([O-])=O", "zwitterion", "zwitterion"),
        # dinitro group
        ("O=[N+]([O-])C1=CC=CC([N+]([O-])=O)=C1", "uncharged", "uncharged"),
        # carboxylic acid drawn netural
        ("CCC(O)=O", "uncharged", "acid"),
        # carboxylic acid drawn charged
        ("CCC([O-])=O", "acid", "acid"),
        # glycine ester
        ("NCC(OC)=O", "uncharged", "base"),
    ],
)
def test_get_ion_class(smiles, expected_as_drawn, expected_dominant_microspecies):
    oemol = oemol_from_smiles(smiles)
    assert get_ion_class(oemol) == expected_as_drawn
    assert get_major_microspecies_ion_class(oemol) == expected_dominant_microspecies
