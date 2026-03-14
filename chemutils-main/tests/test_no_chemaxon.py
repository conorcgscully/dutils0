import pytest

from chemutils.chemaxon import ChemaxonNotInstalledError, as_cxmol
from chemutils.molecule import oemol_from_smiles
from chemutils.props import LogP, OEXLogP, calculate_properties_for_smiles


@pytest.mark.no_chemaxon
def test_calculate_property_not_chemaxon():
    calculate_properties_for_smiles("CCCC(=O)O", [OEXLogP])


@pytest.mark.no_chemaxon
def test_calculate_property_chemaxon():
    # Chemaxon error is wrapped behind ValueError
    with pytest.raises(ValueError):
        calculate_properties_for_smiles("CCCC(=O)O", [LogP])


@pytest.mark.no_chemaxon
def test_as_cxmol():
    oemol = oemol_from_smiles("CCCC(=O)O")
    with pytest.raises(ChemaxonNotInstalledError):
        _ = as_cxmol(oemol)
