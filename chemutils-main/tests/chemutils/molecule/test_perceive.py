import pytest
from openeye import oechem

from chemutils.molecule import oemol_from_smiles
from chemutils.molecule.perceive import perceive_atom_properties


@pytest.mark.parametrize("smiles", ["CC(=O)O", "O=C(O)C"])
def test_perceive_atom_properties(smiles):
    oemol = oemol_from_smiles(smiles)

    assert all(atom.GetHyb() == oechem.OEHybridization_Unknown for atom in oemol.GetAtoms())

    perceive_atom_properties(oemol)

    assert all(atom.GetHyb() != oechem.OEHybridization_Unknown for atom in oemol.GetAtoms())
