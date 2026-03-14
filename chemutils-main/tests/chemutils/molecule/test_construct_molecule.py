import numpy as np
import pytest

from chemutils.molecule.construct_molecule import construct_molecule

from .test_standardise import STANDARDISE_CASES


@pytest.mark.parametrize(
    ["smiles", "coords_shape"],
    [
        ("C", (2, 3)),
        ("CCC", (2, 3)),
        ("CCC", (4, 2)),
        ("CC", (2, 4)),
        ("CC", (2, 1)),
        ("CC", (2, 2)),
        ("CC", (2, 0)),
    ],
)
def test_construct_molecule_wrong_coords_size(smiles, coords_shape):
    coordinates = np.zeros(coords_shape)
    with pytest.raises(ValueError):
        _ = construct_molecule(smiles=smiles, coordinates=coordinates)


@pytest.mark.parametrize(["smiles", "num_atoms", "expected_smiles"], STANDARDISE_CASES)
def test_construct_molecule(smiles, num_atoms, expected_smiles):
    pos = np.random.randn(num_atoms, 3).astype(np.float32)
    mol = construct_molecule(smiles=smiles, coordinates=pos)
    assigned_pos = np.array(list(mol.GetCoords().values()))
    assert np.all(pos == assigned_pos), "Assigned atom positions are different from intended"
    assert mol.GetDimension() == 3
