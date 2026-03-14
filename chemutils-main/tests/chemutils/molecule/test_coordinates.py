import numpy as np
import pytest
from openeye import oechem

from chemutils.molecule import get_coordinates, oemol_from_smiles, set_coordinates
from chemutils.molecule.embed import embed_molecule_3d


@pytest.mark.parametrize(
    ["dim", "coords"],
    [
        (1, np.array([[0], [1], [2]])),
        (2, np.array([[0, 0], [0, 1], [2, -1]])),
        (3, np.array([[0, 0, 1], [1, 2, 0], [-2, 1, 0]])),
    ],
)
@pytest.mark.parametrize(
    "copy",
    (True, False),
)
def test_set_coordinates_dimensionality(dim, coords, copy):
    smiles = "CC=O"
    oemol = oemol_from_smiles(smiles)
    assert oemol.GetDimension() == 0
    oemol2 = set_coordinates(oemol, coords, copy_mol=copy)
    if copy:
        assert oemol.GetDimension() == 0
    else:
        assert oemol2 is None
        oemol2 = oemol
    assert oemol2.GetDimension() == dim
    np.testing.assert_allclose(get_coordinates(oemol2), coords)


def test_get_coordinates_none():
    smiles = "CC=O"
    oemol = oemol_from_smiles(smiles)
    with pytest.raises(ValueError):
        _ = get_coordinates(oemol)


@pytest.mark.parametrize(
    ["dim", "coords"],
    [
        (1, np.array([[[0], [1], [2]], [[3], [-3], [0]]])),
        (2, np.array([[[0, 0], [0, 1], [2, -1]], [[1, -1], [2, -2], [3, -3]]])),
        (3, np.array([[[0, 0, 1], [1, 2, 0], [-2, 1, 0]], [[1, 2, 3], [-1, -2, -3], [0, 0, 0]]])),
    ],
)
def test_set_coordinates_conformers(dim, coords):
    smiles = "CC=O"
    oemol = oemol_from_smiles(smiles)
    assert oemol.GetDimension() == 0
    oemol2 = set_coordinates(oemol, coords, copy_mol=True)
    assert isinstance(oemol2, oechem.OEMol)
    assert oemol2.GetDimension() == dim
    np.testing.assert_allclose(get_coordinates(oemol2), coords)


@pytest.mark.parametrize(
    ["dim", "coords"],
    [
        (1, np.array([[[0], [1], [2]], [[3], [-3], [0]]])),
        (2, np.array([[[0, 0], [0, 1], [2, -1]], [[1, -1], [2, -2], [3, -3]]])),
        (3, np.array([[[0, 0, 1], [1, 2, 0], [-2, 1, 0]], [[1, 2, 3], [-1, -2, -3], [0, 0, 0]]])),
    ],
)
def test_set_coordinates_conformers_nocopy(dim, coords):
    smiles = "CC=O"
    oemol = oemol_from_smiles(smiles)
    embed_molecule_3d(oemol, seed=0)
    oemol = oechem.OEMol(oemol)
    assert get_coordinates(oemol).shape == (1, 3, 3)
    set_coordinates(oemol, coords)
    assert oemol.GetDimension() == dim
    np.testing.assert_allclose(get_coordinates(oemol), coords)


@pytest.mark.parametrize(
    ["dim", "coords"],
    [
        (1, np.array([[[0], [1], [2]], [[3], [-3], [0]]])),
        (2, np.array([[[0, 0], [0, 1], [2, -1]], [[1, -1], [2, -2], [3, -3]]])),
        (3, np.array([[[0, 0, 1], [1, 2, 0], [-2, 1, 0]], [[1, 2, 3], [-1, -2, -3], [0, 0, 0]]])),
    ],
)
def test_set_coordinates_conformers_error_oegraphmol(dim, coords):
    smiles = "CC=O"
    oemol = oemol_from_smiles(smiles)
    with pytest.raises(ValueError):
        set_coordinates(oemol, coords)
