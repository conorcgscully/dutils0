import numpy as np
import pytest

from chemutils.molecule import get_adjacency_matrix, std_oemol_from_smiles

from ...cip_validation import CIP_VALIDATION_OEMOLS


@pytest.mark.parametrize(
    ["smiles", "matrix"], [("CC", [[0, 1], [1, 0]]), ("CCC", [[0, 1, 0], [1, 0, 1], [0, 1, 0]])]
)
def test_get_adjaceny_matrix(smiles, matrix):
    oemol = std_oemol_from_smiles(smiles)
    np.testing.assert_equal(get_adjacency_matrix(oemol), np.array(matrix, dtype=np.bool_))


@pytest.mark.parametrize("oemol", CIP_VALIDATION_OEMOLS)
def test_get_adjaceny_matrix_properties(oemol):
    matrix = get_adjacency_matrix(oemol)
    assert matrix.dtype == np.bool_
    assert matrix.shape == (oemol.NumAtoms(), oemol.NumAtoms())
    # The leading diagonal should be False
    for i in range(oemol.NumAtoms()):
        assert not matrix[i, i]
    # Adjacency matrix must be symmetric (equal to transpose)
    np.testing.assert_equal(matrix, matrix.T)
    # The number of True's will be 2 * number of bonds (a to b and b to a).
    assert matrix.sum() == oemol.NumBonds() * 2
