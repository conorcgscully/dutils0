import numpy as np
import pytest

from chemutils.math import get_signed_dihedral_angle
from chemutils.molecule.tetrahedron import IDEALIZED_TETRAHEDRON


@pytest.mark.parametrize(
    "a, b, c, d, expected",
    [
        (
            np.array([1, 0, 0]),
            np.array([0, 0, 0]),
            np.array([0, 1, 0]),
            np.array([0, 0, 1]),
            -np.pi / 2,
        ),
        # Add more test cases as needed
    ],
)
def test_get_tetrahedral_dihedral(a, b, c, d, expected):
    result = get_signed_dihedral_angle(a, b, c, d)
    np.testing.assert_almost_equal(result, expected)


def test_dihedral_inversion():
    a1 = get_signed_dihedral_angle(
        IDEALIZED_TETRAHEDRON[0],
        IDEALIZED_TETRAHEDRON[1],
        IDEALIZED_TETRAHEDRON[2],
        IDEALIZED_TETRAHEDRON[3],
    )
    a2 = get_signed_dihedral_angle(
        IDEALIZED_TETRAHEDRON[0],
        IDEALIZED_TETRAHEDRON[1],
        IDEALIZED_TETRAHEDRON[3],
        IDEALIZED_TETRAHEDRON[2],
    )
    assert a1 == -1 * a2
