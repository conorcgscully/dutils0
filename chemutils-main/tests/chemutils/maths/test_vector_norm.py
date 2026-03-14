import numpy as np
import pytest
from openeye import oechem

from chemutils.maths import clamp_magnitude, distance_between, distance_matrix, length, normalized
from chemutils.molecule import create_atom, create_bond


@pytest.mark.parametrize(
    ["arg", "expected"],
    [
        ([1, 1, 1], 1.732),
        ([0.5, 1.5, 2.5], 2.958),
        ([[1, 1, 1], [0.5, 1.5, 2.5]], [1.732, 2.958]),
    ],
)
def test_length(arg, expected):
    assert length(arg) == pytest.approx(expected, abs=1e-3)


def test_length_bond():
    mol = oechem.OEMol()
    atom1 = create_atom(mol, atomic_number=6, coordinates=(-1, 0, 2))
    atom2 = create_atom(mol, atomic_number=6, coordinates=(0, 1, 1))
    bond = create_bond(mol, atom1=atom1, atom2=atom2)

    assert length(bond) == pytest.approx(1.732, abs=1e-3)


def test_normalized():
    vector = np.array([1, -2, 2])

    assert length(vector) == pytest.approx(3, abs=1e-3)

    normed = normalized(vector)

    assert length(normed) == pytest.approx(1, abs=1e-3)
    np.testing.assert_allclose(normed, [1 / 3, -2 / 3, 2 / 3], atol=1e-3)


def test_distance_between_vectors():
    vec1 = np.array([-1, 2, 0])
    vec2 = np.array([1, 1, 2])

    assert distance_between(vec1, vec2) == pytest.approx(3.0, abs=1e-3)


def test_distance_between_vector_and_vector_list():
    vec1 = np.array([[-1, 2, 0], [0, 2, -1]])
    vec2 = np.array([1, 3, -1])

    assert distance_between(vec1, vec2) == pytest.approx([2.450, 1.414], abs=1e-3)


def test_distance_between_vector_lists():
    vec1 = np.array([[-1, 2, 0], [0, 2, -1]])
    vec2 = np.array([[1, 3, -1], [1, 1, 2]])

    np.testing.assert_allclose(distance_between(vec1, vec2), [2.450, 3.317], atol=1e-3)


def test_distance_between_vector_lists_outer():
    vec1 = np.array([[-1, 2, 0], [0, 2, -1]])
    vec2 = np.array([[1, 3, -1], [1, 1, 2]])

    distances = distance_between(vec1, vec2, mode="outer")

    np.testing.assert_allclose(
        distance_between(vec1, vec2, mode="outer"), [[2.450, 3.0], [1.414, 3.317]], atol=1e-3
    )

    assert distances[0, 0] == pytest.approx(distance_between(vec1[0], vec2[0]))
    assert distances[1, 0] == pytest.approx(distance_between(vec1[1], vec2[0]))


def test_distance_between_atoms():
    mol = oechem.OEMol()
    atom1 = create_atom(mol, atomic_number=6, coordinates=(-1, 0, 1))
    atom2 = create_atom(mol, atomic_number=6, coordinates=(5, -3, 2))

    assert distance_between(atom1, atom2) == pytest.approx(6.782, abs=1e-3)


def test_distance_matrix_vector():
    vec1 = np.array([[-1, 2, 0], [0, 2, -1]])

    np.testing.assert_allclose(distance_matrix(vec1), [[0, 1.414], [1.414, 0]], atol=1e-3)


def test_distance_matrix_vector_lists():
    vec1 = np.array([[-1, 2, 0], [0, 2, -1]])
    vec2 = np.array([[1, 3, -1], [1, 1, 2]])

    np.testing.assert_allclose(
        distance_matrix(vec1, vec2), [[2.450, 3.0], [1.414, 3.317]], atol=1e-3
    )


def test_clamp_magnitude():
    vec1 = np.array([[2, 1, 2]])

    assert length(vec1) == pytest.approx(3)

    vec1_clamped = clamp_magnitude(vec1, max_magnitude=2)

    assert length(vec1_clamped) == pytest.approx(2)
    np.testing.assert_allclose(vec1_clamped, [[4 / 3, 2 / 3, 4 / 3]], atol=1e-3)

    vec2 = np.array([[2, 1, 2]])

    vec2_clamped = clamp_magnitude(vec2, max_magnitude=5)
    np.testing.assert_allclose(vec2_clamped, vec2, atol=1e-3)


def test_clamp_magnitude_vector_list():
    vecs = np.array([[2, 1, 2], [-2, 3, 6]])

    vecs_clamped = clamp_magnitude(vecs, max_magnitude=5)

    np.testing.assert_allclose(
        vecs_clamped, [[2, 1, 2], [-2 * 5 / 7, 3 * 5 / 7, 6 * 5 / 7]], atol=1e-3
    )
