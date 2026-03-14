import numpy as np
import pytest
from openeye import oechem

from chemutils.molecule import (
    create_atom,
    get_coordinates,
    oemol_from_smiles,
    rotate_molecule,
    transform_molecule,
    translate_molecule,
)
from chemutils.molecule.transform import get_identity_transformation

TRANSLATION = [1, 4, -2]

ROTATION = [[0, -1, 0], [1, 0, 0], [0, 0, 1]]

ORIGIN = [0, 1, 0]

COORDS_ORIG = [[1, 0, 0], [0, 1, 0]]

COORDS_TRANSLATED = [[2, 4, -2], [1, 5, -2]]

COORDS_ROTATED = [[0, 1, 0], [-1, 0, 0]]

COORDS_ROTATED_ORIGIN = [[1, 2, 0], [0, 1, 0]]


@pytest.fixture()
def mol():
    oemol = oechem.OEGraphMol()
    create_atom(oemol, atomic_number=1, coordinates=COORDS_ORIG[0])
    create_atom(oemol, atomic_number=17, coordinates=COORDS_ORIG[1])
    oemol.SetDimension(3)
    return oemol


def test_translate_mol(mol):
    newmol = translate_molecule(mol, translation=TRANSLATION)
    np.testing.assert_allclose(get_coordinates(mol), COORDS_ORIG)
    np.testing.assert_allclose(get_coordinates(newmol), COORDS_TRANSLATED)


def test_translate_copy_mol(mol):
    newmol = translate_molecule(mol, translation=TRANSLATION, copy_mol=False)
    assert newmol is mol
    np.testing.assert_allclose(get_coordinates(mol), COORDS_TRANSLATED)


def test_rotate_mol(mol):
    newmol = rotate_molecule(mol, rotation=ROTATION)
    np.testing.assert_allclose(get_coordinates(mol), COORDS_ORIG)
    np.testing.assert_allclose(get_coordinates(newmol), COORDS_ROTATED)


def test_rotate_copy_mol(mol):
    newmol = rotate_molecule(mol, rotation=ROTATION, copy_mol=False)
    assert newmol is mol
    np.testing.assert_allclose(get_coordinates(mol), COORDS_ROTATED)


def test_rotate_origin(mol):
    newmol = rotate_molecule(mol, rotation=ROTATION, origin=ORIGIN)
    np.testing.assert_allclose(get_coordinates(mol), COORDS_ORIG)
    np.testing.assert_allclose(get_coordinates(newmol), COORDS_ROTATED_ORIGIN)


def test_transform_molecule_no_3d():
    oemol = oemol_from_smiles("CCC")

    with pytest.raises(ValueError) as exc:
        _ = transform_molecule(oemol, transformation=get_identity_transformation)

    assert str(exc.value) == "Molecule must have 3D coordinates to be transformed."


def test_transform_molecule_not_rigid(mol):
    with pytest.raises(ValueError) as exc:
        _ = transform_molecule(
            mol, transformation=[[1, 0, 0, 0], [0, 2, 0, 0], [0, 0, 3, 0], [0, 0, 0, 1]]
        )

    assert str(exc.value) == "Transformation matrix must be a rigid transformation."


def test_transform_molecule_not_affine(mol):
    with pytest.raises(ValueError) as exc:
        _ = transform_molecule(
            mol, transformation=[[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 2]]
        )

    assert (
        str(exc.value)
        == "Transformation matrix must be an affine transformation: expected bottom row of `[0 0 0 1]`, but got `[0 0 0 2]`."
    )


def test_translate_molecule_no_3d():
    oemol = oemol_from_smiles("CCC")

    with pytest.raises(ValueError) as exc:
        _ = translate_molecule(oemol, translation=[1, 2, 3])

    assert str(exc.value) == "Molecule must have 3D coordinates to be transformed."


def test_translate_molecule_invalid_offset(mol):
    with pytest.raises(ValueError) as exc:
        _ = translate_molecule(mol, translation=[1, 2])

    assert (
        str(exc.value) == "Translation must be a 3D vector: expected shape `(3, )`, but got `(2,)`."
    )


def test_rotate_molecule_no_3d():
    oemol = oemol_from_smiles("CCC")

    with pytest.raises(ValueError) as exc:
        _ = rotate_molecule(oemol, rotation=[[1, 0, 0], [0, 1, 0], [0, 0, 1]])

    assert str(exc.value) == "Molecule must have 3D coordinates to be transformed."


def test_rotate_molecule_rotation_not_orthogonal(mol):
    with pytest.raises(ValueError) as exc:
        _ = rotate_molecule(mol, rotation=[[1, 0, 0], [0, 1, 0], [0, 0, 2]])

    assert str(exc.value) == "Rotation matrix must be orthogonal."


def test_rotate_molecule_rotation_not_3x3(mol):
    with pytest.raises(ValueError) as exc:
        _ = rotate_molecule(mol, rotation=[[1, 0], [0, 1]])

    assert (
        str(exc.value)
        == "Rotation must be a 3x3 matrix: expected shape `(3, 3)`, but got `(2, 2)`."
    )


def test_rotate_molecule_origin_not_3d(mol):
    with pytest.raises(ValueError) as exc:
        _ = rotate_molecule(mol, rotation=[[1, 0, 0], [0, 1, 0], [0, 0, 1]], origin=[1, 2])

    assert str(exc.value) == "Origin must be a 3D vector: expected shape `(3, )`, but got `(2,)`."
