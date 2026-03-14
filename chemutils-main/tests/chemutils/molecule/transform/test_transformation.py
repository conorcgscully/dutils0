import numpy as np
import pytest

from chemutils.molecule.transform import (
    get_identity_transformation,
    get_transformation_from_rotation_matrix,
    get_transformation_from_translation,
    transform_points,
)


def test_get_identify_transformation():
    np.testing.assert_allclose(
        get_identity_transformation(),
        [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]],
    )


@pytest.mark.parametrize("offset", [[0, 0, 0], (1, 4, -2), np.array([-2, 0, 3])])
def test_get_transformation_from_translation(offset):
    np.testing.assert_allclose(
        get_transformation_from_translation(offset),
        [[1, 0, 0, offset[0]], [0, 1, 0, offset[1]], [0, 0, 1, offset[2]], [0, 0, 0, 1]],
    )


@pytest.mark.parametrize("offset", [[0, 0], (1, 4, -2, 1), np.array([[-2, 0, 3]])])
def test_get_transformation_from_translation_invalid_offset(offset):
    with pytest.raises(ValueError) as exc:
        _ = get_transformation_from_translation(offset)
    assert (
        str(exc.value)
        == f"Translation must be a 3D vector: expected shape `(3, )`, but got `{np.asarray(offset).shape}`."
    )


@pytest.mark.parametrize(
    "rotation", [[[1, 0, 0], [0, 1, 0], [0, 0, 1]], [[0, 1, 0], [-1, 0, 0], [0, 0, 1]]]
)
def test_get_transformation_from_rotation_matrix(rotation):
    np.testing.assert_allclose(
        get_transformation_from_rotation_matrix(rotation),
        [
            [rotation[0][0], rotation[0][1], rotation[0][2], 0],
            [rotation[1][0], rotation[1][1], rotation[1][2], 0],
            [rotation[2][0], rotation[2][1], rotation[2][2], 0],
            [0, 0, 0, 1],
        ],
    )


@pytest.mark.parametrize(
    "rotation", [[[1, 0, 0], [0, 1, 0], [0, 0, 1]], [[0, 1, 0], [-1, 0, 0], [0, 0, 1]]]
)
def test_get_transformation_from_rotation_matrix_origin(rotation):
    origin = [2, -1, 3]
    affine = get_transformation_from_rotation_matrix(rotation, origin=origin)

    # Check that rotation matrix is orthogonal
    np.testing.assert_allclose(affine[:3, :3].T @ affine[:3, :3], np.eye(3))

    # Check that the origin point is preserved
    np.testing.assert_allclose(transform_points([origin], transformation=affine), [origin])


@pytest.mark.parametrize("rotation", [[[0, 1, 0], [1, 0, 0]], np.array([[1, 1], [0, 1], [0, 1]])])
def test_get_transformation_from_rotation_matrix_rotation_invalid_shape(rotation):
    with pytest.raises(ValueError) as exc:
        _ = get_transformation_from_rotation_matrix(rotation)
    assert (
        str(exc.value)
        == f"Rotation must be a 3x3 matrix: expected shape `(3, 3)`, but got `{np.asarray(rotation).shape}`."
    )


@pytest.mark.parametrize(
    "origin",
    [
        [0, 0],
        (1, 2, 4, -2),
        np.array([[-2, 0, 3]]),
    ],
)
def test_get_transformation_from_rotation_matrix_origin_invalid_shape(origin):
    with pytest.raises(ValueError) as exc:
        _ = get_transformation_from_rotation_matrix(
            [[1, 0, 0], [0, 1, 0], [0, 0, 1]], origin=origin
        )
    assert (
        str(exc.value)
        == f"Origin must be a 3D vector: expected shape `(3, )`, but got `{np.asarray(origin).shape}`."
    )


def test_get_transformation_from_rotation_matrix_rotation_not_orthogonal():
    rotation = [[1, 0, 0], [0, 2, 0], [0, 0, -2]]

    with pytest.raises(ValueError) as exc:
        _ = get_transformation_from_rotation_matrix(rotation)
    assert str(exc.value) == "Rotation matrix must be orthogonal."


@pytest.mark.parametrize(
    "points, transformation, expected",
    [
        ([[1, 0, 3], [0, -1, 2]], get_identity_transformation(), [[1, 0, 3], [0, -1, 2]]),
        (
            [[1, 0, 3], [0, -1, 2]],
            get_transformation_from_translation([1, 2, 3]),
            [[2, 2, 6], [1, 1, 5]],
        ),
        (
            [[1, 0, 3], [0, -1, 2]],
            get_transformation_from_rotation_matrix([[0, 1, 0], [1, 0, 0], [0, 0, 1]]),
            [[0, 1, 3], [-1, 0, 2]],
        ),
    ],
)
def test_transform_points(points, transformation, expected):
    np.testing.assert_allclose(transform_points(points, transformation=transformation), expected)
