from collections.abc import Sequence

import numpy as np
import numpy.typing as npt

AffineTransformationLike = npt.NDArray[np.float32] | Sequence[Sequence[float]]
"""Any type that can be converted to a NumPy array of shape `(4, 4)`."""

AffineTransformation = npt.NDArray[np.float32]
"""NumPy array of shape `(4, 4)`."""

Vector3DLike = npt.NDArray[np.float32] | Sequence[float]
"""Any type that can be converted to a NumPy array of shape `(3,)`."""

Vector3DListLike = npt.NDArray[np.float32] | Sequence[Sequence[float]]
"""Any type that can be converted to a NumPy array of shape `(N, 3)`."""

Vector3DList = npt.NDArray[np.float32]
"""A NumPy array of shape `(N, 3)`."""

RotationMatrixLike = npt.NDArray[np.float32] | Sequence[Sequence[float]]
"""Any type that can be converted to a NumPy array of shape `(3, 3)`."""


def get_identity_transformation() -> AffineTransformation:
    """Get the identity affine transformation matrix."""
    return np.eye(4, dtype=np.float32)


def get_transformation_from_translation(translation: Vector3DLike, /) -> AffineTransformation:
    """
    Get the affine transformation matrix for a 3D translation.

    Args:
        translation: Translation, as a NumPy array or an array-like object that can be converted
            to a NumPy array of length 3.

    Returns:
        Affine transformation matrix representing the translation, as a NumPy array of shape `(4, 4)`.

    Raises:
        ValueError: If the translation is not a 3D vector.
    """
    translation = np.asarray(translation)
    if translation.shape != (3,):
        raise ValueError(
            f"Translation must be a 3D vector: expected shape `(3, )`, but got `{translation.shape}`."
        )

    affine = get_identity_transformation()
    affine[:3, 3] = translation
    return affine


def get_transformation_from_rotation_matrix(
    rotation: RotationMatrixLike, /, *, origin: Vector3DLike | None = None
) -> AffineTransformation:
    """
    Get the affine transformation matrix for a 3D rotation from a rotation matrix.

    The rotation matrix may be proper (maintains handedness) or improper (reverses handedness,
    involving a reflection).

    Args:
        rotation: Rotation matrix, as a NumPy array or an array-like object that can be converted to
            a NumPy array of shape `(3, 3)`.
        origin: Origin of the rotation, as a NumPy array or an array-liek object that can be converted
            to a NumPy array of length `3`. When not provided, `(0, 0, 0)` will be used.

    Returns:
        The affine transformation matrix for the rotation.
    """
    if origin is not None:
        origin = np.asarray(origin)
        if origin.shape != (3,):
            raise ValueError(
                f"Origin must be a 3D vector: expected shape `(3, )`, but got `{origin.shape}`."
            )
        return (
            get_transformation_from_translation(origin)
            @ get_transformation_from_rotation_matrix(rotation)
            @ get_transformation_from_translation(-origin)
        )

    rotation = np.asarray(rotation)
    if rotation.shape != (3, 3):
        raise ValueError(
            f"Rotation must be a 3x3 matrix: expected shape `(3, 3)`, but got `{rotation.shape}`."
        )

    if not np.allclose(rotation.T @ rotation, np.eye(3)):
        raise ValueError("Rotation matrix must be orthogonal.")

    affine = get_identity_transformation()
    affine[:3, :3] = rotation
    return affine


def transform_points(
    points: Vector3DListLike, /, *, transformation: AffineTransformationLike
) -> Vector3DList:
    """
    Transform a list of 3D points using an affine transformation matrix.

    The transformation can be any affine transformation, including translations, rotations,
    reflections, scalings, shears and any combination thereof. It is guaranteed to preserve lines
    and parallelism, but not distances or angles.

    Args:
        points: List of 3D points to transform. Accepts NumPy arrays or array-like objects such as
            list or tuples that can be converted to NumPy arrays.
        transformation: Affine transformation matrix, in augmented matrix form. Accepts a NumPy array
            or an array-like object that can be converted to a 4x4 NumPy array.

    Raises:
        ValueError: If the points are not convertible to a NumPy array with shape `(N, 3)`.
        ValueError: If the transformation matrix is not convertible to a NumPy array with shape `(4, 4)`.
        ValueError: If the last row of the transformation matrix is not `[0, 0, 0, 1]`, and hence it does
            does not represent an affine transformation.
    """
    points = np.asarray(points)
    if len(points.shape) != 2 or points.shape[1] != 3:
        raise ValueError(
            f"Points must be a list of 3D vectors: expected shape `(N, 3)`, but got `{points.shape}`"
        )
    transformation = np.asarray(transformation)
    if transformation.shape != (4, 4):
        raise ValueError(
            f"Transformation must be an affine transformation matrix: expected shape `(4, 4)`, but got `{transformation.shape}`"
        )
    if not np.allclose(transformation[3, :], [0, 0, 0, 1]):
        raise ValueError(
            f"Transformation must be an affine transformation matrix: expected bottom row of `[0, 0, 0, 1]`, but got `{transformation[3, :]}`."
        )

    points = np.asarray(points).T
    points = np.append(points, np.ones((1, points.shape[1])), axis=0)
    points = transformation @ points
    return points[:3, :].T
