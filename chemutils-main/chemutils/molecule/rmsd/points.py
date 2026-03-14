import numpy as np
from openeye import oechem

from chemutils.molecule.transform import (
    AffineTransformation,
    Vector3DListLike,
    get_transformation_from_rotation_matrix,
    get_transformation_from_translation,
)


def get_rmsd_between_points(src: Vector3DListLike, dest: Vector3DListLike, /) -> float:
    """
    Calculate the RMSD between two sets of points, without applying any transformation.

    No transformation is applied. For the minimized RMSD obtained by superimposing the points,
    use `get_optimal_rmsd_and_transformation_between_points`.

    This function will return the same result regardless of the order of the two arguments.

    Args:
        src: Source points.
        dest: Destination points.

    Returns:
        RMSD between the two sets of points.

    Raises:
        ValueError: If the points are not lists of 3D vectors.
        ValueError: If the two sets of points are not the same length.
    """
    src = np.asarray(src)
    dest = np.asarray(dest)

    if len(src.shape) != 2 or src.shape[1] != 3:
        raise ValueError(
            f"RMSD argument must be list of 3D vectors: expected shape `(N, 3)`, but got `{src.shape}`"
        )
    if len(dest.shape) != 2 or dest.shape[1] != 3:
        raise ValueError(
            f"RMSD argument must be list of 3D vectors: expected shape `(N, 3)`, but got `{dest.shape}`"
        )
    if src.shape != dest.shape:
        raise ValueError(
            f"RMSD arguments must be the same length: lhs has shape `{src.shape}`, rhs has shape `{dest.shape}"
        )

    size = src.shape[0]
    src = oechem.OEDoubleArray(src.flatten())
    dest = oechem.OEDoubleArray(dest.flatten())

    rmsd: float = oechem.OERMSD(dest, src, size)
    if rmsd < 0:
        raise ValueError("RMSD calculation failed.")
    return rmsd


def get_optimal_rmsd_and_transformation_between_points(
    *, src: Vector3DListLike, dest: Vector3DListLike
) -> tuple[float, AffineTransformation]:
    """
    Calculate the RMSD between two sets of points and the transformation that must be applied to obtain it.

    Args:
        src: Source points.
        dest: Destination points.

    Returns:
        RMSD between the two sets of points, and the affine transformation matrix to be applied to the source
        points to obtain the optimal alignment.

    Raises:
        ValueError: If the points are not lists of 3D vectors.
        ValueError: If the two sets of points are not the same length.
    """
    src = np.asarray(src)
    dest = np.asarray(dest)

    if len(src.shape) != 2 or src.shape[1] != 3:
        raise ValueError(
            f"RMSD argument must be list of 3D vectors: expected shape `(N, 3)`, but got `{src.shape}`"
        )
    if len(dest.shape) != 2 or dest.shape[1] != 3:
        raise ValueError(
            f"RMSD argument must be list of 3D vectors: expected shape `(N, 3)`, but got `{dest.shape}`"
        )
    if src.shape != dest.shape:
        raise ValueError(
            f"RMSD arguments must be the same length: lhs has shape `{src.shape}`, rhs has shape `{dest.shape}"
        )

    size = src.shape[0]
    src = oechem.OEDoubleArray(src.flatten())
    dest = oechem.OEDoubleArray(dest.flatten())
    oe_rotation_matrix = oechem.OEDoubleArray(9)
    oe_translation = oechem.OEDoubleArray(3)

    rmsd = oechem.OERMSD(dest, src, size, True, oe_rotation_matrix, oe_translation)
    if rmsd < 0:
        raise ValueError("RMSD calculation failed.")

    rotation = get_transformation_from_rotation_matrix(np.reshape(oe_rotation_matrix, (3, 3)))
    translation = get_transformation_from_translation(np.reshape(oe_translation, (3,)))
    return rmsd, translation @ rotation
