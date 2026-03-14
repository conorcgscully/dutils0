import numpy as np
from openeye import oechem

from .transformation import (
    AffineTransformationLike,
    RotationMatrixLike,
    Vector3DLike,
    get_transformation_from_rotation_matrix,
    get_transformation_from_translation,
)


def transform_molecule(
    mol: oechem.OEMolBase, /, *, transformation: AffineTransformationLike, copy_mol: bool = True
) -> oechem.OEMolBase:
    """
    Transform the 3D coordinates of the given molecule using an affine transformation matrix.

    The affine transformation must be a rigid transformation, namely it should be distance-preserving
    and consist of purely a composition of translations, rotations and reflections.

    Args:
        mol: Molecule to transform.
        transformation: Affine transformation matrix, in augmented matrix form. Accepts a NumPy array
            or an array-like object that can be converted to a 4x4 NumPy array. An exception will be
            raised if the matrix does not represent a rigid transformation.
        copy_mol: If `False`, the molecule is modified in place. When `True` (the default), a copy
            of the molecule is returned, leaving the original with its coordinates unaltered.

    Returns:
        The transformed molecule. If `copy_mol` was `False`, this will be the same object as `mol`.

    Raises:
        ValueError: If the molecule does not have 3D coordinates.
        ValueError: If the transformation matrix is not a rigid transformation.
        ValueError: If the transformation matrix is not an affine transformation.
    """
    if not mol.GetDimension() == 3:
        raise ValueError("Molecule must have 3D coordinates to be transformed.")

    if copy_mol:
        mol = mol.CreateCopy()

    transformation = np.asarray(transformation)

    rotation_matrix = transformation[:3, :3]

    if not np.allclose(rotation_matrix.T @ rotation_matrix, np.eye(3), atol=0.001):
        raise ValueError("Transformation matrix must be a rigid transformation.")

    translation = transformation[:3, 3]

    bottom = transformation[3, :]
    if not np.allclose(bottom, [0, 0, 0, 1]):
        raise ValueError(
            f"Transformation matrix must be an affine transformation: expected bottom row of `[0 0 0 1]`, but got `{bottom}`."
        )

    oe_rotation_matrix = oechem.OEDoubleArray(rotation_matrix.flatten().tolist())
    oe_translation = oechem.OEDoubleArray(translation.flatten().tolist())

    oechem.OERotate(mol, oe_rotation_matrix)
    oechem.OETranslate(mol, oe_translation)

    return mol


def translate_molecule(
    mol: oechem.OEMolBase, /, *, translation: Vector3DLike, copy_mol: bool = True
) -> oechem.OEMolBase:
    """
    Translate a molecule in 3D space by a translation vector.

    Args:
        mol: Molecule to translate.
        translation: Translation as a NumPy array, or as a array-like object that can be converted
            to a NumPy array of length 3.
        copy_mol: If `False`, the molecule is modified in place. When `True` (the default), a copy
            of the molecule is returned, leaving the original with its coordinates unaltered.

    Returns:
        The translated molecule. If `copy_mol` was `False`, this will be the same object as `mol`.
    """
    transformation = get_transformation_from_translation(translation)
    return transform_molecule(mol, transformation=transformation, copy_mol=copy_mol)


def rotate_molecule(
    mol: oechem.OEMolBase,
    /,
    *,
    rotation: RotationMatrixLike,
    origin: Vector3DLike | None = None,
    copy_mol: bool = True,
) -> oechem.OEMolBase:
    """
    Rotate the given molecule by the given rotation matrix.

    Args:
        mol: Molecule to rotate.
        rotation: Rotation matrix. Accepts a NumPy array or an array-like object that can be
            converted to a 3x3 NumPy array.
        origin: Origin about which to perform the rotation. When not provided, `(0, 0, 0)` will be used.
        copy_mol: If `False`, the molecule is modified in place. When `True` (the default), a copy
            of the molecule is returned, leaving the original with its coordinates unaltered.

    Returns:
        The rotated molecule. If `copy_mol` was `False`, this will be the same object as `mol`.
    """
    transformation = get_transformation_from_rotation_matrix(rotation, origin=origin)
    return transform_molecule(mol, transformation=transformation, copy_mol=copy_mol)
