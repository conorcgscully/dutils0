import numpy as np
from openeye import oechem

from chemutils.molecule.transform import (
    AffineTransformation,
    get_transformation_from_rotation_matrix,
    get_transformation_from_translation,
)


def get_rmsd_for_match(match: oechem.OEMatchBase, /) -> float:
    """
    Calculate the RMSD between two sets of atoms, given a match between them.

    No transformation is applied. For the minimized RMSD obtained by superimposing the atoms,
    use `get_optimal_rmsd_and_transformation_for_match`.

    Args:
        match: OpenEye match describing a correspondance between the atoms.
    """
    src = next(iter(match.GetPatternAtoms())).GetParent()
    dest = next(iter(match.GetTargetAtoms())).GetParent()

    rmsd: float = oechem.OERMSD(src, dest, match)
    if rmsd < 0:
        raise ValueError(f"RMSD calculation failed for match `{match}`.")

    return rmsd


def get_optimal_rmsd_and_transformation_for_match(
    match: oechem.OEMatchBase, /
) -> tuple[float, AffineTransformation]:
    """
    Calculate the minimised RMSD between two sets of atoms and the transformation that must be applied to obtain it.

    This uses the Kabsch algorithm to align the source atoms to the destination atoms.

    Args:
        match: OpenEye match describing a correspondance between the atoms.

    Returns:
        Tuple of the RMSD between the two sets of atoms, and the affine transformation matrix to
        be applied to the source atoms to obtain the optimal alignment.
    """
    src = next(iter(match.GetPatternAtoms())).GetParent()
    dest = next(iter(match.GetTargetAtoms())).GetParent()

    oe_rotation_matrix = oechem.OEDoubleArray(9)
    oe_translation = oechem.OEDoubleArray(3)

    # When using matches, `OERMSD` treats the query as the 'destination' and the target as the 'source'
    # We do the other way, and invert the transformation
    rmsd: float = oechem.OERMSD(src, dest, match, True, oe_rotation_matrix, oe_translation)
    if rmsd < 0:
        raise ValueError(f"RMSD calculation failed for match `{match}`.")

    rotation = get_transformation_from_rotation_matrix(np.reshape(oe_rotation_matrix, (3, 3)))
    translation = get_transformation_from_translation(np.reshape(oe_translation, (3,)))
    return rmsd, np.linalg.inv(translation @ rotation)
