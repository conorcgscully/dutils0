import numpy as np
from openeye import oechem

from chemutils.molecule import smiles_from_oemol
from chemutils.molecule.transform import (
    AffineTransformation,
    get_transformation_from_rotation_matrix,
    get_transformation_from_translation,
)


def get_rmsd_between_molecules(src: oechem.OEMolBase, dest: oechem.OEMolBase, /) -> float:
    """
    Calculate the RMSD between two molecules, given their current coordinates.

    No transformation is applied. For the minimized RMSD obtained by superimposing the molecules,
    use `get_optimal_rmsd_and_transformation_between_molecules`.

    This function will return the same result regardless of the order of the two molecules.

    Args:
        src: Source molecule.
        dest: Destination molecule.

    Returns:
        RMSD between the two molecules, in Angstroms.

    Raises:
        ValueError: If the molecules do not represent the same molecule.
    """
    rmsd: float = oechem.OERMSD(dest, src)
    if rmsd < 0:
        raise ValueError(
            f"RMSD calculation failed: failed to map `{smiles_from_oemol(src)}` to `{smiles_from_oemol(dest)}`."
        )
    return rmsd


def get_optimal_rmsd_and_transformation_between_molecules(
    *, src: oechem.OEMolBase, dest: oechem.OEMolBase
) -> tuple[float, AffineTransformation]:
    """
    Calculate the minimised RMSD between two molecules and the transformation that must be applied to obtain it.

    This uses the Kabsch algorithm to align the source molecule to the destination molecule.

    Args:
        src: Source molecule.
        dest: Destination molecule.

    Returns:
        Tuple of the RMSD between the two molecules, in Angstroms, and the affine transformation matrix to
        be applied to the source molecule to obtain the optimal alignment.

    Raises:
        ValueError: If the molecules do not represent the same molecule.
    """
    oe_rotation_matrix = oechem.OEDoubleArray(9)
    oe_translation = oechem.OEDoubleArray(3)

    rmsd = oechem.OERMSD(dest, src, True, True, True, oe_rotation_matrix, oe_translation)
    if rmsd < 0:
        raise ValueError(
            f"RMSD calculation failed: failed to map `{smiles_from_oemol(src)}` to `{smiles_from_oemol(dest)}`."
        )

    rotation = get_transformation_from_rotation_matrix(np.reshape(oe_rotation_matrix, (3, 3)))
    translation = get_transformation_from_translation(np.reshape(oe_translation, (3,)))
    return rmsd, translation @ rotation
