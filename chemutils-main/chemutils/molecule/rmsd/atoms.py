import numpy as np
from openeye import oechem

from chemutils.molecule.transform import AffineTransformation, Vector3DList

from .points import get_optimal_rmsd_and_transformation_between_points, get_rmsd_between_points


def get_rmsd_for_atom_mapping(mapping: dict[oechem.OEAtomBase, oechem.OEAtomBase], /) -> float:
    """
    Calculate the RMSD between two sets of atoms, given a mapping between them.

    No transformation is applied. For the minimized RMSD obtained by superimposing the atoms,
    use `get_optimal_rmsd_and_transformation_for_atom_mapping`.

    Args:
        mapping: Mapping between two sets of atoms.

    Returns:
        RMSD between the two sets of atoms.
    """
    return get_rmsd_between_atoms(list(mapping.keys()), list(mapping.values()))


def get_optimal_rmsd_and_transformation_for_atom_mapping(
    mapping: dict[oechem.OEAtomBase, oechem.OEAtomBase], /
) -> tuple[float, AffineTransformation]:
    """
    Calculate the minimised RMSD between two sets of atoms and the transformation that must be applied to obtain it.

    This uses the Kabsch algorithm to align the source atoms to the destination atoms.

    Args:
        mapping: Mapping between two sets of atoms.

    Returns:
        Tuple of the RMSD between the two sets of atoms, and the affine transformation matrix to
        be applied to the source atoms to obtain the optimal alignment.
    """
    return get_optimal_rmsd_and_transformation_between_atoms(
        src=list(mapping.keys()), dest=list(mapping.values())
    )


def get_rmsd_between_atoms(src: list[oechem.OEAtomBase], dest: list[oechem.OEAtomBase]) -> float:
    """
    Calculate the RMSD between two sets of atoms, given their current coordinates.

    No transformation is applied. For the minimized RMSD obtained by superimposing the atoms,
    use `get_optimal_rmsd_and_transformation_between_atoms`.

    This function will return the same result regardless of the order of the two sets of atoms.

    Args:
        src: Source atoms.
        dest: Destination atoms.

    Returns:
        RMSD between the two sets of atoms, in Angstroms.
    """
    src_coords = get_coordinates_for_atoms(src)
    dest_coords = get_coordinates_for_atoms(dest)
    return get_rmsd_between_points(src_coords, dest_coords)


def get_optimal_rmsd_and_transformation_between_atoms(
    *, src: list[oechem.OEAtomBase], dest: list[oechem.OEAtomBase]
) -> tuple[float, AffineTransformation]:
    """
    Calculate the minimised RMSD between two sets of atoms and the transformation that must be applied to obtain it.

    This uses the Kabsch algorithm to align the source atoms to the destination atoms.

    Args:
        src: Source atoms.
        dest: Destination atoms.

    Returns:
        Tuple of the RMSD between the two sets of atoms, in Angstroms, and the affine transformation matrix to
        be applied to the source atoms to obtain the optimal alignment.
    """
    src_coords = get_coordinates_for_atoms(src)
    dest_coords = get_coordinates_for_atoms(dest)
    return get_optimal_rmsd_and_transformation_between_points(src=src_coords, dest=dest_coords)


def get_coordinates_for_atoms(atoms: list[oechem.OEAtomBase], /) -> Vector3DList:
    """Get the coordinates of a list of atoms."""
    parent: oechem.OEMolBase = atoms[0].GetParent()
    return np.array([parent.GetCoords(atom) for atom in atoms])
