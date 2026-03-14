from functools import lru_cache
from itertools import permutations
from typing import Literal

import numpy as np
import numpy.typing as npt
from openeye import oechem

from chemutils.math import get_signed_dihedral_angle
from chemutils.molecule.cip.assignment import (
    TetrahedralCIPAssignment,
    get_tetrahedral_cip_from_ordered_neighbours,
)
from chemutils.molecule.cip.priority import neighbours_ordered_by_cip_priority
from chemutils.molecule.coordinates import get_coordinates
from chemutils.molecule.stereo import check_tetrahedral_coordinates

# Taken from RFAA SI - a tetrahedron on the unit sphere
IDEALIZED_TETRAHEDRON: npt.NDArray[np.float32] = np.array(
    [
        [(8 / 9) ** 0.5, 0, -1 / 3],
        [-((2 / 9) ** 0.5), (2 / 3) ** 0.5, -1 / 3],
        [-((2 / 9) ** 0.5), -((2 / 3) ** 0.5), -1 / 3],
        [0.0, 0.0, 1.0],
    ]
)

NEIGHBOUR_COUNT = Literal[3, 4]


def check_tetrahedral_neighbour_count(*, neighbour_count: NEIGHBOUR_COUNT) -> None:
    if neighbour_count not in [3, 4]:
        raise ValueError(f"Neighbour count must be 3 or 4, not {neighbour_count}")


@lru_cache
def get_idealized_tetrahedron_coords_from_cip_label(
    *, stereo_label: TetrahedralCIPAssignment | None, neighbour_count: NEIGHBOUR_COUNT
) -> npt.NDArray[np.float32]:
    check_tetrahedral_neighbour_count(neighbour_count=neighbour_count)
    ordered_indices = np.arange(4)
    if stereo_label is TetrahedralCIPAssignment.S:
        ordered_indices = ordered_indices[[0, 1, 3, 2]] if neighbour_count == 4 else ordered_indices
    elif stereo_label is TetrahedralCIPAssignment.R:
        ordered_indices = ordered_indices[[0, 1, 3, 2]] if neighbour_count == 3 else ordered_indices
    else:
        raise ValueError("Idealized tetrahedron computation only defined for R/S stereo labels")
    idealized_coordinates: npt.NDArray[np.float32] = IDEALIZED_TETRAHEDRON[ordered_indices]

    # Set first point to center of idealized tetrahedron if we are considering a chiral centers with only 3 explicitly represented neighbours
    if neighbour_count == 3:
        return np.concatenate(
            [np.zeros((1, 3), dtype=idealized_coordinates.dtype), idealized_coordinates[:-1]],
            axis=0,
        )
    return idealized_coordinates


@lru_cache
def get_idealized_tetrahedron_angles_from_cip_label(
    *, stereo_label: TetrahedralCIPAssignment | None, neighbour_count: NEIGHBOUR_COUNT
) -> npt.NDArray[np.float32]:
    """
    Compute the angles between faces of an idealized tetrahedron; either on the unit sphere or the origin as the top point.

    With the sign depending on the chiral label.
    Since only 4 unique arguments are allowed, we can lru_cache this.

    Args:
        stereo_label: oechem.OEAtomStereo_LeftHanded or oechem.OEAtomStereo_RightHanded
        neighbour_count: 3 (implicit hydrogen; amine) or 4

    Returns:
        npt.NDArray[np.float32]: Shape (24,) angles fo all permutations between two faces; sign depends on chiral label.
    """
    check_tetrahedral_neighbour_count(neighbour_count=neighbour_count)
    tetrahedon_coords = get_idealized_tetrahedron_coords_from_cip_label(
        stereo_label=stereo_label, neighbour_count=neighbour_count
    )
    all_perms_angles = []
    for perm in permutations([0, 1, 2, 3]):
        all_perms_angles.append(
            get_signed_dihedral_angle(
                tetrahedon_coords[perm[0]],
                tetrahedon_coords[perm[1]],
                tetrahedon_coords[perm[2]],
                tetrahedon_coords[perm[3]],
            )
        )
    return np.stack(all_perms_angles)


def get_tetrahedron_angles_from_coords_and_cip_label(
    *,
    atom_coordinates: npt.NDArray[np.float32],
    cip_ordered_neighbour_coordinates: npt.NDArray[np.float32],
) -> npt.NDArray[np.float32]:
    """
    Get all pairwise angles (upper triangular portion) between the normals of the planes.

    Through (neighbour_coordinates[0], atom_coordinates, neighbour_coordinates[1]),
    (neighbour_coordinates[1], atom_coordinates, neighbour_coordinates[2]), ..., (neighbour_coordinates[-1], atom_coordinates, neighbour_coordinates[0]).

    Args:
        atom_coordinates: Coordinates of the central atom, as a numpy array of shape `(3, )`
        cip_ordered_neighbour_coordinates: Coordinates of the neighbouring atoms, as a numpy array of shape
            `(3, 3)` or `(4, 3)`. Need to be sorted by CIP priority.

    Returns:
        shape 6 or 12 numpy array of the upper triangular portion of the pairwise angle matrix (depending on the number of chiral neighbours)

    Raises:
        ValueError: Coordinates are wrong shape.
    """
    check_tetrahedral_coordinates(
        atom_coordinates=atom_coordinates, neighbour_coordinates=cip_ordered_neighbour_coordinates
    )
    relative_neighbour_coordinates = cip_ordered_neighbour_coordinates - atom_coordinates[None, :]
    num_neighbours = cip_ordered_neighbour_coordinates.shape[0]

    if num_neighbours == 3:
        # If we have three neighbours, we need to add the center of the tetrahedron
        relative_neighbour_coordinates = np.concatenate(
            [np.zeros((1, 3)), relative_neighbour_coordinates], axis=0
        )
    all_perms_angles = []
    for perm in permutations([0, 1, 2, 3]):
        all_perms_angles.append(
            get_signed_dihedral_angle(
                relative_neighbour_coordinates[perm[0]],
                relative_neighbour_coordinates[perm[1]],
                relative_neighbour_coordinates[perm[2]],
                relative_neighbour_coordinates[perm[3]],
            )
        )
    return np.stack(all_perms_angles).T


def safe_get_stereo_label_from_cip_ordered_neighbours(
    atom: oechem.OEAtomBase,
) -> TetrahedralCIPAssignment | None:
    if atom.IsChiral() and atom.HasStereoSpecified(oechem.OEAtomStereo_Tetrahedral):
        try:
            cip_ordered_neighbours = neighbours_ordered_by_cip_priority(atom)
            return get_tetrahedral_cip_from_ordered_neighbours(
                atom, cip_ordered_neighbours=cip_ordered_neighbours
            )
        except ValueError:
            return None
    return None


def get_idealized_tetrahedron_angles_from_oemol(*, oemol: oechem.OEMol) -> npt.NDArray[np.float32]:
    idealized_angles = []
    for atom in oemol.GetAtoms():
        stereo_label = safe_get_stereo_label_from_cip_ordered_neighbours(atom)
        neigh_count = atom.GetExplicitDegree()
        if neigh_count not in (3, 4) or stereo_label is None:
            idealized_angles.append(np.zeros((24,), dtype=np.float32))
            continue
        idealized_angles.append(
            get_idealized_tetrahedron_angles_from_cip_label(
                stereo_label=stereo_label, neighbour_count=neigh_count
            )
        )
    return np.stack(idealized_angles)


def get_tetrahedron_angles_from_oemol_coords(*, oemol: oechem.OEMol) -> npt.NDArray[np.float32]:
    angles = []
    coordinates = get_coordinates(oemol)
    if coordinates.ndim == 3:
        if coordinates.shape[0] != 1:
            raise ValueError("Expecting single-conformer molecule")
        coordinates = coordinates[0]
    for atom_coords, atom in zip(coordinates, oemol.GetAtoms(), strict=True):
        stereo_label = safe_get_stereo_label_from_cip_ordered_neighbours(atom)
        neigh_count = atom.GetExplicitDegree()
        if neigh_count not in (3, 4) or stereo_label is None:
            angles.append(np.zeros((24,), dtype=np.float32))
            continue
        # Get sorted neighbour idxs
        neighbours = np.array(
            [atom.GetIdx() for atom in neighbours_ordered_by_cip_priority(atom)], dtype=np.int32
        )
        angles.append(
            get_tetrahedron_angles_from_coords_and_cip_label(
                atom_coordinates=atom_coords,
                cip_ordered_neighbour_coordinates=coordinates[neighbours],
            )
        )
    return np.stack(angles)
