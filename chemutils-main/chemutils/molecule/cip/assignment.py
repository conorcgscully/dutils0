from enum import Enum

from openeye import oechem


class TetrahedralCIPAssignment(str, Enum):
    R = "R"
    S = "S"


def get_tetrahedral_cip_from_ordered_neighbours(
    atom: oechem.OEAtomBase, *, cip_ordered_neighbours: list[oechem.OEAtomBase]
) -> TetrahedralCIPAssignment | None:
    """
    Get the CIP assignment for a given tetrahedral stereocenter, given neighbouring atoms ordered by CIP priority.

    Args:
        atom: OpenEye atom to calculate CIP assignment for.
        cip_ordered_neighbours: Neighbours ordered from low priority to high priority, based on the CIP classification
            rules.

    Returns:
        'R' or 'S' if an assignment is made, or None if the atom is not a stereocenter or
        stereochemistry is not specified for that atom.
    """
    # CIP stereochemistry is defined based on clockwise/anticlockwise arrangement of the top three neighbours
    # by priority, starting at the highest priority.
    relevant_neighbours = list(reversed(cip_ordered_neighbours))[:3]
    stereo_tetra = atom.GetStereo(relevant_neighbours, oechem.OEAtomStereo_Tetrahedral)
    if stereo_tetra == oechem.OEAtomStereo_Right:
        return TetrahedralCIPAssignment.R
    if stereo_tetra == oechem.OEAtomStereo_Left:
        return TetrahedralCIPAssignment.S
    return None
