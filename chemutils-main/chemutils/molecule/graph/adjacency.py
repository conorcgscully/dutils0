import numpy as np
import numpy.typing as npt
from openeye import oechem


def get_adjacency_matrix(oemol: oechem.OEMolBase, /) -> npt.NDArray[np.bool_]:
    """
    Get the adjacency matrix for an OpenEye molecule.

    The order of atoms is the order of atoms of the `oemol`, as returned by `GetAtoms`.

    Args:
        oemol: OpenEye molecule.

    Returns:
        `NxN` boolean NumPy matrix, where True at `[i, j]` indicates a bond between atoms with
        indices `i` and `j`.
    """
    return np.array(
        [
            [oemol.GetBond(atom1, atom2) is not None for atom1 in oemol.GetAtoms()]
            for atom2 in oemol.GetAtoms()
        ],
        dtype=np.bool_,
    )
