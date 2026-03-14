from collections.abc import Callable

import numpy as np
import numpy.typing as npt
from openeye import oechem


def get_atom_mask(
    oemol: oechem.OEMolBase,
    predicate: oechem.OEUnaryAtomPred | Callable[[oechem.OEAtomBase], bool],
    /,
) -> npt.NDArray[np.bool_]:
    """
    Get a boolean mask over the atoms in a molecule, given a predicate.

    The order of atoms is the order of atoms of the `oemol`, as returned by `GetAtoms`.

    Args:
        oemol: OpenEye molecule.
        predicate: Either an OpenEye `OEUnaryAtomPred`, or any function which takes in an `OEAtomBase`
            and returns a boolean.

    Returns:
        NumPy boolean array of shape `(n_atoms, )`.
    """
    return np.array([predicate(atom) for atom in oemol.GetAtoms()], dtype=np.bool_)
