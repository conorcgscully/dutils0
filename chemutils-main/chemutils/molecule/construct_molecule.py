import numpy as np
import numpy.typing as npt
from openeye import oechem

from .smiles import std_oemol_from_smiles


def construct_molecule(
    *,
    smiles: str,
    coordinates: npt.NDArray[np.float32],
) -> oechem.OEGraphMol:
    """
    Construct molecule from SMILES string and atom coordinates.

    The atom coordinates are assumed to be in the same order as the OEMol created by standardising
    the provided SMILES string using `standardise_oemol`.

    Args:
        smiles: String representation of the molecule in SMILES specification.
        coordinates: Atomic postions as a NumPy array of shape (N_atoms, 3).

    Returns:
        An `OEGraphMol` object with assigned coordinates
    """
    mol = std_oemol_from_smiles(smiles)
    mol.SetCoords(np.asarray(coordinates).reshape(-1))
    mol.SetDimension(3)
    return mol
