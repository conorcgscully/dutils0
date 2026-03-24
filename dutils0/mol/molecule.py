"""Molecule coordinate helpers."""

from __future__ import annotations

import numpy as np
from rdkit import Chem


def get_molecule(mol: Chem.Mol) -> np.ndarray:
    """
    Return atomic coordinates from an RDKit molecule as a NumPy array.

    Parameters
    ----------
    mol
        RDKit Mol. Must have exactly one conformer.

    Returns
    -------
    np.ndarray
        Shape ``(n_atoms, 3)``, dtype float.

    Raises
    ------
    ValueError
        If mol is None, does not have exactly one conformer, or has no atoms.
    """
    if mol is None:
        raise ValueError("mol is None")
    if mol.GetNumConformers() != 1:
        raise ValueError(
            f"mol must have exactly 1 conformer, got {mol.GetNumConformers()}"
        )
    if mol.GetNumAtoms() == 0:
        raise ValueError("mol has no atoms")

    conf = mol.GetConformer()
    return np.array(conf.GetPositions(), dtype=float)
