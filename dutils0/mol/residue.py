"""Residue-level helpers for RDKit molecules."""

from __future__ import annotations

from typing import Dict, List, Tuple

from rdkit import Chem

ResidueKey = Tuple[str, int, str, str]


def get_residue(mol: Chem.Mol) -> Dict[ResidueKey, List[int]]:
    """
    Group atom indices by residue.

    Parameters
    ----------
    mol
        RDKit Mol with AtomPDBResidueInfo on atoms.

    Returns
    -------
    dict
        Keys are (chain_id, residue_number, insertion_code, residue_name).
        Values are lists of atom indices belonging to that residue.

    Raises
    ------
    ValueError
        If no atoms with AtomPDBResidueInfo are found.
    """
    if mol is None:
        raise ValueError("mol is None")

    residue_to_atom_indices: Dict[ResidueKey, List[int]] = {}

    for atom in mol.GetAtoms():
        pdb_info = atom.GetPDBResidueInfo()
        if pdb_info is None:
            continue
        key: ResidueKey = (
            pdb_info.GetChainId().strip(),
            int(pdb_info.GetResidueNumber()),
            pdb_info.GetInsertionCode().strip(),
            pdb_info.GetResidueName().strip(),
        )
        residue_to_atom_indices.setdefault(key, []).append(atom.GetIdx())

    if not residue_to_atom_indices:
        raise ValueError("No atoms with AtomPDBResidueInfo were found in mol")

    return residue_to_atom_indices
