"""Protein-ligand proximity helpers."""

from __future__ import annotations

import numpy as np
from rdkit import Chem

from .molecule import get_molecule
from .residue import get_residue, ResidueKey


def flag_protein_residues_by_ligand_proximity(
    protein_mol: Chem.Mol,
    ligand_mol: Chem.Mol,
    cutoff: float = 4.5,
) -> dict:
    """
    Flag each protein residue as 1/0 depending on whether any protein atom in that
    residue lies within `cutoff` Å of any ligand atom.

    Parameters
    ----------
    protein_mol
        RDKit Mol for the protein. Must have exactly one conformer and protein atoms
        should carry AtomPDBResidueInfo.
    ligand_mol
        RDKit Mol for the ligand. Must have exactly one conformer.
    cutoff
        Distance cutoff in Å.

    Returns
    -------
    dict
        Keys are strings of the form ``{resname}{resnum}.{chain}{inscode}``.
        Values are 1 if any atom in that residue is within cutoff of any ligand atom,
        0 otherwise.
    """
    if cutoff <= 0:
        raise ValueError("cutoff must be > 0")

    lig_coords = get_molecule(ligand_mol)
    prot_coords = get_molecule(protein_mol)
    residue_to_atom_indices = get_residue(protein_mol)

    cutoff_sq = cutoff * cutoff
    flags: dict[ResidueKey, int] = {}

    for residue_key, atom_indices in residue_to_atom_indices.items():
        flagged = 0
        for atom_idx in atom_indices:
            d2 = np.sum((lig_coords - prot_coords[atom_idx]) ** 2, axis=1)
            if np.any(d2 <= cutoff_sq):
                flagged = 1
                break
        flags[residue_key] = flagged

    return {
        f"{resname}{resnum}.{chain}{inscode if inscode else ''}": flag
        for (chain, resnum, inscode, resname), flag in flags.items()
    }
