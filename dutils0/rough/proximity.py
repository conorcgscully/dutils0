def flag_protein_residues_by_ligand_proximity(
    protein_mol: Mol,
    ligand_mol: Mol,
    cutoff: float = 4.5,
) -> Dict[ResidueKey, int]:
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
        Keys are residue tuples:
            (chain_id, residue_number, insertion_code, residue_name)
        Values are:
            1 if any atom in that residue is within cutoff of any ligand atom
            0 otherwise
    """
    if protein_mol is None:
        raise ValueError("protein_mol is None")
    if ligand_mol is None:
        raise ValueError("ligand_mol is None")
    if cutoff <= 0:
        raise ValueError("cutoff must be > 0")

    if protein_mol.GetNumConformers() != 1:
        raise ValueError(
            f"protein_mol must have exactly 1 conformer, got {protein_mol.GetNumConformers()}"
        )
    if ligand_mol.GetNumConformers() != 1:
        raise ValueError(
            f"ligand_mol must have exactly 1 conformer, got {ligand_mol.GetNumConformers()}"
        )

    prot_conf = protein_mol.GetConformer()
    lig_conf = ligand_mol.GetConformer()

    # Ligand coordinates: shape (n_lig_atoms, 3)
    lig_coords = np.array(
        [
            [
                lig_conf.GetAtomPosition(i).x,
                lig_conf.GetAtomPosition(i).y,
                lig_conf.GetAtomPosition(i).z,
            ]
            for i in range(ligand_mol.GetNumAtoms())
        ],
        dtype=float,
    )

    if lig_coords.size == 0:
        raise ValueError("ligand_mol has no atoms")

    cutoff_sq = cutoff * cutoff

    # Collect protein atoms by residue
    residue_to_atom_indices: Dict[ResidueKey, list[int]] = {}

    for atom in protein_mol.GetAtoms():
        pdb_info = atom.GetPDBResidueInfo()
        if pdb_info is None:
            continue

        residue_key: ResidueKey = (
            pdb_info.GetChainId().strip(),
            int(pdb_info.GetResidueNumber()),
            pdb_info.GetInsertionCode().strip(),
            pdb_info.GetResidueName().strip(),
        )
        residue_to_atom_indices.setdefault(residue_key, []).append(atom.GetIdx())

    if not residue_to_atom_indices:
        raise ValueError(
            "No protein atoms with AtomPDBResidueInfo were found in protein_mol"
        )

    flags: Dict[ResidueKey, int] = {}

    for residue_key, atom_indices in residue_to_atom_indices.items():
        flagged = 0

        for atom_idx in atom_indices:
            p = prot_conf.GetAtomPosition(atom_idx)
            p_xyz = np.array([p.x, p.y, p.z], dtype=float)

            # squared distances from this protein atom to all ligand atoms
            d2 = np.sum((lig_coords - p_xyz) ** 2, axis=1)

            if np.any(d2 <= cutoff_sq):
                flagged = 1
                break

        flags[residue_key] = flagged

    return {
        f"{resname}{resnum}.{chain}{inscode if inscode else ''}": flag
        for (chain, resnum, inscode, resname), flag in flags.items()
    }