from openeye import oechem


def get_alpha_carbon_mapping_by_residue_number(
    src: oechem.OEMol, dest: oechem.OEMol, /
) -> dict[oechem.OEAtomBase, oechem.OEAtomBase]:
    """Assumes that src and dest proteins have the same residue numbers."""
    src_residues = {}
    for atom in src.GetAtoms():
        if atom.GetAtomicNum() == 6 and atom.GetName().strip() == "CA":
            residue_num = oechem.OEAtomGetResidue(atom).GetResidueNumber()
            if residue_num in src_residues:
                raise ValueError(
                    f"Multiple alpha carbons with the same residue number: {residue_num} in `src` protein"
                )
            src_residues[residue_num] = atom

    dest_residues = {}
    for atom in dest.GetAtoms():
        if atom.GetAtomicNum() == 6 and atom.GetName().strip() == "CA":
            residue_num = oechem.OEAtomGetResidue(atom).GetResidueNumber()
            if residue_num in dest_residues:
                raise ValueError(
                    f"Multiple alpha carbons with the same residue number: {residue_num} in `dest` protein"
                )
            dest_residues[residue_num] = atom

    alpha_carbon_mapping = {}
    for residue_num in src_residues:
        if residue_num in dest_residues:
            alpha_carbon_mapping[src_residues[residue_num]] = dest_residues[residue_num]

    return alpha_carbon_mapping
