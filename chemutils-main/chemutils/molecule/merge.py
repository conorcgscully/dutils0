from openeye import oechem


def merge_molecules(
    *oemols: oechem.OEGraphMol,
    atom_mapping: dict[oechem.OEAtomBase, oechem.OEAtomBase] | None = None,
    bond_mapping: dict[oechem.OEBondBase, oechem.OEBondBase] | None = None,
) -> oechem.OEGraphMol:
    mol = oechem.OEGraphMol()
    for oemol in oemols:
        merge_atom_mapping, merge_bond_mapping = oechem.OEAddMols(mol, oemol)
        if atom_mapping is not None:
            for src_atom in oemol.GetAtoms():
                atom_mapping[src_atom] = merge_atom_mapping[src_atom.GetIdx()]
        if bond_mapping is not None:
            for src_bond in oemol.GetBonds():
                bond_mapping[src_bond] = merge_bond_mapping[src_bond.GetIdx()]
    return mol
