from openeye import oechem


def canonicalize_molecule(mol: oechem.OEMolBase, /) -> None:
    """
    Canonicalize the order of the atoms and bonds in the molecule.

    This orders the atom and bond objects in the molecule such that they will be in the same order,
    independent of their initial order.

    No actual modifications of the molecule will occur, it purely reindexes the atoms and bonds.

    Args:
        mol: Molecule to canonicalize.
    """
    # Order the atoms according to the canonical isomeric SMILES ordering
    mol.OrderAtoms(
        oechem.OEAtomVector(list(oechem.OEGetSmiStringOrder(mol, oechem.OESMILESFlag_ISOMERIC)))
    )
    # Order the bonds in the same order
    oechem.OECanonicalOrderBonds(mol)
    # Reindex atoms and bonds again
    mol.Sweep()
    for bond in mol.GetBonds():
        if bond.GetBgnIdx() > bond.GetEndIdx():
            bond.SwapEnds()
