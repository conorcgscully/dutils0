from openeye import oechem


def perceive_atom_properties(mol: oechem.OEGraphMol, /) -> None:
    """Perceive standard properties, including aromaticity, chirality and hybridisation."""
    oechem.OEAssignAromaticFlags(mol, oechem.OEAroModel_OpenEye)

    oechem.OEClearChiralPerception(mol)  # This will clear chiral perception with carbonOnly=True
    oechem.OEPerceiveChiral(
        mol, carbonOnly=False
    )  # If carbonOnly=True, then nitrogen chiral centers will be ignored

    # Need to clear perception, otherwise OE3DToInternalStereo does nothing
    mol.SetPerceived(oechem.OEPerceived_AtomStereo, False)
    mol.SetPerceived(oechem.OEPerceived_BondStereo, False)
    oechem.OE3DToInternalStereo(mol)
    oechem.OEAssignHybridization(mol)
