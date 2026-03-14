from openeye import oechem

from chemutils.molecule.stereo import reflect_stereocenter

from .smiles import smiles_from_oemol


def molecule_with_flipped_stereocenters(
    *, oemol: oechem.OEMolBase, atom_idxs: list[int]
) -> oechem.OEMolBase:
    """
    Return a copy of a molecule with the specified stereocenters flipped.

    Args:
        oemol: OpenEye molecule to flip;
        atom_idxs: Indices of the atoms who's tetrahedral stereochemistry should be flipped.

    Returns:
        New OpenEye molecule with the given tetrahedral centers flipped.
    """
    oemol = oemol.CreateCopy()
    for atom_idx in atom_idxs:
        reflect_stereocenter(oemol.GetAtom(oechem.OEHasAtomIdx(atom_idx)))

    # Reperceive MDL properties based on new stereochemistry
    oechem.OEMDLPerceiveParity(oemol)
    oechem.OEMDLPerceiveBondStereo(oemol)
    return oemol


def convert_enhanced_stereo_to_unspecified(oemol: oechem.OEMolBase, /) -> oechem.OEMolBase:
    """
    Convert any AND/OR enhanced stereochemistry to unspecified stereochemistry and remove enhanced stereo information.

    This allows a CXSMILES to be converted into a standard SMILES, which may be required for certain applications.

    If the provided molecule does not have any enhanced stereochemistry information, this is a no-op.

    Args:
        oemol: Molecule with potential enhanced stereochemistry information.

    Returns:
        None
    """
    atom_indices: set[int] = set()
    for group in oemol.GetGroups(oechem.OEIsMDLStereoGroup()):
        match group.GetGroupType():
            case oechem.OEGroupType_MDLAndStereo | oechem.OEGroupType_MDLOrStereo:
                atom_indices |= {
                    atom.GetIdx()
                    for atom in group.GetAtoms()
                    if atom.HasStereoSpecified(oechem.OEAtomStereo_Tetrahedral)
                }

    for atom_idx in atom_indices:
        atom = oemol.GetAtom(oechem.OEHasAtomIdx(atom_idx))
        if atom.HasStereoSpecified(oechem.OEAtomStereo_Tetrahedral):
            atom.SetStereo([], oechem.OEAtomStereo_Tetrahedral, oechem.OEAtomStereo_Undefined)

    for group in list(oemol.GetGroups(oechem.OEIsMDLStereoGroup())):
        oemol.DeleteGroup(group)

    return oemol


def get_stereoisomer_mixtures_from_oemol(
    oemol: oechem.OEMolBase,
    /,
) -> list[tuple[oechem.OEMolBase, ...]]:
    """
    Get the stereoisomers in mixtures described by enhanced stereochemistry.

    Args:
        oemol: OpenEye molecule containing enhanced stereochemistry information.

    Returns:
        List of lists of OpenEye molecules, where the outer list is the unique mixtures
        and the inner lists are the stereoisomers in each mixture.

    Raises:
        ValueError: Nonchiral atoms marked with enhanced stereo labels.
    """
    # Ensure we have chirality perceived
    oechem.OEPerceiveChiral(oemol)

    # Identify and- and or- stereo groups
    and_stereo_atoms = []
    or_stereo_atoms = []
    for group in oemol.GetGroups(oechem.OEIsMDLStereoGroup()):
        match group.GetGroupType():
            case oechem.OEGroupType_MDLAndStereo:
                atoms = [
                    atom.GetIdx()
                    for atom in group.GetAtoms()
                    if atom.HasStereoSpecified(oechem.OEAtomStereo_Tetrahedral)
                ]
                and_stereo_atoms.append(atoms)
            case oechem.OEGroupType_MDLOrStereo:
                atoms = [
                    atom.GetIdx()
                    for atom in group.GetAtoms()
                    if atom.HasStereoSpecified(oechem.OEAtomStereo_Tetrahedral)
                ]
                or_stereo_atoms.append([atom.GetIdx() for atom in group.GetAtoms()])

    # Generate unique mixture compositions, based on or groups
    mixture_start_points: list[oechem.OEMolBase] = [oemol]
    for or_group in or_stereo_atoms:
        for mol in list(mixture_start_points):
            mixture_start_points.append(
                molecule_with_flipped_stereocenters(oemol=mol, atom_idxs=or_group)
            )

    mixtures: list[list[oechem.OEMolBase]] = []
    for mixture_start_point in mixture_start_points:
        mixture = [mixture_start_point]
        for and_group in and_stereo_atoms:
            for mol in list(mixture):
                mixture.append(molecule_with_flipped_stereocenters(oemol=mol, atom_idxs=and_group))
        mixtures.append(mixture)

    unique_mixtures: list[tuple[oechem.OEMolBase, ...]] = []

    # Remove duplicates
    for mixture in mixtures:
        smiles_in_mixture = {}
        for mol in list(mixture):
            smiles = oechem.OEMolToSmiles(mol)
            if smiles not in smiles_in_mixture:
                smiles_in_mixture[smiles] = mol
        unique_mixtures.append(
            tuple(smiles_in_mixture[smiles] for smiles in sorted(smiles_in_mixture.keys()))
        )

    return sorted(
        unique_mixtures, key=lambda mixture: tuple(oechem.OEMolToSmiles(mol) for mol in mixture)
    )


def remove_redundant_enhanced_stereo(oemol: oechem.OEMolBase, /) -> None:
    """
    Remove redundant stereo labels from an OpenEye molecule.

    When combined with a following call to get_canonical_stereoisomer, this gives a
    molecule whose CXSMILES will be reduced to a simple canonical form.

    Args:
        oemol: OpenEye molecule.
    """
    for group in list(oemol.GetGroups(oechem.OEIsMDLStereoGroup())):
        match group.GetGroupType():
            case oechem.OEGroupType_MDLAbsStereo:
                # Absolute stereo groups are already implicit
                oemol.DeleteGroup(group)
            case oechem.OEGroupType_MDLAndStereo | oechem.OEGroupType_MDLOrStereo:
                atoms = [
                    atom.GetIdx()
                    for atom in group.GetAtoms()
                    if atom.HasStereoSpecified(oechem.OEAtomStereo_Tetrahedral)
                ]
                epimer = molecule_with_flipped_stereocenters(oemol=oemol, atom_idxs=atoms)
                # If flipping this center gives the sames SMILES, then we have a redundant label
                if smiles_from_oemol(oemol) == smiles_from_oemol(epimer):
                    oemol.DeleteGroup(group)


def standardise_stereo_group_labels(oemol: oechem.OEMolBase, /) -> None:
    """
    Standardise the stereo group order in an OpenEye molecule.

    OpenEye automatically renumbers the groups from 1 when writing to CXSMILES. However, the order
    of the groups depends on the order they were in the input. This reorders them by looking at the
    atoms in each group, and where they would fall in a canonical isomeric SMILES. It reorders the
    groups such that groups with lower atom indices appear before groups with higher atom indices.

    Args:
        oemol: OpenEye molecule with enhanced stereochemistry.
    """
    canonical_atom_order = list(oechem.OEGetSmiStringOrder(oemol, oechem.OESMILESFlag_ISOMERIC))
    groups = list(oemol.GetGroups(oechem.OEIsMDLStereoGroup()))
    # Order groups by the index of their atoms in the canonical SMILES
    groups = sorted(
        groups,
        key=lambda group: tuple(canonical_atom_order.index(atom) for atom in group.GetAtoms()),
    )
    # Delete and recreate groups to reorder them
    for group in groups:
        oemol.DeleteGroup(group)
        oemol._NewGroup_Atoms(group.GetGroupType(), list(group.GetAtoms()))


def get_canonical_stereoisomer(oemol: oechem.OEMolBase, /) -> oechem.OEMolBase:
    """
    Get the canonical stereoisomer with enhanced stereochemistry for a give molecule.

    This picks the first stereoisomer alphabetically by isomeric canonical smiles, and renormalises
    the enhanced stereo groups to be ordered by the atoms they refer to.

    Args:
        oemol: OpenEye molecule with enhanced stereochemistry.

    Returns:
        New OpenEye molecule which will be a stereoisomer with enhanced stereochemistry, that represents
        the same mixture as the original molecule.
    """
    oemol = get_stereoisomer_mixtures_from_oemol(oemol)[0][0]
    standardise_stereo_group_labels(oemol)
    return oemol
