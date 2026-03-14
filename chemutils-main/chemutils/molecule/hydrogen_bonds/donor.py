from openeye import oechem


def can_be_hydrogen_donor(atom: oechem.OEAtomBase, /) -> bool:
    """
    Could the given atom be a hydrogen donor?

    Defined to be any Nitrogen, Oxygen, Fluorine or Sulfur.

    Args:
        atom: OpenEye atom.

    Returns:
        True if the given atom could be a hydrogen donor if it has hydrogens.
    """
    return atom.GetAtomicNum() in [7, 8, 9, 16]


def is_hydrogen_donor(atom: oechem.OEAtomBase, /) -> bool:
    """
    Is the given atom a hydrogen donor?

    A hydrogen donor is any Nitrogen, Oxygen, Fluorine or Sulfur with at least one hydrogen.

    Args:
        atom: OpenEye atom.

    Returns:
        True if the given atom could be a hydrogen donor.
    """
    return can_be_hydrogen_donor(atom) and atom.GetTotalHCount() > 0


def get_num_hydrogen_donors(mol: oechem.OEMolBase, /) -> int:
    """
    Calculate the number of hydrogen donating atoms.

    Each atom counts once, such that an amine -NH2 counts as a single hydrogen donor.

    Args:
        mol: OpenEye molecule.

    Returns:
        Number of atoms that could be hydrogen donors.
    """
    return sum(is_hydrogen_donor(atom) for atom in mol.GetAtoms())


def get_num_donor_hydrogens(mol: oechem.OEMolBase, /) -> int:
    """
    Calculate the number of hydrogens that can be donated.

    Each hydrogens counts separately, such that an amine -NH2 has two donor hydrogens.

    Args:
        mol: OpenEye molecule.

    Returns:
        Number of hydrogens that could be donated.
    """
    return sum(atom.GetTotalHCount() for atom in mol.GetAtoms() if is_hydrogen_donor(atom))
