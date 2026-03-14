from openeye import oechem

from .electrons import get_num_sp_lone_pairs_on_atom


def can_be_hydrogen_acceptor(atom: oechem.OEAtomBase, /) -> bool:
    """
    Could the given atom be a hydrogen acceptor?

    Defined to be any Nitrogen, Oxygen or non-aromatic SP2-hybridized Sulfur.

    Args:
        atom: OpenEye atom.

    Returns:
        True if the given atom could be a hydrogen acceptor if it has lone pairs.
    """
    if atom.GetAtomicNum() == 16:
        return (
            oechem.OEGetHybridization(atom) == oechem.OEHybridization_sp2 and not atom.IsAromatic()
        )
    return atom.GetAtomicNum() in [7, 8]


def is_hydrogen_acceptor(atom: oechem.OEAtomBase, /) -> bool:
    """
    Is the given atom a hydrogen acceptor?

    A hydrogen acceptor is any Nitrogen, Oxygen or non-aromatic SP2-hybridized Sulfur with at least
    one lone pair in an s or sp hybrid orbital.

    Args:
        atom: OpenEye atom.

    Returns:
        True if the given atom could be a hydrogen acceptor.
    """
    return can_be_hydrogen_acceptor(atom) and get_num_sp_lone_pairs_on_atom(atom) > 0


def get_num_hydrogen_acceptors(mol: oechem.OEMolBase, /) -> int:
    """
    Calculate the number of hydrogen accepting atoms.

    Each atom counts once, such that an carbonyl =O counts as a single hydrogen acceptor.

    Args:
        mol: OpenEye molecule.

    Returns:
        Number of atoms that could be hydrogen acceptors.
    """
    return sum(is_hydrogen_acceptor(atom) for atom in mol.GetAtoms())


def get_num_acceptor_lone_pairs(mol: oechem.OEMolBase, /) -> int:
    """
    Calculate the number of lone pairs in s or sp orbitals at in a molecule.

    Each lone pair counts separately, such that a carbonyl =O has two lone pairs.

    Args:
        mol: OpenEye molecule.

    Returns:
        Number of lone pairs that could accept a hydrogen.
    """
    return sum(
        get_num_sp_lone_pairs_on_atom(atom) for atom in mol.GetAtoms() if is_hydrogen_acceptor(atom)
    )
