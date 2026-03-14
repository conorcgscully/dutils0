from openeye import oechem, oemolprop, oeomega

from chemutils.props.property import MolecularProperty


def get_num_rotatable_bonds(mol: oechem.OEMolBase, /) -> int:
    """
    Get the number of rotatable bonds in a molecule.

    Rotatable bonds in OpenEye are defined as any single non-ring bond, bounded to nonterminal nonhydrogen
    atoms. In addition, the carbon-carbon triple bond in acetylene is rotatable. It excludes single bonds
    observed in a certain groups, such as sulfonamides, esters, and amidine.
    """
    return oemolprop.OEGetRotatableBondCount(mol)  # type: ignore


def get_num_aromatic_rings(mol: oechem.OEMolBase, /) -> int:
    """Get the number of aromatic rings in a molecule."""
    return oemolprop.OEGetAromaticRingCount(mol)  # type: ignore


def get_largest_ring_size(mol: oechem.OEMolBase, /) -> int:
    """Get the size of the largest ring in a molecule."""
    largest_ring_size = 0
    for atom in mol.GetAtoms():
        ring_size = oechem.OEAtomGetSmallestRingSize(atom)
        largest_ring_size = max(largest_ring_size, ring_size)
    return largest_ring_size


def is_macrocycle(mol: oechem.OEMolBase, /) -> bool:
    """Does the molecule contain a macrocycle of size 10 or more."""
    return oeomega.OEIsMacrocycle(mol)  # type: ignore


NumberRotatableBonds = MolecularProperty(
    id="num_rotatable_bonds", oemol_func=get_num_rotatable_bonds
)

NumberAromaticRings = MolecularProperty(id="num_aromatic_rings", oemol_func=get_num_aromatic_rings)

LargestRingSize = MolecularProperty(id="largest_ring_size", oemol_func=get_largest_ring_size)

IsMacrocycle = MolecularProperty(id="is_macrocycle", oemol_func=is_macrocycle)
