from openeye import oechem

# Matches `oequacpac.OERemoveFormalCharge()` behaviour, determined by extensive testing
HANDLED_ATOMIC_NUMBERS = [
    # Carbon
    6,
    # Nitrogen
    7,
    # Oxygen
    8,
    # Fluorine
    9,
    # Silicon
    14,
    # Phosphorus
    15,
    # Sulfur
    16,
    # Chlorine
    17,
    # Arsenic
    33,
    # Selenium
    34,
    # Bromine
    35,
    # Iodine
    53,
]


def remove_atom_formal_charge(atom: oechem.OEAtomBase, /) -> None:
    """
    Adjust the atom's charge by adding and removing hydrogens to achieve 0 charge.

    The built-in `oequacpac.OERemoveFormalCharge()` fails when there is any one atom it cannot handle.
    For example, if there is any ions such as `[Ca+2]`, the function will fail to fix ions that it could
    have fixed if the offending ions were not present.

    Args:
        atom: OpenEye atom to adjust.
    """
    if atom.GetAtomicNum() not in HANDLED_ATOMIC_NUMBERS:
        return
    charge = atom.GetFormalCharge()
    implicit_h_count = atom.GetImplicitHCount()
    explicit_h_count = atom.GetExplicitHCount()
    h_count = implicit_h_count + explicit_h_count
    if charge > 0 and h_count >= charge:
        atom.SetFormalCharge(0)
        total_h_to_remove = charge
        if implicit_h_count >= total_h_to_remove:
            atom.SetImplicitHCount(atom.GetImplicitHCount() - total_h_to_remove)
        else:
            atom.SetImplicitHCount(0)
            total_h_to_remove -= implicit_h_count
            nbrs = [nbr for nbr in atom.GetAtoms() if nbr.GetAtomicNum() == 1]
            for to_remove in nbrs[:total_h_to_remove]:
                atom.GetParent().DeleteAtom(to_remove)
    elif charge < 0:
        atom.SetFormalCharge(0)
        atom.SetImplicitHCount(atom.GetImplicitHCount() + abs(charge))


def remove_formal_charges(mol: oechem.OEMolBase, /) -> None:
    """
    Remove formal charges from atoms.

    This adds and removes hydrogens to change an atom's charge to `0`. It is only applied to the elements
    C, N, O, F, Si, P, S, Cl, As, Se, Br and I.

    This function mimics the behaviour of `oequacpac.OERemoveFormalCharge()`. That means if any nonstandard
    element has a charge, no changes are made.

    Args:
        mol: OpenEye molecule.
    """
    # Copy OpenEye's behaviour of not doing anything if any atoms fail
    if any(
        atom.GetAtomicNum() not in HANDLED_ATOMIC_NUMBERS and atom.GetFormalCharge() != 0
        for atom in mol.GetAtoms()
    ):
        return
    for atom in list(mol.GetAtoms()):
        remove_atom_formal_charge(atom)
