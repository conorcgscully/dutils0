from typing import Literal

from openeye import oechem


def remove_atom(
    atom: oechem.OEAtomBase,
    *,
    bond_handling: Literal["dangling", "hydrogens", "rgroups"],
    adjust_rings: bool = True,
    adjust_aromaticity: bool = True,
) -> None:
    """
    Remove an OpenEye atom from a molecule.

    By default, this recalculates ring membership and aromaticity if the atom is in a ring or is aromatic.

    Args:
        atom: Atom to remove.
        bond_handling: When breaking bonds, whether to just remove them (`dangling`), add hydrogens to the
            neighbouring atoms (`hydrogens`), or add R-groups to the neighbouring atoms (`rgroups`).
        adjust_rings: Whether to recalculate ring membership.
        adjust_aromaticity: Whether to recalculate aromaticity.

    Raises:
        ValueError: If the atom has explicit hydrogens
        ValueError: If the atom has neighbours with explicit hydrogens.
    """
    # Aromaticity perception implicitly adjusts rings, so we can't adjust rings without adjusting aromaticity.
    if not adjust_rings and adjust_aromaticity:
        raise ValueError("Cannot set `adjust_rings` to `False` if `adjust_aromaticity` is `True`.")

    mol: oechem.OEMolBase = atom.GetParent()

    if atom.GetExplicitHCount() > 0:
        raise ValueError("Cannot remove an atom with explicit hydrogens.")

    is_aromatic = atom.IsAromatic()
    is_in_ring = atom.IsInRing()
    nbrs_with_bond_orders: dict[oechem.OEAtomBase, int] = {
        bond.GetNbr(atom): bond.GetOrder() for bond in atom.GetBonds()
    }
    if bond_handling == "hydrogens" and any(
        nbr.GetExplicitHCount() > 0 for nbr in nbrs_with_bond_orders
    ):
        raise ValueError("Cannot remove an atom with neighbours that have explicit hydrogens.")
    mol.DeleteAtom(atom)

    if bond_handling == "hydrogens":
        for nbr, bond_order in nbrs_with_bond_orders.items():
            nbr.SetImplicitHCount(nbr.GetImplicitHCount() + bond_order)
    elif bond_handling == "rgroups":
        rgroup_idx = 1
        for nbr, bond_order in nbrs_with_bond_orders.items():
            rgroup = mol.NewAtom(0)
            rgroup.SetMapIdx(rgroup_idx)
            mol.NewBond(nbr, rgroup, bond_order)
            rgroup_idx += 1

    if adjust_rings and is_in_ring:
        oechem.OEFindRingAtomsAndBonds(mol)
    if adjust_aromaticity and is_aromatic:
        oechem.OEAssignAromaticFlags(mol)


def remove_bond(
    bond: oechem.OEBondBase,
    *,
    bond_handling: Literal["dangling", "hydrogens", "electrons", "rgroups"],
    adjust_rings: bool = True,
    adjust_aromaticity: bool = True,
) -> None:
    """
    Remove an OpenEye bond from a molecule.

    By default, this recalculates ring membership and aromaticity if the atom is in a ring or is aromatic.

    Args:
        bond: Bond to remove.
        bond_handling: Whether to leave the connected atoms as if (`dangling`), add hydrogens to them (`hydrogens`),
        decrease their charge (`electrons`) or add R-groups (`rgroups`).
        adjust_rings: Whether to recalculate ring membership.
        adjust_aromaticity: Whether to recalculate aromaticity.

    Raises:
        ValueError: If the bond has an adjacent atom with explicit hydrogens.
    """
    # Aromaticity perception implicitly adjusts rings, so we can't adjust rings without adjusting aromaticity.
    if not adjust_rings and adjust_aromaticity:
        raise ValueError("Cannot set `adjust_rings` to `False` if `adjust_aromaticity` is `True`.")

    mol = bond.GetParent()

    if bond.GetBgn().GetExplicitHCount() > 0:
        raise ValueError("Cannot remove a bond adjacent to atom with explicit hydrogens.")
    if bond.GetEnd().GetExplicitHCount() > 0:
        raise ValueError("Cannot remove a bond adjacent to atom with explicit hydrogens.")

    is_aromatic = bond.IsAromatic()
    is_in_ring = bond.IsInRing()

    bond_order = bond.GetOrder()
    nbrs: list[oechem.OEAtomBase] = [bond.GetBgn(), bond.GetEnd()]

    mol.DeleteBond(bond)

    if bond_handling == "hydrogens":
        for nbr in nbrs:
            nbr.SetImplicitHCount(nbr.GetImplicitHCount() + bond_order)
    elif bond_handling == "rgroups":
        rgroup_idx = 1
        for nbr in nbrs:
            rgroup = mol.NewAtom(0)
            rgroup.SetMapIdx(rgroup_idx)
            mol.NewBond(nbr, rgroup, bond_order)
            rgroup_idx += 1
    elif bond_handling == "electrons":
        for nbr in nbrs:
            nbr.SetFormalCharge(nbr.GetFormalCharge() - bond_order)

    if adjust_rings and is_in_ring:
        oechem.OEFindRingAtomsAndBonds(mol)
    if adjust_aromaticity and is_aromatic:
        oechem.OEAssignAromaticFlags(mol)
