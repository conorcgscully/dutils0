from openeye import oechem


def make_hydrogens_implicit(mol: oechem.OEMolBase, /, *, remove_isotopic_hydrogens: bool) -> None:
    """
    Convert all explicit hydrogens to implicit hydrogens.

    Isotopic hydrogen such as deutrium are a special case, as you will lose
    isotope information by making them implicit. Set `remove_isotopic_hydrogens`
    if you are prepared to lose this information.

    Args:
        mol: Molecule to convert hydrogens for.
        remove_isotopic_hydrogens: Should isotopic hydrogens be removed.

    As all hydrogens are converted to being just a property of their neighbouring heavy atom,
    ksotopic information about hydrogens will be lost.
    """
    oechem.OESuppressHydrogens(mol, False, False, not remove_isotopic_hydrogens)


class ProtonationError(ValueError):
    pass


class DeprotonationError(ValueError):
    pass


def protonate(atom: oechem.OEAtomBase, /) -> None:
    """
    Protonate an atom by adding a hydrogen to it and increasing its formal charge by one.

    If the molecule has explicit hydrogens, the hydrogen will be added as an explicit hydrogen. If
    the molecule also has coordinates, then the hydrogen atom will be placed at a sensible position.

    Args:
        atom: Atom to protonate.

    Raises:
        ProtonationError: If the atom is a hydrogen atom with atomic number 1.
        ProtonationError: If the atom has both implicit and explicit hydrogens.
    """
    if atom.GetAtomicNum() == 1:
        raise ProtonationError("Cannot protonate a hydrogen atom.")

    mol: oechem.OEMolBase = atom.GetParent()
    num_implicits = sum(atom.GetImplicitHCount() for atom in mol.GetAtoms(oechem.OEIsHeavy()))
    num_explicits = sum(atom.GetExplicitHCount() for atom in mol.GetAtoms(oechem.OEIsHeavy()))
    if num_implicits > 0 and num_explicits > 0:
        raise ProtonationError(
            "Cannot protonate an atom with both implicit and explicit hydrogens."
        )
    has_coords = mol.GetDimension() == 3
    use_explicit = num_explicits > 0 or (num_implicits == 0 and has_coords)

    atom.SetImplicitHCount(atom.GetImplicitHCount() + 1)
    atom.SetFormalCharge(atom.GetFormalCharge() + 1)
    if use_explicit:
        oechem.OEAddExplicitHydrogens(mol, False, has_coords)


def deprotonate(atom: oechem.OEAtomBase, /) -> None:
    """
    Deprotonate an atom by removing a hydrogen from it and decreasing its formal charge by one.

    If a hydrogen atom is passed, then this is the hydrogen that will be removed. Otherwise, a hydrogen
    will be removed from the given atom.

    Args:
        atom: Atom to deprotonate.

    Raises:
        DeprotonationError: If the atom is a hydrogen atom with two neighbours.
        DeprotonationError: If the atom is not a hydrogen and has no hydrogens.
        DeprotonationError: If the atom has more than one explicit hydrogens.
        DeprotonationError: If the atom has both implicit and explicit hydrogens.
    """
    if atom.GetAtomicNum() == 1:
        neighs = list(atom.GetAtoms())
        if len(neighs) > 1:
            raise DeprotonationError(
                "Cannot deprotonate by removing a hydrogen atom with two neighbours."
            )
        if len(neighs) == 0:
            raise DeprotonationError("Cannot deprotonate a hydrogen atom with no neighbours.")
        neigh: oechem.OEAtomBase = neighs[0]
        mol: oechem.OEMolBase = atom.GetParent()
        mol.DeleteAtom(atom)
        neigh.SetFormalCharge(neigh.GetFormalCharge() - 1)
    else:
        num_implicit = atom.GetImplicitHCount()
        num_explicit = atom.GetExplicitHCount()
        match num_implicit, num_explicit:
            case 0, 0:
                raise DeprotonationError("Cannot deprotonate an atom with no hydrogens.")
            case 0, 1:
                explicit_hydrogen = next(iter(atom.GetAtoms(oechem.OEIsHydrogen())))
                deprotonate(explicit_hydrogen)
                return
            case 0, _:
                raise DeprotonationError(
                    "Cannot deprotonate an atom with more than one explicit hydrogens."
                )
            case _, 0:
                atom.SetImplicitHCount(atom.GetImplicitHCount() - 1)
                atom.SetFormalCharge(atom.GetFormalCharge() - 1)
                return
            case _, _:
                raise DeprotonationError(
                    "Cannot deprotonate an atom with both implicit and explicit hydrogens."
                )
