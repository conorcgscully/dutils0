from openeye import oechem

from .subset import get_subset
from .types import AtomsSelection


class NoMatchingAtomError(IndexError):
    pass


class MultipleMatchingAtomsError(IndexError):
    pass


def get_atom(oemol: oechem.OEMolBase, /, *, selection: AtomsSelection) -> oechem.OEAtomBase:
    """
    Get a single atom from an OpenEye molecule that matches a selection.

    If the selection matches multiple atoms, an IndexError is raised.

    Args:
        oemol: OpenEye molecule to search.
        selection: Selection to match atoms in the molecule.

    Returns:
        Atom that matches the selection.

    Raises:
        NoMatchingAtomError: If the selection matches no atoms.
        MultipleMatchingAtomsError: If the selection matches multiple atoms.
    """
    subset = get_subset(oemol, selection=selection)
    if subset.NumAtoms() == 0:
        raise NoMatchingAtomError(f"No atom matched selection `{selection}`.")
    if subset.NumAtoms() > 1:
        raise MultipleMatchingAtomsError(f"Multiple atoms matched selection `{selection}`.")
    return next(iter(subset.GetAtoms()))


def get_first_atom(oemol: oechem.OEMolBase, /, *, selection: AtomsSelection) -> oechem.OEAtomBase:
    """
    Get the first atom from an OpenEye molecule that matches a selection.

    Args:
        oemol: OpenEye molecule to search.
        selection: Selection to match atoms in the molecule.

    Returns:
        Atom that matches the selection.

    Raises:
        NoMatchingAtomError: If the selection matches no atoms.
    """
    selection = get_subset(oemol, selection=selection)
    if selection.NumAtoms() == 0:
        raise NoMatchingAtomError(f"No atom matched selection `{selection}`.")
    return next(iter(selection.GetAtoms()))


def get_atoms(oemol: oechem.OEMolBase, /, *, selection: AtomsSelection) -> list[oechem.OEAtomBase]:
    """
    Get a list of atoms from an OpenEye molecule that matches a selection.

    Args:
        oemol: OpenEye molecule to search.
        selection: Selection to match atoms in the molecule.

    Returns:
        List of atoms that match the selection.
    """
    selection = get_subset(oemol, selection=selection)
    return list(selection.GetAtoms())
