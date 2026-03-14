from typing import Literal

from openeye import oechem

from .types import (
    AtomsSelection,
    _selection_as_atombondset,
    _selection_as_atombondset_or_unaryatompred,
)


def add_internal_bonds_to_selection(set: oechem.OEAtomBondSet, /) -> oechem.OEAtomBondSet:
    new_set = set.CreateCopy()
    for atom in set.GetAtoms():
        for bond in atom.GetBonds():
            if set.HasAtom(bond.GetNbr(atom)):
                new_set.AddBond(bond)
    return new_set


def get_subset(
    mol: oechem.OEGraphMol, /, *, selection: AtomsSelection, include_internal_bonds: bool = True
) -> oechem.OEAtomBondSet:
    """
    Get an `oechem.OEAtomBondSet` representing the subset of the molecule that matches the given selection.

    Args:
        mol: OpenEye molecule for which to get the subset.
        selection: Query to match atoms in the molecule. Can be one of the following:
            * An `oechem.OEUnaryAtomPred`, such as `oechem.OEIsHeavy()` or `oechem.OEHasAtomicNum(4)`
            * A function or lambda function that takes an atom and returns a boolean indicating whether the atom should be included in the subset.
            * An `oechem.OEAtomBondSet` containing the atoms and bonds to include in the subset.
            * A list of `oechem.OEAtomBase`.
        include_internal_bonds: Whether to automatically include bonds in the subset if both connected atoms are in the subset. Defaults to `True`.

    Returns:
        An `oechem.OEAtomBondSet` contains the lists of atoms and bonds chosen for this subset.
    """
    selection = _selection_as_atombondset(selection, mol=mol)
    if include_internal_bonds:
        selection = add_internal_bonds_to_selection(selection)
    return selection


class SubsetError(ValueError):
    pass


def subset_molecule(
    mol: oechem.OEGraphMol,
    /,
    *,
    selection: AtomsSelection,
    bond_handling: Literal["dangling", "hydrogens", "rgroups"],
    include_internal_bonds: bool = True,
    atom_mapping: dict[oechem.OEAtomBase, oechem.OEAtomBase] | None = None,
    bond_mapping: dict[oechem.OEBondBase, oechem.OEBondBase] | None = None,
) -> oechem.OEGraphMol:
    """
    Get an new OpenEye molecule of a specific subset of the provided molecule that matches the given selection.

    Args:
        mol: OpenEye molecule for which to get the subset.
        selection: Selection to match atoms in the molecule. Can be one of the following:
            * An `oechem.OEUnaryAtomPred`, such as `oechem.OEIsHeavy()` or `oechem.OEHasAtomicNum(4)`
            * A function or lambda function that takes an atom and returns a boolean indicating whether the atom should be included in the subset.
            * An `oechem.OEAtomBondSet` containing the atoms and bonds to include in the subset.
            * A list of `oechem.OEAtomBase`.
        bond_handling: How to handle bonds broken between atoms in the subset and atoms outside the subset. Can be one of the following:
            * `dangling`: Break the bond and leave the atoms in the with an unsatisfied valence.
            * `hydrogens`: Insert implicit hydrogens when bonds are broken, to maintain the same valence.
            * `rgroups`: Insert R-Groups when bonds are broken, to maintain the same valence.
        include_internal_bonds: Whether to automatically include bonds in the subset if both connected atoms are in the subset. Defaults to `True`.
        atom_mapping: If provided, populate a dictionary mapping atoms in the original molecule to atoms in the subset molecule.
        bond_mapping: If provided, populate a dictionary mapping bonds in the original molecule to bonds in the subset molecule.

    Returns:
        A new OpenEye molecule containing the subset of the original molecule that matches the query.
    """
    selection = _selection_as_atombondset_or_unaryatompred(selection, mol=mol)

    # Subset behaviour is to include internal bonds for atom predicates, but not for atom-bond sets
    # Make this behaviour consistent based on include_internal_bonds
    if isinstance(selection, oechem.OEAtomBondSet) and include_internal_bonds:
        selection = add_internal_bonds_to_selection(selection)
    elif not isinstance(selection, oechem.OEAtomBondSet) and not include_internal_bonds:
        selection = get_subset(mol, selection=selection, include_internal_bonds=False)

    new_mol = oechem.OEGraphMol()
    adjust_h_counts = bond_handling == "hydrogens"
    r_group = bond_handling == "rgroups"

    atommap = oechem.OEAtomArray(mol.GetMaxAtomIdx())
    bondmap = oechem.OEBondArray(mol.GetMaxBondIdx())
    if not oechem.OESubsetMol(new_mol, mol, selection, adjust_h_counts, r_group, atommap, bondmap):
        raise SubsetError("Failed to create a subset of the molecule.")
    if atom_mapping is not None:
        atom_mapping.update(
            {
                mol.GetAtom(oechem.OEHasAtomIdx(i)): atommap[i]
                for i in range(mol.GetMaxAtomIdx())
                if atommap[i] is not None
            }
        )
    if bond_mapping is not None:
        bond_mapping.update(
            {
                mol.GetBond(oechem.OEHasBondIdx(i)): bondmap[i]
                for i in range(mol.GetMaxBondIdx())
                if bondmap[i] is not None
            }
        )
    return new_mol
