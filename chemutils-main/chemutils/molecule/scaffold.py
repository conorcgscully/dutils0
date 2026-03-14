from collections import defaultdict
from typing import Literal

from openeye import oechem

from chemutils.molecule.remove import remove_atom
from chemutils.molecule.selection import subset_molecule

VALENCIES = {
    1: [1],  # Hydrogen
    5: [3],  # Boron
    6: [4],  # Carbon
    7: [3],  # Nitrogen
    8: [2],  # Oxygen
    9: [1],  # Fluorine
    14: [4],  # Silicon
    15: [3, 5],  # Phosphorus
    16: [2, 4, 6],  # Sulfur
    17: [1],  # Chlorine
    34: [2],  # Selenium
    35: [1],  # Bromine
    53: [1],  # Iodine
}


def get_ring_linker_subset(oemol: oechem.OEMolBase, /) -> oechem.OEAtomBondSet:
    """
    Get the subset of rings and linkers of a molecule, as defined by Bemis-Murcko.

    This includes all ring systems, as well as any direct paths between them. This
    is obtained by iteratively removing atoms of degree 1.

    Args:
        oemol: Molecule to obtain the rings and linkers subset of.

    Returns:
        AtomBondSet of the atoms and bonds that make up the rings and linkers.
    """
    # New way
    atombondset = oechem.OEAtomBondSet()

    atom_idx_to_atom = {}
    bond_idx_to_bond = {}

    potential_sidechain_atom_idxs = set()
    fragment_atom_idxs = set()
    fragment_bond_idxs = set()

    adjacent_atom_idxs = defaultdict(set)
    adjacent_bond_idxs = defaultdict(set)

    for atom in oemol.GetAtoms():
        idx = atom.GetIdx()
        atom_idx_to_atom[idx] = atom
        if atom.GetHvyDegree() == 1:
            potential_sidechain_atom_idxs.add(idx)
        fragment_atom_idxs.add(idx)

    for bond in oemol.GetBonds():
        idx = bond.GetIdx()
        bond_idx_to_bond[idx] = bond
        atom1_idx = bond.GetBgn().GetIdx()
        atom2_idx = bond.GetEnd().GetIdx()
        fragment_bond_idxs.add(idx)
        adjacent_atom_idxs[atom1_idx].add(atom2_idx)
        adjacent_atom_idxs[atom2_idx].add(atom1_idx)
        adjacent_bond_idxs[atom1_idx].add(idx)
        adjacent_bond_idxs[atom2_idx].add(idx)

    # Strictly decreasing, starting at terminal atoms and iterating inwards
    while potential_sidechain_atom_idxs:
        new_potential_sidechain_atom_idxs = set()
        for potential_sidechain_atom_idx in potential_sidechain_atom_idxs:
            if potential_sidechain_atom_idx not in fragment_atom_idxs:
                continue
            fragment_nbr_idxs = [
                nbr2
                for nbr2 in adjacent_atom_idxs[potential_sidechain_atom_idx]
                if nbr2 in fragment_atom_idxs
            ]
            if len(fragment_nbr_idxs) == 1:
                fragment_atom_idxs.remove(potential_sidechain_atom_idx)
                new_potential_sidechain_atom_idxs.add(fragment_nbr_idxs[0])
                fragment_bond_idxs -= adjacent_bond_idxs[potential_sidechain_atom_idx]
        potential_sidechain_atom_idxs = new_potential_sidechain_atom_idxs

    for atom_idx in fragment_atom_idxs:
        atombondset.AddAtom(atom_idx_to_atom[atom_idx])
    for bond_idx in fragment_bond_idxs:
        atombondset.AddBond(bond_idx_to_bond[bond_idx])
    return atombondset


def get_bemis_murcko_scaffold(
    oemol: oechem.OEGraphMol,
    /,
    external_double_bond_behaviour: Literal["placeholder", "keep", "remove"] = "placeholder",
) -> oechem.OEGraphMol:
    """
    Get the Bemis-Murcko scaffold of a molecule.

    By default, this will behave as described in the original paper, where double bonds connected to a ring
    or linker are maintained, but replaced with a generic R group (in the original paper, they were represented
    as pairs of electrons).

    This differs from the behaviour of RDKit's MurckoScaffold, which includes the atom at the end of these double
    bonds. To replicate this behaviour, set `external_double_bond_behaviour` to `"keep"`.

    These double bonds can also be removed completely by specifying `external_double_bond_behaviour` as `"remove"`.

    See https://github.com/rdkit/rdkit/discussions/6844 for the motivation behind this implementation.

    Args:
        oemol: Molecule to obtain the Bemis-Murcko scaffold of.
        external_double_bond_behaviour: Behaviour of double bonds connected to the scaffold. Can be `"placeholder"`,
            `"keep"`, or `"remove"`.
            When `"placeholder"`, double bonds connecting a scaffold atom to another atom are replaced with a double
            bond to a placeholder `*` atom.
            When `"keep"`, the double bond and it's attached atom are included in the scaffold.
            When `"remove"`, the double bond and it's attached atom are excluded from the scaffold.

    Returns:
        Molecule of the Bemis-Murcko scaffold.
    """
    ring_linker = get_ring_linker_subset(oemol)

    match external_double_bond_behaviour:
        case "keep":
            for atom in ring_linker.GetAtoms():
                for bond in atom.GetBonds():
                    if bond.GetOrder() > 1:
                        ring_linker.AddBond(bond)
                        ring_linker.AddAtom(bond.GetNbr(atom))
            submol = subset_molecule(
                oemol,
                selection=ring_linker,
                bond_handling="hydrogens",
                include_internal_bonds=False,
            )
        case "remove":
            submol = subset_molecule(
                oemol,
                selection=ring_linker,
                bond_handling="hydrogens",
                include_internal_bonds=False,
            )
            for atom in submol.GetAtoms():
                if (atomic_num := atom.GetAtomicNum()) in [15, 16]:
                    # Correct valency
                    degree = atom.GetHvyDegree()
                    desired_valency = min(val for val in VALENCIES[atomic_num] if val >= degree)
                    atom.SetImplicitHCount(desired_valency - degree)
        case "placeholder":
            submol = subset_molecule(
                oemol, selection=ring_linker, bond_handling="rgroups", include_internal_bonds=False
            )
            for atom in submol.GetAtoms():
                if atom.GetAtomicNum() == 0:
                    if atom.GetValence() == 1:
                        remove_atom(atom, bond_handling="hydrogens")
                    else:
                        # Relabel R-groups without numbers, as otherwise the SMILES will contain `[R<n>]`.
                        atom.SetMapIdx(0)
        case _:
            raise ValueError(
                f"Invalid external double bond behaviour: {external_double_bond_behaviour}"
            )

    return submol


def get_generic_bemis_murcko_scaffold(
    oemol: oechem.OEGraphMol,
    /,
) -> oechem.OEGraphMol:
    """
    Get the generic Bemis-Murcko scaffold of a molecule.

    The generic scaffold (also called the graph framework or cylic skeleton) is the scaffold where all bonds are
    replaced with aliphatic single bonds and all atoms are replaced with carbons. This can be obtained by setting
    `generic = True`.

    Args:
        oemol: Molecule to obtain the Bemis-Murcko scaffold of.

    Returns:
        Molecule of the generic Bemis-Murcko scaffold.
    """
    ring_linker = get_ring_linker_subset(oemol)

    submol = subset_molecule(
        oemol, selection=ring_linker, bond_handling="hydrogens", include_internal_bonds=False
    )
    for bond in submol.GetBonds():
        bond.SetOrder(1)
        bond.SetAromatic(False)
    for atom in submol.GetAtoms():
        atom.SetAtomicNum(6)
        atom.SetImplicitHCount(4 - atom.GetHvyDegree())
    return submol
