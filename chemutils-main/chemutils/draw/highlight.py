from collections import defaultdict
from collections.abc import Callable, Collection
from dataclasses import dataclass
from typing import TypeAlias

from openeye import oechem
from rdkit import Chem

from .color import Color, parse_color

HighlightSelection = list[oechem.OEAtomBase] | list[oechem.OEBondBase]

AtomPredicate: TypeAlias = (
    int
    | oechem.OEAtomBase
    | list[oechem.OEAtomBase]
    | list[int]
    | oechem.OEAtomIter
    | oechem.OEUnaryAtomBoolFunc
    | Callable[[oechem.OEAtomBase], bool]
    | oechem.OEAtomBondSet
    | oechem.OEMatchBase
)
BondPredicate: TypeAlias = (
    int
    | oechem.OEBondBase
    | list[oechem.OEBondBase]
    | list[int]
    | oechem.OEBondIter
    | oechem.OEUnaryBondBoolFunc
    | Callable[[oechem.OEBondBase], bool]
    | oechem.OEAtomBondSet
    | oechem.OEMatchBase
)
AtomBondPredicate: TypeAlias = (
    list[oechem.OEAtomBase | oechem.OEBondBase] | oechem.OEAtomBondSet | oechem.OEMatchBase
)

DEFAULT_HIGHLIGHT = "gold"


@dataclass(kw_only=True, frozen=True)
class Highlight:
    """
    Highlight atoms and bonds in a molecule.

    Either a rule that applies to both atoms and bonds must be applied (using `atom_and_bonds`), or
    a combination of `atoms` and `bonds` rules can be used.

    The `atoms` argument can be one of:
        * An atom index
        * An `OEAtomBase`
        * A list of atom indices.
        * A list of `OEAtomBase`
        * An `OEAtomIter`
        * An `OEUnaryAtomBoolFunc`
        * A function that takes an `OEAtomBase` and returns a boolean.
        * An `OEAtomBondSet`
        * An `OEMatchBase`

    The `bonds` argument can be one of:
        * A bond index
        * An `OEBondBase`
        * A list of bond indices.
        * A list of `OEBondBase`
        * An `OEBondIter`
        * An `OEUnaryBondBoolFunc`
        * A function that takes an `OEBondBase` and returns a boolean.
        * An `OEAtomBondSet`
        * An `OEMatchBase`

    The `atoms_and_bonds` argument can be one of:
        * A list of `OEAtomBase`` and `OEBondBase``
        * An `OEAtomBondSet`
        * An `OEMatchBase`

    Args:
        atoms: Rule to select which atoms to highlight.
        bonds: Rule to select which bonds to highlight.
        atoms_and_bonds: Rule to select which atoms and bonds to highlight.
        color: Color to highlight with. Can be either a hex code or a color name, such as
            `#ffcc88` or `white`.
        connected: Should bonds between highlighted atoms be highlighted?
        overwrite: Should this highlight overwrite other highlights, rather than be split?

    Returns:
        Highlight object to be passed to `draw_molecule`.
    """

    color: Color
    atoms: AtomPredicate | None = None
    bonds: BondPredicate | None = None
    atoms_and_bonds: AtomBondPredicate | None = None
    connected: bool = False
    overwrite: bool = False


def get_atoms_and_bonds(
    *,
    oeatoms: Collection[oechem.OEAtomBase],
    oebonds: Collection[oechem.OEBondBase],
    atoms: AtomPredicate | None = None,
    bonds: BondPredicate | None = None,
    atoms_and_bonds: AtomBondPredicate | None = None,
) -> tuple[list[oechem.OEAtomBase], list[oechem.OEBondBase]]:
    if atoms_and_bonds is None and atoms is None and bonds is None:
        raise ValueError(
            "Molecule highlighting error: must provide either `atoms` and/or `bonds`, or `atoms_and_bonds`."
        )
    if atoms_and_bonds is not None and (atoms is not None and bonds is not None):
        raise ValueError(
            "Molecule highlighting error: must provide either `atoms` and/or `bonds`, or `atoms_and_bonds`."
        )
    if atoms_and_bonds:
        match atoms_and_bonds:
            case list():
                return [item for item in atoms_and_bonds if isinstance(item, oechem.OEAtomBase)], [
                    item for item in atoms_and_bonds if isinstance(item, oechem.OEBondBase)
                ]
            case oechem.OEAtomBondSet() | oechem.OEMatchBase():
                return get_atoms_and_bonds(
                    oeatoms=oeatoms, oebonds=oebonds, atoms=atoms_and_bonds, bonds=atoms_and_bonds
                )
            case _:
                raise ValueError(
                    f"Molecule highlighting failed: unknown type {type(atoms_and_bonds)} for `atoms_and_bonds`."
                )
    else:
        atom_list = []
        bond_list = []
        match atoms:
            case None:
                pass
            case int():
                atom_list = [atom for atom in oeatoms if atom.GetIdx() == atoms]
            case oechem.OEAtomBase():
                atom_list = [atom for atom in oeatoms if atom == atom]
            case list():
                if len(atoms) == 0:
                    atom_list = []
                else:
                    if isinstance(atoms[0], int):
                        atom_list = [atom for atom in oeatoms if atom.GetIdx() in atoms]
                    else:
                        atom_list = atoms
            case oechem.OEAtomIter():
                atom_list = list(atoms)
            case oechem.OEUnaryAtomBoolFunc():
                atom_list = [atom for atom in oeatoms if atoms(atom)]
            case _ if callable(atoms):
                atom_list = [atom for atom in oeatoms if atoms(atom)]
            case oechem.OEAtomBondSet():
                atom_list = list(atoms.GetAtoms())
            case oechem.OEMatchBase():
                atom_list = list(atoms.GetPatternAtoms()) + list(atoms.GetTargetAtoms())
            case _:
                raise ValueError(
                    f"Molecule highlighting failed: unknown type {type(atoms)} for `atoms`."
                )
        match bonds:
            case None:
                pass
            case int():
                bond_list = [bond for bond in oebonds if bond.GetIdx() == bonds]
            case oechem.OEBondBase():
                bond_list = [bonds]
            case list():
                if len(bonds) == 0:
                    bond_list = []
                else:
                    if isinstance(bonds[0], int):
                        bond_list = [bond for bond in oebonds if bond.GetIdx() in bonds]
                    else:
                        bond_list = bonds
            case oechem.OEBondIter():
                bond_list = list(bonds)
            case oechem.OEUnaryBondBoolFunc():
                bond_list = [bond for bond in oebonds if bonds(bond)]
            case _ if callable(bonds):
                bond_list = [bond for bond in oebonds if bonds(bond)]
            case oechem.OEAtomBondSet():
                bond_list = list(bonds.GetBonds())
            case oechem.OEMatchBase():
                bond_list = list(bonds.GetPatternBonds()) + list(bonds.GetTargetBonds())
            case _:
                raise ValueError(
                    f"Molecule highlighting failed: unknown type {type(bonds)} for `bonds`."
                )
        return atom_list, bond_list


def get_highlights(
    *highlights: Highlight,
    rdmol: Chem.Mol,
    oeatom_to_rdatom: dict[oechem.OEAtomBase, Chem.Atom],
    oebond_to_rdbond: dict[oechem.OEBondBase, Chem.Bond],
) -> tuple[
    dict[int, list[tuple[float, float, float]]], dict[int, list[tuple[float, float, float]]]
]:
    atom_highlights = defaultdict(list)
    bond_highlights = defaultdict(list)

    for highlight in highlights:
        color = parse_color(highlight.color)

        atoms, bonds = get_atoms_and_bonds(
            oeatoms=oeatom_to_rdatom,
            oebonds=oebond_to_rdbond,
            atoms=highlight.atoms,
            bonds=highlight.bonds,
            atoms_and_bonds=highlight.atoms_and_bonds,
        )

        rdatom_idxs = {oeatom_to_rdatom[atom].GetIdx() for atom in atoms}
        rdbond_idxs = {oebond_to_rdbond[bond].GetIdx() for bond in bonds}

        for rdatom_idx in rdatom_idxs:
            if highlight.overwrite:
                atom_highlights[rdatom_idx] = [color]
            else:
                atom_highlights[rdatom_idx].append(color)

        for rdbond_idx in rdbond_idxs:
            if highlight.overwrite:
                bond_highlights[rdbond_idx] = [color]
            else:
                bond_highlights[rdbond_idx].append(color)

        if highlight.connected:
            for bond in rdmol.GetBonds():
                if bond.GetBeginAtomIdx() in rdatom_idxs and bond.GetEndAtomIdx() in rdatom_idxs:
                    if highlight.overwrite:
                        bond_highlights[bond.GetIdx()] = [color]
                    else:
                        bond_highlights[bond.GetIdx()].append(color)

    return dict(atom_highlights), dict(bond_highlights)
