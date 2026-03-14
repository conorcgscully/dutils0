import pytest

from chemutils.molecule import (
    get_atom,
    oemol_from_smiles,
    remove_atom,
    remove_bond,
    smiles_from_oemol,
    subset_molecule,
)

CASES = [
    ("CCC", 0, "hydrogens", "CC"),
    ("CCC", 0, "dangling", "C[CH2]"),
    ("CCC", 0, "rgroups", "[R1]CC"),
    ("CCC", 1, "hydrogens", "C.C"),
    ("CCC", 1, "dangling", "[CH3].[CH3]"),
    ("CCC", 1, "rgroups", "[R1]C.[R2]C"),
    ("C1CCCCC1", 3, "hydrogens", "CCCCC"),
    ("C1CCCCC1", 3, "dangling", "[CH2]CCC[CH2]"),
    ("C1CCCCC1", 3, "rgroups", "[R1]CCCCC[R2]"),
    ("c1ccccc1", 4, "hydrogens", "CC=CC=C"),
    ("c1ccccc1", 4, "dangling", "[CH]=CC=C[CH]"),
    ("c1ccccc1", 4, "rgroups", "[R1]=CC=CC=C[R2]"),
]


@pytest.mark.parametrize(
    ["smiles", "index", "bond_handling", "new_smiles"],
    CASES,
)
def test_remove_atom(smiles, index, bond_handling, new_smiles):
    oemol = oemol_from_smiles(smiles)
    atom = get_atom(oemol, selection=index)
    remove_atom(atom, bond_handling=bond_handling)
    assert smiles_from_oemol(oemol) == new_smiles


@pytest.mark.parametrize(
    ["smiles", "index", "bond_handling", "new_smiles"],
    CASES,
)
def test_remove_atom_same_behaviour_as_subset_molecule(smiles, index, bond_handling, new_smiles):
    oemol = oemol_from_smiles(smiles)
    subset = subset_molecule(oemol, selection=f"not index {index}", bond_handling=bond_handling)
    assert smiles_from_oemol(subset) == new_smiles


def test_remove_atom_dont_adjust_ring_aromaticity():
    oemol = oemol_from_smiles("c1ccccc1")

    atom = get_atom(oemol, selection=0)

    assert all(atom.IsInRing() for atom in oemol.GetAtoms())
    assert all(atom.IsAromatic() for atom in oemol.GetAtoms())

    remove_atom(atom, bond_handling="hydrogens", adjust_rings=False, adjust_aromaticity=False)

    assert all(atom.IsInRing() for atom in oemol.GetAtoms())
    assert all(atom.IsAromatic() for atom in oemol.GetAtoms())


def test_remove_atom_dont_adjust_aromaticity():
    oemol = oemol_from_smiles("c1ccccc1")

    atom = get_atom(oemol, selection=0)

    assert all(atom.IsInRing() for atom in oemol.GetAtoms())
    assert all(atom.IsAromatic() for atom in oemol.GetAtoms())

    remove_atom(atom, bond_handling="hydrogens", adjust_rings=True, adjust_aromaticity=False)

    assert all(not atom.IsInRing() for atom in oemol.GetAtoms())
    assert all(atom.IsAromatic() for atom in oemol.GetAtoms())


def test_remove_atom_dont_adjust_ring():
    oemol = oemol_from_smiles("c1ccccc1")

    atom = get_atom(oemol, selection=0)

    assert all(atom.IsInRing() for atom in oemol.GetAtoms())
    assert all(atom.IsAromatic() for atom in oemol.GetAtoms())

    with pytest.raises(ValueError):
        remove_atom(atom, bond_handling="hydrogens", adjust_rings=False, adjust_aromaticity=True)


@pytest.mark.parametrize(
    ["smiles", "index", "bond_handling", "new_smiles"],
    [
        ("CCCC", 0, "hydrogens", "C.CCC"),
        ("CCCC", 0, "dangling", "[CH3].CC[CH2]"),
        ("CCCC", 0, "electrons", "[CH3-].CC[CH2-]"),
        ("CCCC", 0, "rgroups", "[R1]C.[R2]CCC"),
        ("CCCC", 1, "hydrogens", "CC.CC"),
        ("CCCC", 1, "dangling", "C[CH2].C[CH2]"),
        ("CCCC", 1, "electrons", "C[CH2-].C[CH2-]"),
        ("CCCC", 1, "rgroups", "[R1]CC.[R2]CC"),
        ("CC=CC", 1, "hydrogens", "CC.CC"),
        ("CC=CC", 1, "dangling", "C[CH].C[CH]"),
        ("CC=CC", 1, "electrons", "C[CH-2].C[CH-2]"),
        ("CC=CC", 1, "rgroups", "[R1]=CC.[R2]=CC"),
        ("C1CCCCC1", 3, "hydrogens", "CCCCCC"),
        ("C1CCCCC1", 3, "dangling", "[CH2]CCCC[CH2]"),
        ("C1CCCCC1", 3, "electrons", "[CH2-]CCCC[CH2-]"),
        ("C1CCCCC1", 3, "rgroups", "[R1]CCCCCC[R2]"),
        ("c1ccccc1", 3, "hydrogens", "C=CC=CC=C"),
        ("c1ccccc1", 3, "dangling", "[CH]=CC=CC=[CH]"),
        ("c1ccccc1", 3, "electrons", "[CH-]=CC=CC=[CH-]"),
        ("c1ccccc1", 3, "rgroups", "[R1]C=CC=CC=C[R2]"),
    ],
)
def test_remove_bond(smiles, index, bond_handling, new_smiles):
    oemol = oemol_from_smiles(smiles)
    bond = list(oemol.GetBonds())[index]
    remove_bond(bond, bond_handling=bond_handling)
    assert smiles_from_oemol(oemol) == new_smiles
