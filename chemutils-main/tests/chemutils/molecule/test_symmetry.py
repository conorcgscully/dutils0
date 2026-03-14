import numpy as np
import pytest

from chemutils.molecule.smiles import oemol_from_smiles
from chemutils.molecule.symmetry import get_symmetry_classes, order_symmetry_classes


@pytest.mark.parametrize(
    ["smiles", "expected_symmetry"],
    [
        ("C", [1]),
        ("CC", [1, 1]),
        ("C=O", [1, 2]),
        ("CCC", [1, 2, 1]),
        ("C=CC", [1, 2, 3]),
        ("C=CC=C", [1, 2, 2, 1]),
        ("c1ccccc1", [1, 1, 1, 1, 1, 1]),
        ("C12C3C4C1C5C2C3C45", [1, 1, 1, 1, 1, 1, 1, 1]),
        ("C12C3C1C4C2C34", [1, 1, 1, 1, 1, 1]),
        ("CN1C=NC2=C1C(=O)N(C(=O)N2C)C", 1 + np.arange(14)),
    ],
)
def test_get_symmetry_classes_valid(smiles, expected_symmetry):
    mol = oemol_from_smiles(smiles)
    assert np.array_equal(get_symmetry_classes(mol), expected_symmetry)


@pytest.mark.parametrize(
    ["array", "expected"],
    [
        ([1, 2, 3], [1, 2, 3]),
        ([1, 3, 1], [1, 2, 1]),
        ([3, 2, 1], [1, 2, 3]),
        ([3, 2, 3, 2, 1], [1, 2, 1, 2, 3]),
        ([1], [1]),
        ([3], [1]),
        ([], []),
    ],
)
def test_order_symmetry_classes(array, expected):
    assert np.array_equal(np.array(order_symmetry_classes(array)), expected)


@pytest.mark.parametrize(
    ["smiles", "num_atoms", "num_symmetries"],
    [
        # ibuprofen, 2 symm methyls, 2 symm oxygens in carboxylic acid, 2x2 symm carbons in ring
        ("CC(C)Cc1ccc(cc1)C(C)C(=O)O", 15, 11),
        # AMP, 3 symmetric oxygens
        ("O=P(O)(O)OC[C@H]3O[C@@H](n2cnc1c(ncnc12)N)[C@H](O)[C@@H]3O", 23, 21),
        # nitrobenzene with N substitution to break ring symmetry. 2 symmetric oxygens
        ("c1cnc(cc1)[N+](=O)[O-]", 9, 8),
        # 2-Thienylboronic acid, 2 symmetric oxygens
        ("c1sc(cc1)[B](-O)[O-]", 8, 7),
    ],
)
def test_delocalized_symmetry(smiles, num_atoms, num_symmetries):
    mol = oemol_from_smiles(smiles)
    symms = get_symmetry_classes(mol)
    assert (mol.NumAtoms() - np.unique(symms).shape[0]) == (num_atoms - num_symmetries)
