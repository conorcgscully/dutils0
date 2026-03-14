import pytest

from chemutils.molecule import oemol_from_smiles
from chemutils.props.bonds_rings import (
    get_largest_ring_size,
    get_num_aromatic_rings,
    get_num_rotatable_bonds,
    is_macrocycle,
)


@pytest.mark.parametrize(
    ["smiles", "expected"],
    [
        ("CC", 0),
        ("CCC", 0),
        ("CCCC", 1),
        ("C=CC", 0),
        ("C=CCC", 1),
        ("CC=CC", 0),
        ("CC(O)OCC", 2),
        # An ester is not considered rotatable
        ("CC(=O)OCC", 1),
        ("C1CCCCC1", 0),
        ("c1ccccc1", 0),
    ],
)
def test_rotatable_bonds(smiles, expected):
    mol = oemol_from_smiles(smiles)
    assert get_num_rotatable_bonds(mol) == expected


@pytest.mark.parametrize(
    ["smiles", "expected"],
    [
        ("C1CCCCC1", 0),  # Cyclohexane
        ("c1ccccc1", 1),  # Benzene
        ("C12=C(C=CN2)C=CC=C1", 2),  # Indole
        ("c1ccoc1", 1),  # Furan
        ("c1cnc[nH]1", 1),  # Imidazole
        ("c1c2c(nc[nH]2)ncn1", 2),  # Purine
        ("c1c2ccccc2ccc1", 2),  # Naphthalene
        ("c1ccc2cc3ccccc3cc2c1", 3),  # Anthracene
        ("c1ccncc1", 1),  # Pyridine
        ("c1ccc2c(c1)ccc3c2ccc4c3cccc4", 4),  # Chrysene
    ],
)
def test_aromatic_rings(smiles, expected):
    mol = oemol_from_smiles(smiles)
    assert get_num_aromatic_rings(mol) == expected


@pytest.mark.parametrize(
    ["smiles", "expected"],
    [
        ("CCCCC", 0),  # Pentane, no rings
        ("C1CCCCC1", 6),  # 6-member ring
        ("C1CCCCCCC1", 8),  # 8-member ring
        ("C1CCCCCCCCC1", 10),  # 10-member ring
        ("C1CCCCCCCCCCC1", 12),  # 12-member ring
        ("C1CNCCCNCCC1", 10),  # 10-member ring, heteroatoms
        ("C12=C(C=CN2)C=CC=C1", 6),  # Indole, largest ring is 6
        (r"O=C1/C=C(\C=C/C=C1/O)C(C)C", 7),  # Hinokitiol, 7-member
        ("O1CCOCCOCCOCC1", 12),  # 12-crown-4 crown ether
        ("C[C@@H]1CCCCCCCCCCCCC(=O)C1", 15),  # Muscone, 15-member
        ("O1CCOCCOCCOCCOCCOCC1", 18),  # 18-crown-6 crown ether
        ("c1cc2cc(c1)-c3ccnc(n3)Nc4ccc(c(c4)COC/C=C/COC2)OCCN5CCCC5", 18),  # Pacritinib, 18-member
    ],
)
def test_largest_ring_size(smiles, expected):
    mol = oemol_from_smiles(smiles)
    assert get_largest_ring_size(mol) == expected


@pytest.mark.parametrize(
    ["smiles", "expected"],
    [
        ("C1CCCCC1", False),  # 6-member ring
        ("C1CCCCCCC1", False),  # 8-member ring
        ("C1CCCCCCCCC1", True),  # 10-member ring
        ("C1CCCCCCCCCCC1", True),  # 12-member ring
        ("C1CNCCCNCCC1", True),  # 10-member ring, heteroatoms
    ],
)
def test_is_macrocycle(smiles, expected):
    mol = oemol_from_smiles(smiles)
    assert is_macrocycle(mol) == expected
