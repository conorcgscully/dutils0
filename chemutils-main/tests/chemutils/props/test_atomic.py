import pytest
from openeye import oechem
from rdkit import Chem

from chemutils.molecule.smiles import oemol_from_smiles
from chemutils.props.atomic import (
    get_fraction_carbon_sp3,
    get_molecular_weight,
    get_num_chiral_atoms,
    get_num_heavy_atoms,
    get_num_spiro_atoms,
    get_num_unspecified_chiral_atoms,
    get_spacial_score,
)


@pytest.mark.parametrize(
    ["smiles", "expected"],
    [
        ("C", 1),
        ("[CH4]", 1),
        ("[C]", 1),
        ("[C]([H])([H])([H])[H]", 1),
        ("CO", 2),
        ("[C@H](O)(N)C", 4),
        ("O[C@@H](N)C", 4),
    ],
)
def test_get_num_heavy_atoms(smiles, expected):
    mol = oemol_from_smiles(smiles)
    assert get_num_heavy_atoms(mol) == expected


CHIRAL_ATOM_CASES = [
    ("F/C=C/C[C@@H](Cl)Br", 1, 1),
    ("F/C=C/C[C@H](Cl)Br", 1, 1),
    ("F/C=C/C[CH](Cl)Br", 1, 1),
    ("ClC(Br)C1CCCC(C1)/C=C/Cl", 3, 3),
    ("[C@H](O)(N)C", 1, 1),
    ("c1ccccc1", 0, 0),
    ("c1ccccc1C2=NCC(=O)N(C)c3ccc(Cl)cc23", 0, 0),
    ("CCN(C)CCC", 1, 0),
    ("CN1CC(C)(C)N(CC)CC1", 2, 0),
    ("CN1CC(C)(C)N(CC)C(C)C1", 3, 1),
]


@pytest.mark.parametrize(
    ["smiles", "num_chiral"],
    [
        (smiles, num_chiral_with_nitrogens)
        for smiles, num_chiral_with_nitrogens, _ in CHIRAL_ATOM_CASES
    ],
)
def test_get_num_chiral_atoms_with_nitrogens(smiles, num_chiral):
    mol = oemol_from_smiles(smiles)
    assert get_num_chiral_atoms(mol) == num_chiral


@pytest.mark.parametrize(
    ["smiles", "num_chiral"],
    [(smiles, num_chiral_no_nitrogens) for smiles, _, num_chiral_no_nitrogens in CHIRAL_ATOM_CASES],
)
def test_get_num_chiral_atoms_no_nitrogens(smiles, num_chiral):
    mol = oemol_from_smiles(smiles)
    assert get_num_chiral_atoms(mol, include_nitrogens=False) == num_chiral


@pytest.mark.parametrize(
    ["smiles", "expected"],
    [
        ("F/C=C/C[C@@H](Cl)Br", 0),
        ("F/C=C/C[C@H](Cl)Br", 0),
        ("F/C=C/C[CH](Cl)Br", 1),
        ("ClC(Br)C1CCCC(C1)/C=C/Cl", 3),
        ("[CH](O)(N)C", 1),
        ("[C@H](O)(N)C", 0),
        ("c1ccccc1", 0),
        ("c1ccccc1C2=NCC(=O)N(C)c3ccc(Cl)cc23", 0),
    ],
)
def test_get_num_unspecified_chiral_atoms(smiles, expected):
    mol = oemol_from_smiles(smiles)
    assert get_num_unspecified_chiral_atoms(mol) == expected


@pytest.mark.parametrize(
    ["smiles", "expected"],
    [
        ("c1ccccc1", 0),
        ("c1ccccc1C2=NCC(=O)N(C)c3ccc(Cl)cc23", 0),
        ("C1CCCCC12CCCCC2", 1),
        ("C12CCCCC1CCCC2", 0),
    ],
)
@pytest.mark.parametrize("rdkit", [True, False])
def test_get_num_spiro_atoms(smiles, expected, rdkit):
    mol = Chem.MolFromSmiles(smiles) if rdkit else oemol_from_smiles(smiles)
    assert get_num_spiro_atoms(mol) == expected


@pytest.mark.parametrize(
    ["smiles", "expected"],
    [
        # Examples from the paper https://pubs.acs.org/doi/10.1021/acs.jmedchem.3c00689
        # Table 1 Example
        ("c1ccccc1C2CCCC2OC=C", 26.0),
        # Table 2
        ("C=CCCBr", 7.4),  # Example 1
        ("CCCCC", 8.4),  # Example 2
        ("CCCCBr", 8.4),  # Example 3
        ("CC(C)CBr", 9.6),  # Example 4
        ("CC=CCBr", 10.0),  # Example 5
        ("CC(C)(C)Br", 12.0),  # Example 6
        ("CCC(C)Br", 15.0),  # Example 7
        ("c1ccccc1", 8.0),  # Example 8
        ("C1=CC=CC=N1", 8.0),  # Example 9
        ("C1=CC=CCC1", 18.7),  # Example 10
        ("C1=CCCCC1", 21.4),  # Example 11
        ("C1CCCCC1", 24.0),  # Example 12
        ("C1CCCC1C", 25.5),  # Example 13
        ("C1CCC1(C)C", 29.0),  # Example 14
        ("O=C1COCCN1c2ccc(cc2)N3C[C@@H](OC3=O)CNC(=O)c4ccc(s4)Cl", 19.4),  # Example 15, rivaroxaban
        (
            "C[C@H](CCCC(C)C)[C@H]1CC[C@@H]2[C@@]1(CC[C@H]3[C@H]2CC=C4[C@@]3(CC[C@@H](C4)O)C)C",
            46.5,
        ),  # Example 16, cholesterol
    ],
)
@pytest.mark.parametrize("rdkit", [True, False])
@pytest.mark.parametrize("explicit_hs", [True, False])
def test_get_spacial_score(smiles, expected, rdkit, explicit_hs):
    if rdkit:
        mol = Chem.MolFromSmiles(smiles)
        if explicit_hs:
            mol = Chem.AddHs(mol)
    else:
        mol = oemol_from_smiles(smiles)
        if explicit_hs:
            oechem.OEAddExplicitHydrogens(mol)
    assert get_spacial_score(mol) == pytest.approx(expected, abs=0.1)


@pytest.mark.parametrize(
    ["smiles", "expected"],
    [
        ("C", 16.04),
        ("[CH4]", 16.04),
        ("[C]", 12.01),
        ("[C]([H])([H])([H])[H]", 16.04),
    ],
)
def test_get_molecular_weight(smiles, expected):
    mol = oemol_from_smiles(smiles)
    assert get_molecular_weight(mol) == pytest.approx(expected, abs=0.01)


@pytest.mark.parametrize(
    ["smiles", "expected"],
    [
        ("C", 1 / 1),
        ("CC", 2 / 2),
        ("C(C)(C)(C)C", 5 / 5),
        ("CCO", 2 / 2),
        ("C=C", 0 / 2),
        ("CC=C", 1 / 3),
        ("C1CCCCC1", 6 / 6),
        ("c1ccccc1", 0 / 6),
        ("C#CC", 1 / 3),
        ("C=O", 0 / 1),
        ("CO", 1 / 1),
    ],
)
def test_fsp3(smiles, expected):
    mol = oemol_from_smiles(smiles)
    assert get_fraction_carbon_sp3(mol) == pytest.approx(expected, abs=0.01)
