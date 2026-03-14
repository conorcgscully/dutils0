import pytest
from openeye import oechem

from chemutils.molecule.equal import assert_mol_equal
from chemutils.molecule.smiles import (
    canonicalize_smiles,
    cxsmiles_from_oemol,
    oemol_from_smiles,
    smiles_from_oemol,
)

from ..cip_validation import CIP_VALIDATION_OEMOLS
from ..props.drug_smiles import DRUG_SMILES

TEST_SMILES = [
    "C",  # Single carbon
    "CC",  # Ethane
    "C=C",  # Ethene, double bond
    "N/C=C/O",  # Trans bond
    r"N\C=C/O",  # Cis bond
    "C#C",  # Ethyne, triple bond,
    "C1CCCCC1",  # Cyclohexane, single ring
    "c1ccccc1",  # Benzene, single aromatic ring,
    "[NH4+]",  # Ammonium ion,
    "[C@](F)(Cl)(I)Br",  # Chiral molecule
    "[C@@](F)(Cl)(I)Br",  # Chiral molecule
    r"C1/C=C\C1",  # Ring of size 4 with double bond
    r"C1C/C=C\C1",  # Ring of size 5 with double bond
    r"C1C/C=C\CC1",  # Ring of size 6 with double bond
    r"C1CC/C=C\CC1",  # Ring of size 7 with double bond
    "C(C)C",  # SMILES in non canonical order
]


# Test that oemol_from_smiles (using OEParseSmiles) gives the same result as OESmilesToMol.
@pytest.mark.parametrize(
    "smiles",
    TEST_SMILES
    + list(DRUG_SMILES.values())
    + [oechem.OEMolToSmiles(mol) for mol in CIP_VALIDATION_OEMOLS],
)
def test_oemol_from_smiles_equivalent_OESmilesToMol(smiles):
    mol_OESmilesToMol = oechem.OEGraphMol()
    assert oechem.OESmilesToMol(mol_OESmilesToMol, smiles)
    mol_oemol_from_smiles = oemol_from_smiles(smiles)
    assert_mol_equal(mol_OESmilesToMol, mol_oemol_from_smiles)


# Test that smiles_from_oemol gives the same result as OEMolToSmiles.
@pytest.mark.parametrize(
    "smiles",
    TEST_SMILES
    + list(DRUG_SMILES.values())
    + [oechem.OEMolToSmiles(mol) for mol in CIP_VALIDATION_OEMOLS],
)
def test_oemol_from_smile(smiles):
    mol = oemol_from_smiles(smiles)
    assert smiles_from_oemol(mol) == oechem.OEMolToSmiles(mol)


@pytest.mark.parametrize(
    "smiles",
    TEST_SMILES
    + list(DRUG_SMILES.values())
    + [oechem.OEMolToSmiles(mol) for mol in CIP_VALIDATION_OEMOLS]
    + ["[C@H](F)(Cl)I |&1:1|"],  # Enhanced stereo
)
def test_canonicalize_smiles(smiles):
    mol = oemol_from_smiles(smiles)
    assert smiles_from_oemol(mol) == canonicalize_smiles(smiles)


def test_smiles_from_oemol_atom_order():
    mol = oemol_from_smiles("C(C)C")
    atoms = list(mol.GetAtoms())
    order = []
    assert smiles_from_oemol(mol, atom_order=order) == "CCC"
    assert order == [atoms[1], atoms[0], atoms[2]]


@pytest.mark.parametrize(
    ["cxsmiles", "expected_cxsmiles", "expected_groups"],
    [
        (
            "C[C@H]([C@@H](N)Cl)O |&1:1,&2:2|",
            "C[C@H]([C@@H](N)Cl)O |&1:1,&2:2|",
            [oechem.OEGroupType_MDLAndStereo, oechem.OEGroupType_MDLAndStereo],
        ),
        (
            "C[C@H]([C@@H](N)Cl)O |o1:1,&2:2|",
            "C[C@H]([C@@H](N)Cl)O |o1:1,&1:2|",
            [oechem.OEGroupType_MDLOrStereo, oechem.OEGroupType_MDLAndStereo],
        ),
        (
            "C[C@H]([C@@H](N)Cl)O |o1:1,a:2|",
            "C[C@H]([C@@H](N)Cl)O |a:2,o1:1|",
            [oechem.OEGroupType_MDLOrStereo, oechem.OEGroupType_MDLAbsStereo],
        ),
        # With extra options
        (
            "C[C@H]([C@@H](N)Cl)O |o1:1,K:4,&1:2|",
            "C[C@H]([C@@H](N)Cl)O |o1:1,&1:2|",
            [oechem.OEGroupType_MDLOrStereo, oechem.OEGroupType_MDLAndStereo],
        ),
        (
            "C[C@H]([C@@H](N)Cl)O |o1:1,a:2,c:3|",
            "C[C@H]([C@@H](N)Cl)O |a:2,o1:1|",
            [oechem.OEGroupType_MDLOrStereo, oechem.OEGroupType_MDLAbsStereo],
        ),
    ],
)
def test_oemol_from_smiles_parses_cxsmiles(cxsmiles, expected_cxsmiles, expected_groups):
    oemol = oemol_from_smiles(cxsmiles)
    assert oemol.GetTitle() == ""
    groups = [group.GetGroupType() for group in oemol.GetGroups(oechem.OEIsMDLStereoGroup())]
    assert groups == expected_groups
    assert cxsmiles_from_oemol(oemol) == expected_cxsmiles


def test_cxsmiles_unknown_tags_regression():
    smiles = "Cl[C@H](F)I |&1:1,c:2|"

    # Raw OpenEye
    oemol_raw = oechem.OEGraphMol()
    oechem.OESmilesToMol(oemol_raw, smiles)

    # Title is set to the tags
    assert oemol_raw.GetTitle() == "|&1:1,c:2|"
    # No enhanced stereo is parsed
    assert oemol_raw.NumGroups() == 0
    # SMILES is rearranged, but tags are left unchanged, so now it's invalid 😡
    assert cxsmiles_from_oemol(oemol_raw) == "[C@H](F)(Cl)I |&1:1,c:2|"

    # Corrected chemutils parsing
    oemol = oemol_from_smiles(smiles)

    # Title is cleared
    assert oemol.GetTitle() == ""
    # Enhanced stereo is parsed
    assert oemol.NumGroups() == 1
    # SMILES is rearranged, and so are tags 😃
    assert cxsmiles_from_oemol(oemol) == "[C@H](F)(Cl)I |&1:0|"


@pytest.mark.parametrize(["smiles", "num_atoms"], [("C", 1), ("CC", 2), ("CCC", 3)])
def test_smiles_atoms(smiles, num_atoms):
    mol = oemol_from_smiles(smiles)
    assert mol.NumAtoms() == num_atoms


@pytest.mark.parametrize(
    ["smiles", "num_bonds"], [("C", 0), ("CC", 1), ("C(C)(C)C", 3), ("c1ccccc1", 6)]
)
def test_bonds(smiles, num_bonds):
    mol = oemol_from_smiles(smiles)
    assert mol.NumBonds() == num_bonds


@pytest.mark.parametrize(["smiles", "num_double_bonds"], [("CC", 0), ("C=C", 1), ("C=CC=C", 2)])
def test_double_bond(smiles, num_double_bonds):
    mol = oemol_from_smiles(smiles)
    assert len([bond for bond in mol.GetBonds() if bond.GetOrder() == 2]) == num_double_bonds


@pytest.mark.parametrize(
    ["smiles", "num_triple_bonds"],
    [
        ("CC", 0),
        ("C#C", 1),
    ],
)
def test_triple_bond(smiles, num_triple_bonds):
    mol = oemol_from_smiles(smiles)
    assert len([bond for bond in mol.GetBonds() if bond.GetOrder() == 3]) == num_triple_bonds


@pytest.mark.parametrize(
    ["smiles", "num_aromatic_atoms", "num_aromatic_bonds"], [("C1CCCCC1", 0, 0), ("c1ccccc1", 6, 6)]
)
def test_aromatic(smiles, num_aromatic_atoms, num_aromatic_bonds):
    mol = oemol_from_smiles(smiles)
    assert len([atom for atom in mol.GetAtoms() if atom.IsAromatic()]) == num_aromatic_atoms
    assert len([bond for bond in mol.GetBonds() if bond.IsAromatic()]) == num_aromatic_bonds


@pytest.mark.parametrize("smiles", ["X", "C£C", 'C"C'])
def test_invalid(smiles):
    with pytest.raises(ValueError):
        _ = oemol_from_smiles(smiles)


@pytest.mark.parametrize(
    "smiles", ["N", "[NH3]", "[NH4+]", "C", "C=O", "C#C", "CCC", "C(C)C", "c1ccccc1", "C1CCCCC1"]
)
def test_equivalence_to_oemol(smiles):
    mol_theirs = oechem.OEGraphMol()
    oechem.OESmilesToMol(mol_theirs, smiles)
    mol_ours = oemol_from_smiles(smiles)
    assert_mol_equal(mol_ours, mol_theirs)
