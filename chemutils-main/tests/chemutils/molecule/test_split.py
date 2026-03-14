import pytest
from openeye import oechem

from chemutils.molecule.equal import assert_mol_equal
from chemutils.molecule.smiles import oemol_from_smiles
from chemutils.molecule.split import split_molecule

TEST_CASES = [
    (
        "C1=CNC=C1CO",
        "CCCCCC",
    ),
    (
        "CN1CCOC(CO)(CO)C1=O",
        "CCCCCC",
    ),
    (
        "C1=CNC=C1CO",
        "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O",
    ),
    (
        "C1=CNC=C1CO",
        "CC1=C(C=C(C=C1)NC(=O)C2=CC=C(C=C2)CN3CCN(CC3)C)NC4=NC=CC(=N4)C5=CN=CC=C5",
    ),
]


@pytest.mark.parametrize(
    ["first_smiles_str", "second_smiles_str"],
    TEST_CASES,
)
def test_split_molecule(first_smiles_str, second_smiles_str):
    # Load the OEMol object containing two separate molecules

    expected_first_oemol = oemol_from_smiles(first_smiles_str)

    expected_second_oemol = oemol_from_smiles(second_smiles_str)

    both_smiles_str = f"{first_smiles_str}.{second_smiles_str}"

    expected_both_oemol = oemol_from_smiles(both_smiles_str)

    # Test splitting the object into two separate objects
    actual_first_oemol, actual_second_oemol = list(split_molecule(expected_both_oemol))

    # Add the two separate molecule objects back into a new OEMol.
    actual_both_oemol = actual_first_oemol.CreateCopy()
    oechem.OEAddMols(actual_both_oemol, actual_second_oemol)

    # Assert that the original and re-stitched OEMol are the same.
    assert_mol_equal(actual_both_oemol, expected_both_oemol)

    # Assert that the original and re-stitched first OEMol are the same.
    assert_mol_equal(actual_first_oemol, expected_first_oemol)

    # Assert that the original and re-stitched second OEMol are the same.
    assert_mol_equal(actual_second_oemol, expected_second_oemol)
