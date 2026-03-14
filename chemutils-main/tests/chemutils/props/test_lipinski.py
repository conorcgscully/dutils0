import pytest

from chemutils.molecule.smiles import oemol_from_smiles
from chemutils.props.hydrogen_bonds import (
    get_lipinski_hydrogen_acceptors,
    get_lipinski_hydrogen_donors,
)

from .drug_smiles import DRUG_SMILES


@pytest.mark.parametrize(
    ["drug_name", "expected_donors", "expected_acceptors"],
    [
        # Table 1, Lipisnki (1997), with donors modified to be per heavy atom rather than per hydrogen.
        ("aciclovir", 3, 8),
        ("alprazolam", 0, 4),
        ("aspirin", 1, 4),
        ("atenolol", 3, 5),
        ("azithromycin", 5, 14),
        ("azidothymidine", 2, 9),
        ("benzylpenicillin", 2, 6),
        ("caffeine", 0, 6),
        ("candoxatril", 2, 8),
        ("captopril", 1, 4),
        ("carbamazepine", 1, 3),
        ("chloramphenicol", 3, 7),
        ("cimetidine", 3, 6),
        ("clonidine", 2, 3),
        ("cyclosporine", 5, 23),
    ],
)
def test_drug(drug_name, expected_donors, expected_acceptors):
    smiles = DRUG_SMILES[drug_name]
    oemol = oemol_from_smiles(smiles)
    assert get_lipinski_hydrogen_donors(oemol) == expected_donors
    assert get_lipinski_hydrogen_acceptors(oemol) == expected_acceptors
