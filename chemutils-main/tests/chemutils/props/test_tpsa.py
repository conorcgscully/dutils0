import pytest

from chemutils.molecule.smiles import oemol_from_smiles
from chemutils.props.tpsa import get_2d_tpsa

from .drug_smiles import DRUG_SMILES


@pytest.mark.parametrize(
    ["drug_name", "twod_tpsa"],
    [
        # Test example drugs from Table 3, Ertl, Rohde and Selzer (2000).
        ("metoprolol", 50.7),
        ("nordiazepam", 41.5),
        ("diazepam", 32.7),
        ("oxprenolol", 50.7),
        ("phenazone", 26.9),
        ("oxazepam", 61.7),
        ("alprenolol", 41.5),  # Typo in paper (was 41.9), corrected by summing manually
        ("practolol", 70.6),
        ("pindolol", 57.3),
        ("ciprofloxacin", 74.6),
        ("metolazone", 92.5),
        ("tranexamic acid", 63.3),
        ("atenolol", 84.6),
        ("sulpiride", 101.7),
        ("mannitol", 121.4),
        ("foscarnet", 94.8),
        ("sulfasalazine", 141.3),
        ("olsalazine", 139.8),
        ("lactulose", 189.5),  # Typo in paper (was 197.4), corrected by summing manually
        ("raffinose", 268.7),
    ],
)
def test_drug(drug_name, twod_tpsa):
    smiles = DRUG_SMILES[drug_name]
    mol = oemol_from_smiles(smiles)
    assert get_2d_tpsa(mol) == pytest.approx(twod_tpsa, abs=0.05)


# Test each fragment individually, from Table 1, Ertl, Rohde and Selzer (2000).
@pytest.mark.parametrize(
    ["smiles", "twod_tpsa"],
    [
        ("CN(C)C", 3.24),
        ("C=NC", 12.36),
        ("C#N", 23.79),
        ("C=N(C)=C", 11.68),
        ("C#N=C", 13.60),
        ("CN1CC1", 3.01),
        ("CNC", 12.03),
        ("C1NC1", 21.94),
        ("C=N", 23.85),
        ("CN", 26.02),
        ("C[N+](C)(C)C", 0),
        ("C=[N+](C)C", 3.01),
        ("C#[N+]C", 4.36),
        ("C[NH+](C)C", 4.44),
        ("C[NH+]=C", 13.97),
        ("C[NH2+]C", 16.61),
        ("C=[NH2+]", 25.59),
        ("C[NH3+]", 27.64),
        ("c1ccncc1", 12.89),
        ("c12n3cccc3ccc1cccc2", 4.41),
        ("c1n(C)ccc1", 4.93),
        ("c1cccc2c1c[nH]c2", 15.79),
        ("c1cc[nH+]cc1", 14.14),
        ("COC", 9.23),
        ("C1OC1", 12.53),
        ("C=O", 17.07),
        ("CO", 20.23),
        ("C[O-]", 23.06),
        ("c1ccoc1", 13.14),
    ],
)
def test_fragment(smiles, twod_tpsa):
    mol = oemol_from_smiles(smiles)
    assert get_2d_tpsa(mol) == pytest.approx(twod_tpsa, abs=0.005)
