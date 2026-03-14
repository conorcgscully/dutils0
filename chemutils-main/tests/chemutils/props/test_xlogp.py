import pytest

from chemutils.molecule.smiles import oemol_from_smiles
from chemutils.props.logp import get_oexlogp

from .drug_smiles import DRUG_SMILES


@pytest.mark.parametrize(
    ["drug_name", "oexlogp"],
    [
        # 19 drugs from original XLogP paper: J. Chem. Inf. Comput. Sci. 1997, 37, 3, 615-621
        ("atropine", 2.41),
        ("chloramphenicol", 0.02),
        ("chlorothiazide", 0.73),
        ("chlorpromazine", 4.58),
        ("cimetidine", -0.63),
        ("diazepam", 3.21),
        ("diltiazem", 2.95),
        ("diphenhydramine", 3.7),
        ("flufenamic acid", 3.08),
        ("haloperidol", 4.26),
        ("imipramine", 3.77),
        ("lidocaine", 2.58),
        ("phenobarbital", 1.3),
        ("phenytoin", 2.19),
        ("procainamide", 1.27),
        ("propranolol", 2.77),
        ("tetracaine", 2.82),
        ("trimethoprim", 1.08),
        ("verapamil", 5.71),
    ],
)
def test_drug(drug_name, oexlogp):
    smiles = DRUG_SMILES[drug_name]
    mol = oemol_from_smiles(smiles)
    assert get_oexlogp(mol) == pytest.approx(oexlogp, abs=0.005)


@pytest.mark.parametrize(
    "smiles",
    [
        "Oc1ccc(cc1[N+](=O)[O-])[As](O)(O)=O",
    ],
)
def test_invalid(smiles):
    mol = oemol_from_smiles(smiles)
    assert get_oexlogp(mol) is None
