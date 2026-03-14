import pytest
from rdkit import Chem

from chemutils.molecule.smiles import oemol_from_smiles
from chemutils.props.logp import get_wildman_crippen_logp

from .drug_smiles import DRUG_SMILES


@pytest.mark.parametrize(
    ["drug_name", "expected_logp"],
    [
        # No example values in paper, so same drug list as XLogP paper.
        # Values calculated through RDKit.
        ("atropine", 1.93),
        ("chloramphenicol", 0.91),
        ("chlorothiazide", 0.13),
        ("chlorpromazine", 4.89),
        ("cimetidine", 0.60),
        ("diazepam", 3.15),
        ("diltiazem", 3.37),
        ("diphenhydramine", 3.35),
        ("flufenamic acid", 4.15),
        ("haloperidol", 4.43),
        ("imipramine", 3.88),
        ("lidocaine", 2.58),
        ("phenobarbital", 0.70),
        ("phenytoin", 1.77),
        ("procainamide", 1.34),
        ("propranolol", 2.58),
        ("tetracaine", 2.62),
        ("trimethoprim", 1.26),
        ("verapamil", 5.09),
    ],
)
@pytest.mark.parametrize("rdkit", [True, False])
def test_drug(drug_name, expected_logp, rdkit):
    smiles = DRUG_SMILES[drug_name]
    mol = Chem.MolFromSmiles(smiles) if rdkit else oemol_from_smiles(smiles)
    assert get_wildman_crippen_logp(mol) == pytest.approx(expected_logp, abs=0.005)


@pytest.mark.parametrize(
    "smiles",
    [*DRUG_SMILES.values()],
)
@pytest.mark.parametrize("rdkit", [True, False])
def test_agrees_rdkit(smiles, rdkit):
    rdkit_logp = Chem.Crippen.MolLogP(Chem.MolFromSmiles(smiles))
    mol = Chem.MolFromSmiles(smiles) if rdkit else oemol_from_smiles(smiles)
    assert get_wildman_crippen_logp(mol) == pytest.approx(rdkit_logp)
