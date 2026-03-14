import pytest
from openeye import oechem
from rdkit import Chem
from rdkit.Chem import QED

from chemutils.molecule.smiles import oemol_from_smiles
from chemutils.props.druglike import get_qed

from .drug_smiles import DRUG_SMILES


# Test cases that would fool standard RDKit QED.
@pytest.mark.parametrize(
    "smiles",
    [r"CC12C3([Ir]1456(C2(C4(C53C)C)C)[N]7=CC=C(C=C7C(=O)[N]6(CCCc8ccc(cc8)S(=O)(=O)N)Cl)O)C"],
)
def test_unusual_smiles(smiles):
    mol = oemol_from_smiles(smiles)
    get_qed(mol)


@pytest.mark.parametrize(
    ["drug_name", "expected_qed"],
    [
        # Example drugs from
        # Values are slightly different due to differences between RDKit and paper
        ("paroxetine", 0.934),  # Was 0.935 in paper
        ("leflunomide", 0.911),  # Was 0.929 in paper
        ("granisetron", 0.926),
        ("pergolide", 0.909),  # Was 0.923 in paper
        ("molindone", 0.918),  # Was 0.920 in paper
    ],
)
@pytest.mark.parametrize("rdkit", [True, False])
@pytest.mark.parametrize("explicit_hs", [True, False])
def test_drug(drug_name, expected_qed, rdkit, explicit_hs):
    smiles = DRUG_SMILES[drug_name]
    if rdkit:
        mol = Chem.MolFromSmiles(smiles)
        if explicit_hs:
            mol = Chem.AddHs(mol)
    else:
        mol = oemol_from_smiles(smiles)
        if explicit_hs:
            oechem.OEAddExplicitHydrogens(mol)
    assert get_qed(mol) == pytest.approx(expected_qed, abs=0.0005)


@pytest.mark.parametrize(
    "smiles",
    [*DRUG_SMILES.values()],
)
def test_agrees_rdkit(smiles):
    rdkit_qed = QED.qed(Chem.MolFromSmiles(smiles))
    oemol = oemol_from_smiles(smiles)
    assert get_qed(oemol) == pytest.approx(rdkit_qed)
