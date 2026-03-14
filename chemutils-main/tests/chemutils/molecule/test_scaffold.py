import pytest
from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold

from chemutils.molecule import (
    canonicalize_smiles,
    get_bemis_murcko_scaffold,
    get_generic_bemis_murcko_scaffold,
    oemol_from_smiles,
    smiles_from_oemol,
)
from chemutils.rdmol_from_oemol import rdmol_from_oemol

CASES = [
    # Benzene
    {
        "smiles": "c1ccccc1",
        "placeholder": "c1ccccc1",
        "remove": "c1ccccc1",
        "keep": "c1ccccc1",
        "generic": "C1CCCCC1",
    },
    # Quinone
    {
        "smiles": "C1=CC(=O)C=CC1=O",
        "placeholder": "*=C1C=CC(=*)C=C1",
        "remove": "C1C=CCC=C1",
        "keep": "C1=CC(=O)C=CC1=O",
        "generic": "C1CCCCC1",
    },
    # Quinone methide
    {
        "smiles": "C1=CC(=O)C=CC1=C",
        "placeholder": "*=C1C=CC(=*)C=C1",
        "remove": "C1C=CCC=C1",
        "keep": "C=C1C=CC(=O)C=C1",
        "generic": "C1CCCCC1",
    },
    # Thioridazine
    {
        "smiles": "CN1CCCCC1CCN2c3ccccc3Sc4c2cc(cc4)SC",
        "placeholder": "c1ccc2c(c1)N(c3ccccc3S2)CCC4CCCCN4",
        "remove": "c1ccc2c(c1)N(c3ccccc3S2)CCC4CCCCN4",
        "keep": "c1ccc2c(c1)N(c3ccccc3S2)CCC4CCCCN4",
        "generic": "C1CCC(CC1)CCC2C3CCCCC3CC4C2CCCC4",
    },
    # Caffeine
    {
        "smiles": "Cn1cnc2c1c(=O)n(c(=O)n2C)C",
        "placeholder": "*=c1c2c([nH]c(=*)[nH]1)nc[nH]2",
        "remove": "c1[nH]c2c(n1)NCNC2",
        "keep": "c1[nH]c2c(n1)[nH]c(=O)[nH]c2=O",
        "generic": "C1CCC2CCCC2C1",
    },
    # Example case from https://github.com/rdkit/rdkit/discussions/6844
    {
        "smiles": "Fc1ccc(OC(=O)C2CS(=O)(=NC)C2)cc1",
        "placeholder": "*=C(C1CS(=*)(=*)C1)Oc2ccccc2",
        "remove": "c1ccc(cc1)OCC2CSC2",
        "keep": "c1ccc(cc1)OC(=O)C2CS(=N)(=O)C2",
        "generic": "C1CCC(CC1)CCC2CCC2",
    },
    # Phosphate test
    {
        "smiles": "C1CC1OP(=O)(OC1CC1)OC1CC1",
        "placeholder": "*=P(OC1CC1)(OC2CC2)OC3CC3",
        "remove": "C1CC1OP(OC2CC2)OC3CC3",
        "keep": "C1CC1OP(=O)(OC2CC2)OC3CC3",
        "generic": "C1CC1CC(CC2CC2)CC3CC3",
    },
]


@pytest.mark.parametrize(
    ["smiles", "scaffold"], [(case["smiles"], case["placeholder"]) for case in CASES]
)
def test_get_bemis_murcko_scaffold(smiles, scaffold):
    oemol = oemol_from_smiles(smiles)
    scaffold_oemol = get_bemis_murcko_scaffold(oemol)
    assert smiles_from_oemol(scaffold_oemol) == scaffold


@pytest.mark.parametrize(
    ["smiles", "scaffold"], [(case["smiles"], case["remove"]) for case in CASES]
)
def test_get_bemis_murcko_scaffold_remove(smiles, scaffold):
    oemol = oemol_from_smiles(smiles)
    scaffold_oemol = get_bemis_murcko_scaffold(oemol, external_double_bond_behaviour="remove")
    assert smiles_from_oemol(scaffold_oemol) == scaffold


@pytest.mark.parametrize(["smiles", "scaffold"], [(case["smiles"], case["keep"]) for case in CASES])
def test_get_bemis_murcko_scaffold_keep(smiles, scaffold):
    oemol = oemol_from_smiles(smiles)
    scaffold_oemol = get_bemis_murcko_scaffold(oemol, external_double_bond_behaviour="keep")
    assert smiles_from_oemol(scaffold_oemol) == scaffold


@pytest.mark.parametrize(
    ["smiles", "scaffold"], [(case["smiles"], case["generic"]) for case in CASES]
)
def test_get_generic_bemis_murcko_scaffold(smiles, scaffold):
    oemol = oemol_from_smiles(smiles)
    scaffold_oemol = get_generic_bemis_murcko_scaffold(oemol)
    assert smiles_from_oemol(scaffold_oemol) == scaffold


@pytest.mark.parametrize("smiles", [case["smiles"] for case in CASES])
def test_get_bemis_murcko_scaffold_keep_matches_rdkit(smiles):
    oemol = oemol_from_smiles(smiles)
    scaffold_oemol = get_bemis_murcko_scaffold(oemol, external_double_bond_behaviour="keep")
    chemutils_smiles = smiles_from_oemol(scaffold_oemol)
    rdmol = rdmol_from_oemol(scaffold_oemol)
    scaffold = MurckoScaffold.GetScaffoldForMol(rdmol)
    rdkit_smiles = canonicalize_smiles(Chem.MolToSmiles(scaffold))
    assert chemutils_smiles == rdkit_smiles


@pytest.mark.parametrize("smiles", [case["smiles"] for case in CASES])
def test_get_generic_bemis_murcko_scaffold_matches_rdkit(smiles):
    oemol = oemol_from_smiles(smiles)
    scaffold_oemol = get_generic_bemis_murcko_scaffold(oemol)
    chemutils_smiles = smiles_from_oemol(scaffold_oemol)
    rdmol = rdmol_from_oemol(scaffold_oemol)
    scaffold = MurckoScaffold.GetScaffoldForMol(rdmol)
    scaffold = MurckoScaffold.MakeScaffoldGeneric(scaffold)
    rdkit_smiles = canonicalize_smiles(Chem.MolToSmiles(scaffold))
    assert chemutils_smiles == rdkit_smiles
