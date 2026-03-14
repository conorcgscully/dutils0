import pytest
from rdkit import Chem

from chemutils.molecule.smiles import oemol_from_smiles
from chemutils.props.synthetic_accessibility import get_synthetic_accessibility_score


@pytest.mark.parametrize(
    ["smiles", "expected_sascore"],
    [
        # From supplementary information, Ertl & Schuffenhauer (2009)
        (r"COc4ccc3nc(NC(=O)CSc2nnc(c1ccccc1C)[nH]2)sc3c4", 2.4935),
        (
            r"OC8Cc7c(O)c(C2C(O)C(c1ccc(O)c(O)c1)Oc3cc(O)ccc23)c(O)c(C5C(O)C(c4ccc(O)c(O)c4)Oc6cc(O)ccc56)c7OC8c9ccc(O)c(O)c9",
            5.3203,
        ),
        (r"NC(=O)Nc1nsnc1C(=O)Nc2ccccc2", 2.3245),
        (r"C=CCn5c(=O)[nH]c(=O)c(=C4CC(c2ccc1OCOc1c2)N(c3ccccc3)N4)c5=O", 3.2639),
        (r"Oc1c(Cl)cc(Cl)cc1CNc2cccc3cn[nH]c23", 2.4007),
        (r"CC(C)C(C)C=CC(C)C1CCC3C1(C)CCC4C2(C)CCC(O)CC25CCC34OO5", 6.4954),
        (r"CC45CC(O)C1C(CC=C2CC3(CCC12)OCCO3)C4CCC56OCCO6", 5.5019),
        (r"CCc2ccc(c1ccccc1)cc2", 1.1697),
        (r"CC5C4C(CC3C2CC=C1CC(OC(C)=O)C(O)C(O)C1(C)C2CC(O)C34C)OC56CCC(=C)CO6", 5.9761),
        (r"CSc2ncnc3cn(C1OC(CO)C(O)C1O)nc23", 4.1693),
        (
            r"CCc1c(C)c2cc5nc(nc4[nH]c(cc3nc(cc1[nH]2)C(=O)C3(C)CC)c(CCC(=O)OC)c4C)C(C)(O)C5(O)CCC(=O)OC",
            6.5846,
        ),
        (r"CN(COC(C)=O)c1nc(N(C)COC(C)=O)nc(N(C)COC(C)=O)n1", 3.3497),
        (r"CSc2ccc(OCC(=O)Nc1ccc(C(C)C)cc1)cc2", 1.7294),
        (r"Cc2ccc(C(=O)Nc1ccccc1)cc2", 1.1391),
        (r"CC5CC(C)C(O)(CC4CC3OC2(CCC1(OC(C=CCCC(O)=O)CC=C1)O2)C(C)CC3O4)OC5C(Br)=C", 6.7016),
        (
            r"COc8ccc(C27C(CC1C5C(CC=C1C2c3cc(OC)ccc3O)C(=O)N(c4cccc(C(O)=O)c4)C5=O)C(=O)N(Nc6ccc(Cl)cc6Cl)C7=O)cc8",
            4.4622,
        ),
        (
            r"CC=CC(O)CC=CCC(C)C(O)CC(=O)NCC(O)C(C)C(=O)NCCCC2OC1(CCCC(CCC(C)C=C(C)C(C)O)O1)CCC2C",
            6.0551,
        ),
        (r"CCC(C)=CC(=O)OC1C(C)CC3OC1(O)C(O)C2(C)CCC(O2)C(C)(C)C=CC(C)C3=O", 6.8287),
        (r"CCC(CO)NC(=O)c2cccc(S(=O)(=O)N1CCCCCC1)c2", 2.4430),
        (r"CCCCCC1OC(=O)CCCCCCCCC=CC1=O", 3.4577),
        (r"COc1ccc(Cl)cc1", 1.0490),
        (r"CC(C)(C)C(Br)C(=O)NC(C)(C)C1CCC(C)(NC(=O)C(Br)C(C)(C)C)CC1", 4.5539),
        (r"COc2cc(CNc1ccccc1)ccc2OCC(=O)Nc3ccc(Cl)cc3", 1.7354),
        (
            r"COC4C=C(C)CC(C=CC=CC#CC1CC1Cl)OC(=O)CC3(O)CC(OC2OC(C)C(O)C(C)(O)C2OC)C(C)C(O3)C4C",
            7.3266,
        ),
        (r"CCc2ccc(OC(=O)c1ccccc1Cl)cc2", 1.5282),
        (r"COc1ccccc1c2ccccc2", 1.1772),
        (
            r"CCCC(NC(=O)C1CC2CN1C(=O)C(C(C)(C)C)NC(=O)Cc3cccc(OCCCO2)c3)C(=O)C(=O)NCC(=O)NC(C(O)=O)c4ccc(NS(N)(=O)=O)cc4",
            6.2989,
        ),
        (
            r"COC4C(O)C(C)OC(OCC3C=CC=CC(=O)C(C)CC(C)C(OC2OC(C)CC1(OC(=O)OC1C)C2O)C(C)C=CC(=O)OC3C)C4OC",
            6.9473,
        ),
        (r"CC(C)(C)c4ccc(C(=O)Nc3nc2C(CC(=O)NCC#C)C1(C)CCC(O)C(C)(CO)C1Cc2s3)cc4", 4.2412),
        (
            r"CCC7(C4OC(C3OC2(COC(c1ccc(OC)cc1)O2)C(C)CC3C)CC4C)CCC(C6(C)CCC5(CC(OCC=C)C(C)C(C(C)C(OC)C(C)C(O)=O)O5)O6)O7",
            6.9575,
        ),
        (r"O=C(OCc1ccccc1)c2ccccc2", 1.2062),
        (
            r"CC(C)CC(NC(=O)C(CC(=O)NC2OC(CO)C(OC1OC(CO)C(O)C(O)C1NC(C)=O)C(O)C2NC(C)=O)NC(=O)c3ccccc3)C(=O)NC(C(C)O)C(N)=O",
            5.3024,
        ),
        (
            r"CCCC5OC(=O)C(C)C(=O)C(C)C(OC1OC(C)CC(N(C)C)C1O)C(C)(OCC=Cc3cnc2ccc(OC)cc2c3)CC(C)C4=NCCN6C(C4C)C5(C)OC6=O",
            6.9177,
        ),
        (r"COC(=O)c1ccccc1NC(=O)CC(c2ccccc2)c3ccccc3", 1.7863),
        (r"Cc4onc5c1ncc(Cl)cc1n(C3CCCC(CNC(=O)OCc2ccccc2)C3)c(=O)c45", 3.5166),
        (r"CC(C)OCCCNC(=O)c3cc2c(=O)n(C)c1ccccc1c2n3C", 2.4205),
        (r"COC(=O)N4CCCC(N3CCC(n1c(=O)n(S(C)(=O)=O)c2ccccc12)CC3)C4", 3.0924),
        (r"Cc5c(C=NN3C(=O)C2C1CC(C=C1)C2C3=O)c4ccccc4n5Cc6ccc([N+](=O)[O-])cc6", 4.4654),
        (
            r"CCC5OC(=O)C(C)C(=O)C(C)C(OC1OC(C)CC(N(C)C)C1O)C(C)(OCC#Cc4cc(c3ccc2ccccc2n3)no4)CC(C)C(=O)C(C)C6NC(=O)OC56C",
            6.2462,
        ),
        (r"CC(=O)Nc1ccccc1NC(=O)COc2ccccc2", 1.5617),
    ],
)
@pytest.mark.parametrize("rdkit", [True, False])
def test_paper(smiles, expected_sascore, rdkit):
    mol = Chem.MolFromSmiles(smiles) if rdkit else oemol_from_smiles(smiles)
    assert get_synthetic_accessibility_score(mol) == pytest.approx(expected_sascore, abs=0.00005)
