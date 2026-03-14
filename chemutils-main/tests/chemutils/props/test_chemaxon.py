import pytest

from chemutils.chemaxon.molecule import cxmol_from_smiles
from chemutils.molecule import oemol_from_smiles
from chemutils.props import LogD, LogP, PKa, PKaConjugateAcid

CASES = [("C[C@@H](C(=O)O)N", 2.475, 9.477, -2.841, -2.844)]


@pytest.mark.parametrize(
    ["smiles", "pka", "pka_conjugate_acid", "logp", "logd"],
    CASES,
)
def test_property_smiles(smiles, pka, pka_conjugate_acid, logp, logd):
    assert PKa.get(smiles) == pytest.approx(pka, rel=1e-2)
    assert PKaConjugateAcid.get(smiles) == pytest.approx(pka_conjugate_acid, rel=1e-2)
    assert LogP.get(smiles) == pytest.approx(logp, rel=1e-2)
    assert LogD.get(smiles) == pytest.approx(logd, rel=1e-2)


@pytest.mark.parametrize(
    ["smiles", "pka", "pka_conjugate_acid", "logp", "logd"],
    CASES,
)
def test_property_oemol(smiles, pka, pka_conjugate_acid, logp, logd):
    oemol = oemol_from_smiles(smiles)
    assert PKa.get(oemol) == pytest.approx(pka, rel=1e-2)
    assert PKaConjugateAcid.get(oemol) == pytest.approx(pka_conjugate_acid, rel=1e-2)
    assert LogP.get(oemol) == pytest.approx(logp, rel=1e-2)
    assert LogD.get(oemol) == pytest.approx(logd, rel=1e-2)


@pytest.mark.parametrize(
    ["smiles", "pka", "pka_conjugate_acid", "logp", "logd"],
    CASES,
)
def test_property_cxmol(smiles, pka, pka_conjugate_acid, logp, logd):
    cxmol = cxmol_from_smiles(smiles)
    assert PKa.get(cxmol) == pytest.approx(pka, rel=1e-2)
    assert PKaConjugateAcid.get(cxmol) == pytest.approx(pka_conjugate_acid, rel=1e-2)
    assert LogP.get(cxmol) == pytest.approx(logp, rel=1e-2)
    assert LogD.get(cxmol) == pytest.approx(logd, rel=1e-2)
