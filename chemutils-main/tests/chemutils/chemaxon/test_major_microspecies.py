from dataclasses import dataclass

import pytest

from chemutils.chemaxon import get_major_microspecies
from chemutils.molecule import oemol_from_smiles, smiles_from_oemol
from chemutils.molecule.equal import AllAtomProperties, AtomProperty, assert_mol_equal


@dataclass
class Case:
    smiles: str
    acidic: str
    neutral: str
    basic: str


CASES = [
    Case(smiles="CC(=O)O", acidic="CC(=O)O", neutral="CC(=O)[O-]", basic="CC(=O)[O-]"),
    Case(smiles="CC(=O)[O-]", acidic="CC(=O)O", neutral="CC(=O)[O-]", basic="CC(=O)[O-]"),
    Case(smiles="CCN", acidic="CC[NH3+]", neutral="CC[NH3+]", basic="CCN"),
    Case(smiles="CC[NH3+]", acidic="CC[NH3+]", neutral="CC[NH3+]", basic="CCN"),
    Case(
        smiles="C[C@@H](C(=O)O)N",
        acidic="C[C@@H](C(=O)O)[NH3+]",
        neutral="C[C@@H](C(=O)[O-])[NH3+]",
        basic="C[C@@H](C(=O)[O-])N",
    ),
]


@pytest.mark.parametrize(["smiles", "expected"], [(case.smiles, case.neutral) for case in CASES])
def test_get_major_microspecies_neutral(smiles, expected):
    assert get_major_microspecies(smiles) == expected


@pytest.mark.parametrize(["smiles", "expected"], [(case.smiles, case.acidic) for case in CASES])
def test_get_major_microspecies_acidic(smiles, expected):
    assert get_major_microspecies(smiles, pH=0) == expected


@pytest.mark.parametrize(["smiles", "expected"], [(case.smiles, case.basic) for case in CASES])
def test_get_major_microspecies_basic(smiles, expected):
    assert get_major_microspecies(smiles, pH=14) == expected


@pytest.mark.parametrize(["smiles", "expected"], [(case.smiles, case.neutral) for case in CASES])
def test_get_major_microspecies_neutral_oemol(smiles, expected):
    oemol = oemol_from_smiles(smiles)
    oemol_microspecies = get_major_microspecies(oemol, pH=7.4)
    assert smiles_from_oemol(oemol_microspecies) == expected
    assert_mol_equal(
        oemol,
        oemol_microspecies,
        atom_properties=AllAtomProperties
        & ~AtomProperty.Degree
        & ~AtomProperty.ImplicitHCount
        & ~AtomProperty.FormalCharge
        & ~AtomProperty.Valence,
    )


@pytest.mark.parametrize(["smiles", "expected"], [(case.smiles, case.acidic) for case in CASES])
def test_get_major_microspecies_acidic_oemol(smiles, expected):
    oemol = oemol_from_smiles(smiles)
    oemol_microspecies = get_major_microspecies(oemol, pH=0)
    assert smiles_from_oemol(oemol_microspecies) == expected
    assert_mol_equal(
        oemol,
        oemol_microspecies,
        atom_properties=AllAtomProperties
        & ~AtomProperty.Degree
        & ~AtomProperty.ImplicitHCount
        & ~AtomProperty.FormalCharge
        & ~AtomProperty.Valence,
    )


@pytest.mark.parametrize(["smiles", "expected"], [(case.smiles, case.basic) for case in CASES])
def test_get_major_microspecies_basic_oemol(smiles, expected):
    oemol = oemol_from_smiles(smiles)
    oemol_microspecies = get_major_microspecies(oemol, pH=14)
    assert smiles_from_oemol(oemol_microspecies) == expected
    assert_mol_equal(
        oemol,
        oemol_microspecies,
        atom_properties=AllAtomProperties
        & ~AtomProperty.Degree
        & ~AtomProperty.ImplicitHCount
        & ~AtomProperty.FormalCharge
        & ~AtomProperty.Valence,
    )


@pytest.mark.parametrize("smiles", [case.smiles for case in CASES])
@pytest.mark.parametrize("pH", [0, 7.4, 14])
def test_get_major_microspecies_atomic(smiles, pH):
    major_microspecies = get_major_microspecies(smiles, pH=pH)
    assert major_microspecies == get_major_microspecies(major_microspecies, pH=pH)


@pytest.mark.parametrize("smiles", [case.smiles for case in CASES])
def test_get_major_microspecies_preserves_metadata(smiles):
    oemol = oemol_from_smiles(smiles)
    list(oemol.GetAtoms())[1].SetName("Test")
    oemol_microspecies = get_major_microspecies(oemol, pH=7.4)
    assert list(oemol_microspecies.GetAtoms())[1].GetName() == "Test"
