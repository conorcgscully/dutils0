import pytest

from chemutils.molecule import (
    AtomProperty,
    AutomorphAtomProperties,
    AutomorphBondProperties,
    BondProperty,
    has_homomorphism,
    iterate_homomorphism_matches,
    oemol_from_smiles,
)


@pytest.mark.parametrize(
    ["smiles1", "smiles2", "num_matches"],
    [
        ["CC", "CCC", 0],
        ["CCC", "CCC", 2],
        ["c1ccccc1C", "c1ccccc1C", 2],
        ["c1ccccc1", "c1ccccc1", 12],
        ["C1CCCC1", "C1CCCC1", 10],
    ],
)
def test_homomorphism(smiles1, smiles2, num_matches):
    oemol1 = oemol_from_smiles(smiles1)
    oemol2 = oemol_from_smiles(smiles2)

    matches = list(iterate_homomorphism_matches(pattern=oemol1, target=oemol2))
    assert len(matches) == num_matches
    assert has_homomorphism(pattern=oemol1, target=oemol2) == (num_matches > 0)


def test_homomorphism_no_props():
    oemol1 = oemol_from_smiles("c1ccccc1CC(C)C")
    oemol2 = oemol_from_smiles("c1ncccc1CC(=O)O")

    matches = list(iterate_homomorphism_matches(pattern=oemol1, target=oemol2))
    assert len(matches) == 0

    matches_noprops = list(
        iterate_homomorphism_matches(
            pattern=oemol1,
            target=oemol2,
            atom_properties=AtomProperty.Nothing,
            bond_properties=BondProperty.Nothing,
        )
    )
    assert len(matches_noprops) == 4


def test_homomorphism_atomic_number():
    oemol1 = oemol_from_smiles("CCC")
    oemol2 = oemol_from_smiles("CCN")

    matches_noatomnum = list(
        iterate_homomorphism_matches(
            pattern=oemol1,
            target=oemol2,
            atom_properties=AutomorphAtomProperties & ~AtomProperty.AtomicNumber,
        )
    )
    assert len(matches_noatomnum) == 2

    matches_atomicnumber = list(iterate_homomorphism_matches(pattern=oemol1, target=oemol2))
    assert len(matches_atomicnumber) == 0


def test_homomorphism_aromatic():
    oemol1 = oemol_from_smiles("c1ccccc1")
    oemol2 = oemol_from_smiles("C1CCCCC1")

    matches_noatomnum = list(
        iterate_homomorphism_matches(
            pattern=oemol1,
            target=oemol2,
            atom_properties=AutomorphAtomProperties & ~AtomProperty.Aromatic,
            bond_properties=AutomorphBondProperties & ~BondProperty.Aromatic,
        )
    )
    assert len(matches_noatomnum) == 12

    matches_atomicnumber = list(iterate_homomorphism_matches(pattern=oemol1, target=oemol2))
    assert len(matches_atomicnumber) == 0


def test_homomorphism_charge():
    oemol1 = oemol_from_smiles("CCO")
    oemol2 = oemol_from_smiles("CC[O-]")

    matches_charge = list(
        iterate_homomorphism_matches(
            pattern=oemol1,
            target=oemol2,
        )
    )
    assert len(matches_charge) == 1

    matches_nocharge = list(
        iterate_homomorphism_matches(
            pattern=oemol1,
            target=oemol2,
            atom_properties=AutomorphAtomProperties | AtomProperty.FormalCharge,
        )
    )

    assert len(matches_nocharge) == 0
