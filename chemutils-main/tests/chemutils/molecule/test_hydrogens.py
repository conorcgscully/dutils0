import pytest
from openeye import oechem

from chemutils.molecule import oemol_from_smiles, smiles_from_oemol
from chemutils.molecule.conformer import generate_single_conformer
from chemutils.molecule.hydrogens import (
    DeprotonationError,
    ProtonationError,
    deprotonate,
    make_hydrogens_implicit,
    protonate,
)


@pytest.mark.parametrize(
    ["smiles", "num_hydrogens_isotopic"],
    [
        ("C", 0),
        ("[CH4]", 0),
        ("[CH3][H]", 0),
        ("[CH3][2H]", 1),
    ],
)
def test_make_hydrogens_implicit(smiles, num_hydrogens_isotopic):
    oemol = oemol_from_smiles(smiles)
    make_hydrogens_implicit(oemol, remove_isotopic_hydrogens=True)
    assert len([atom for atom in oemol.GetAtoms() if atom.GetAtomicNum() == 1]) == 0
    oemol = oemol_from_smiles(smiles)
    make_hydrogens_implicit(oemol, remove_isotopic_hydrogens=False)
    assert (
        len([atom for atom in oemol.GetAtoms() if atom.GetAtomicNum() == 1])
        == num_hydrogens_isotopic
    )


def test_protonate_implicit():
    oemol = oemol_from_smiles("N")
    nitrogen: oechem.OEAtomBase = oemol.GetAtom(oechem.OEHasAtomicNum(7))

    assert oemol.NumAtoms() == 1
    assert nitrogen.GetImplicitHCount() == 3
    assert nitrogen.GetFormalCharge() == 0
    assert smiles_from_oemol(oemol) == "N"

    protonate(nitrogen)

    assert oemol.NumAtoms() == 1
    assert nitrogen.GetImplicitHCount() == 4
    assert nitrogen.GetFormalCharge() == 1
    assert smiles_from_oemol(oemol) == "[NH4+]"


def test_protonate_explicit_no_coordinates():
    oemol = oemol_from_smiles("N")
    oechem.OEAddExplicitHydrogens(oemol, False, False)
    nitrogen: oechem.OEAtomBase = oemol.GetAtom(oechem.OEHasAtomicNum(7))

    assert oemol.NumAtoms() == 4
    assert nitrogen.GetImplicitHCount() == 0
    assert nitrogen.GetExplicitHCount() == 3
    assert nitrogen.GetFormalCharge() == 0
    assert smiles_from_oemol(oemol) == "N"

    protonate(nitrogen)

    assert oemol.NumAtoms() == 5
    assert nitrogen.GetImplicitHCount() == 0
    assert nitrogen.GetExplicitHCount() == 4
    assert nitrogen.GetFormalCharge() == 1
    assert smiles_from_oemol(oemol) == "[NH4+]"
    assert oemol.GetDimension() == 0
    for atom in oemol.GetAtoms():
        assert oemol.GetCoords(atom) == (0, 0, 0)


def test_protonate_explicit_coordinates():
    oemol = oemol_from_smiles("N")
    oechem.OEAddExplicitHydrogens(oemol, False, False)
    oemol = generate_single_conformer(oemol)
    nitrogen: oechem.OEAtomBase = oemol.GetAtom(oechem.OEHasAtomicNum(7))

    assert oemol.NumAtoms() == 4
    assert oemol.NumBonds() == 3
    assert nitrogen.GetImplicitHCount() == 0
    assert nitrogen.GetExplicitHCount() == 3
    assert nitrogen.GetFormalCharge() == 0
    assert smiles_from_oemol(oemol) == "N"

    protonate(nitrogen)

    assert oemol.NumAtoms() == 5
    assert oemol.NumBonds() == 4
    assert nitrogen.GetImplicitHCount() == 0
    assert nitrogen.GetExplicitHCount() == 4
    assert nitrogen.GetFormalCharge() == 1
    assert smiles_from_oemol(oemol) == "[NH4+]"
    assert oemol.GetDimension() == 3
    for atom in oemol.GetAtoms():
        assert oemol.GetCoords(atom) != (0, 0, 0)


def test_protonate_hydrogen():
    oemol = oemol_from_smiles("N")
    oechem.OEAddExplicitHydrogens(oemol, False, False)
    oemol = generate_single_conformer(oemol)
    hydrogen: oechem.OEAtomBase = oemol.GetAtom(oechem.OEHasAtomicNum(1))

    with pytest.raises(ProtonationError) as exc:
        protonate(hydrogen)

    assert str(exc.value) == "Cannot protonate a hydrogen atom."


def test_protonate_mixed_implicit_explicit():
    oemol = oemol_from_smiles("[NH2]")
    oechem.OEAddExplicitHydrogens(oemol, False, False)
    nitrogen: oechem.OEAtomBase = oemol.GetAtom(oechem.OEHasAtomicNum(7))
    nitrogen.SetImplicitHCount(1)

    assert smiles_from_oemol(oemol) == "N"

    with pytest.raises(ProtonationError) as exc:
        protonate(nitrogen)

    assert str(exc.value) == "Cannot protonate an atom with both implicit and explicit hydrogens."


def test_deprotonate_explicit_h():
    oemol = oemol_from_smiles("O")
    oechem.OEAddExplicitHydrogens(oemol, False, False)
    oxygen: oechem.OEAtomBase = oemol.GetAtom(oechem.OEHasAtomicNum(8))
    hydrogen: oechem.OEAtomBase = oemol.GetAtom(oechem.OEHasAtomicNum(1))

    assert oemol.NumAtoms() == 3
    assert oemol.NumBonds() == 2
    assert oxygen.GetImplicitHCount() == 0
    assert oxygen.GetExplicitHCount() == 2
    assert oxygen.GetFormalCharge() == 0
    assert smiles_from_oemol(oemol) == "O"

    deprotonate(hydrogen)

    assert oemol.NumAtoms() == 2
    assert oemol.NumBonds() == 1
    assert oxygen.GetImplicitHCount() == 0
    assert oxygen.GetExplicitHCount() == 1
    assert oxygen.GetFormalCharge() == -1
    assert smiles_from_oemol(oemol) == "[OH-]"


def test_deprotonate_explicit_h_multiple_neighbours():
    oemol = oemol_from_smiles("[O][H][O]")
    oechem.OEAddExplicitHydrogens(oemol, False, False)
    hydrogen: oechem.OEAtomBase = oemol.GetAtom(oechem.OEHasAtomicNum(1))

    with pytest.raises(DeprotonationError) as exc:
        deprotonate(hydrogen)

    assert str(exc.value) == "Cannot deprotonate by removing a hydrogen atom with two neighbours."


def test_deprotonate_no_hydrogens():
    oemol = oemol_from_smiles("[O]")
    oxygen: oechem.OEAtomBase = oemol.GetAtom(oechem.OEHasAtomicNum(8))

    with pytest.raises(DeprotonationError) as exc:
        deprotonate(oxygen)

    assert str(exc.value) == "Cannot deprotonate an atom with no hydrogens."


def test_deprotonate_one_explicit():
    oemol = oemol_from_smiles("Cl")
    oechem.OEAddExplicitHydrogens(oemol, False, False)
    chlorine: oechem.OEAtomBase = oemol.GetAtom(oechem.OEHasAtomicNum(17))

    assert oemol.NumAtoms() == 2
    assert oemol.NumBonds() == 1
    assert chlorine.GetImplicitHCount() == 0
    assert chlorine.GetExplicitHCount() == 1
    assert chlorine.GetFormalCharge() == 0

    deprotonate(chlorine)

    assert oemol.NumAtoms() == 1
    assert oemol.NumBonds() == 0
    assert chlorine.GetImplicitHCount() == 0
    assert chlorine.GetExplicitHCount() == 0
    assert chlorine.GetFormalCharge() == -1


def test_deprotonate_multiple_explicit():
    oemol = oemol_from_smiles("O")
    oechem.OEAddExplicitHydrogens(oemol, False, False)
    oxygen: oechem.OEAtomBase = oemol.GetAtom(oechem.OEHasAtomicNum(8))

    with pytest.raises(DeprotonationError) as exc:
        deprotonate(oxygen)

    assert str(exc.value) == "Cannot deprotonate an atom with more than one explicit hydrogens."


def test_deprotonate_multiple_implicit():
    oemol = oemol_from_smiles("O")
    oxygen: oechem.OEAtomBase = oemol.GetAtom(oechem.OEHasAtomicNum(8))

    assert oemol.NumAtoms() == 1
    assert oemol.NumBonds() == 0
    assert oxygen.GetImplicitHCount() == 2
    assert oxygen.GetExplicitHCount() == 0
    assert oxygen.GetFormalCharge() == 0

    deprotonate(oxygen)

    assert oemol.NumAtoms() == 1
    assert oemol.NumBonds() == 0
    assert oxygen.GetImplicitHCount() == 1
    assert oxygen.GetExplicitHCount() == 0
    assert oxygen.GetFormalCharge() == -1


def test_deprotonate_mixed_implicit_explicit():
    oemol = oemol_from_smiles("[NH2]")
    oechem.OEAddExplicitHydrogens(oemol, False, False)
    nitrogen: oechem.OEAtomBase = oemol.GetAtom(oechem.OEHasAtomicNum(7))
    nitrogen.SetImplicitHCount(1)

    assert smiles_from_oemol(oemol) == "N"

    with pytest.raises(DeprotonationError) as exc:
        deprotonate(nitrogen)

    assert str(exc.value) == "Cannot deprotonate an atom with both implicit and explicit hydrogens."
