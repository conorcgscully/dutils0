import numpy as np
import pytest
from openeye import oechem

from chemutils.molecule import oemol_from_smiles
from chemutils.molecule.coordinates import get_coordinates, set_coordinates
from chemutils.molecule.embed import embed_molecule_3d


def get_bond_length(bond: oechem.OEBondBase):
    x1, y1, z1 = bond.GetParent().GetCoords(bond.GetBgn())
    x2, y2, z2 = bond.GetParent().GetCoords(bond.GetEnd())
    return np.linalg.norm([x2 - x1, y2 - y1, z2 - z1])


@pytest.mark.parametrize(
    ["smiles", "bond_length"],
    [
        ("CC", 1.54),
        ("C=C", 1.33),
        ("C#C", 1.21),
        ("c1ccccc1", 1.39),
    ],
)
@pytest.mark.parametrize("forcefield", [None, "UFF", "MMFF94"])
def test_embed_molecule_3d(smiles, bond_length, forcefield):
    mol = oemol_from_smiles(smiles)
    embed_molecule_3d(mol, seed=0, optimize=forcefield is not None, forcefield=forcefield)
    bond = next(iter(mol.GetBonds()))
    assert get_bond_length(bond) == pytest.approx(bond_length, abs=0.05)


@pytest.mark.parametrize(
    ["smiles", "timeout", "success"],
    [
        ("CC", 1, True),
        (
            "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC",
            1,
            False,
        ),
    ],
)
def test_embed_molecule_3d_timeout(smiles, timeout, success):
    if not success:
        with pytest.raises(ValueError):
            mol = oemol_from_smiles(smiles)
            embed_molecule_3d(mol, seed=0, timeout=timeout)
    else:
        mol = oemol_from_smiles(smiles)
        embed_molecule_3d(mol, seed=0, timeout=timeout)
        assert get_coordinates(mol) is not None


@pytest.mark.parametrize("smiles", ["C", "CC=C", "c1ccccc1"])
def test_embed_molecule_3d_reproducible(smiles):
    mol = oemol_from_smiles(smiles)
    embed_molecule_3d(mol, seed=0)
    mol2 = oemol_from_smiles(smiles)
    embed_molecule_3d(mol2, seed=0)
    np.testing.assert_allclose(get_coordinates(mol), get_coordinates(mol2), atol=0.001)


# Test constraining a whole molecule
@pytest.mark.parametrize("forcefield", ["UFF", "MMFF94"])
def test_embed_molecule_3d_constraint(forcefield):
    mol = oemol_from_smiles("CC=C")
    embed_molecule_3d(mol, seed=0)

    coordinates = get_coordinates(mol)
    coordinates[:, 0] += 20.0
    set_coordinates(mol, coordinates=coordinates)

    mol2 = oemol_from_smiles("CC=C")
    embed_molecule_3d(mol2, seed=0, constraint=mol, forcefield=forcefield)

    np.testing.assert_allclose(get_coordinates(mol), get_coordinates(mol2), atol=0.1)


# Test constraining a subset of a molecule
@pytest.mark.parametrize("forcefield", ["UFF", "MMFF94"])
def test_embed_molecule_3d_constraint_subset(forcefield):
    mol = oemol_from_smiles("CC=C")
    embed_molecule_3d(mol, seed=0, optimize=True)

    coordinates = get_coordinates(mol)
    coordinates[:, 0] += 20.0
    set_coordinates(mol, coordinates=coordinates)
    mol2 = oemol_from_smiles("CC=CO")
    embed_molecule_3d(mol2, seed=0, constraint=mol, forcefield=forcefield)

    # Check shared body is constrained
    np.testing.assert_allclose(get_coordinates(mol)[:3], get_coordinates(mol2)[:3], atol=0.1)
    # Check that the unconstrained CO bond is reasonable
    assert get_bond_length(list(mol2.GetBonds())[-1]) == pytest.approx(1.39, abs=0.1)


# Test trying to constrain to a molecule that is not a subset
@pytest.mark.parametrize("forcefield", ["UFF", "MMFF94"])
def test_embed_molecule_3d_constraint_not_subset(forcefield):
    mol = oemol_from_smiles("CC=CC")
    embed_molecule_3d(mol, seed=0)

    mol2 = oemol_from_smiles("C=O")
    with pytest.raises(ValueError):
        embed_molecule_3d(mol2, seed=0, constraint=mol, forcefield=forcefield)


def test_embed_molecule_3d_can_change_chirality():
    # Check that [C@] is conserved under embedding
    mol = oemol_from_smiles("[C@](F)(Cl)(Br)I")
    embed_molecule_3d(mol, seed=0)

    mol_copy = mol.CreateCopy()
    assert oechem.OE3DToInternalStereo(mol_copy)
    assert oechem.OEMolToSmiles(mol_copy) == "[C@](F)(Cl)(Br)I"

    # Check that [C@@] is conserved under embedding
    mol2 = oemol_from_smiles("[C@@](F)(Cl)(Br)I")
    embed_molecule_3d(mol2, seed=0)

    mol2_copy = mol2.CreateCopy()
    assert oechem.OE3DToInternalStereo(mol2_copy)
    assert oechem.OEMolToSmiles(mol2_copy) == "[C@@](F)(Cl)(Br)I"

    # Check that [C@@] can be constrained to become [C@]
    mol3 = oemol_from_smiles("[C@@](F)(Cl)(Br)I")
    embed_molecule_3d(mol3, seed=0, constraint=mol)

    mol3_copy = mol3.CreateCopy()
    assert oechem.OE3DToInternalStereo(mol3_copy)
    assert oechem.OEMolToSmiles(mol3_copy) == "[C@](F)(Cl)(Br)I"

    # Check that [C@] can be constrained to become [C@@]
    mol4 = oemol_from_smiles("[C@](F)(Cl)(Br)I")
    embed_molecule_3d(mol4, seed=0, constraint=mol2)

    mol4_copy = mol4.CreateCopy()
    assert oechem.OE3DToInternalStereo(mol4_copy)
    assert oechem.OEMolToSmiles(mol4_copy) == "[C@@](F)(Cl)(Br)I"


def test_embed_molecule_3d_multiple():
    smiles = "c1ccccc1.c1ccccc1"
    oemol = oemol_from_smiles(smiles)
    embed_molecule_3d(oemol, seed=0, optimize=True)

    for atom_idx1 in range(oemol.NumAtoms()):
        for atom_idx2 in range(atom_idx1):
            atom_distance = oechem.OEGetDistance(
                oemol,
                oemol.GetAtom(oechem.OEHasAtomIdx(atom_idx1)),
                oemol.GetAtom(oechem.OEHasAtomIdx(atom_idx2)),
            )
            assert atom_distance >= 0.5 and atom_distance < 20
