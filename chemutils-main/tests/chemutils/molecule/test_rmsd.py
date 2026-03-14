import numpy as np
import pytest
from openeye import oechem

from chemutils.molecule import get_coordinates, oemol_from_smiles, set_coordinates
from chemutils.molecule.rmsd import (
    get_optimal_rmsd_and_transformation_between_atoms,
    get_optimal_rmsd_and_transformation_between_molecules,
    get_optimal_rmsd_and_transformation_between_points,
    get_optimal_rmsd_and_transformation_for_atom_mapping,
    get_optimal_rmsd_and_transformation_for_match,
    get_rmsd_between_atoms,
    get_rmsd_between_molecules,
    get_rmsd_between_points,
    get_rmsd_for_atom_mapping,
    get_rmsd_for_match,
)
from chemutils.molecule.transform import transform_molecule, transform_points


def test_get_rmsd_between_points():
    points1 = [[0, 0, 0], [1, 0, 0]]
    points2 = [[1, 0, 0], [1, 1, 0]]

    assert get_rmsd_between_points(points1, points2) == 1.0
    assert get_rmsd_between_points(points2, points1) == 1.0

    # Transform from points1 to points2
    rmsd1, transform1 = get_optimal_rmsd_and_transformation_between_points(
        src=points1, dest=points2
    )
    assert rmsd1 == 0.0
    np.testing.assert_allclose(
        transform1, np.array([[0, -1, 0, 1], [1, 0, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]]), atol=0.001
    )
    # Check applying transformation to points1 gives rmsd
    assert get_rmsd_between_points(
        transform_points(points1, transformation=transform1), points2
    ) == pytest.approx(rmsd1)

    # Transform from points2 to points1
    rmsd2, transform2 = get_optimal_rmsd_and_transformation_between_points(
        src=points2, dest=points1
    )
    assert rmsd2 == 0.0
    np.testing.assert_allclose(
        transform2, np.array([[0, 1, 0, 0], [1, 0, 0, -1], [0, 0, -1, 0], [0, 0, 0, 1]]), atol=0.001
    )

    # Check applying transformation to points1 gives rmsd
    assert get_rmsd_between_points(
        transform_points(points2, transformation=transform2), points1
    ) == pytest.approx(rmsd2)


def test_get_rmsd_between_molecules_same_order():
    mol1 = oemol_from_smiles("CCC(=O)O")
    set_coordinates(
        mol1,
        np.loadtxt(
            "tests/chemutils/molecule/get_rmsd_between_molecules_same_order_mol1_coords.csv",
            delimiter=",",
        ),
    )
    mol2 = oemol_from_smiles("CCC(=O)O")
    set_coordinates(
        mol2,
        np.loadtxt(
            "tests/chemutils/molecule/get_rmsd_between_molecules_same_order_mol2_coords.csv",
            delimiter=",",
        ),
    )

    assert get_rmsd_between_molecules(mol1, mol2) == pytest.approx(0.908, abs=0.001)
    # Same RMSD as doing points directly
    assert get_rmsd_between_points(get_coordinates(mol1), get_coordinates(mol2)) == pytest.approx(
        0.908, abs=0.001
    )

    rmsd1, transform1 = get_optimal_rmsd_and_transformation_between_molecules(src=mol1, dest=mol2)
    assert rmsd1 == pytest.approx(0.482, abs=0.001)
    assert get_rmsd_between_molecules(
        transform_molecule(mol1, transformation=transform1), mol2
    ) == pytest.approx(rmsd1)

    rmsd2, transform2 = get_optimal_rmsd_and_transformation_between_molecules(src=mol2, dest=mol1)
    assert rmsd2 == pytest.approx(0.482, abs=0.001)
    assert get_rmsd_between_molecules(
        transform_molecule(mol2, transformation=transform2), mol1
    ) == pytest.approx(rmsd2)


def test_get_rmsd_between_molecules_different_order():
    mol1 = oemol_from_smiles("CCC(=O)O")
    set_coordinates(
        mol1,
        np.loadtxt(
            "tests/chemutils/molecule/get_rmsd_between_molecules_different_order_mol1_coords.csv",
            delimiter=",",
        ),
    )
    mol2 = oemol_from_smiles("O=C(O)CC")
    set_coordinates(
        mol2,
        np.loadtxt(
            "tests/chemutils/molecule/get_rmsd_between_molecules_different_order_mol2_coords.csv",
            delimiter=",",
        ),
    )

    assert get_rmsd_between_molecules(mol1, mol2) == pytest.approx(0.496, abs=0.001)
    # Molecule is in a different order, so pure RMSD of coordinates is wrong
    assert get_rmsd_between_points(get_coordinates(mol1), get_coordinates(mol2)) == pytest.approx(
        2.462, abs=0.001
    )


def test_get_rmsd_between_atoms():
    mol1 = oemol_from_smiles("CCC(=O)O")
    set_coordinates(
        mol1,
        np.loadtxt(
            "tests/chemutils/molecule/get_rmsd_between_atoms_mol1_coords.csv", delimiter=","
        ),
    )
    mol2 = oemol_from_smiles("CCC(=O)O")
    set_coordinates(
        mol2,
        np.loadtxt(
            "tests/chemutils/molecule/get_rmsd_between_atoms_mol2_coords.csv", delimiter=","
        ),
    )

    # Atom lists for carboyxlic acid
    atoms1 = list(mol1.GetAtoms())[2:5]
    atoms2 = list(mol2.GetAtoms())[2:5]

    # Atom mapping for carboxylic acid
    mapping = dict(zip(atoms1, atoms2, strict=False))

    # Match for carboxylic acid
    match = oechem.OEMatch()
    for atom1, atom2 in mapping.items():
        match.AddPair(atom1, atom2)

    # All three methods should return the same result
    assert get_rmsd_between_atoms(atoms1, atoms2) == pytest.approx(1.451, abs=0.001)
    assert get_rmsd_for_atom_mapping(mapping) == pytest.approx(1.451, abs=0.001)
    assert get_rmsd_for_match(match) == pytest.approx(1.451, abs=0.001)

    expected_transform = [
        [0.749, 0.354, -0.560, 0.170],
        [0.485, 0.283, 0.828, -0.311],
        [0.451, -0.892, 0.040, -0.261],
        [0.000, 0.000, 0.000, 1.000],
    ]

    rmsd_atoms, transform_atoms = get_optimal_rmsd_and_transformation_between_atoms(
        src=atoms1, dest=atoms2
    )
    assert rmsd_atoms == pytest.approx(0.049, abs=0.001)
    np.testing.assert_allclose(transform_atoms, expected_transform, atol=0.001)

    rmsd_mapping, transform_mapping = get_optimal_rmsd_and_transformation_for_atom_mapping(mapping)
    assert rmsd_mapping == pytest.approx(0.049, abs=0.001)
    np.testing.assert_allclose(transform_mapping, expected_transform, atol=0.001)

    rmsd_match, transform_match = get_optimal_rmsd_and_transformation_for_match(match)
    assert rmsd_match == pytest.approx(0.049, abs=0.001)
    np.testing.assert_allclose(transform_match, expected_transform, atol=0.001)
