import numpy as np
from openeye import oechem

from chemutils.maths.convert import as_coordinates, as_vectors
from chemutils.molecule import create_atom


def test_as_coordinates_np_array():
    np.testing.assert_allclose(as_coordinates(np.array([1, -0.5, 2.5])), np.array([1, -0.5, 2.5]))


def test_as_coordinates_tuple():
    np.testing.assert_allclose(as_coordinates((1, -0.5, 2.5)), np.array([1, -0.5, 2.5]))


def test_as_coordinates_nested():
    np.testing.assert_allclose(
        as_coordinates([(1, -0.5, 2.5), (-0.5, -2.1, 1), (0.5, 1.5, 0.5)]),
        np.array([[1, -0.5, 2.5], [-0.5, -2.1, 1], [0.5, 1.5, 0.5]]),
    )


def test_as_coordinates_atom():
    mol = oechem.OEGraphMol()
    atom = create_atom(mol, atomic_number=6, coordinates=(1, -0.5, 2.5))

    np.testing.assert_allclose(as_coordinates(atom), np.array([1, -0.5, 2.5]))


def test_as_coordinates_mol():
    mol = oechem.OEGraphMol()
    create_atom(mol, atomic_number=6, coordinates=(1, -0.5, 2.5))
    create_atom(mol, atomic_number=6, coordinates=(-0.5, -2.1, 1))
    create_atom(mol, atomic_number=8, coordinates=(0.5, 1.5, 0.5))

    np.testing.assert_allclose(
        as_coordinates(mol), np.array([[1, -0.5, 2.5], [-0.5, -2.1, 1], [0.5, 1.5, 0.5]])
    )


def test_as_coordinates_mol_empty():
    mol = oechem.OEGraphMol()

    np.testing.assert_allclose(as_coordinates(mol), np.array([]))


def test_as_coordinates_atom_iter():
    mol = oechem.OEGraphMol()
    create_atom(mol, atomic_number=6, coordinates=(1, -0.5, 2.5))
    create_atom(mol, atomic_number=6, coordinates=(-0.5, -2.1, 1))
    create_atom(mol, atomic_number=8, coordinates=(0.5, 1.5, 0.5))

    atom_iter = mol.GetAtoms(oechem.OEIsCarbon())
    np.testing.assert_allclose(
        as_coordinates(atom_iter), np.array([[1, -0.5, 2.5], [-0.5, -2.1, 1]])
    )


def test_as_coordinates_atom_list():
    mol = oechem.OEGraphMol()
    atom1 = create_atom(mol, atomic_number=6, coordinates=(1, -0.5, 2.5))
    create_atom(mol, atomic_number=6, coordinates=(-0.5, -2.1, 1))
    atom2 = create_atom(mol, atomic_number=8, coordinates=(0.5, 1.5, 0.5))

    atoms = [atom1, atom2]
    np.testing.assert_allclose(as_coordinates(atoms), np.array([[1, -0.5, 2.5], [0.5, 1.5, 0.5]]))


def test_as_coordinates_atom_list_empty():
    np.testing.assert_allclose(as_coordinates([]), np.array([]))


def test_as_coordinates_atom_bond_set():
    mol = oechem.OEGraphMol()
    atom1 = create_atom(mol, atomic_number=6, coordinates=(1, -0.5, 2.5))
    create_atom(mol, atomic_number=6, coordinates=(-0.5, -2.1, 1))
    atom3 = create_atom(mol, atomic_number=8, coordinates=(0.5, 1.5, 0.5))

    atombondset = oechem.OEAtomBondSet()
    atombondset.AddAtom(atom1)
    atombondset.AddAtom(atom3)

    np.testing.assert_allclose(
        as_coordinates(atombondset), np.array([[1, -0.5, 2.5], [0.5, 1.5, 0.5]])
    )


def test_as_coordinates_atom_bond_set_empty():
    atombondset = oechem.OEAtomBondSet()

    np.testing.assert_allclose(as_coordinates(atombondset), np.array([]))


def test_as_vectors_np_array():
    np.testing.assert_allclose(as_vectors(np.array([1, -0.5, 2.5])), np.array([1, -0.5, 2.5]))


def test_as_vectors_tuple():
    np.testing.assert_allclose(as_vectors((1, -0.5, 2.5)), np.array([1, -0.5, 2.5]))


def test_as_vectors_nested():
    np.testing.assert_allclose(
        as_vectors([(1, -0.5, 2.5), (-0.5, -2.1, 1), (0.5, 1.5, 0.5)]),
        np.array([[1, -0.5, 2.5], [-0.5, -2.1, 1], [0.5, 1.5, 0.5]]),
    )


def test_as_vectors_bond():
    mol = oechem.OEGraphMol()
    atom1 = create_atom(mol, atomic_number=6, coordinates=(1, -0.5, 2.5))
    atom2 = create_atom(mol, atomic_number=6, coordinates=(-0.5, -2.1, 1))
    bond = mol.NewBond(atom1, atom2)

    np.testing.assert_allclose(as_vectors(bond), np.array([-1.5, -1.6, -1.5]))


def test_as_vectors_bond_iter():
    mol = oechem.OEGraphMol()
    atom1 = create_atom(mol, atomic_number=6, coordinates=(1, -0.5, 2.5))
    atom2 = create_atom(mol, atomic_number=6, coordinates=(-0.5, -2.1, 1))
    mol.NewBond(atom1, atom2)
    atom3 = create_atom(mol, atomic_number=8, coordinates=(0.5, 1.5, 0.5))
    mol.NewBond(atom2, atom3)

    bond_iter = mol.GetBonds()
    np.testing.assert_allclose(
        as_vectors(bond_iter), np.array([[-1.5, -1.6, -1.5], [1.0, 3.6, -0.5]])
    )


def test_as_vectors_bond_iter_empty():
    mol = oechem.OEGraphMol()

    bond_iter = mol.GetBonds()
    np.testing.assert_allclose(as_vectors(bond_iter), np.array([]))


def test_as_vectors_atom_bond_set():
    mol = oechem.OEGraphMol()
    atom1 = create_atom(mol, atomic_number=6, coordinates=(1, -0.5, 2.5))
    atom2 = create_atom(mol, atomic_number=6, coordinates=(-0.5, -2.1, 1))
    bond1 = mol.NewBond(atom1, atom2)
    atom3 = create_atom(mol, atomic_number=8, coordinates=(0.5, 1.5, 0.5))
    bond2 = mol.NewBond(atom2, atom3)

    atombondset = oechem.OEAtomBondSet()
    atombondset.AddBond(bond1)
    atombondset.AddBond(bond2)

    np.testing.assert_allclose(
        as_vectors(atombondset), np.array([[-1.5, -1.6, -1.5], [1.0, 3.6, -0.5]])
    )


def test_as_vectors_atom_bond_set_empty():
    atombondset = oechem.OEAtomBondSet()

    np.testing.assert_allclose(as_vectors(atombondset), np.array([]))
