import numpy as np
import numpy.testing as npt
import pytest
from openeye import oechem

from chemutils.molecule import standardise_oemol, std_oemol_from_smiles
from chemutils.molecule.cip.assignment import TetrahedralCIPAssignment
from chemutils.molecule.conformer import generate_single_conformer
from chemutils.molecule.coordinates import get_coordinates, set_coordinates
from chemutils.molecule.tetrahedron import (
    IDEALIZED_TETRAHEDRON,
    get_idealized_tetrahedron_angles_from_cip_label,
    get_idealized_tetrahedron_angles_from_oemol,
    get_idealized_tetrahedron_coords_from_cip_label,
    get_tetrahedron_angles_from_coords_and_cip_label,
    get_tetrahedron_angles_from_oemol_coords,
)


@pytest.mark.parametrize(
    "stereo_label, neighbour_count, expected",
    [
        (
            TetrahedralCIPAssignment.R,
            3,
            np.array(
                [
                    [0.0, 0.0, 0.0],
                    [(8 / 9) ** 0.5, 0, -1 / 3],
                    [-((2 / 9) ** 0.5), (2 / 3) ** 0.5, -1 / 3],
                    [0.0, 0.0, 1.0],
                ]
            ),
        ),
        (
            TetrahedralCIPAssignment.S,
            4,
            np.array(
                [
                    [(8 / 9) ** 0.5, 0, -1 / 3],
                    [-((2 / 9) ** 0.5), (2 / 3) ** 0.5, -1 / 3],
                    [0.0, 0.0, 1.0],
                    [-((2 / 9) ** 0.5), -((2 / 3) ** 0.5), -1 / 3],
                ]
            ),
        ),
    ],
)
def test_get_idealized_tetrahedron_coords_from_r_s_label(stereo_label, neighbour_count, expected):
    result = get_idealized_tetrahedron_coords_from_cip_label(
        stereo_label=stereo_label, neighbour_count=neighbour_count
    )
    npt.assert_almost_equal(result, expected)
    with pytest.raises(ValueError):
        _ = get_idealized_tetrahedron_coords_from_cip_label(stereo_label=None, neighbour_count=3)


@pytest.mark.parametrize(
    "stereo_label, neighbour_count, expected",
    [
        (TetrahedralCIPAssignment.R, 3, 0.6154797086703875),
        (TetrahedralCIPAssignment.S, 3, -0.6154797086703875),
        (TetrahedralCIPAssignment.R, 4, 1.2309594173407747),
        (TetrahedralCIPAssignment.S, 4, -1.2309594173407747),
    ],
)
def test_get_idealized_tetrahedron_angles_from_r_s_label(stereo_label, neighbour_count, expected):
    result = get_idealized_tetrahedron_angles_from_cip_label(
        stereo_label=stereo_label, neighbour_count=neighbour_count
    )
    # Check only the first permutation, as that is what we have precomputed
    npt.assert_almost_equal(result[0], expected)


@pytest.mark.parametrize(
    "atom_coordinates, neighbour_coordinates, stereo_label",
    [
        (np.array([0, 0, 0]), IDEALIZED_TETRAHEDRON[:3], TetrahedralCIPAssignment.S),
        (np.array([0, 0, 0]), IDEALIZED_TETRAHEDRON[[0, 1, 3]], TetrahedralCIPAssignment.R),
        (np.array([0, 0, 0]), IDEALIZED_TETRAHEDRON, TetrahedralCIPAssignment.R),
        (np.array([0, 0, 0]), IDEALIZED_TETRAHEDRON[[0, 1, 3, 2]], TetrahedralCIPAssignment.S),
    ],
)
def test_get_tetrahedron_angles_from_coords_and_r_s_label(
    atom_coordinates, neighbour_coordinates, stereo_label
):
    result = get_tetrahedron_angles_from_coords_and_cip_label(
        atom_coordinates=atom_coordinates,
        cip_ordered_neighbour_coordinates=neighbour_coordinates,
    )
    ideal_result = get_idealized_tetrahedron_angles_from_cip_label(
        stereo_label=stereo_label,
        neighbour_count=neighbour_coordinates.shape[0],
    )
    npt.assert_almost_equal(result, ideal_result)


def test_get_oemol_idealized_tetrahedron_angles():
    achiral_smiles = "CCCCCCC"
    achiral_mol = std_oemol_from_smiles(achiral_smiles)
    achiral_angles = get_idealized_tetrahedron_angles_from_oemol(oemol=achiral_mol)
    assert achiral_angles.shape == (achiral_mol.NumAtoms(), 24)
    assert np.all(achiral_angles == 0.0)

    chiral_smiles = "C[C@@H](C(=O)O)N"
    chiral_mol = std_oemol_from_smiles(chiral_smiles)
    chiral_angles = get_idealized_tetrahedron_angles_from_oemol(oemol=chiral_mol)
    assert chiral_angles.shape == (chiral_mol.NumAtoms(), 24)
    assert not np.all(chiral_angles == 0.0)
    # One chiral center
    assert np.any(chiral_angles != 0.0, axis=1).sum() == 1

    flipped_chiral_smiles = "C[C@H](C(=O)O)N"
    flipped_chiral_mol = std_oemol_from_smiles(flipped_chiral_smiles)
    flipped_chiral_angles = get_idealized_tetrahedron_angles_from_oemol(oemol=flipped_chiral_mol)
    assert flipped_chiral_angles.shape == (flipped_chiral_mol.NumAtoms(), 24)
    assert not np.all(flipped_chiral_angles == 0.0)
    # One chiral center
    assert np.any(flipped_chiral_angles != 0.0, axis=1).sum() == 1
    assert not np.allclose(chiral_angles, flipped_chiral_angles)
    assert np.allclose(chiral_angles, -flipped_chiral_angles)


def get_conformer_from_smiles(smiles: str) -> oechem.OEMol:
    mol = std_oemol_from_smiles(smiles)
    conf = generate_single_conformer(mol)
    assert conf.NumAtoms() != mol.NumAtoms()
    standardise_oemol(conf)
    assert conf.NumAtoms() == mol.NumAtoms()
    return conf


def test_get_tetrahedron_angles_from_oemol_coords():
    achiral_conf = get_conformer_from_smiles("CCCCCCC")
    achiral_angles = get_tetrahedron_angles_from_oemol_coords(oemol=achiral_conf)
    assert achiral_angles.shape == (achiral_conf.NumAtoms(), 24)
    assert np.all(achiral_angles == 0.0)

    chiral_conf = get_conformer_from_smiles("C[C@@H](C(=O)O)N")
    chiral_angles = get_tetrahedron_angles_from_oemol_coords(oemol=chiral_conf)
    assert chiral_angles.shape == (chiral_conf.NumAtoms(), 24)
    assert not np.all(chiral_angles == 0.0)
    # One chiral center
    assert np.any(chiral_angles != 0.0, axis=1).sum() == 1

    # Check that mirroring the coordinates whilst keeping the chirality flips the angles
    flipped_chiral_conf = oechem.OEMol(chiral_conf)
    set_coordinates(flipped_chiral_conf, -get_coordinates(flipped_chiral_conf))
    flipped_chiral_angles = get_tetrahedron_angles_from_oemol_coords(oemol=flipped_chiral_conf)
    assert flipped_chiral_angles.shape == chiral_angles.shape
    assert np.allclose(chiral_angles, -flipped_chiral_angles)

    # Check that flipping the chirality has a similar effect
    other_chiral_conf = get_conformer_from_smiles("C[C@H](C(=O)O)N")
    other_chiral_angles = get_tetrahedron_angles_from_oemol_coords(oemol=other_chiral_conf)
    assert other_chiral_angles.shape == chiral_angles.shape
    assert np.allclose(other_chiral_angles, flipped_chiral_angles, atol=5e-2)

    # Check that the distance to the true angles is similar for both enantiomers
    true_angles = get_idealized_tetrahedron_angles_from_oemol(oemol=chiral_conf)
    other_true_angles = get_idealized_tetrahedron_angles_from_oemol(oemol=other_chiral_conf)
    assert np.allclose(
        chiral_angles[1] - true_angles[1],
        -1 * (other_chiral_angles[1] - other_true_angles[1]),
        atol=5e-2,
    )
