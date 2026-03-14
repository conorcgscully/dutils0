from dataclasses import dataclass
from pathlib import Path

import fsutils as fs
import pytest
from openeye import oechem

from chemutils.molecule import oemol_from_smiles
from chemutils.molecule.hydrogen_bonds import (
    get_hybridization,
    get_num_acceptor_lone_pairs,
    get_num_donor_hydrogens,
    get_num_hydrogen_acceptors,
    get_num_hydrogen_donors,
    is_hydrogen_acceptor,
    is_hydrogen_donor,
)


@dataclass
class HBondCase:
    name: str
    smiles: str
    num_hydrogen_donors: int
    num_donor_hydrogens: int
    hydrogen_donor_indices: list[int]
    num_hydrogen_acceptors: int
    num_acceptor_lone_pairs: int
    hydrogen_acceptor_indices: list[int]


# Use this default to allow fields to be omitted from the yaml when there's no acceptors/donors
DEFAULT = {
    "num_hydrogen_donors": 0,
    "num_donor_hydrogens": 0,
    "hydrogen_donor_indices": [],
    "num_hydrogen_acceptors": 0,
    "num_acceptor_lone_pairs": 0,
    "hydrogen_acceptor_indices": [],
}

CASES = [
    HBondCase(**(DEFAULT | case)) for case in fs.read_yaml(Path(__file__).parent / "hbonds.yaml")
]


@pytest.mark.parametrize(
    ["smiles", "num_hydrogen_donors"], [(case.smiles, case.num_hydrogen_donors) for case in CASES]
)
def test_get_num_hydrogen_donors(smiles, num_hydrogen_donors):
    assert get_num_hydrogen_donors(oemol_from_smiles(smiles)) == num_hydrogen_donors


@pytest.mark.parametrize(
    ["smiles", "num_donor_hydrogens"], [(case.smiles, case.num_donor_hydrogens) for case in CASES]
)
def test_get_num_donor_hydrogens(smiles, num_donor_hydrogens):
    assert get_num_donor_hydrogens(oemol_from_smiles(smiles)) == num_donor_hydrogens


@pytest.mark.parametrize(
    ["smiles", "hydrogen_donor_indices"],
    [(case.smiles, case.hydrogen_donor_indices) for case in CASES],
)
def test_is_hydrogen_donor(smiles, hydrogen_donor_indices):
    oemol = oemol_from_smiles(smiles)
    assert [
        i for i, atom in enumerate(oemol.GetAtoms()) if is_hydrogen_donor(atom)
    ] == hydrogen_donor_indices


@pytest.mark.parametrize(
    ["smiles", "num_hydrogen_acceptors"],
    [(case.smiles, case.num_hydrogen_acceptors) for case in CASES],
)
def test_get_num_hydrogen_acceptors(smiles, num_hydrogen_acceptors):
    assert get_num_hydrogen_acceptors(oemol_from_smiles(smiles)) == num_hydrogen_acceptors


@pytest.mark.parametrize(
    ["smiles", "num_acceptor_lone_pairs"],
    [(case.smiles, case.num_acceptor_lone_pairs) for case in CASES],
)
def test_get_num_acceptor_lone_pairs(smiles, num_acceptor_lone_pairs):
    assert get_num_acceptor_lone_pairs(oemol_from_smiles(smiles)) == num_acceptor_lone_pairs


@pytest.mark.parametrize(
    ["smiles", "hydrogen_acceptor_indices"],
    [(case.smiles, case.hydrogen_acceptor_indices) for case in CASES],
)
def test_is_hydrogen_acceptor(smiles, hydrogen_acceptor_indices):
    oemol = oemol_from_smiles(smiles)
    assert [
        i for i, atom in enumerate(oemol.GetAtoms()) if is_hydrogen_acceptor(atom)
    ] == hydrogen_acceptor_indices


@pytest.mark.parametrize(
    ["smiles", "hydroxyl_index", "hybridisation"],
    [
        ("CC[O-]", 2, oechem.OEHybridization_sp3),
        ("C(=O)[O-]", 2, oechem.OEHybridization_sp2),
        ("Cl(=O)(=O)(=O)[O-]", 4, oechem.OEHybridization_sp3),
        ("S(=O)(=O)([O-])[O-]", 3, oechem.OEHybridization_sp3),
        ("P(=O)([O-])([O-])[O-]", 2, oechem.OEHybridization_sp3),
    ],
)
def test_hydroxyl(smiles, hydroxyl_index, hybridisation):
    oemol = oemol_from_smiles(smiles)
    assert get_hybridization(oemol.GetAtom(oechem.OEHasAtomIdx(hydroxyl_index))) == hybridisation
