import re

import numpy as np
import pytest
from openeye import oechem

from chemutils.maths.metric.gdt_ts import _unaligned_gdt, get_unaligned_gdt_ts
from chemutils.molecule.create import create_atom


def create_protein_with_ca_atoms(ca_coords: list[tuple[float, float, float]]) -> oechem.OEMol:
    mol = oechem.OEMol()

    for i, (x, y, z) in enumerate(ca_coords):
        create_atom(
            mol,
            atomic_number=6,
            name="CA",
            coordinates=(x, y, z),
            residue_number=i,
            chain_id="A",
        )

    return mol


@pytest.mark.parametrize(
    ["ca_coords1", "ca_coords2", "expected_gdt"],
    [
        (
            [(0.0, 0.0, 0.0), (1.5, 0.0, 0.0), (3.0, 0.0, 0.0)],
            [(0.0, 0.0, 0.0), (1.5, 0.0, 0.0), (3.0, 0.0, 0.0)],
            1.0,
        ),
        (
            [(0.0, 0.0, 0.0), (1.5, 0.0, 0.0), (3.0, 0.0, 0.0)],
            [(0.0, 0.0, 0.0), (1.5, 1.5, 0.0), (3.0, 3.0, 0.0)],
            0.75,
        ),
        (
            [(0.0, 0.0, 0.0), (1.5, 1.5, 0.0), (3.0, 3.0, 0.0)],
            [(0.0, 0.0, 0.0), (1.5, 0.0, 0.0), (3.0, 0.0, 0.0)],
            0.75,
        ),
        (
            [(0.0, 0.0, 0.0), (1.5, 0.0, 0.0), (3.0, 0.0, 0.0)],
            [(0.0, 0.0, 0.0), (1.5, 0.0, 0.0)],
            1.0,
        ),
    ],
    ids=[
        "same protein",
        "different proteins",
        "different proteins in swapped order",
        "one protein with less CA atoms than the other",
    ],
)
def test_get_unaligned_gdt_ts(ca_coords1, ca_coords2, expected_gdt):
    protein1 = create_protein_with_ca_atoms(ca_coords1)
    protein2 = create_protein_with_ca_atoms(ca_coords2)

    assert get_unaligned_gdt_ts(protein1, protein2) == pytest.approx(expected_gdt)


def test_get_unaligned_gdt_ts_including_non_ca_atoms():
    ca_coords = [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (2.0, 0.0, 0.0)]
    protein1 = create_protein_with_ca_atoms(ca_coords)
    create_atom(
        protein1,
        atomic_number=8,
        name="O",
        coordinates=(3.0, 0.0, 0.0),
        residue_number=4,
        chain_id="A",
    )

    protein2 = create_protein_with_ca_atoms(ca_coords)
    create_atom(
        protein2,
        atomic_number=7,
        name="N",
        coordinates=(3.0, 0.0, 0.0),
        residue_number=4,
        chain_id="A",
    )

    assert get_unaligned_gdt_ts(protein1, protein2) == pytest.approx(1.0)


def test_get_unaligned_gdt_ts_multiple_alpha_carbons_with_same_residue_number_src_protein():
    ca_coords = [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (2.0, 0.0, 0.0)]
    protein1 = create_protein_with_ca_atoms(ca_coords)
    protein2 = create_protein_with_ca_atoms(ca_coords)
    create_atom(
        protein1,
        atomic_number=6,
        name="CA",
        coordinates=(0.0, 0.0, 0.0),
        residue_number=1,
        chain_id="A",
    )

    with pytest.raises(
        ValueError, match="Multiple alpha carbons with the same residue number: 1 in `src` protein"
    ):
        get_unaligned_gdt_ts(protein1, protein2)


def test_get_unaligned_gdt_ts_multiple_alpha_carbons_with_same_residue_number_dest_protein():
    ca_coords = [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (2.0, 0.0, 0.0)]
    protein1 = create_protein_with_ca_atoms(ca_coords)
    protein2 = create_protein_with_ca_atoms(ca_coords)
    create_atom(
        protein2,
        atomic_number=6,
        name="CA",
        coordinates=(0.0, 0.0, 0.0),
        residue_number=1,
        chain_id="A",
    )

    with pytest.raises(
        ValueError, match="Multiple alpha carbons with the same residue number: 1 in `dest` protein"
    ):
        get_unaligned_gdt_ts(protein1, protein2)


def test_get_unaligned_gdt_ts_multiple_alpha_carbons_with_same_residue_number_multiple_chains_src_protein():
    ca_coords = [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (2.0, 0.0, 0.0)]
    protein1 = create_protein_with_ca_atoms(ca_coords)
    protein2 = create_protein_with_ca_atoms(ca_coords)
    create_atom(
        protein1,
        atomic_number=6,
        name="CA",
        coordinates=(0.0, 0.0, 0.0),
        residue_number=1,
        chain_id="B",
    )

    with pytest.raises(
        ValueError, match="Multiple alpha carbons with the same residue number: 1 in `src` protein"
    ):
        get_unaligned_gdt_ts(protein1, protein2)


def test_get_unaligned_gdt_ts_multiple_alpha_carbons_with_same_residue_number_multiple_chains_dest_protein():
    ca_coords = [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (2.0, 0.0, 0.0)]
    protein1 = create_protein_with_ca_atoms(ca_coords)
    protein2 = create_protein_with_ca_atoms(ca_coords)
    create_atom(
        protein2,
        atomic_number=6,
        name="CA",
        coordinates=(0.0, 0.0, 0.0),
        residue_number=1,
        chain_id="B",
    )

    with pytest.raises(
        ValueError, match="Multiple alpha carbons with the same residue number: 1 in `dest` protein"
    ):
        get_unaligned_gdt_ts(protein1, protein2)


def test_unaligned_gdt_shapes_mismatch():
    coords1 = np.array([(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (2.0, 0.0, 0.0)])
    coords2 = np.array([(0.0, 0.0, 0.0), (1.0, 0.0, 0.0)])
    with pytest.raises(ValueError, match=re.escape("Shape mismatch: (3, 3) != (2, 3)")):
        _unaligned_gdt(coords1, coords2, cutoffs=(1.0, 2.0, 4.0, 8.0))
