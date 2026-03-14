import pytest
from openeye import oechem

from chemutils.molecule.cip.assignment import (
    get_tetrahedral_cip_from_ordered_neighbours,
)
from chemutils.molecule.cip.priority import neighbours_ordered_by_cip_priority

from ...cip_validation import CIP_VALIDATION_TETRA_STEREO_CASES


@pytest.mark.parametrize(
    ["oemol", "atom_idx", "validation_label"], CIP_VALIDATION_TETRA_STEREO_CASES
)
def test_cip_validation_suite(oemol, atom_idx, validation_label):
    atom = oemol.GetAtom(oechem.OEHasAtomIdx(atom_idx))

    ordered_neighbours = neighbours_ordered_by_cip_priority(atom)
    chemutils_label = get_tetrahedral_cip_from_ordered_neighbours(
        atom, cip_ordered_neighbours=ordered_neighbours
    )
    assert chemutils_label == validation_label
