import pickle

import pytest

from chemutils.molecule import assert_mol_equal

from .cip_validation import CIP_VALIDATION_OEMOLS_3D


@pytest.mark.parametrize("oemol", CIP_VALIDATION_OEMOLS_3D)
def test_pickle(oemol):
    # Test to make sure that the fix to OEGraphMol means that OEGraphMol now pickles correctly
    bytes = pickle.dumps(oemol)
    oemol2 = pickle.loads(bytes)
    assert_mol_equal(oemol, oemol2)
