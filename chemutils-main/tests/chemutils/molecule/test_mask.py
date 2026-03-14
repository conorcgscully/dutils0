import numpy as np
import pytest
from openeye import oechem

from chemutils.molecule import oemol_from_smiles
from chemutils.molecule.mask import get_atom_mask


@pytest.mark.parametrize(
    ["smiles", "predicate", "mask"],
    [
        ("C=CC", oechem.OEHasDoubleBond, [1, 1, 0]),
        ("CC(=O)C", oechem.OEHasDoubleBondO, [0, 1, 0, 0]),
        ("NCO", oechem.OEHasAtomicNum(7), [1, 0, 0]),
        ("CC=C", lambda atom: atom.GetImplicitHCount() == 2, [0, 0, 1]),
    ],
)
def test_get_atom_mask(smiles, predicate, mask):
    oemol = oemol_from_smiles(smiles)
    result = get_atom_mask(oemol, predicate)
    assert result.shape == (oemol.NumAtoms(),)
    assert result.dtype == np.bool_
    np.testing.assert_equal(result, np.array(mask, dtype=np.bool_))
