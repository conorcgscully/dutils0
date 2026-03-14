import pytest

from chemutils.molecule import mol_equal, oemol_from_smiles
from chemutils.molecule.canonicalize import canonicalize_molecule


@pytest.mark.parametrize("smiles_set", [["CC(=O)O", "O=C(O)C", "OC(=O)C", "C(C)(=O)O"]])
def test_canonicalize_molecule(smiles_set):
    oemols = [oemol_from_smiles(smiles) for smiles in smiles_set]
    for oemol in oemols:
        canonicalize_molecule(oemol)

    for i in range(1, len(oemols)):
        assert mol_equal(oemols[0], oemols[i])
