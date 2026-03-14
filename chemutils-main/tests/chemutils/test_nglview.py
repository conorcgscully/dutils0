import numpy as np
import pytest

from chemutils.molecule import oemol_from_smiles, set_coordinates
from chemutils.molecule.embed import embed_molecule_3d
from chemutils.nglview.widget import (
    OEMolStructure,
    OEMolTrajectory,
    get_nglview_structure,
)


def test_get_nglview_structure_no_coords():
    oemol = oemol_from_smiles("c1ccccc1")
    with pytest.raises(RuntimeError):
        _ = get_nglview_structure(oemol)


def test_get_nglview_structure_structure():
    oemol = oemol_from_smiles("c1ccccc1")
    embed_molecule_3d(oemol, seed=0)
    widget = get_nglview_structure(oemol)
    assert isinstance(widget, OEMolStructure)


def test_get_nglview_structure_coords():
    oemol = oemol_from_smiles("c1ccccc1")
    coords = np.random.random((3, 6, 3))
    oemol = set_coordinates(oemol, coordinates=coords, copy_mol=True)
    widget = get_nglview_structure(oemol)
    assert isinstance(widget, OEMolTrajectory)
    assert widget.n_frames == 3
    np.testing.assert_allclose(widget.get_coordinates(1), coords[1])
