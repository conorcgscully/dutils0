from typing import Any

import nglview
import numpy as np
from openeye import oechem

from chemutils import fs
from chemutils.molecule import write_molecule_str
from chemutils.molecule.coordinates import get_coordinates


class OEMolStructure(nglview.Structure):  # type: ignore[misc]
    """Wrapper around an OpenEye molecule for use with nglview."""

    def __init__(self, oemol: oechem.OEMolBase):
        super().__init__()
        self.ext = "pdb"
        self.oemol = oemol

    def get_structure_string(self) -> str:
        return write_molecule_str(self.oemol, format="pdb")


class OEMolTrajectory(nglview.Structure, nglview.Trajectory):  # type: ignore[misc]
    """Wrapper around an OpenEye molecule for use with nglview."""

    def __init__(self, oemol: oechem.OEMCMolBase):
        super().__init__()
        self.ext = "pdb"
        self.oemol = oemol

    def get_structure_string(self) -> str:
        return write_molecule_str(self.oemol, format="pdb")

    def get_coordinates(self, index: int) -> np.typing.NDArray[np.float32]:
        return get_coordinates(self.oemol)[index]  # type: ignore[no-any-return]

    @property
    def n_frames(self) -> int:
        return self.oemol.NumConfs()  # type: ignore[no-any-return]


def get_nglview_structure(obj: Any) -> nglview.NGLWidget:
    """Show an arbitrary object using nglview."""
    if isinstance(obj, str):
        obj = fs.read_molecule(obj)

    if isinstance(obj, oechem.OEMCMolBase) and obj.NumConfs() > 1:
        return OEMolTrajectory(obj)

    if isinstance(obj, oechem.OEMolBase):
        if obj.GetDimension() != 3:
            raise RuntimeError("OpenEye molecule does not have coordinates.")
        return OEMolStructure(obj)

    raise ValueError(f"Cannot work out how to show object {obj} using nglview")
