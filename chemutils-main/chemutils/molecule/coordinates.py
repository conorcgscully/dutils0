from typing import Literal, overload

import numpy as np
from numpy.typing import NDArray
from openeye import oechem


def get_coordinates(oemol: oechem.OEMolBase, /) -> NDArray[np.float32]:
    """
    Get the coordinates from an OpenEye molecule, which can represent either a single structure or multiple conformers.

    If the molecule is an instance of `oechem.OEMCMolBase` and can represent one or more conformers, then an array
    of shape (number of conformers, number of atoms, dimension) is returned. Otherwise, an array of shape
    (number of atoms, dimension) is returned.

    Args:
        oemol: OpenEye molecule with coordinates.

    Returns:
        NumPy array of shape (NumberOfAtoms, Dimensions).

    Raises:
        ValueError: The molecule has dimension 0, and hence no coordinates.
    """
    dim = oemol.GetDimension()
    if dim == 0:
        raise ValueError("Molecule does not have coordinates.")
    if isinstance(oemol, oechem.OEMCMolBase):
        return np.array([list(conf.GetCoords().values()) for conf in oemol.GetConfs()])[:, :, :dim]
    return np.array(list(oemol.GetCoords().values()))[:, :dim]


@overload
def set_coordinates(
    oemol: oechem.OEGraphMol, /, coordinates: NDArray[np.float32], *, copy_mol: Literal[True]
) -> oechem.OEGraphMol:
    pass


@overload
def set_coordinates(
    oemol: oechem.OEGraphMol, /, coordinates: NDArray[np.float32], *, copy_mol: Literal[False] = ...
) -> None:
    pass


@overload
def set_coordinates(
    oemol: oechem.OEGraphMol, /, coordinates: NDArray[np.float32], *, copy_mol: bool = ...
) -> oechem.OEGraphMol | None:
    pass


def set_coordinates(
    oemol: oechem.OEGraphMol, /, coordinates: NDArray[np.float32], *, copy_mol: bool = False
) -> oechem.OEGraphMol | oechem.OEMCMolBase | None:
    """
    Set an OpenEye molecule's coordinates.

    This modifies the molecule in-place by default, but returns a copy of the molecule
    if `copy_mol` is `True`.

    Args:
        oemol: Original OpenEye molecule.
        coordinates: Coordinates as a NumPy array of shape (NumberOfAtoms, Dimensions).
        copy_mol: Should the molecule be copied before modifying the coordinates.

    Returns:
        Copy of `oemol` with specified coordinates if `copy_mol` is `True`, else `None`.

    Raises:
        ValueError: Coordinates array is invalid shape.
    """
    if copy_mol:
        oemol = oemol.CreateCopy()
    if (
        len(coordinates.shape) not in (2, 3)
        or coordinates.shape[-2] != oemol.NumAtoms()
        or coordinates.shape[-1] not in {1, 2, 3}
    ):
        raise ValueError(
            "Coordinates must be of shape (NumberOfAtoms, Dimensions) or (NumberOfConformers, NumberOfAtoms, Dimensions) with 1 <= Dimensions <= 3"
        )
    if len(coordinates.shape) == 3:
        if not copy_mol and isinstance(oemol, oechem.OEGraphMol):
            raise ValueError(
                "Cannot set multiple conformers on a non-multi-conformer molecule. Use `copy_mol=True`."
            )
        if not isinstance(oemol, oechem.OEMCMolBase):
            oemol = oechem.OEMol(oemol)
        oemol.DeleteConfs()
        for coords in coordinates:
            # Pad coordinates to be 3xN with 0's, and then convert to a flat array.
            coords = np.pad(coords, ((0, 0), (0, 3 - coords.shape[-1]))).reshape(-1)
            oemol.NewConf(oechem.OEFloatArray(coords))
        oemol.SetDimension(coordinates.shape[-1])
    else:
        oemol.SetDimension(coordinates.shape[-1])
        oemol.SetCoords(np.pad(coordinates, ((0, 0), (0, 3 - coordinates.shape[-1]))).reshape(-1))

    if copy_mol:
        return oemol
    return None
