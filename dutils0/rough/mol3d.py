"""Small RDKit 3D geometry helpers."""

from __future__ import annotations

from typing import cast

import numpy as np
from rdkit import Chem


def get_coord_array(mol: Chem.Mol, conf_id: int | None = None) -> np.ndarray:
    """
    Return atomic coordinates from an RDKit conformer as a NumPy array.

    Parameters
    ----------
    mol
        RDKit molecule (ROMol / Chem.Mol).
    conf_id
        Optional conformer ID. If provided, that conformer is used. If omitted
        (None), the molecule's *default* conformer is used via ``mol.GetConformer()``
        (RDKit's "first"/default conformer; typically the only conformer, or the
        one with ID 0 when conformers are embedded).

    Returns
    -------
    np.ndarray
        Array of shape ``(mol.GetNumAtoms(), 3)`` with dtype float, in the same
        atom order as ``mol``.

    Raises
    ------
    ValueError
        If the molecule has no conformers, if the requested conformer ID does
        not exist on the molecule, or if the selected conformer is not 3D.
    """
    if mol.GetNumConformers() == 0:
        raise ValueError("molecule has no conformers; cannot extract coordinates")

    if conf_id is None:
        # RDKit's default selection: the "first" conformer (usually ID 0).
        conf = mol.GetConformer()
    else:
        conf_ids = [c.GetId() for c in mol.GetConformers()]
        if conf_id not in conf_ids:
            raise ValueError(
                f"conformer id {conf_id} does not exist on molecule; available ids: {conf_ids}"
            )
        conf = mol.GetConformer(int(conf_id))

    if not conf.Is3D():
        raise ValueError(
            "selected conformer is not marked as 3D (2D coordinates are not accepted)"
        )

    coords = np.asarray(conf.GetPositions(), dtype=float)

    expected_shape = (mol.GetNumAtoms(), 3)
    if coords.shape != expected_shape:
        raise ValueError(f"expected coordinate array of shape {expected_shape}, got {coords.shape}")

    # NumPy typing: ensure we return a floating array.
    coords = cast(np.ndarray, coords)
    if not np.issubdtype(coords.dtype, np.floating):
        raise ValueError(f"expected float dtype for coordinates, got {coords.dtype}")

    return coords


def _validate_point3d(arr: np.ndarray, name: str = "point") -> None:
    """Raise ValueError if *arr* is not a valid 3D point array."""
    arr = np.asarray(arr)
    if arr.shape != (3,):
        raise ValueError(f"{name} must have shape (3,), got {arr.shape}")
    if not np.issubdtype(arr.dtype, np.number):
        raise ValueError(f"{name} must have a numeric dtype, got {arr.dtype}")
    if not np.all(np.isfinite(arr)):
        raise ValueError(f"{name} contains non-finite values (NaN or inf)")


def _validate_coords(coords: np.ndarray, name: str = "coords") -> None:
    """Raise ValueError if *coords* is not a valid (N, 3) coordinate array."""
    coords = np.asarray(coords)
    if coords.ndim != 2 or coords.shape[1] != 3:
        raise ValueError(f"{name} must have shape (N, 3), got {coords.shape}")
    if coords.shape[0] == 0:
        raise ValueError(f"{name} must contain at least one point")
    if not np.issubdtype(coords.dtype, np.number):
        raise ValueError(f"{name} must have a numeric dtype, got {coords.dtype}")
    if not np.all(np.isfinite(coords)):
        raise ValueError(f"{name} contains non-finite values (NaN or inf)")


def distance(a: np.ndarray, b: np.ndarray) -> float:
    """
    Return the Euclidean distance between two 3D points.

    Parameters
    ----------
    a, b
        Arrays of shape ``(3,)`` representing 3D coordinates.

    Returns
    -------
    float
        Euclidean distance between ``a`` and ``b``.

    Raises
    ------
    ValueError
        If either point is not a finite numeric array of shape ``(3,)``.
    """
    _validate_point3d(a, "a")
    _validate_point3d(b, "b")
    return float(np.linalg.norm(b - a))


def centroid(coords: np.ndarray) -> np.ndarray:
    """
    Return the centroid of a set of 3D coordinates.

    Parameters
    ----------
    coords
        Array of shape ``(N, 3)`` containing atomic coordinates.

    Returns
    -------
    np.ndarray
        Array of shape ``(3,)`` representing the centroid.

    Raises
    ------
    ValueError
        If *coords* is not a finite numeric array of shape ``(N, 3)`` with at
        least one point.
    """
    _validate_coords(coords)
    return np.mean(coords, axis=0)
