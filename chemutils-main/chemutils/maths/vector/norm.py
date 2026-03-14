from typing import Literal

import numpy as np

from chemutils.maths.convert import (
    CoordinatesLike,
    VectorsLike,
    as_coordinates,
    as_vectors,
)
from chemutils.maths.types import (
    Scalars,
    Vectors,
)


def length(arg: VectorsLike, /) -> Scalars:
    """
    Calculate the length of a N-dimensional vector or set of vectors.

    This function works an N-dimensional vector, or a (nested) list of such vectors.
    """
    arg = as_vectors(arg)
    return np.linalg.norm(arg, axis=-1)  # type: ignore


def normalized(v: VectorsLike) -> Vectors:
    v = as_vectors(v)
    return v / np.linalg.norm(v, axis=-1, keepdims=True)  # type: ignore


def distance_between(
    x: CoordinatesLike, y: CoordinatesLike, /, mode: Literal["pairwise", "outer"] = "pairwise"
) -> Scalars:
    """Calculate the distance between two sets of N-dimensional coordinates."""
    x = as_coordinates(x)
    y = as_coordinates(y)
    if mode == "outer":
        return np.linalg.norm(x[..., :, None, :] - y[..., None, :, :], axis=-1)  # type: ignore
    else:
        return np.linalg.norm(x - y, axis=-1)  # type: ignore


def distance_matrix(x: CoordinatesLike, y: CoordinatesLike | None = None, /) -> Scalars:
    """
    Calculate the distance matrix, containing the pairwise distances between the arguments.

    If only one argument is provided, a *self* distance matrix will be calculated.
    """
    return distance_between(x, y if y is not None else x, mode="outer")


def clamp_magnitude(v: VectorsLike, /, *, max_magnitude: float) -> Vectors:
    """Clamp the magnitude of a vector or collection of vectors such that it is not larger than a given value."""
    v = as_vectors(v)
    return v * np.minimum(max_magnitude / length(v), 1)[..., None]  # type: ignore
