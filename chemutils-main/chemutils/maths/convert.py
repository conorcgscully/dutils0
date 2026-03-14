from typing import TypeAlias

import numpy as np
import numpy.typing as npt
from openeye import oechem

from chemutils.molecule import get_coordinates

from .types import Vectors

CoordinatesLike: TypeAlias = (
    npt.ArrayLike
    | list[oechem.OEAtomBase]
    | oechem.OEAtomBase
    | oechem.OEMolBase
    | oechem.OEAtomIter
    | oechem.OEAtomBondSet
)
"""Any object that can be converted to a (nested) array of coordinates."""
VectorsLike: TypeAlias = (
    npt.ArrayLike
    | list[oechem.OEBondBase]
    | oechem.OEBondBase
    | oechem.OEBondIter
    | oechem.OEAtomBondSet
)
"""Any object that can be converted to a (nested) array of vectors."""


def as_coordinates(arg: CoordinatesLike, /) -> Vectors:
    """
    Interpret an argument as a NumPy array representing coordinates.

    If the argument can be converted to a NumPy array, such as a nested combination of lists or tuples,
    then it will be converted to a NumPy array.

    If the argument is an OpenEye molecule, then the coordinates of the molecule will be returned as a
    NumPy array of shape `(N, d)`, where `N` is the number of atoms in the molecule and `d` is the dimensionality
    of the coordinates (usually 3).

    If the argument is an OpenEye atom, then the coordinates of the atom will be returned, of shape `(d, )`,
    where `d` is the dimensionality of the coordinates (usually 3).
    """
    match arg:
        case oechem.OEMolBase():
            if arg.NumAtoms() == 0:
                return np.array([], dtype=float)
            return get_coordinates(arg)
        case oechem.OEAtomBase():
            return np.array(arg.GetParent().GetCoords(arg))
        case list() if len(arg) == 0:
            return np.array([], dtype=float)
        case list() if isinstance(arg[0], oechem.OEAtomBase):
            return np.array([atom.GetParent().GetCoords(atom) for atom in arg])
        case oechem.OEAtomBondSet():
            return np.array([atom.GetParent().GetCoords(atom) for atom in arg.GetAtoms()])
        case oechem.OEAtomIter():
            return np.array([atom.GetParent().GetCoords(atom) for atom in arg])
        case _:
            return np.asarray(arg)


def as_vectors(arg: VectorsLike, /) -> Vectors:
    """
    Convert an argument to a NumPy array representing vectors.

    As well as handling standard conversion of lists, tuples etc. to NumPy arrays,
    this also handles OpenEye bonds and bond iterators.
    """
    match arg:
        case oechem.OEBondBase():
            return _get_vector_from_bond(arg)
        case oechem.OEBondIter():
            return np.array([_get_vector_from_bond(bond) for bond in arg])
        case oechem.OEAtomBondSet():
            return np.array([_get_vector_from_bond(bond) for bond in arg.GetBonds()])
        case list() if len(arg) == 0:
            return np.array([], dtype=float)
        case list() if isinstance(arg[0], oechem.OEBondBase):
            return np.array([_get_vector_from_bond(bond) for bond in arg])
        case _:
            return np.asarray(arg)


def _get_vector_from_bond(bond: oechem.OEBondBase) -> npt.NDArray[np.float32]:
    return np.array(bond.GetParent().GetCoords(bond.GetEnd())) - np.array(  # type: ignore
        bond.GetParent().GetCoords(bond.GetBgn())
    )
