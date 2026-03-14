from typing import Any

import numpy as np
from openeye import oechem
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.DataStructs import cDataStructs

from chemutils.props.property import MolecularProperty
from chemutils.rdmol_from_oemol import as_rdmol, clean_rdmol_hydrogens


def get_morgan_fingerprint(
    mol: Chem.Mol | oechem.OEMolBase,
    radius: int = 3,
    num_bits: int = 2048,
    useChirality: bool = False,
) -> Any:
    """
    Get the morgan fingerprint of a molecule.

    The Morgan fingerprint measures the environment of each atom.

    This returns a 2048 bit array, packed into 256 NumPy array of uint8. A full array of 2048 items
    can be reproduced by using `numpy.unpackbits`.

    Arguments:
        mol: RDKit molecule.
        radius: Distance to measure a Morgan Fingerprint
        num_bits: Number of bits to encode the fingerprint in.
        useChirality: Whether to use chirality in the fingerprint.

    Returns:
        Morgan fingerprint as a NumPy array of num_bits bits packed into 256 uint8's.
    """
    mol = as_rdmol(mol)
    mol = clean_rdmol_hydrogens(mol)

    fingerprint = rdMolDescriptors.GetMorganFingerprintAsBitVect(
        mol, radius=radius, nBits=num_bits, useChirality=useChirality
    )
    array = np.zeros((0,), dtype=np.uint8)
    cDataStructs.ConvertToNumpyArray(fingerprint, array)
    array = np.packbits(array)
    return array


MorganFingerprint = MolecularProperty(id="morgan_fingerprint", rdkit_func=get_morgan_fingerprint)
