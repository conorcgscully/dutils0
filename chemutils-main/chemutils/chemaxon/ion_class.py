from enum import Enum

import numpy as np
from openeye import oechem

from .major_microspecies import get_major_microspecies


class IonClass(str, Enum):
    Zwitterion = "zwitterion"
    Base = "base"
    Acid = "acid"
    Uncharged = "uncharged"


def neutralise_adjacent_charges(
    mol: oechem.OEMolBase,
) -> oechem.OEMolBase:
    """
    Neutralise charge seperated charges such as nitro groups.

    This prepares a molecule for ion classification by remove charges that are directly bonded,
    such as the `N+-O-` bond in nitro or nitrones which would otherwise falsely classified zwitterionic.

    Args:
        mol: OpenEye molecule.

    Returns:
        oemol with adjacent charges neutralised.
    """
    copy = mol.CreateCopy()
    for atom in copy.GetAtoms():
        chargeA = atom.GetFormalCharge()
        if chargeA != 0:
            nbrs = atom.GetAtoms()
            for nbr in nbrs:
                chargeB = nbr.GetFormalCharge()
                if chargeB == -chargeA:  # check charge is equal and opposite
                    atom.SetFormalCharge(0)
                    nbr.SetFormalCharge(0)
    return copy


def get_ion_class(
    mol: oechem.OEMolBase,
) -> IonClass:
    """
    Determine the ion class of a molecule, which can be one of 'zwitterion', 'base', 'acid', or 'uncharged'.

    This checks to remove charges that are directly bonded, i.e. N+-O- bond in nitro or nitrones.

    Args:
        mol: OpenEye molecule.

    Returns:
        String determining the ion class - one of 'zwitterion', 'base', 'acid', or 'neutral'.
    """
    edited_mol = neutralise_adjacent_charges(mol)

    count_positive = 0.0
    count_negative = 0.0
    for atom in edited_mol.GetAtoms():
        atom_charge = atom.GetFormalCharge()
        if atom_charge > 0.0:
            count_positive += atom_charge
        elif atom_charge < 0.0:
            count_negative += atom_charge

    # If the molecule has equal positive and negative charges, it's a zwitterion
    # if the molecule has two basic groups and one acidic group, it should count as basic
    if count_positive == np.abs(count_negative) and count_positive > 0:
        return IonClass.Zwitterion
    elif count_positive > np.abs(count_negative):
        return IonClass.Base
    elif np.abs(count_negative) > count_positive:
        return IonClass.Acid
    else:
        return IonClass.Uncharged


def get_major_microspecies_ion_class(
    mol: oechem.OEMolBase,
    pH: float = 7.4,
) -> IonClass:
    """
    Determines the ion class of a molecule based on the dominant microspecies.

    See `get_ion_class` for more details.

    Args:
        mol: OpenEye molecule.
        pH: pH to calculate the major microspecies at.

    Returns:
        String determining the ion class - one of 'zwitterion', 'base', 'acid', or 'neutral'.
    """
    return get_ion_class(get_major_microspecies(mol, pH=pH))
