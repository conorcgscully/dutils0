from openeye import oechem, oemolprop
from rdkit import Chem
from rdkit.Chem import rdmolops
from rdkit.Chem.rdMolDescriptors import CalcNumSpiroAtoms
from rdkit.Chem.SpacialScore import SPS

from chemutils.rdmol_from_oemol import as_rdmol

from .property import MolecularProperty


def get_num_heavy_atoms(mol: oechem.OEGraphMol) -> int:
    """
    Get the number of non-hydrogen atoms in a molecule.

    Arguments:
        mol: Molecule to compute number of heavy atoms for.

    Returns:
        Number of heavy atoms in query molecule.
    """
    """"""
    return sum(atom.GetAtomicNum() != 1 for atom in mol.GetAtoms())


def get_num_chiral_atoms(mol: oechem.OEGraphMol, include_nitrogens: bool = True) -> int:
    """
    Get the number of chiral atoms (tetrahedral stereocenters).

    Arguments:
        mol: Molecule to compute number of stereocenters.
        include_nitrogens: Wether nitrogen stereocenters count. default = True

    Returns:
        Number of chiral atoms in query molecule.
    """
    if include_nitrogens:
        return sum(atom.IsChiral() for atom in mol.GetAtoms())
    return sum(atom.IsChiral() for atom in mol.GetAtoms() if not atom.IsNitrogen())


def get_num_unspecified_chiral_atoms(mol: oechem.OEGraphMol, include_nitrogens: bool = True) -> int:
    """
    Get the number of undefined chiral atoms (tetrahedral stereocenters).

    Arguments:
        mol: Molecule to compute number of stereocenters.
        include_nitrogens: Wether nitrogen stereocenters count. default = True

    Returns:
        Number of undefined chiral atoms in query molecule.
    """
    if include_nitrogens:
        return sum(not atom.HasStereoSpecified() for atom in mol.GetAtoms() if atom.IsChiral())
    return sum(
        not atom.HasStereoSpecified()
        for atom in mol.GetAtoms()
        if atom.IsChiral() and not atom.IsNitrogen()
    )


def get_num_spiro_atoms(mol: Chem.Mol | oechem.OEMolBase) -> int:
    """
    Get the number of spiro atoms (tetrahedral stereocenters).

    Arguments:
        mol: Molecule to compute number of stereocenters.

    Returns:
        Number of spiro atoms in query molecule.
    """
    mol = as_rdmol(mol)
    return CalcNumSpiroAtoms(mol)  # type: ignore


def get_spacial_score(mol: Chem.Mol | oechem.OEMolBase) -> int:
    """
    Calculate the normalized Spacial Score, which is a measure of complexity.

    Make sure to standardise OEMols before using this function as it is sensitive to hybridization perception.
    Paper: https://pubs.acs.org/doi/10.1021/acs.jmedchem.3c00689.

    Arguments:
        mol: Molecule to compute spacial score for.

    Returns:
        Normalized spacial score for query molecule.
    """
    mol = as_rdmol(mol)
    # Explicit hydrogens cause a massive overestimation of the spacial score
    # Annoying that RDKit doesn't account for this
    mol = rdmolops.RemoveAllHs(mol)
    return SPS(mol)  # type: ignore


def get_molecular_weight(mol: oechem.OEGraphMol) -> float:
    """
    Get the molecular weight of a molecule in daltons.

    Arguments:
        mol: Molecule to compute the molecular weight.

    Returns:
        Molecular weight of query molecule in Daltons.
    """
    return oechem.OECalculateMolecularWeight(mol)  # type: ignore


def get_fraction_carbon_sp3(mol: oechem.OEMolBase) -> float:
    """Get the fraction of carbons in the molecule that have sp3 hybridisation, as described by Lovering et al. (2009)."""
    return oemolprop.OEGetFractionCsp3(mol)  # type: ignore


NumberHeavyAtoms = MolecularProperty(id="num_heavy_atoms", oemol_func=get_num_heavy_atoms)

NumberSpiroAtoms = MolecularProperty(id="num_spiro_atoms", rdkit_func=get_num_spiro_atoms)

SpacialScore = MolecularProperty(id="spacial_score", rdkit_func=get_spacial_score)

NumberChiralAtoms = MolecularProperty(id="num_chiral_atoms", oemol_func=get_num_chiral_atoms)

MolecularWeight = MolecularProperty(id="molecular_weight", oemol_func=get_molecular_weight)

FractionCarbonSP3 = MolecularProperty(id="fsp3", oemol_func=get_fraction_carbon_sp3)
"""Fraction of carbons that are SP3 hybridised, as described by Lovering et al. (2009)."""
