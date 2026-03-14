"""
Code for determining a representative version of a molecule.

Converting a molecule to a 'representative' form aims to standardise many features of a molecule,
whilst being more generic than `standardise_oemol` (which is aimed at DLF).

Namely, it supports both molecules with multiple components, neutralises molecules with a more effective
algorithm, and standardises various enhanced stereochemistry labels.
"""

from openeye import oechem

from .canonicalize import canonicalize_molecule
from .enhanced_stereo import (
    convert_enhanced_stereo_to_unspecified,
    get_canonical_stereoisomer,
    remove_redundant_enhanced_stereo,
)
from .hydrogens import make_hydrogens_implicit
from .neutralize import neutralize_molecule
from .perceive import perceive_atom_properties
from .smiles import cxsmiles_from_oemol, oemol_from_smiles, smiles_from_oemol


def representative_oemol_from_oemol(oemol: oechem.OEMolBase, /) -> oechem.OEMolBase:
    """
    Convert an OpenEye molecule to a representative molecule that can be used for comparisons and databases.

    The steps taken are:
    * Neutralising the molecule using `neutralize_molecule`. This fixes valencies and charges, whilst setting
      certain groups such as nitros to have charge separation.
    * Make hydrogens implicit. Unlike `standardise_oemol`, this does not remove isotopic hydrogens.
    * Percieve atom properties such as hybridization, aromaticity, and stereochemistry from 3D coordinates.
    * Canonicalize the molecule using `canonicalize_molecule`. This orders the atoms based on the order they appear in
      the canonical smiles.
    * If enhanced sterochemistry is present, standardise the order & remove redunant labels such as absolute stereo labels
      and labels for meso compounds.

    Args:
      oemol: OpenEye molecule to convert to a representative form.
    """
    oemol = oemol.CreateCopy()

    neutralize_molecule(oemol)
    make_hydrogens_implicit(oemol, remove_isotopic_hydrogens=False)
    perceive_atom_properties(oemol)
    canonicalize_molecule(oemol)

    # Enhanced stereo cleanup
    if oemol.NumGroups() > 0:
        remove_redundant_enhanced_stereo(oemol)
        oemol = get_canonical_stereoisomer(oemol)

    return oemol


def representative_oemol_from_smiles(smiles: str, /) -> oechem.OEGraphMol:
    """
    Convert a SMILES string to a representative OpenEye molecule.

    See `convert_molecule_to_representative` for more details.

    Args:
      smiles: SMILES string.

    Returns:
      OpenEye molecule in a representative form.
    """
    oemol = oemol_from_smiles(smiles)
    oemol = representative_oemol_from_oemol(oemol)
    return oemol


def representative_smiles_from_smiles(smiles: str, /) -> str:
    """
    Convert a SMILES string to a representative form.

    See `convert_molecule_to_representative` for more details.

    If a CXSMILES is provided, the enhanced stereochemistry will be converted to unspecified stereochemistry.

    Args:
      smiles: SMILES string.

    Returns:
      SMILES string of representative molecule.
    """
    oemol = representative_oemol_from_smiles(smiles)
    if oemol.NumGroups() > 0:
        convert_enhanced_stereo_to_unspecified(oemol)
    return smiles_from_oemol(oemol)


def representative_cxsmiles_from_smiles(smiles: str, /) -> str:
    """
    Convert a CXSMILES string to a representative form.

    See `convert_molecule_to_representative` for more details.

    Args:
      smiles: CXSMILES string.

    Returns:
      CXSMILES string of representative molecule.
    """
    oemol = representative_oemol_from_smiles(smiles)
    return cxsmiles_from_oemol(oemol)


def representative_smiles_from_oemol(
    oemol: oechem.OEMolBase, /, *, already_representative: bool = False
) -> str:
    """
    Get the representative SMILES of a molecule.

    Args:
      oemol: OpenEye molecule.
      already_representative: Whether the molecule is already in a representative form. If this is True, this is
        simply `smiles_from_oemol`, but any enhanced stereo is converted to unspecified.

    Returns:
      Representative SMILES.
    """
    if not already_representative:
        oemol = representative_oemol_from_oemol(oemol)
    if oemol.NumGroups() > 0:
        oemol = oemol.CreateCopy()
        convert_enhanced_stereo_to_unspecified(oemol)
    return smiles_from_oemol(oemol)


def representative_cxsmiles_from_oemol(
    oemol: oechem.OEMolBase, /, *, already_representative: bool = False
) -> str:
    """
    Get the representative CXSMILES of a molecule.

    Args:
      oemol: OpenEye molecule.
      already_representative: Whether the molecule is already in a representative form. If this is True, this is
        simply `cxsmiles_from_oemol`.

    Returns:
      Representative CXSMILES.
    """
    if not already_representative:
        oemol = representative_oemol_from_oemol(oemol)
    return cxsmiles_from_oemol(oemol)
