from openeye import oechem

from .canonicalize import canonicalize_molecule
from .charge import remove_formal_charges
from .hydrogens import make_hydrogens_implicit
from .perceive import perceive_atom_properties


def standardise_oemol(mol: oechem.OEGraphMol, /) -> None:
    """
    Standardise an OEGraphMol in-place.

    This does the following:
    * Remove all formal charges.
    * Removes all hydrogens, include isotopic hydrogen.
    * Assigns aromaticity based on the OpenEye molecule.
    * Perceives stereochemistry based on three-dimensional structure.
    * Order the atoms based on the order they appear in a canonical isomeric SMILES.
    * Canonically order the bonds based on the atoms.
    * Reindex both atoms and bonds sequentially.
    * Ensure all bonds are ordered so their start index is less than their end index.

    Args:
        mol: Molecule to process (in-place).
    """
    oechem.OEDeleteEverythingExceptTheFirstLargestComponent(mol)

    remove_formal_charges(mol)
    make_hydrogens_implicit(mol, remove_isotopic_hydrogens=True)
    perceive_atom_properties(mol)
    canonicalize_molecule(mol)
