from collections.abc import Generator

from openeye import oechem

from ..equal import AtomProperty, BondProperty
from .properties import AutomorphAtomProperties, AutomorphBondProperties
from .substructure import iterate_substructure_matches


def iterate_homomorphism_matches(
    *,
    pattern: oechem.OEMolBase,
    target: oechem.OEMolBase,
    atom_properties: AtomProperty = AutomorphAtomProperties,
    bond_properties: BondProperty = AutomorphBondProperties,
    max_matches: int | None = None,
) -> Generator[oechem.OEMatchBase, None, None]:
    """
    Iterate over all homomorphisms between a pattern molecule and a target molecule.

    This is a special case of substructure searching, where the source molecule is a complete substructure
    of the target molecule.

    By default, atoms are considered equal if their elements, aromaticity, ring membership and heavy degree
    match. Bonds are considered equal if their aromaticity match.

    The `atom_properties` and `bond_properties` can be used to decide what properties are considered.

    Args:
        pattern: Source molecule.
        target: Target molecule.
        atom_properties: Atom properties to consider when matching.
        bond_properties: Bond properties to consider when matching.
        max_matches: Maximum number of matches to return.

    Yields:
        OpenEye matches between the pattern and target molecules.
    """
    if pattern.NumAtoms() != target.NumAtoms():
        return
    if pattern.NumBonds() != target.NumBonds():
        return
    yield from iterate_substructure_matches(
        pattern=pattern,
        target=target,
        atom_properties=atom_properties,
        bond_properties=bond_properties,
        unique=False,
        max_matches=max_matches,
    )


def has_homomorphism(
    *,
    pattern: oechem.OEMolBase,
    target: oechem.OEMolBase,
    atom_properties: AtomProperty = AutomorphAtomProperties,
    bond_properties: BondProperty = AutomorphBondProperties,
) -> bool:
    return any(
        iterate_homomorphism_matches(
            pattern=pattern,
            target=target,
            atom_properties=atom_properties,
            bond_properties=bond_properties,
            max_matches=1,
        )
    )
