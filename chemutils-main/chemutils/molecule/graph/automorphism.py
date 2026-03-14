from collections.abc import Generator

from openeye import oechem

from ..equal import AtomProperty, BondProperty
from .properties import AutomorphAtomProperties, AutomorphBondProperties
from .substructure import iterate_substructure_matches


def iterate_automorphism_matches(
    oemol: oechem.OEMolBase,
    /,
    *,
    atom_properties: AtomProperty = AutomorphAtomProperties,
    bond_properties: BondProperty = AutomorphBondProperties,
    max_matches: int | None = None,
) -> Generator[oechem.OEMatchBase, None, None]:
    """
    Iterate over all automorphisms where a molecule maps to itself.

    Unlike `OEGetAutomorphs`, this function is not arbitarily capped at 1024 automorphisms.

    This is a special case of substructure searching, where the molecules is considered a substructure
    of itself.

    By default, atoms are considered equal if their elements, aromaticity, ring membership and heavy degree
    match. Bonds are considered equal if their aromaticity match.

    The `atom_properties` and `bond_properties` can be used to decide what properties are considered.

    Args:
        oemol: Molecule to search for automorphisms.
        atom_properties: Atom properties to consider when matching.
        bond_properties: Bond properties to consider when matching.
        unique: Whether to only return unique matches.
        max_matches: Maximum number of matches to return.

    Yields:
        OpenEye matches corresponding to automorphisms of the molecule.
    """
    yield from iterate_substructure_matches(
        pattern=oemol,
        target=oemol,
        atom_properties=atom_properties,
        bond_properties=bond_properties,
        unique=False,
        max_matches=max_matches,
    )


def has_nontrivial_automorphism(
    oemol: oechem.OEMolBase,
    /,
    *,
    atom_properties: AtomProperty = AutomorphAtomProperties,
    bond_properties: BondProperty = AutomorphBondProperties,
) -> bool:
    matches = list(
        iterate_automorphism_matches(
            oemol, atom_properties=atom_properties, bond_properties=bond_properties, max_matches=2
        )
    )
    return len(matches) > 1
