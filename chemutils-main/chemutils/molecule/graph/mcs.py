from collections.abc import Generator

from openeye import oechem

from ..equal import AtomProperty, BondProperty
from ..smiles import oemol_from_smiles
from .properties import SubstructureAtomProperties, SubstructureBondProperties
from .query import get_query_from_oemol
from .utils import InvalidSMARTSError


def iterate_maximum_common_substructure_matches(
    *,
    pattern: oechem.OEMolBase | str,
    target: oechem.OEMolBase | str,
    atom_properties: AtomProperty | None = None,
    bond_properties: BondProperty | None = None,
    unique: bool = False,
    max_matches: int | None = None,
    min_atoms: int | None = None,
) -> Generator[oechem.OEMatchBase, None, None]:
    """
    Iterate over all maximum common substructure matches between a pattern molecule or SMARTS and a target molecule.

    By default, atoms are considered equal if their elements, aromaticity and formal charge match. Bonds
    are considered equal if their orders and aromaticity match.

    The `atom_properties` and `bond_properties` can be used to decide what properties are considered.

    Args:
        pattern: Pattern molecule or SMARTS, treated as the pattern in the match.
        target: Target molecule, treated as the target in the match.
        atom_properties: Atom properties to consider when matching. Cannot be used with a SMARTS string.
        bond_properties: Bond properties to consider when matching. Cannot be used with a SMARTS string.
        unique: If true, only unique matches are returned.
        max_matches: Maximum number of matches to return.
        min_atoms: Minimum number of matched atoms in each match.

    Yields:
        OpenEye matches between the pattern and target molecules.
    """
    if isinstance(pattern, str):
        if atom_properties is not None or bond_properties is not None:
            raise ValueError(
                "Cannot use `atom_properties` or `bond_properties` with SMARTS patterns."
            )
        mcss = oechem.OEMCSSearch(pattern)
        if not mcss.IsValid():
            raise InvalidSMARTSError(f"Invalid SMARTS pattern: {pattern}")
    else:
        query = get_query_from_oemol(
            pattern,
            atom_properties=atom_properties
            if atom_properties is not None
            else SubstructureAtomProperties,
            bond_properties=bond_properties
            if bond_properties is not None
            else SubstructureBondProperties,
        )
        query_idx_to_pattern_atom = {
            query_atom.GetIdx(): pattern_atom
            for query_atom, pattern_atom in zip(query.GetAtoms(), pattern.GetAtoms(), strict=True)
        }
        query_idx_to_pattern_bond = {
            query_bond.GetIdx(): pattern_bond
            for query_bond, pattern_bond in zip(query.GetBonds(), pattern.GetBonds(), strict=True)
        }
        mcss = oechem.OEMCSSearch(query)

    if isinstance(target, str):
        target = oemol_from_smiles(target)

    if max_matches is not None:
        mcss.SetMaxMatches(max_matches)
    if min_atoms is not None:
        mcss.SetMinAtoms(min_atoms)

    for match in mcss.Match(target, unique):
        if isinstance(pattern, oechem.OEMolBase):
            new_match = oechem.OEMatch()
            for pair in match.GetAtoms():
                new_match.AddPair(query_idx_to_pattern_atom[pair.pattern.GetIdx()], pair.target)
            for pair in match.GetBonds():
                new_match.AddPair(query_idx_to_pattern_bond[pair.pattern.GetIdx()], pair.target)
            match = new_match
        yield match


def has_maximum_common_substructure(
    *,
    pattern: oechem.OEMolBase | str,
    target: oechem.OEMolBase | str,
    atom_properties: AtomProperty | None = None,
    bond_properties: BondProperty | None = None,
    min_atoms: int | None = None,
) -> bool:
    return any(
        iterate_maximum_common_substructure_matches(
            pattern=pattern,
            target=target,
            atom_properties=atom_properties,
            bond_properties=bond_properties,
            unique=False,
            max_matches=1,
            min_atoms=min_atoms,
        )
    )
