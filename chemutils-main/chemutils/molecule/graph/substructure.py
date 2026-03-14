from collections.abc import Generator

import numpy as np
import numpy.typing as npt
from openeye import oechem

from ..equal import AtomProperty, BondProperty
from ..smiles import oemol_from_smiles
from .properties import SubstructureAtomProperties, SubstructureBondProperties
from .query import get_query_from_oemol
from .utils import InvalidSMARTSError, get_atom_mapping_from_match, get_target_mask_from_match


def iterate_substructure_matches(
    *,
    pattern: oechem.OEMolBase | str,
    target: oechem.OEMolBase | str,
    atom_properties: AtomProperty | None = None,
    bond_properties: BondProperty | None = None,
    unique: bool = False,
    max_matches: int | None = None,
) -> Generator[oechem.OEMatchBase, None, None]:
    """
    Iterate over all substructure matches where a pattern molecule or SMARTS is part of a target molecules.

    By default, atoms are considered equal if their elements, aromaticity and formal charge match. Bonds
    are considered equal if their orders and aromaticity match.

    The `atom_properties` and `bond_properties` can be used to decide what properties are considered.

    Args:
        pattern: Pattern molecule or SMARTS, treated as the pattern in the match.
        target: Target molecule, treated as the target in the match.
        atom_properties: Atom properties to consider when matching. Cannot be used with a SMARTS string.
        bond_properties: Bond properties to consider when matching. Cannot be used with a SMARTS string.
        unique: If true, only unique matches are returned. For example, searching for a benzene ring inside
            a phenol, if `unique=True` there is only one ring. If `unique=False`, then technically there are
            12 different ways of mapping the same benzene ring due to rotations and reflections.
        max_matches: Maximum number of matches to return.

    Yields:
        OpenEye matches where the pattern molecule is a substructure of the target molecules.
    """
    if isinstance(pattern, str):
        if atom_properties is not None or bond_properties is not None:
            raise ValueError(
                "Cannot use `atom_properties` or `bond_properties` with SMARTS patterns."
            )

        subsearch = oechem.OESubSearch(pattern)
        if not subsearch.IsValid():
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
        subsearch = oechem.OESubSearch(query)

    if isinstance(target, str):
        target = oemol_from_smiles(target)

    if max_matches:
        subsearch.SetMaxMatches(max_matches)
    else:
        subsearch.SetMaxMatches(0)

    for match in subsearch.Match(target, unique):
        if isinstance(pattern, oechem.OEMolBase):
            new_match = oechem.OEMatch()
            for pair in match.GetAtoms():
                new_match.AddPair(query_idx_to_pattern_atom[pair.pattern.GetIdx()], pair.target)
            for pair in match.GetBonds():
                new_match.AddPair(query_idx_to_pattern_bond[pair.pattern.GetIdx()], pair.target)
            match = new_match
        yield match


def has_substructure(
    pattern: oechem.OEMolBase | str,
    target: oechem.OEMolBase | str,
    atom_properties: AtomProperty | None = None,
    bond_properties: BondProperty | None = None,
) -> bool:
    return any(
        iterate_substructure_matches(
            pattern=pattern,
            target=target,
            atom_properties=atom_properties,
            bond_properties=bond_properties,
            unique=False,
            max_matches=1,
        )
    )


def get_substructure_indices(
    *,
    pattern: oechem.OEMolBase | str,
    target: oechem.OEMolBase | str,
    atom_properties: AtomProperty | None = None,
    bond_properties: BondProperty | None = None,
    max_matches: int | None = None,
) -> Generator[int, None, None]:
    """
    Iterate atom indices in a query (target) molecule that match a reference (pattern) molecule or SMARTS string.

    Args:
        pattern: Pattern molecule or SMARTS, treated as the pattern in the match.
        target: Target molecule, treated as the target in the match.
        atom_properties: Atom properties to consider when matching. Cannot be used with a SMARTS string.
        bond_properties: Bond properties to consider when matching. Cannot be used with a SMARTS string.
        max_matches: Maximum number of matches to return.

    Yields:
        Indices of atoms in the molecule which are matched by one or more substructures.
    """
    for match in iterate_substructure_matches(
        pattern=pattern,
        target=target,
        atom_properties=atom_properties,
        bond_properties=bond_properties,
        unique=True,
        max_matches=max_matches,
    ):
        for atom_base in get_atom_mapping_from_match(match).values():
            yield atom_base.GetIdx()


def get_substructure_counts(
    *,
    pattern: oechem.OEMolBase | str,
    target: oechem.OEMolBase | str,
    atom_properties: AtomProperty | None = None,
    bond_properties: BondProperty | None = None,
    unique: bool = False,
    max_matches: int | None = None,
) -> npt.NDArray[np.int_]:
    """
    Get an array of counts of substructure matches in a molecule.

    Args:
        pattern: Pattern molecule or SMARTS, treated as the pattern in the match.
        target: Target molecule, treated as the target in the match.
        atom_properties: Atom properties to consider when matching. Cannot be used with a SMARTS string.
        bond_properties: Bond properties to consider when matching. Cannot be used with a SMARTS string.
        unique: If true, only unique matches are returned. For example, searching for a benzene ring inside
            a phenol, if `unique=True` there is only one ring. If `unique=False`, then technically there are
            12 different ways of mapping the same benzene ring due to rotations and reflections.
        max_matches: Maximum number of matches to return.

    Returns:
        Array of counts of substructure matches in the molecule.
    """
    if isinstance(target, str):
        target = oemol_from_smiles(target)

    output_counts = np.zeros(target.NumAtoms(), dtype=int)
    for match in iterate_substructure_matches(
        pattern=pattern,
        target=target,
        atom_properties=atom_properties,
        bond_properties=bond_properties,
        unique=unique,
        max_matches=max_matches,
    ):
        target_mask = get_target_mask_from_match(match)
        output_counts[target_mask] += 1
    return output_counts


def get_substructure_mask(
    *,
    pattern: oechem.OEMolBase | str,
    target: oechem.OEMolBase | str,
    atom_properties: AtomProperty | None = None,
    bond_properties: BondProperty | None = None,
    unique: bool = False,
    max_matches: int | None = None,
) -> npt.NDArray[np.bool_]:
    """
    Get a boolean mask of atoms in a molecule that match a given SMARTS pattern.

    Args:
        pattern: Pattern molecule or SMARTS, treated as the pattern in the match.
        target: Target molecule, treated as the target in the match.
        atom_properties: Atom properties to consider when matching. Cannot be used with a SMARTS string.
        bond_properties: Bond properties to consider when matching. Cannot be used with a SMARTS string.
        unique: If true, only unique matches are returned. For example, searching for a benzene ring inside
            a phenol, if `unique=True` there is only one ring. If `unique=False`, then technically there are
            12 different ways of mapping the same benzene ring due to rotations and reflections.
        max_matches: Maximum number of matches to return.

    Returns:
        Boolean mask of atoms in the molecule that match the SMARTS pattern.
    """
    count_mask = get_substructure_counts(
        pattern=pattern,
        target=target,
        atom_properties=atom_properties,
        bond_properties=bond_properties,
        unique=unique,
        max_matches=max_matches,
    )
    return count_mask > 0
