from __future__ import annotations

from collections.abc import Generator

from openeye import oechem
from typing_extensions import deprecated

from .graph.substructure import (
    get_substructure_counts,
    get_substructure_indices,
    get_substructure_mask,
)

__all__ = [
    "get_substructure_counts",
    "get_substructure_indices",
    "get_substructure_mask",
]


class InvalidSMARTSError(ValueError):
    """Error thrown when an invalid SMARTS string is parsed."""

    pass


def _construct_subsearch(*, smarts: str) -> oechem.OESubSearch:
    sub_search = oechem.OESubSearch()
    sub_search.SetMaxMatches(0)  # Set to zero to disable thresholding of matches
    if not sub_search.Init(smarts):
        raise InvalidSMARTSError(f"Invalid SMARTS string: {smarts}")
    return sub_search


@deprecated(
    "`has_substructure_match` is deprecated. Use `iterate_substructure_matches` and check it returns any matches."
)
def has_substructure_match(*, smarts: str, oemol: oechem.OEGraphMol) -> bool:
    """
    Check if a given SMARTS string is matched by a given molecule.

    Args:
        smarts: SMARTS string.
        oemol: OpenEye molecule.

    Returns:
        True if the SMARTS pattern is matched in `oemol`, and False otherwise.

    Raises:
        InvalidSMARTSError: Provided SMARTS string is invalid.
    """
    sub_search = _construct_subsearch(smarts=smarts)
    oechem.OEPrepareSearch(oemol, sub_search)
    return sub_search.SingleMatch(oemol)  # type: ignore


@deprecated("`get_substructure_matches` is deprecated. Use `iterate_substructure_matches`.")
def get_substructure_matches(
    *, smarts: str, oemol: oechem.OEGraphMol
) -> Generator[oechem.OEMatchBase, None, None]:
    """
    Iterate over OpenEye matches for a given substructure (specified as a SMARTS).

    Args:
        smarts: SMARTS string.
        oemol: OpenEye molecule.

    Yields:
        OpenEye `OEMatchBase` objects.

    Raises:
        InvalidSMARTSError: Provided SMARTS string is invalid.
    """
    sub_search = _construct_subsearch(smarts=smarts)
    oechem.OEPrepareSearch(oemol, sub_search)
    match_iter = sub_search.Match(oemol, True)
    yield from match_iter
