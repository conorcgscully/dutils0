from collections.abc import Callable
from typing import overload

HEADER_RECORD_IN_ORDER: list[str | tuple[str, ...]] = [
    "HEADER",
    "OBSLTE",
    "TITLE ",
    "SPLIT ",
    "CAVEAT",
    "COMPND",
    "SOURCE",
    "KEYWDS",
    "EXPDTA",
    "NUMMDL",
    "MDLTYP",
    "AUTHOR",
    "REVDAT",
    "SPRSDE",
    "JRNL  ",
    # Remarks are ordered from "REMARK   1" to "REMARK 999"
    *(f"REMARK {i!s:>3}" for i in range(1, 1000)),
    "DBREF ",
    # Don't mess with ordering of DBREF1 and DBREF2
    ("DBREF1", "DBREF2"),
    "SEQADV",
    "SEQRES",
    "MODRES",
    "HET   ",
    "HETNAM",
    "HETSYN",
    "FORMUL",
    "HELIX ",
    "SHEET ",
    "SSBOND",
    "LINK  ",
    "CISPEP",
    "SITE  ",
    "CRYST1",
    # Likewise, don't mess with ordering of these three sets
    ("ORIGX1", "ORIGX2", "ORIGX3"),
    ("SCALE1", "SCALE2", "SCALE3"),
    ("MTRIX1", "MTRIX2", "MTRIX3"),
]
"""
List of PDB header records by order they should appear, with equal priority records in tuples.

The order of records in a PDB is defined here: https://www.wwpdb.org/documentation/file-format-content/format33/sect1.html.
"""

HEADER_RECORD_ORDER = {
    record: order
    for order, records in enumerate(HEADER_RECORD_IN_ORDER)
    for record in (records if isinstance(records, tuple) else (records,))
}
"""Mapping of PDB record types to their order (lower being earlier in the file).)"""

HEADER_RECORDS = set(HEADER_RECORD_ORDER.keys())
"""List of records considered to be in the header."""

BODY_RECORDS = {
    "MODEL ",
    "ATOM  ",
    "ANISOU",
    "TER   ",
    "HETATM",
    "ENDMDL",
    "CONECT",
    "MASTER",
    "END   ",
}
"""List of records that are not considered to be in the header."""


def get_record(line: str, /) -> str:
    """
    Get the record type (such as `DBREF `) from a line.

    For REMARK records, also returns the mark number, such as `REMARK  20`.

    Args:
        line: Line from a PDB file.

    Returns:
        Record type as a string.
    """
    if len(line) < 6:
        line = line.ljust(6)

    if line[:6] == "REMARK":
        return line[:10]
    return line[:6]


def get_line_order(line: str, /) -> int:
    """Get the order that a given line should appear in a PDB header."""
    return HEADER_RECORD_ORDER[get_record(line)]


@overload
def insert_lines_into_pdb(*, pdb: str, lines: list[str]) -> str: ...


@overload
def insert_lines_into_pdb(*, pdb: bytes, lines: list[str]) -> bytes: ...


@overload
def insert_lines_into_pdb(*, pdb: list[str], lines: list[str]) -> list[str]: ...


def insert_lines_into_pdb(
    *, pdb: str | bytes | list[str], lines: list[str]
) -> str | bytes | list[str]:
    """
    Insert one or more lines in a PDB header in the correct order.

    The order of records in a PDB header are defined here:
    https://www.wwpdb.org/documentation/file-format-content/format33/sect1.html

    Lines are inserted after existing lines of the same record type.

    If multiple lines are provided of the same record type (i.e. `DBREF `), their order is maintained.

    Args:
        pdb: PDB, as either `str` or `bytes` contents, or as a list of lines.
        lines: A list of lines to insert into the PDB.

    Returns:
        PDB with the lines inserted, in the same format as the input.
    """
    if len(lines) == 0:
        return pdb

    match pdb:
        case str():
            existing_lines = pdb.splitlines()
        case bytes():
            existing_lines = pdb.decode().splitlines()
        case list():
            existing_lines = pdb

    # Check all inserted lines are header lines (and hence have an order)
    for line in lines:
        if (record := get_record(line)) not in HEADER_RECORDS:
            raise ValueError(f"Cannot insert invalid PDB record `{record}`.")

    # Sort lines by their order
    # `sorted` is a stable sort, so lines with the same record order will be kept in the order they are provided
    lines = sorted(lines, key=lambda line: get_line_order(line))

    new_lines = []
    next_line_order: int | None = get_line_order(lines[0])
    in_header = True
    for existing_line in existing_lines:
        record = get_record(existing_line)
        if record in BODY_RECORDS:
            # If we get to the body, then we're done with the header
            # Insert all of our lines here
            in_header = False
            new_lines.extend(lines)
            lines = []
        elif record in HEADER_RECORDS:
            # If we come back to a header after seeing the body, something is wrong
            if not in_header:
                raise ValueError(
                    f"PDB has header record `{record}` following body records. Cannot insert."
                )
            # Get all our new lines that have lower priority than this header line, and insert them
            while next_line_order is not None and next_line_order < get_line_order(existing_line):
                new_lines.append(lines.pop(0))
                next_line_order = get_line_order(lines[0]) if len(lines) > 0 else None
        else:
            raise ValueError(f"PDB has unknown record `{record}`.")
        new_lines.append(existing_line)

    # If we reach the end and there's still lines left (and hence the PDB does not have a body)
    # Then add the rest of the lines to the end
    if len(lines) > 0:
        new_lines.extend(lines)

    match pdb:
        case str():
            return "\n".join(new_lines)
        case bytes():
            return "\n".join(new_lines).encode()
        case list():
            return new_lines


@overload
def filter_lines_in_pdb(*, pdb: str, filter: Callable[[str], bool]) -> str: ...


@overload
def filter_lines_in_pdb(*, pdb: bytes, filter: Callable[[str], bool]) -> bytes: ...


@overload
def filter_lines_in_pdb(*, pdb: list[str], filter: Callable[[str], bool]) -> list[str]: ...


def filter_lines_in_pdb(
    *, pdb: str | bytes | list[str], filter: Callable[[str], bool]
) -> str | bytes | list[str]:
    """
    Filter the lines in a PDB using a function.

    Args:
        pdb: PDB, as either `str` or `bytes` contents, or as a list of lines.
        filter: Function that takes a single line, and returned true if that line should be kept.

    Returns:
        PDB with only lines that pass the filter remaining.
    """
    match pdb:
        case str():
            existing_lines = pdb.splitlines()
        case bytes():
            existing_lines = pdb.decode().splitlines()
        case list():
            existing_lines = pdb

    lines = [line for line in existing_lines if filter(line)]

    match pdb:
        case str():
            return "\n".join(lines)
        case bytes():
            return "\n".join(lines).encode()
        case list():
            return lines
