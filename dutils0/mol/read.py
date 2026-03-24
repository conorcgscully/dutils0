"""Read PDB/mmCIF into a Gemmi Structure."""

from __future__ import annotations

import gemmi
from os import PathLike
from pathlib import Path


def read_protein_bytes(source: bytes | str | PathLike[str]) -> gemmi.Structure:
    """
    Read PDB or mmCIF data from a path or raw bytes and return a Gemmi Structure.

    Path-like input: path to a single file (.pdb, .cif, .ent; .gz allowed). Format is
    auto-detected. Bytes input: raw UTF-8 file content; format is inferred from content.

    Args:
        source: Path (str or PathLike) to a file, or raw file bytes.

    Returns:
        gemmi.Structure with at least one model (model → chain → residue → atom).

    Raises:
        FileNotFoundError: Path does not exist.
        OSError: I/O error when reading the file.
        ValueError: Path is not a regular file; invalid UTF-8 (bytes); parse failure;
            no _atom_site block (bytes CIF); or structure has no model.

    Example:
        st = read_protein_bytes("1abc.pdb")
        st = read_protein_bytes(Path("file.cif"))
        st = read_protein_bytes(open("file.cif", "rb").read())
    """
    if isinstance(source, bytes):
        return _read_from_bytes(source)
    return _read_from_path(source)


def _read_from_path(source: str | PathLike[str]) -> gemmi.Structure:
    path = Path(source).resolve()
    if not path.exists():
        raise FileNotFoundError(path)
    if not path.is_file():
        raise ValueError("not a file")
    try:
        structure = gemmi.read_structure(str(path))
    except RuntimeError as e:
        raise ValueError(f"Failed to parse PDB/CIF: {e}") from e
    _require_non_empty_structure(structure)
    return structure


def _require_non_empty_structure(structure: gemmi.Structure) -> None:
    """Raise ValueError if structure has no model or first model has no atoms."""
    if len(structure) < 1:
        raise ValueError("no model")
    if structure[0].count_atom_sites() == 0:
        raise ValueError("no model")


def _read_from_text(text: str) -> gemmi.Structure:
    """Parse PDB or mmCIF from a string. Format inferred from content."""
    # First non-empty, non-comment line determines format (CIF often has leading #)
    first_line = ""
    for line in text.splitlines():
        stripped_line = line.strip()
        if stripped_line and not stripped_line.startswith("#"):
            first_line = stripped_line
            break
    if first_line.startswith("data_") or first_line.startswith("_"):
        # CIF/mmCIF
        doc = gemmi.cif.read_string(text)
        block = None
        for b in doc:
            if b.find_loop("_atom_site.id"):
                block = b
                break
        if block is None:
            raise ValueError("no _atom_site block")
        try:
            structure = gemmi.make_structure_from_block(block)
        except RuntimeError as e:
            raise ValueError(f"Failed to parse PDB/CIF: {e}") from e
    else:
        # PDB
        try:
            structure = gemmi.read_pdb_string(text)
        except RuntimeError as e:
            raise ValueError(f"Failed to parse PDB/CIF: {e}") from e
    _require_non_empty_structure(structure)
    return structure


def _read_from_bytes(source: bytes) -> gemmi.Structure:
    try:
        text = source.decode("utf-8")
    except UnicodeDecodeError as e:
        raise ValueError(f"Invalid UTF-8: {e}") from e
    return _read_from_text(text)


def read_protein_str(source: str) -> gemmi.Structure:
    """
    Read PDB or mmCIF data from an in-memory string and return a Gemmi Structure.

    Format is inferred from content (first non-empty, non-comment line: data_ or _
    → CIF/mmCIF, else PDB). Same parsing and validation as read_protein_bytes.

    Args:
        source: PDB or mmCIF file content as a string.

    Returns:
        gemmi.Structure with at least one model (model → chain → residue → atom).

    Raises:
        ValueError: Empty or non-parseable content; parse failure; no _atom_site block
            (CIF); or structure has no model.

    Example:
        st = read_protein_str(open("file.cif").read())
    """
    return _read_from_text(source)
