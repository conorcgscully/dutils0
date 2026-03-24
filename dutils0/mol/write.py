"""Write Gemmi Structure to PDB or mmCIF."""

from __future__ import annotations

import gemmi

_ENCODING = "utf-8"


def write_protein_bytes(structure: gemmi.Structure, fmt: str) -> bytes:
    """
    Serialize a Gemmi Structure to PDB or mmCIF as UTF-8 bytes. No file I/O.

    Format (case-insensitive): "pdb" → PDB; "cif" or "mmcif" → mmCIF. Same format
    values as write_protein_str. Unsupported format raises ValueError.

    Args:
        structure: The Gemmi Structure to serialize (model → chain → residue → atom).
        fmt: Output format: "pdb", "cif", or "mmcif".

    Returns:
        bytes: PDB or mmCIF file content as UTF-8-encoded bytes.

    Raises:
        TypeError: structure is not a gemmi.Structure.
        ValueError: Unsupported format value or serialization failure.

    Example:
        raw = write_protein_bytes(st, "pdb")
        raw = write_protein_bytes(st, "cif")
    """
    return write_protein_str(structure, fmt).encode(_ENCODING)


def write_protein_str(structure: gemmi.Structure, fmt: str) -> str:
    """
    Serialize a Gemmi Structure to a PDB or mmCIF string. Format is chosen by the fmt argument.

    Format (case-insensitive): "pdb" → PDB text; "cif" or "mmcif" → mmCIF text. No file I/O.

    Args:
        structure: The Gemmi Structure to serialize (model → chain → residue → atom).
        fmt: Output format: "pdb", "cif", or "mmcif".

    Returns:
        str: PDB or mmCIF file content as a string.

    Raises:
        TypeError: structure is not a gemmi.Structure.
        ValueError: Unsupported format value or serialization failure.

    Example:
        s = write_protein_str(st, "pdb")
        s = write_protein_str(st, "cif")
    """
    if not isinstance(structure, gemmi.Structure):
        raise TypeError(
            f"structure must be a gemmi.Structure, got {type(structure).__name__}"
        )
    normalized = fmt.strip().lower()
    if normalized not in ("pdb", "cif", "mmcif"):
        raise ValueError(
            f"Unsupported format: {normalized!r}. Use 'pdb', 'cif', or 'mmcif'."
        )
    try:
        if normalized == "pdb":
            return structure.make_pdb_string()
        return structure.make_mmcif_document().as_string()
    except RuntimeError as e:
        raise ValueError(f"Failed to serialize structure: {e}") from e
