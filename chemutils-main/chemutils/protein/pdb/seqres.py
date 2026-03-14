import math
import warnings
from collections.abc import Iterable

SEQRES_SUBSEQUENT_RECORDS = (
    "HET",
    "FORMUL",
    "HETNAM",
    "HETSYN",
    "HELIX",
    "SHEET",
    "SSBOND",
    "LINK",
    "CISPEP",
    "SITE",
    "CRYST1",
    "MTRIX",
    "ORIGX",
    "SCALE",
    "MODEL",
    "ATOM",
    "ANISOU",
    "TER",
    "HETATM",
    "ENDMDL",
    "CONECT",
    "MASTER",
    "END",
)
"""
PDB Record names that follow SEQRES, according to
https://www.wwpdb.org/documentation/file-format-content/format33/v3.3.html
"""


def parse_seqres_lines(lines: Iterable[str], /) -> dict[str, list[str]]:
    """
    Read sequences from SEQRES lines, as mappings of chain ID to list of residue names.

    Args:
        lines: List of SEQRES lines from the PDB.

    Returns:
        Dictionary mapping chain IDs to lists of residue names.

    Raises:
        ValueError: Line does not start with SEQRES.
    """
    sequences: dict[str, list[str]] = {}
    for line in lines:
        if not line.startswith("SEQRES"):
            raise ValueError(f"Line `{line}` is not a valid SEQRES line.")
        chain_id = line[11]
        residue_names = line[19:].split()
        sequences[chain_id] = sequences.get(chain_id, []) + residue_names
    return sequences


def read_seqres_from_pdb(pdb: str | bytes, /) -> dict[str, list[str]]:
    """
    Read PDB sequences from SEQRES, as mappings of chain ID to list of residue names.

    Args:
        pdb: PDB contents as either string or bytes.

    Returns:
        Dictionary mapping chain IDs to lists of residue names.
    """
    if isinstance(pdb, bytes):
        pdb = pdb.decode()
    return parse_seqres_lines(line for line in pdb.splitlines() if line.startswith("SEQRES"))


SEQRES_FORMAT_STRING = "SEQRES {:>3} {:>1} {:>4}  {:>3} {:>3} {:>3} {:>3} {:>3} {:>3} {:>3} {:>3} {:>3} {:>3} {:>3} {:>3} {:>3}          "
"""
SEQRES format string, with correct column spacings and justifications.
See https://www.wwpdb.org/documentation/file-format-content/format33/sect3.html for more details.
"""


def generate_seqres(*, sequence_residue_names: dict[str, list[str]]) -> list[str]:
    """
    Generate SEQRES records for a set of sequences.

    Args:
        sequence_residue_names: Sequences as a mapping of chain IDs to a list of residues.

    Returns:
        List of lines representing the corresponding SEQRES records.
    """
    lines = []
    for chain_id, residue_names in sequence_residue_names.items():
        sequence_length = len(residue_names)
        # Pad out to multiple of 13
        padded_residue_names = residue_names + [""] * (
            13 * math.ceil(sequence_length / 13) - sequence_length
        )
        for line_index in range(len(padded_residue_names) // 13):
            line = SEQRES_FORMAT_STRING.format(
                line_index + 1,
                chain_id,
                sequence_length,
                *padded_residue_names[line_index * 13 : (line_index + 1) * 13],
            )
            lines.append(line)
    return lines


def add_seqres_to_pdb(*, pdb: str, sequence_residue_names: dict[str, list[str]]) -> str:
    """
    Add or update the SEQRES records in a PDB.

    If SEQRES records already exist, they are overwritten.

    Args:
        pdb: PDB contents as a string.
        sequence_residue_names: Sequences as a mapping of chain IDs to list of residue names.

    Returns:
        PDB with the SEQRES records inserted.
    """
    warnings.warn(
        DeprecationWarning(
            "Use `pdb = filter_lines_in_pdb(pdb=pdb, filter=lambda line: not line.startswith('SEQRES'))` and `pdb = insert_lines_into_pdb(pdb=pdb, lines=generate_seqres(sequence_residue_names=...)` instead."
        ),
        stacklevel=2,
    )
    new_lines = []
    added = False
    for line in pdb.splitlines(keepends=True):
        if line.startswith("SEQRES"):
            continue
        if not added and line.startswith(SEQRES_SUBSEQUENT_RECORDS):
            added = True
            for seqres_line in generate_seqres(sequence_residue_names=sequence_residue_names):
                new_lines.append(seqres_line + "\n")
        new_lines.append(line)
    return "".join(new_lines)
