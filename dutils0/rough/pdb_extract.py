"""
Extract a subset of atoms from a PDB or CIF file and return as a PDB-format string.

Uses BioPython's Bio.PDB to parse (PDBParser or MMCIFParser), filter via custom
Select subclasses, and write to a string via PDBIO and StringIO.
"""

from __future__ import annotations

from io import StringIO

from Bio.PDB import MMCIFParser, PDBIO, PDBParser, Select


class ResidueSelect(Select):
    """Select only atoms belonging to a single residue in a given chain."""

    def __init__(self, chain_id: str, res_num: int) -> None:
        self.chain_id = chain_id
        self.res_num = res_num

    def accept_chain(self, chain: object) -> bool:
        return chain.get_id() == self.chain_id

    def accept_residue(self, residue: object) -> bool:
        return (
            residue.get_parent().get_id() == self.chain_id
            and residue.id[1] == self.res_num
        )

    def accept_atom(self, atom: object) -> bool:
        return True


class AtomIndexSelect(Select):
    """Select only atoms whose PDB serial numbers are in the given set."""

    def __init__(self, indices_set: set[int]) -> None:
        self.indices_set = indices_set

    def accept_residue(self, residue: object) -> bool:
        return any(
            atom.get_serial_number() in self.indices_set
            for atom in residue.get_atoms()
        )

    def accept_atom(self, atom: object) -> bool:
        return atom.get_serial_number() in self.indices_set


def _structure_from_file(path: str):
    """Parse a PDB or CIF file and return a Bio.PDB Structure."""
    path_lower = path.lower()
    if path_lower.endswith(".cif") or path_lower.endswith(".cif.gz"):
        parser = MMCIFParser(QUIET=True)
    else:
        parser = PDBParser(QUIET=True)
    structure = parser.get_structure("structure", path)
    return structure


def _parse_residue_selector(selector: str) -> tuple[str, int]:
    """Parse 'CHAIN:RESNUM' into (chain_id, res_num)."""
    parts = [p.strip() for p in selector.split(":")]
    if len(parts) != 2:
        raise ValueError(
            f"Residue selector must be in form CHAIN:RESNUM (e.g. A:42), got {selector!r}"
        )
    chain_id, res_num_str = parts[0], parts[1]
    if len(chain_id) != 1:
        raise ValueError(
            f"Chain ID must be a single character, got {chain_id!r}"
        )
    try:
        res_num = int(res_num_str)
    except ValueError:
        raise ValueError(
            f"Residue number must be an integer, got {res_num_str!r}"
        ) from None
    return chain_id, res_num


def extract_selected_atoms_as_pdb(
    pdb_path: str,
    residue_selector: str | None = None,
    atom_indices: list[int] | None = None,
) -> str:
    """
    Read a PDB or CIF file and return a PDB-format string containing only selected atoms.

    Exactly one of residue_selector or atom_indices must be provided.
    Atom indices are PDB serial numbers from the original file (columns 7-11 of
    ATOM/HETATM lines). Duplicate serials in the file are treated as one.
    If the file has multiple models, only the first model is used.
    CIF input is detected by a path ending in .cif or .cif.gz.

    Args:
        pdb_path: Path to a PDB or CIF file (.pdb, .cif, .cif.gz, etc.).
        residue_selector: Optional selector in form "CHAIN:RESNUM" (e.g. "A:42").
            Chain ID must be a single character. Insertion codes are not supported;
            if multiple residues share the same sequence number (e.g. 42 and 42A),
            the first encountered in the chain is selected. Keeps only atoms
            belonging to that residue in that chain.
        atom_indices: Optional list of PDB atom serial numbers to keep. All
            requested indices must exist in the structure; if any are missing,
            ValueError is raised.

    Returns:
        PDB-formatted string containing only the selected atoms (first model only).
        Original atom serial numbers are preserved in the output.

    Raises:
        ValueError: If neither or both of residue_selector and atom_indices
            are provided; if the residue selector format is invalid; if the
            requested chain or residue is not found; or if no atoms match or
            any requested atom index is missing.

    Example:
        Residue selection (PDB or CIF)::

            pdb_str = extract_selected_atoms_as_pdb("/path/to/file.pdb", residue_selector="A:42")
            pdb_str = extract_selected_atoms_as_pdb("/path/to/file.cif", residue_selector="A:42")

        Atom-index selection::

            pdb_str = extract_selected_atoms_as_pdb("/path/to/file.pdb", atom_indices=[1, 2, 3, 10])
    """
    if residue_selector is None and atom_indices is None:
        raise ValueError(
            "Must provide exactly one of residue_selector or atom_indices."
        )
    if residue_selector is not None and atom_indices is not None:
        raise ValueError(
            "Cannot provide both residue_selector and atom_indices."
        )

    structure = _structure_from_file(pdb_path)
    model = structure[0]

    if residue_selector is not None:
        chain_id, res_num = _parse_residue_selector(residue_selector)
        chain = None
        for c in model.get_chains():
            if c.get_id() == chain_id:
                chain = c
                break
        if chain is None:
            raise ValueError(
                f"Chain {chain_id!r} not found in structure."
            )
        residue_found = None
        for res in chain.get_residues():
            if res.id[1] == res_num:
                residue_found = res
                break
        if residue_found is None:
            raise ValueError(
                f"Residue {res_num} not found in chain {chain_id!r}."
            )
        select = ResidueSelect(chain_id=chain_id, res_num=res_num)
    else:
        indices_set = set(atom_indices)
        structure_serials = {
            atom.get_serial_number()
            for atom in model.get_atoms()
        }
        if indices_set.isdisjoint(structure_serials):
            raise ValueError(
                "No atoms in structure match the requested atom indices."
            )
        missing = indices_set - structure_serials
        if missing:
            raise ValueError(
                f"Requested atom indices not found in structure: {sorted(missing)}."
            )
        select = AtomIndexSelect(indices_set=indices_set)

    buffer = StringIO()
    io = PDBIO()
    io.set_structure(model)
    io.save(buffer, select=select, preserve_atom_numbering=True)
    return buffer.getvalue()
