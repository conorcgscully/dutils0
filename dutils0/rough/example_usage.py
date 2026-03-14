"""
Example usage of extract_selected_atoms_as_pdb.

Run from repo root with a PDB file path, e.g.:

    python -m rough.example_usage path/to/structure.pdb

Or import and call:

    from rough import extract_selected_atoms_as_pdb

    # Residue selection: keep only chain A, residue 42
    pdb_str = extract_selected_atoms_as_pdb("file.pdb", residue_selector="A:42")

    # Atom-index selection: keep only atoms with PDB serials 1, 2, 3, 10
    pdb_str = extract_selected_atoms_as_pdb("file.pdb", atom_indices=[1, 2, 3, 10])
"""

from __future__ import annotations

import sys

from rough import extract_selected_atoms_as_pdb


def main() -> None:
    if len(sys.argv) < 2:
        print("Usage: python -m rough.example_usage <pdb_path>")
        sys.exit(1)
    pdb_path = sys.argv[1]

    print("Residue selection (A:42):")
    pdb_res = extract_selected_atoms_as_pdb(pdb_path, residue_selector="A:42")
    print(pdb_res[:500] + "..." if len(pdb_res) > 500 else pdb_res)

    print("\nAtom-index selection [1, 2, 3, 10]:")
    pdb_atoms = extract_selected_atoms_as_pdb(pdb_path, atom_indices=[1, 2, 3, 10])
    print(pdb_atoms[:500] + "..." if len(pdb_atoms) > 500 else pdb_atoms)


if __name__ == "__main__":
    main()
