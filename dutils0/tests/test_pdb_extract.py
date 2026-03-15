"""Tests for dutils0.rough.pdb_extract."""

from pathlib import Path

import pytest

from dutils0.rough.pdb_extract import extract_selected_atoms_as_pdb


def test_extract_requires_residue_selector_or_atom_indices():
    """Must provide exactly one of residue_selector or atom_indices."""
    with pytest.raises(ValueError, match="Must provide exactly one"):
        extract_selected_atoms_as_pdb("/nonexistent.pdb")

    with pytest.raises(ValueError, match="Cannot provide both"):
        extract_selected_atoms_as_pdb(
            "/nonexistent.pdb",
            residue_selector="A:1",
            atom_indices=[1, 2],
        )


def test_extract_selected_atoms_as_pdb_accepts_6h41_cif_without_error():
    """Regression: CIF 6H41 loads and is processed without error; returns non-empty PDB string."""
    path = Path(__file__).resolve().parent / "data" / "6H41.cif"
    result = extract_selected_atoms_as_pdb(str(path), residue_selector="A:27")

    assert isinstance(result, str)
    assert len(result) > 0
    assert result.strip().startswith("ATOM ") or result.strip().startswith("HETATM ")
