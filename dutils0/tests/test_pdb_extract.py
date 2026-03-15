"""Tests for dutils0.rough.pdb_extract."""

from pathlib import Path

import pytest

from dutils0.rough.pdb_extract import extract_selected_atoms_as_pdb


def _path_6h41_cif() -> Path:
    """Path to 6H41.cif, cwd-independent (relative to this test file)."""
    return Path(__file__).resolve().parent / "data" / "6H41.cif"


@pytest.fixture
def path_6h41_cif():
    """Fixture: absolute path to tests/data/6H41.cif."""
    return _path_6h41_cif()


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


def _assert_plausible_pdb_string(result: str, context: str = "") -> None:
    """Assert the result is a non-empty string that looks like PDB output."""
    prefix = f"{context}: " if context else ""
    assert isinstance(result, str), f"{prefix}expected str, got {type(result).__name__}"
    stripped = result.strip()
    assert len(stripped) > 0, f"{prefix}output is empty or whitespace-only"
    has_atom = "ATOM " in result or "HETATM " in result
    assert has_atom, (
        f"{prefix}output should contain ATOM or HETATM record markers; "
        f"first 200 chars: {result[:200]!r}"
    )


def test_extract_6h41_cif_residue_selector(path_6h41_cif):
    """6H41.cif: residue selector A:27 runs without error and returns valid-looking PDB."""
    path = path_6h41_cif
    assert path.exists(), f"test data missing: {path}"

    result = extract_selected_atoms_as_pdb(str(path), residue_selector="A:27")

    _assert_plausible_pdb_string(result, "residue_selector A:27")


def test_extract_6h41_cif_atom_indices(path_6h41_cif):
    """6H41.cif: atom_indices [1, 2, 3] runs without error and returns valid-looking PDB."""
    path = path_6h41_cif
    assert path.exists(), f"test data missing: {path}"

    result = extract_selected_atoms_as_pdb(str(path), atom_indices=[1, 2, 3])

    _assert_plausible_pdb_string(result, "atom_indices [1,2,3]")
