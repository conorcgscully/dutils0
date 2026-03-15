"""Tests for dutils0.rough.read_protein_bytes, read_protein_str, gemmi_from_rdkit, gemmi_to_rdkit, gemmi_to_pymol, and protein_to_rdkit."""

from pathlib import Path
from unittest.mock import MagicMock, patch

import gemmi
import pytest
from rdkit import Chem

from dutils0.rough import (
    gemmi_from_rdkit,
    gemmi_to_pymol,
    gemmi_to_rdkit,
    protein_to_rdkit,
    read_protein_bytes,
    read_protein_str,
    write_protein_bytes,
    write_protein_str,
)


def _path_6h41_cif() -> Path:
    """Path to 6H41.cif, cwd-independent (relative to this test file)."""
    return Path(__file__).resolve().parent / "data" / "6H41.cif"


def _path_pdb() -> Path:
    """Path to a .pdb test file if present."""
    return Path(__file__).resolve().parent / "data" / "test.pdb"


@pytest.fixture
def path_6h41_cif():
    """Fixture: absolute path to tests/data/6H41.cif."""
    return _path_6h41_cif()


def _assert_valid_structure(st: gemmi.Structure, context: str = "") -> None:
    """Assert st is a non-empty Gemmi Structure with iterable model/chain/residue/atom."""
    prefix = f"{context}: " if context else ""
    assert isinstance(st, gemmi.Structure), f"{prefix}expected gemmi.Structure, got {type(st)}"
    assert len(st) >= 1, f"{prefix}expected at least one model, got len(structure)={len(st)}"
    model = st[0]
    n = 0
    for chain in model:
        for residue in chain:
            for atom in residue:
                n += 1
                if n >= 1:
                    break
            if n >= 1:
                break
        if n >= 1:
            break
    assert n >= 1, f"{prefix}expected at least one atom in structure"


def test_read_protein_bytes_path_mmcif(path_6h41_cif):
    """Path (mmCIF): read_protein_bytes(path) and read_protein_bytes(Path(path)) return valid Structure."""
    path = path_6h41_cif
    assert path.exists(), f"test data missing: {path}"
    st_str = read_protein_bytes(str(path))
    _assert_valid_structure(st_str, "str path")
    st_path = read_protein_bytes(path)
    _assert_valid_structure(st_path, "Path path")


@pytest.mark.skipif(not _path_pdb().exists(), reason="no .pdb test file")
def test_read_protein_bytes_path_pdb():
    """Path (PDB): same invariants for .pdb file if present."""
    path = _path_pdb()
    st = read_protein_bytes(str(path))
    _assert_valid_structure(st)
    st2 = read_protein_bytes(path)
    _assert_valid_structure(st2)


def test_read_protein_bytes_bytes(path_6h41_cif):
    """Bytes: same file content as bytes returns same invariants as path."""
    path = path_6h41_cif
    assert path.exists(), f"test data missing: {path}"
    raw = path.read_bytes()
    st = read_protein_bytes(raw)
    _assert_valid_structure(st, "bytes")


def test_read_protein_bytes_missing_path():
    """Missing path raises FileNotFoundError."""
    with pytest.raises(FileNotFoundError):
        read_protein_bytes("/nonexistent.pdb")


def test_read_protein_bytes_path_is_directory(tmp_path):
    """Path that is a directory raises ValueError."""
    with pytest.raises(ValueError, match="not a file"):
        read_protein_bytes(tmp_path)


def test_read_protein_bytes_bad_content():
    """Bad or empty content raises ValueError."""
    with pytest.raises(ValueError):
        read_protein_bytes(b"not a structure")
    with pytest.raises(ValueError):
        read_protein_bytes(b"")


def test_read_protein_bytes_invalid_utf8():
    """Invalid UTF-8 bytes raise ValueError."""
    with pytest.raises(ValueError, match="UTF-8"):
        read_protein_bytes(b"\xff\xfe\x00")


def test_read_protein_bytes_cif_no_atom_site():
    """CIF bytes with no _atom_site block raise ValueError."""
    # Minimal CIF: one block, no _atom_site (e.g. dictionary-only block).
    minimal_cif = b"data_dummy\n_entry.id dummy\n"
    with pytest.raises(ValueError, match="no _atom_site block"):
        read_protein_bytes(minimal_cif)


# --- read_protein_str ---


def test_read_protein_str_mmcif(path_6h41_cif):
    """read_protein_str with mmCIF string returns valid Structure."""
    path = path_6h41_cif
    assert path.exists(), f"test data missing: {path}"
    text = path.read_text()
    st = read_protein_str(text)
    _assert_valid_structure(st, "read_protein_str mmCIF")


@pytest.mark.skipif(not _path_pdb().exists(), reason="no .pdb test file")
def test_read_protein_str_pdb():
    """read_protein_str with PDB string returns valid Structure when test.pdb exists."""
    path = _path_pdb()
    text = path.read_text()
    st = read_protein_str(text)
    _assert_valid_structure(st, "read_protein_str PDB")


def test_read_protein_str_bad_content():
    """Bad or empty string raises ValueError."""
    with pytest.raises(ValueError):
        read_protein_str("not a structure")
    with pytest.raises(ValueError):
        read_protein_str("")


def test_read_protein_str_cif_no_atom_site():
    """CIF string with no _atom_site block raises ValueError."""
    minimal_cif = "data_dummy\n_entry.id dummy\n"
    with pytest.raises(ValueError, match="no _atom_site block"):
        read_protein_str(minimal_cif)


def test_read_protein_str_same_result_as_bytes(path_6h41_cif):
    """read_protein_str(text) matches read_protein_bytes(text.encode())."""
    path = path_6h41_cif
    assert path.exists(), f"test data missing: {path}"
    text = path.read_text()
    st_str = read_protein_str(text)
    st_bytes = read_protein_bytes(text.encode("utf-8"))
    _assert_valid_structure(st_str, "str")
    _assert_valid_structure(st_bytes, "bytes")
    assert len(st_str) == len(st_bytes)
    assert st_str[0].count_atom_sites() == st_bytes[0].count_atom_sites()


# --- write_protein_bytes ---


def test_write_protein_bytes_pdb(path_6h41_cif):
    """write_protein_bytes(st, "pdb") returns UTF-8 bytes; round-trip preserves model/atom count."""
    st = read_protein_bytes(path_6h41_cif)
    result = write_protein_bytes(st, "pdb")
    assert isinstance(result, bytes), "return type must be bytes"
    assert len(result) > 0
    assert b"ATOM " in result or b"HETATM " in result
    back = read_protein_bytes(result)
    assert len(back) == len(st), "model count mismatch after round-trip"
    assert back[0].count_atom_sites() == st[0].count_atom_sites(), "atom count mismatch after round-trip"


def test_write_protein_bytes_cif(path_6h41_cif):
    """write_protein_bytes(st, "cif") returns UTF-8 mmCIF bytes; round-trip preserves model/atom count."""
    st = read_protein_bytes(path_6h41_cif)
    result = write_protein_bytes(st, "cif")
    assert isinstance(result, bytes), "return type must be bytes"
    assert len(result) > 0
    text = result.decode("utf-8")
    assert "data_" in text or text.strip().startswith("_") or "_atom_site" in text
    back = read_protein_bytes(result)
    assert len(back) == len(st)
    assert back[0].count_atom_sites() == st[0].count_atom_sites()


def test_write_protein_bytes_mmcif(path_6h41_cif):
    """write_protein_bytes(st, "mmcif") returns same mmCIF bytes as "cif"."""
    st = read_protein_bytes(path_6h41_cif)
    out_cif = write_protein_bytes(st, "cif")
    out_mmcif = write_protein_bytes(st, "mmcif")
    assert out_cif == out_mmcif
    back = read_protein_bytes(out_mmcif)
    _assert_valid_structure(back, "mmcif round-trip")


def test_write_protein_bytes_unsupported_format(path_6h41_cif):
    """Unsupported format raises ValueError with message mentioning format and allowed values."""
    st = read_protein_bytes(path_6h41_cif)
    with pytest.raises(ValueError, match="Unsupported format.*Use 'pdb'"):
        write_protein_bytes(st, "xyz")
    with pytest.raises(ValueError, match="Unsupported format"):
        write_protein_bytes(st, "")


def test_write_protein_bytes_invalid_structure_type():
    """Non-Structure (None or wrong type) raises TypeError."""
    with pytest.raises(TypeError, match="structure must be a gemmi.Structure.*None"):
        write_protein_bytes(None, "pdb")
    with pytest.raises(TypeError, match="structure must be a gemmi.Structure.*str"):
        write_protein_bytes("not a structure", "pdb")


def test_write_protein_bytes_case_insensitive(path_6h41_cif):
    """Format "PDB" and "CIF" are accepted and produce same bytes as lowercase."""
    st = read_protein_bytes(path_6h41_cif)
    assert write_protein_bytes(st, "PDB") == write_protein_bytes(st, "pdb")
    assert write_protein_bytes(st, "CIF") == write_protein_bytes(st, "cif")


# --- write_protein_str ---


def test_write_protein_str_pdb(path_6h41_cif):
    """write_protein_str(st, "pdb") returns non-empty PDB string; round-trip preserves model/atom count."""
    st = read_protein_bytes(path_6h41_cif)
    result = write_protein_str(st, "pdb")
    assert isinstance(result, str)
    assert len(result) > 0
    assert "ATOM " in result or "HETATM " in result
    back = read_protein_str(result)
    assert len(back) == len(st)
    assert back[0].count_atom_sites() == st[0].count_atom_sites()


def test_write_protein_str_cif(path_6h41_cif):
    """write_protein_str(st, "cif") returns non-empty mmCIF string; round-trip preserves model/atom count."""
    st = read_protein_bytes(path_6h41_cif)
    result = write_protein_str(st, "cif")
    assert isinstance(result, str)
    assert len(result) > 0
    assert "data_" in result or result.strip().startswith("_") or "_atom_site" in result
    back = read_protein_str(result)
    assert len(back) == len(st)
    assert back[0].count_atom_sites() == st[0].count_atom_sites()


def test_write_protein_str_mmcif(path_6h41_cif):
    """write_protein_str(st, "mmcif") returns same mmCIF content as "cif"."""
    st = read_protein_bytes(path_6h41_cif)
    out_cif = write_protein_str(st, "cif")
    out_mmcif = write_protein_str(st, "mmcif")
    assert out_cif == out_mmcif
    back = read_protein_str(out_mmcif)
    _assert_valid_structure(back, "mmcif round-trip")


def test_write_protein_str_case_insensitive(path_6h41_cif):
    """Format "PDB" and "CIF" are accepted and produce same content as lowercase."""
    st = read_protein_bytes(path_6h41_cif)
    assert write_protein_str(st, "PDB") == write_protein_str(st, "pdb")
    assert write_protein_str(st, "CIF") == write_protein_str(st, "cif")


def test_write_protein_str_unsupported_format(path_6h41_cif):
    """Unsupported format raises ValueError with message mentioning format and allowed values."""
    st = read_protein_bytes(path_6h41_cif)
    with pytest.raises(ValueError, match="Unsupported format.*Use 'pdb'"):
        write_protein_str(st, "xyz")
    with pytest.raises(ValueError, match="Unsupported format"):
        write_protein_str(st, "")


def test_write_protein_str_invalid_structure_type():
    """Non-Structure (None or wrong type) raises TypeError."""
    with pytest.raises(TypeError, match="structure must be a gemmi.Structure.*None"):
        write_protein_str(None, "pdb")
    with pytest.raises(TypeError, match="structure must be a gemmi.Structure.*str"):
        write_protein_str("not a structure", "pdb")


# --- gemmi_to_rdkit ---


def test_gemmi_to_rdkit_happy_path(path_6h41_cif):
    """gemmi_to_rdkit(st) returns a non-None RDKit Mol; minimal sanity (atom count >= 1)."""
    path = path_6h41_cif
    assert path.exists(), f"test data missing: {path}"
    st = read_protein_bytes(path)
    mol = gemmi_to_rdkit(st)
    assert mol is not None
    assert isinstance(mol, Chem.Mol)
    assert mol.GetNumAtoms() >= 1


def test_gemmi_to_rdkit_invalid_structure_type():
    """gemmi_to_rdkit(None) and gemmi_to_rdkit('not a structure') raise TypeError."""
    with pytest.raises(TypeError, match="gemmi.Structure"):
        gemmi_to_rdkit(None)
    with pytest.raises(TypeError, match="gemmi.Structure"):
        gemmi_to_rdkit("not a structure")


def test_gemmi_to_rdkit_rdkit_returns_none(path_6h41_cif):
    """When MolFromPDBBlock returns None, gemmi_to_rdkit raises ValueError."""
    st = read_protein_bytes(path_6h41_cif)
    with patch("dutils0.rough.convert.Chem.MolFromPDBBlock", return_value=None):
        with pytest.raises(ValueError, match="RDKit failed to parse PDB block"):
            gemmi_to_rdkit(st)


# --- gemmi_from_rdkit ---


def test_gemmi_from_rdkit_happy_path_round_trip(path_6h41_cif):
    """gemmi_from_rdkit(gemmi_to_rdkit(st)) returns valid Structure; round-trip atom count."""
    st = read_protein_bytes(path_6h41_cif)
    mol = gemmi_to_rdkit(st)
    back = gemmi_from_rdkit(mol)
    _assert_valid_structure(back, "gemmi_from_rdkit round-trip")
    # Approximate equality: same number of atoms in first model
    orig_count = st[0].count_atom_sites()
    back_count = back[0].count_atom_sites()
    assert back_count >= 1
    assert back_count == orig_count, f"round-trip atom count {back_count} vs original {orig_count}"


def test_gemmi_from_rdkit_invalid_mol_type():
    """gemmi_from_rdkit(None) and gemmi_from_rdkit('not a mol') raise TypeError."""
    with pytest.raises(TypeError, match="ROMol|Mol"):
        gemmi_from_rdkit(None)
    with pytest.raises(TypeError, match="ROMol|Mol"):
        gemmi_from_rdkit("not a mol")


def test_gemmi_from_rdkit_mol_to_pdb_block_returns_none(path_6h41_cif):
    """When MolToPDBBlock returns None, gemmi_from_rdkit raises ValueError."""
    st = read_protein_bytes(path_6h41_cif)
    mol = gemmi_to_rdkit(st)
    with patch("dutils0.rough.convert.Chem.MolToPDBBlock", return_value=None):
        with pytest.raises(ValueError, match="MolToPDBBlock returned None or empty"):
            gemmi_from_rdkit(mol)


def test_gemmi_from_rdkit_mol_to_pdb_block_returns_empty(path_6h41_cif):
    """When MolToPDBBlock returns empty string, gemmi_from_rdkit raises ValueError."""
    st = read_protein_bytes(path_6h41_cif)
    mol = gemmi_to_rdkit(st)
    with patch("dutils0.rough.convert.Chem.MolToPDBBlock", return_value=""):
        with pytest.raises(ValueError, match="MolToPDBBlock returned None or empty"):
            gemmi_from_rdkit(mol)


def test_gemmi_from_rdkit_parsing_failure_propagates(path_6h41_cif):
    """When read_protein_str raises (e.g. bad PDB), gemmi_from_rdkit propagates the exception."""
    st = read_protein_bytes(path_6h41_cif)
    mol = gemmi_to_rdkit(st)
    with patch("dutils0.rough.convert.read_protein_str", side_effect=ValueError("parse failed")):
        with pytest.raises(ValueError, match="parse failed"):
            gemmi_from_rdkit(mol)


# --- gemmi_to_pymol ---


def test_gemmi_to_pymol_happy_path(path_6h41_cif):
    """gemmi_to_pymol(st, name, cmd=mock) calls read_pdbstr with PDB string and name; returns cmd."""
    st = read_protein_bytes(path_6h41_cif)
    mock_cmd = MagicMock()
    expected_pdb = write_protein_str(st, "pdb")
    result = gemmi_to_pymol(st, "mymol", cmd=mock_cmd)
    mock_cmd.read_pdbstr.assert_called_once()
    call_args = mock_cmd.read_pdbstr.call_args[0]
    assert call_args[0] == expected_pdb
    assert call_args[1] == "mymol"
    assert result is mock_cmd


def test_gemmi_to_pymol_normalizes_object_name(path_6h41_cif):
    """gemmi_to_pymol strips object_name and returns cmd."""
    st = read_protein_bytes(path_6h41_cif)
    mock_cmd = MagicMock()
    result = gemmi_to_pymol(st, "  myobj  ", cmd=mock_cmd)
    assert mock_cmd.read_pdbstr.call_args[0][1] == "myobj"
    assert result is mock_cmd


def test_gemmi_to_pymol_invalid_structure_type():
    """gemmi_to_pymol(None, name) and gemmi_to_pymol('x', name) raise TypeError."""
    mock_cmd = MagicMock()
    with pytest.raises(TypeError, match="gemmi.Structure"):
        gemmi_to_pymol(None, "mymol", cmd=mock_cmd)
    with pytest.raises(TypeError, match="gemmi.Structure"):
        gemmi_to_pymol("not a structure", "mymol", cmd=mock_cmd)
    mock_cmd.read_pdbstr.assert_not_called()


def test_gemmi_to_pymol_invalid_object_name(path_6h41_cif):
    """gemmi_to_pymol(st, '') and gemmi_to_pymol(st, '   ') raise ValueError."""
    st = read_protein_bytes(path_6h41_cif)
    mock_cmd = MagicMock()
    with pytest.raises(ValueError, match="non-empty string"):
        gemmi_to_pymol(st, "", cmd=mock_cmd)
    with pytest.raises(ValueError, match="non-empty string"):
        gemmi_to_pymol(st, "   ", cmd=mock_cmd)
    mock_cmd.read_pdbstr.assert_not_called()


def test_gemmi_to_pymol_serialization_failure_propagates(path_6h41_cif):
    """When write_protein_str raises, gemmi_to_pymol propagates the exception."""
    st = read_protein_bytes(path_6h41_cif)
    mock_cmd = MagicMock()
    with patch("dutils0.rough.convert.write_protein_str", side_effect=ValueError("serialization failed")):
        with pytest.raises(ValueError, match="serialization failed"):
            gemmi_to_pymol(st, "mymol", cmd=mock_cmd)
    mock_cmd.read_pdbstr.assert_not_called()


def test_gemmi_to_pymol_pymol_load_failure(path_6h41_cif):
    """When read_pdbstr raises, gemmi_to_pymol raises ValueError with PyMOL message."""
    st = read_protein_bytes(path_6h41_cif)
    mock_cmd = MagicMock()
    mock_cmd.read_pdbstr.side_effect = RuntimeError("PyMOL internal error")
    with pytest.raises(ValueError, match="PyMOL failed to load PDB"):
        gemmi_to_pymol(st, "mymol", cmd=mock_cmd)
    assert mock_cmd.read_pdbstr.called


# --- protein_to_rdkit ---


def test_protein_to_rdkit_happy_path(path_6h41_cif):
    """protein_to_rdkit(model, cmd=mock) returns non-None RDKit Mol; temp object deleted."""
    st = read_protein_bytes(path_6h41_cif)
    model = st[0]
    mock_cmd = MagicMock()
    # Echo PDB back from get_pdbstr so RDKit gets valid input (same as write_protein_str output).
    pdb_str = write_protein_str(st, "pdb")
    mock_cmd.get_pdbstr.return_value = pdb_str
    mol = protein_to_rdkit(model, cmd=mock_cmd)
    assert mol is not None
    assert isinstance(mol, Chem.Mol)
    assert mol.GetNumAtoms() >= 1
    mock_cmd.read_pdbstr.assert_called_once()
    assert mock_cmd.read_pdbstr.call_args[0][1] == "_dutils0_protein_to_rdkit"
    mock_cmd.get_pdbstr.assert_called_once_with("_dutils0_protein_to_rdkit")
    mock_cmd.delete.assert_called_once_with("_dutils0_protein_to_rdkit")


def test_protein_to_rdkit_accepts_structure(path_6h41_cif):
    """protein_to_rdkit(structure) uses first model and returns same as protein_to_rdkit(st[0])."""
    st = read_protein_bytes(path_6h41_cif)
    mock_cmd = MagicMock()
    pdb_str = write_protein_str(st, "pdb")
    mock_cmd.get_pdbstr.return_value = pdb_str
    mol_from_structure = protein_to_rdkit(st, cmd=mock_cmd)
    assert mol_from_structure is not None
    assert mol_from_structure.GetNumAtoms() >= 1
    mock_cmd.read_pdbstr.assert_called_once()
    mock_cmd.delete.assert_called_once_with("_dutils0_protein_to_rdkit")


def test_protein_to_rdkit_empty_structure_raises():
    """protein_to_rdkit(Structure with no models) raises ValueError."""
    empty_st = gemmi.Structure()
    mock_cmd = MagicMock()
    with pytest.raises(ValueError, match="Structure has no models"):
        protein_to_rdkit(empty_st, cmd=mock_cmd)
    mock_cmd.read_pdbstr.assert_not_called()


def test_protein_to_rdkit_invalid_model_type(path_6h41_cif):
    """protein_to_rdkit(None) and protein_to_rdkit('x') raise TypeError."""
    mock_cmd = MagicMock()
    with pytest.raises(TypeError, match="model_or_structure must be a gemmi.Model or gemmi.Structure.*NoneType"):
        protein_to_rdkit(None, cmd=mock_cmd)
    with pytest.raises(TypeError, match="model_or_structure must be a gemmi.Model or gemmi.Structure.*str"):
        protein_to_rdkit("not a model", cmd=mock_cmd)
    mock_cmd.read_pdbstr.assert_not_called()


def test_protein_to_rdkit_serialization_failure_propagates(path_6h41_cif):
    """When write_protein_str raises, protein_to_rdkit propagates the exception."""
    st = read_protein_bytes(path_6h41_cif)
    model = st[0]
    mock_cmd = MagicMock()
    with patch("dutils0.rough.convert.write_protein_str", side_effect=ValueError("serialization failed")):
        with pytest.raises(ValueError, match="serialization failed"):
            protein_to_rdkit(model, cmd=mock_cmd)
    mock_cmd.read_pdbstr.assert_not_called()


def test_protein_to_rdkit_rdkit_returns_none(path_6h41_cif):
    """When MolFromPDBBlock returns None, protein_to_rdkit raises ValueError; temp object still deleted."""
    st = read_protein_bytes(path_6h41_cif)
    model = st[0]
    mock_cmd = MagicMock()
    mock_cmd.get_pdbstr.return_value = write_protein_str(st, "pdb")
    with patch("dutils0.rough.convert.Chem.MolFromPDBBlock", return_value=None):
        with pytest.raises(ValueError, match="RDKit failed to parse PDB block"):
            protein_to_rdkit(model, cmd=mock_cmd)
    mock_cmd.delete.assert_called_once_with("_dutils0_protein_to_rdkit")


def test_protein_to_rdkit_pymol_load_failure(path_6h41_cif):
    """When read_pdbstr raises, protein_to_rdkit raises ValueError; temp object delete attempted."""
    st = read_protein_bytes(path_6h41_cif)
    model = st[0]
    mock_cmd = MagicMock()
    mock_cmd.read_pdbstr.side_effect = RuntimeError("PyMOL internal error")
    with pytest.raises(ValueError, match="PyMOL failed to load PDB"):
        protein_to_rdkit(model, cmd=mock_cmd)
    mock_cmd.delete.assert_called_once_with("_dutils0_protein_to_rdkit")


def test_protein_to_rdkit_pymol_get_pdbstr_failure(path_6h41_cif):
    """When get_pdbstr raises, protein_to_rdkit raises ValueError; temp object deleted."""
    st = read_protein_bytes(path_6h41_cif)
    model = st[0]
    mock_cmd = MagicMock()
    mock_cmd.get_pdbstr.side_effect = RuntimeError("write failed")
    with pytest.raises(ValueError, match="PyMOL failed to write PDB string"):
        protein_to_rdkit(model, cmd=mock_cmd)
    mock_cmd.delete.assert_called_once_with("_dutils0_protein_to_rdkit")
