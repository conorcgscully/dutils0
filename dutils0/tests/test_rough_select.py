"""Tests for dutils0.rough.pymol_select."""

from unittest.mock import MagicMock

import pytest

from dutils0.rough import pymol_select


def test_pymol_select_happy_path():
    """pymol_select(cmd, selection, object_name) calls count_atoms then create; returns new object name."""
    mock_cmd = MagicMock()
    mock_cmd.count_atoms.return_value = 1
    mock_cmd.create.return_value = None

    result = pymol_select(mock_cmd, "chain A", "newobj")

    mock_cmd.count_atoms.assert_called_once_with("chain A")
    mock_cmd.create.assert_called_once_with("newobj", "chain A")
    assert result == "newobj"


def test_pymol_select_normalizes_name_and_selection():
    """pymol_select strips selection and object_name; create with stripped values; returns new object name."""
    mock_cmd = MagicMock()
    mock_cmd.count_atoms.return_value = 1
    mock_cmd.create.return_value = None

    result = pymol_select(mock_cmd, "  chain A  ", "  myobj  ")

    mock_cmd.count_atoms.assert_called_once_with("chain A")
    mock_cmd.create.assert_called_once_with("myobj", "chain A")
    assert result == "myobj"


def test_pymol_select_invalid_cmd():
    """pymol_select with cmd=None or missing create/count_atoms raises TypeError."""
    with pytest.raises(TypeError, match="cmd must be a PyMOL command object"):
        pymol_select(None, "chain A", "newobj")

    no_create = MagicMock(spec=["count_atoms"])
    no_create.count_atoms.return_value = 1
    with pytest.raises(TypeError, match="cmd must be a PyMOL command object"):
        pymol_select(no_create, "chain A", "newobj")

    no_count = MagicMock(spec=["create"])
    no_count.create.return_value = None
    with pytest.raises(TypeError, match="cmd must be a PyMOL command object"):
        pymol_select(no_count, "chain A", "newobj")


def test_pymol_select_invalid_selection():
    """pymol_select with non-string or empty selection raises ValueError; create not called."""
    mock_cmd = MagicMock()

    with pytest.raises(ValueError, match="selection must be a non-empty string"):
        pymol_select(mock_cmd, 123, "newobj")
    mock_cmd.create.assert_not_called()
    mock_cmd.count_atoms.assert_not_called()

    with pytest.raises(ValueError, match="selection must be a non-empty string"):
        pymol_select(mock_cmd, "", "newobj")
    with pytest.raises(ValueError, match="selection must be a non-empty string"):
        pymol_select(mock_cmd, "   ", "newobj")
    mock_cmd.create.assert_not_called()


def test_pymol_select_invalid_object_name():
    """pymol_select with non-string or empty object_name raises ValueError; create not called."""
    mock_cmd = MagicMock()
    mock_cmd.count_atoms.return_value = 1

    with pytest.raises(ValueError, match="object_name must be a non-empty string"):
        pymol_select(mock_cmd, "chain A", 99)
    mock_cmd.create.assert_not_called()

    with pytest.raises(ValueError, match="object_name must be a non-empty string"):
        pymol_select(mock_cmd, "chain A", "")
    with pytest.raises(ValueError, match="object_name must be a non-empty string"):
        pymol_select(mock_cmd, "chain A", "   ")
    mock_cmd.create.assert_not_called()


def test_pymol_select_empty_selection():
    """When count_atoms returns 0, raise ValueError and do not call create."""
    mock_cmd = MagicMock()
    mock_cmd.count_atoms.return_value = 0

    with pytest.raises(ValueError, match="selection matches no atoms"):
        pymol_select(mock_cmd, "chain A", "newobj")

    mock_cmd.count_atoms.assert_called_once_with("chain A")
    mock_cmd.create.assert_not_called()


def test_pymol_select_creation_failure():
    """When create raises, raise ValueError with chained cause."""
    mock_cmd = MagicMock()
    mock_cmd.count_atoms.return_value = 1
    mock_cmd.create.side_effect = RuntimeError("bad")

    with pytest.raises(ValueError, match="PyMOL selection failed") as exc_info:
        pymol_select(mock_cmd, "chain A", "newobj")

    assert "bad" in str(exc_info.value)
    assert exc_info.value.__cause__ is mock_cmd.create.side_effect
