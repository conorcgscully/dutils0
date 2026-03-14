import pytest
from openeye import oechem

from chemutils.molecule.tags import (
    append_sd_tag,
    delete_sd_tags,
    get_sd_tag,
    get_sd_tags,
    has_sd_tag,
    set_sd_tag,
)


@pytest.fixture
def mol():
    return oechem.OEGraphMol()


def test_get_sd_tags(mol):
    assert get_sd_tags(mol) == {}
    with pytest.raises(LookupError):
        _ = get_sd_tags(mol, name="abc")

    append_sd_tag(mol, name="abc", value="123")

    assert get_sd_tags(mol) == {"abc": "123"}
    assert get_sd_tags(mol, name="abc") == ["123"]

    append_sd_tag(mol, name="def", value="456")

    assert get_sd_tags(mol) == {"abc": "123", "def": "456"}
    assert get_sd_tags(mol, name="abc") == ["123"]
    assert get_sd_tags(mol, name="def") == ["456"]

    append_sd_tag(mol, name="abc", value="789")

    assert get_sd_tags(mol) == {"abc": ["123", "789"], "def": "456"}
    assert get_sd_tags(mol, name="abc") == ["123", "789"]
    assert get_sd_tags(mol, name="def") == ["456"]


def test_has_sd_tag(mol):
    assert not has_sd_tag(mol, name="abc")
    append_sd_tag(mol, name="abc", value="123")
    assert has_sd_tag(mol, name="abc")


def test_get_sd_tag(mol):
    with pytest.raises(LookupError):
        _ = get_sd_tag(mol, name="abc")

    append_sd_tag(mol, name="abc", value="123")

    assert get_sd_tag(mol, name="abc") == "123"

    append_sd_tag(mol, name="abc", value="456")

    with pytest.raises(LookupError):
        _ = get_sd_tag(mol, name="abc")


def test_set_sd_tag(mol):
    set_sd_tag(mol, name="abc", value="123")
    assert get_sd_tags(mol, name="abc") == ["123"]

    set_sd_tag(mol, name="abc", value="456")
    assert get_sd_tags(mol, name="abc") == ["456"]

    append_sd_tag(mol, name="abc", value="123")
    assert get_sd_tags(mol, name="abc") == ["456", "123"]

    with pytest.raises(ValueError):
        set_sd_tag(mol, name="abc", value="789")


def test_append_sd_tag(mol):
    append_sd_tag(mol, name="abc", value="123")
    assert get_sd_tags(mol, name="abc") == ["123"]

    append_sd_tag(mol, name="abc", value="456")
    assert get_sd_tags(mol, name="abc") == ["123", "456"]

    append_sd_tag(mol, name="abc", value="789")
    assert get_sd_tags(mol, name="abc") == ["123", "456", "789"]


def test_delete_sd_tags(mol):
    delete_sd_tags(mol)  # Check we can delete on a molecule with no tags

    append_sd_tag(mol, name="abc", value="123")
    append_sd_tag(mol, name="abc", value="456")
    append_sd_tag(mol, name="def", value="789")
    append_sd_tag(mol, name="ghi", value="XYZ")
    append_sd_tag(mol, name="jkl", value="ABC")

    assert get_sd_tags(mol) == {"abc": ["123", "456"], "def": "789", "ghi": "XYZ", "jkl": "ABC"}

    delete_sd_tags(mol, name="def")
    assert get_sd_tags(mol) == {"abc": ["123", "456"], "ghi": "XYZ", "jkl": "ABC"}

    delete_sd_tags(mol, name="abc")
    assert get_sd_tags(mol) == {"ghi": "XYZ", "jkl": "ABC"}

    delete_sd_tags(mol)
    assert get_sd_tags(mol) == {}
