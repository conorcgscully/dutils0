import pytest
from openeye import oechem

from chemutils.molecule.format import get_format, get_format_and_gzip_from_filename


@pytest.mark.parametrize(
    ["name", "format"],
    [
        ("pdb", oechem.OEFormat_PDB),
        ("PDB", oechem.OEFormat_PDB),
        ("Pdb", oechem.OEFormat_PDB),
        ("mol2", oechem.OEFormat_MOL2),
        ("SDF", oechem.OEFormat_SDF),
    ],
)
def test_valid_format(name, format):
    assert get_format(name) == format


@pytest.mark.parametrize("name", ["not_a_format", "PDB2"])
def test_invalid_format(name):
    with pytest.raises(KeyError):
        _ = get_format(name)


@pytest.mark.parametrize(
    ["filename", "expected_format", "expected_gzip"],
    [
        ("file.pdb", oechem.OEFormat_PDB, False),
        ("file.pdb.gz", oechem.OEFormat_PDB, True),
        ("file.sdf", oechem.OEFormat_SDF, False),
        ("file.sdf.gz", oechem.OEFormat_SDF, True),
        ("file.mol2", oechem.OEFormat_MOL2, False),
        ("file.mol2.gz", oechem.OEFormat_MOL2, True),
        ("file.xyz", oechem.OEFormat_XYZ, False),
        ("file.xyz.gz", oechem.OEFormat_XYZ, True),
        ("some.file.pdb", oechem.OEFormat_PDB, False),
        ("some.file.pdb.gz", oechem.OEFormat_PDB, True),
        ("some.file.sdf", oechem.OEFormat_SDF, False),
        ("some.file.sdf.gz", oechem.OEFormat_SDF, True),
        ("some.file.mol2", oechem.OEFormat_MOL2, False),
        ("some.file.mol2.gz", oechem.OEFormat_MOL2, True),
        ("some.file.xyz", oechem.OEFormat_XYZ, False),
        ("some.file.xyz.gz", oechem.OEFormat_XYZ, True),
        ("http://my-website.com/myfile.sdf?query=1&query2=abc", oechem.OEFormat_SDF, False),
    ],
)
def test_get_format_and_gzip_from_filename(filename, expected_format, expected_gzip):
    format, gzip = get_format_and_gzip_from_filename(filename)
    assert format == expected_format
    assert gzip == expected_gzip


@pytest.mark.parametrize(
    "filename",
    [
        "file",
    ],
)
def test_get_format_and_gzip_from_filename_no_extension(filename):
    with pytest.raises(ValueError) as exc_info:
        _, _ = get_format_and_gzip_from_filename(filename)
    assert str(exc_info.value) == f"Filename {filename} does not have an extension"


@pytest.mark.parametrize(
    "filename",
    ["file.gz"],
)
def test_get_format_and_gzip_from_filename_gz_extension(filename):
    with pytest.raises(ValueError) as exc_info:
        _, _ = get_format_and_gzip_from_filename(filename)
    assert str(exc_info.value) == "Filename only has a .gz extension"


@pytest.mark.parametrize(
    ["filename", "extension"],
    [("file.txt", ".txt"), ("file.blah.gz", ".blah"), ("some.file.pdf", ".pdf")],
)
def test_get_format_and_gzip_from_filename_unknown_extension(filename, extension):
    with pytest.raises(ValueError) as exc_info:
        _, _ = get_format_and_gzip_from_filename(filename)
    assert str(exc_info.value) == f"Unrecognized extension {extension}"
