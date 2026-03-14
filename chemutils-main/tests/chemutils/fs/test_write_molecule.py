import tempfile

import pytest

from chemutils import fs
from chemutils.molecule import oemol_from_smiles, write_molecule_bytes

from .test_read_molecule import FILE_FORMATS, SMILES


@pytest.mark.parametrize("smiles", SMILES)
@pytest.mark.parametrize(
    ["filename", "format", "gzip"],
    FILE_FORMATS,
)
@pytest.mark.parametrize("prefix", ["", "file://"])
# Can't use fake filesystem when dealing with OpenEye input/output
def test_write_molecule_local(smiles, prefix, filename, format, gzip):
    with tempfile.TemporaryDirectory() as tempdir:
        url = f"{prefix}{tempdir}/{filename}"

        oemol = oemol_from_smiles(smiles)

        fs.write_molecule(url, mol=oemol)

        oemol_bytes = write_molecule_bytes(oemol, format=format, gzip=gzip)
        assert fs.read_bytes(url) == oemol_bytes


@pytest.mark.parametrize("smiles", SMILES)
@pytest.mark.parametrize(
    ["filename", "format", "gzip"],
    FILE_FORMATS,
)
def test_write_molecule_s3(s3, s3_bucket, smiles, filename, format, gzip):
    url = f"s3://example-bucket/{filename}"

    oemol = oemol_from_smiles(smiles)

    fs.write_molecule(url, mol=oemol)

    oemol_bytes = write_molecule_bytes(oemol, format=format, gzip=gzip)
    assert fs.read_bytes(url) == oemol_bytes


@pytest.mark.parametrize(
    ["filename", "format", "gzip"],
    FILE_FORMATS,
)
@pytest.mark.parametrize("prefix", ["", "file:///"])
# Can't use fake filesystem when dealing with OpenEye input/output
def test_write_molecules_local(prefix, filename, format, gzip):
    with tempfile.TemporaryDirectory() as tempdir:
        url = f"{prefix}{tempdir}/{filename}"

        oemols = [oemol_from_smiles(smiles) for smiles in SMILES]

        fs.write_molecules(url, mols=oemols)

        oemol_bytes = write_molecule_bytes(*oemols, format=format, gzip=gzip)
        assert fs.read_bytes(url) == oemol_bytes


@pytest.mark.parametrize(
    ["filename", "format", "gzip"],
    FILE_FORMATS,
)
def test_write_molecules_s3(s3, s3_bucket, filename, format, gzip):
    url = f"s3://example-bucket/{filename}"

    oemols = [oemol_from_smiles(smiles) for smiles in SMILES]

    fs.write_molecules(url, mols=oemols)

    oemol_bytes = write_molecule_bytes(*oemols, format=format, gzip=gzip)
    assert fs.read_bytes(url) == oemol_bytes


@pytest.mark.parametrize("smiles", SMILES)
@pytest.mark.parametrize("filename", ["test.test", "test.test.gz"])
def test_write_molecule_no_inferred_format(filename, smiles):
    with tempfile.TemporaryDirectory() as tempdir:
        path = f"{tempdir}/{filename}"
        oemol = oemol_from_smiles(smiles)

        with open(path, "wb") as f:
            f.write(b"test")

            with pytest.raises(ValueError, match="Unrecognized extension .test"):
                fs.write_molecule(path, mol=oemol)


@pytest.mark.parametrize("filename", ["test.test", "test.test.gz"])
def test_write_molecules_no_inferred_format(filename):
    with tempfile.TemporaryDirectory() as tempdir:
        path = f"{tempdir}/{filename}"
        oemols = [oemol_from_smiles(smiles) for smiles in SMILES]

        with open(path, "wb") as f:
            f.write(b"test")

            with pytest.raises(ValueError, match="Unrecognized extension .test"):
                fs.write_molecules(path, mols=oemols)
