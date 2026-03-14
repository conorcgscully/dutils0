import tempfile

import pytest
from openeye import oechem

from chemutils import fs
from chemutils.molecule import (
    assert_mol_equal,
    oemol_from_smiles,
    read_molecule_bytes,
    read_molecules_bytes,
    write_molecule_bytes,
)

SMILES = [
    "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
    "Cn1cnc2c1c(=O)[nH]c(=O)n2C",
    "CCCC1=NN(C2=C1N=C(NC2=O)C3=C(C=CC(=C3)S(=O)(=O)N4CCN(CC4)C)OCC)C",
]

FILE_FORMATS = [
    ("file.mol2", oechem.OEFormat_MOL2, False),
    ("file.sdf", oechem.OEFormat_SDF, False),
    ("file.pdb", oechem.OEFormat_PDB, False),
    ("file.mol2.gz", oechem.OEFormat_MOL2, True),
    ("file.sdf.gz", oechem.OEFormat_SDF, True),
    ("file.pdb.gz", oechem.OEFormat_PDB, True),
]


@pytest.mark.parametrize("smiles", SMILES)
@pytest.mark.parametrize(
    ["filename", "format", "gzip"],
    FILE_FORMATS,
)
@pytest.mark.parametrize("prefix", ["", "file://"])
# Can't use fake filesystem when dealing with OpenEye input/output
def test_read_molecule_local(smiles, prefix, filename, format, gzip):
    with tempfile.TemporaryDirectory() as tempdir:
        url = f"{tempdir}/{filename}"

        oemol = oemol_from_smiles(smiles)

        oemol_bytes = write_molecule_bytes(oemol, format=format, gzip=gzip)

        with open(f"{tempdir}/{filename}", "wb") as f:
            f.write(oemol_bytes)

        mol_from_bytes = read_molecule_bytes(oemol_bytes, format=format, gzip=gzip)

        mol_from_file = fs.read_molecule(f"{prefix}{url}")

        assert_mol_equal(mol_from_bytes, mol_from_file)


@pytest.mark.parametrize("smiles", SMILES)
@pytest.mark.parametrize(
    ["filename", "format", "gzip"],
    FILE_FORMATS,
)
def test_read_molecule_s3(s3, s3_bucket, smiles, filename, format, gzip):
    url = f"s3://example-bucket/{filename}"

    oemol = oemol_from_smiles(smiles)

    oemol_bytes = write_molecule_bytes(oemol, format=format, gzip=gzip)

    s3_bucket.put_object(Key=filename, Body=oemol_bytes)

    mol_from_bytes = read_molecule_bytes(oemol_bytes, format=format, gzip=gzip)

    mol_from_file = fs.read_molecule(url)

    assert_mol_equal(mol_from_bytes, mol_from_file)


@pytest.mark.parametrize("smiles", SMILES)
@pytest.mark.parametrize(
    ["filename", "format", "gzip"],
    FILE_FORMATS,
)
def test_read_molecule_web(mocked_web, smiles, filename, format, gzip):
    url = f"http://website.com/{filename}"

    oemol = oemol_from_smiles(smiles)

    oemol_bytes = write_molecule_bytes(oemol, format=format, gzip=gzip)

    mocked_web.get(
        url,
        body=oemol_bytes,
        status=200,
        content_type="text/plain",
    )

    mol_from_bytes = read_molecule_bytes(oemol_bytes, format=format, gzip=gzip)

    mol_from_file = fs.read_molecule(url)

    assert_mol_equal(mol_from_bytes, mol_from_file)


@pytest.mark.parametrize(
    ["filename", "format", "gzip"],
    FILE_FORMATS,
)
@pytest.mark.parametrize("prefix", ["", "file://"])
# Can't use fake filesystem when dealing with OpenEye input/output
def test_read_molecules_local(prefix, filename, format, gzip):
    with tempfile.TemporaryDirectory() as tempdir:
        url = f"{tempdir}/{filename}"

        oemols = [oemol_from_smiles(smiles) for smiles in SMILES]

        oemol_bytes = write_molecule_bytes(*oemols, format=format, gzip=gzip)

        with open(f"{tempdir}/{filename}", "wb") as f:
            f.write(oemol_bytes)

        mols_from_bytes = list(read_molecules_bytes(oemol_bytes, format=format, gzip=gzip))

        mols_from_file = list(fs.read_molecules(f"{prefix}{url}"))

        assert len(mols_from_bytes) == len(mols_from_file)
        for mol_from_bytes, mol_from_file in zip(mols_from_bytes, mols_from_file, strict=True):
            assert_mol_equal(mol_from_bytes, mol_from_file)


@pytest.mark.parametrize(
    ["filename", "format", "gzip"],
    FILE_FORMATS,
)
def test_read_molecules_s3(s3, s3_bucket, filename, format, gzip):
    url = f"s3://example-bucket/{filename}"

    oemols = [oemol_from_smiles(smiles) for smiles in SMILES]

    oemol_bytes = write_molecule_bytes(*oemols, format=format, gzip=gzip)

    s3_bucket.put_object(Key=filename, Body=oemol_bytes)

    mols_from_bytes = list(read_molecules_bytes(oemol_bytes, format=format, gzip=gzip))

    mols_from_file = list(fs.read_molecules(url))

    assert len(mols_from_bytes) == len(mols_from_file)
    for mol_from_bytes, mol_from_file in zip(mols_from_bytes, mols_from_file, strict=True):
        assert_mol_equal(mol_from_bytes, mol_from_file)


@pytest.mark.parametrize(
    ["filename", "format", "gzip"],
    FILE_FORMATS,
)
def test_read_molecules_web(mocked_web, filename, format, gzip):
    url = f"http://website.com/{filename}"

    oemols = [oemol_from_smiles(smiles) for smiles in SMILES]

    oemol_bytes = write_molecule_bytes(*oemols, format=format, gzip=gzip)

    mocked_web.get(
        url,
        body=oemol_bytes,
        status=200,
        content_type="text/plain",
    )

    mol_from_bytes = read_molecule_bytes(oemol_bytes, format=format, gzip=gzip)

    mol_from_file = fs.read_molecule(url)

    assert_mol_equal(mol_from_bytes, mol_from_file)


@pytest.mark.parametrize("smiles", SMILES)
@pytest.mark.parametrize(
    ["filename", "format", "gzip"],
    FILE_FORMATS,
)
@pytest.mark.parametrize("prefix", ["", "file://"])
# Can't use fake filesystem when dealing with OpenEye input/output
def test_read_molecule_provided_format(smiles, prefix, filename, format, gzip):
    with tempfile.TemporaryDirectory() as tempdir:
        url = f"{tempdir}/{filename}"

        oemol = oemol_from_smiles(smiles)

        oemol_bytes = write_molecule_bytes(oemol, format=format)

        with open(f"{tempdir}/{filename}", "wb") as f:
            f.write(oemol_bytes)

        mol_from_bytes = read_molecule_bytes(oemol_bytes, format=format)

        mol_from_file = fs.read_molecule(f"{prefix}{url}")

        assert_mol_equal(mol_from_bytes, mol_from_file)


@pytest.mark.parametrize("filename", ["test.test", "test.test.gz"])
def test_read_molecule_no_inferred_format(filename):
    with tempfile.TemporaryDirectory() as tempdir:
        path = f"{tempdir}/{filename}"

        with open(path, "wb") as f:
            f.write(b"test")

            with pytest.raises(ValueError, match="Unrecognized extension .test"):
                fs.read_molecule(path, format=None)


@pytest.mark.parametrize("filename", ["test.test", "test.test.gz"])
def test_read_molecules_no_inferred_format(filename):
    with tempfile.TemporaryDirectory() as tempdir:
        path = f"{tempdir}/{filename}"

        with open(path, "wb") as f:
            f.write(b"test")

            with pytest.raises(ValueError, match="Unrecognized extension .test"):
                generator = fs.read_molecules(path, format=None)
                next(generator)
