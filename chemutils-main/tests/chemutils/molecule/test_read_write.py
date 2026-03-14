import pytest
from openeye import oechem

from chemutils.molecule import (
    ReadMoleculeError,
    assert_mol_equal,
    oemol_from_smiles,
    read_molecule_bytes,
    read_molecule_str,
    read_molecules_bytes,
    read_molecules_str,
    write_molecule_bytes,
    write_molecule_str,
)


@pytest.mark.parametrize("smiles", ["CN1C=NC2=C1C(=O)N(C(=O)N2C)C"])
@pytest.mark.parametrize(
    ["format_name", "format"],
    [("mol2", oechem.OEFormat_MOL2), ("SDF", oechem.OEFormat_SDF), ("PDB", oechem.OEFormat_PDB)],
)
def test_read_write_str(smiles, format_name, format):
    orig_mol = oemol_from_smiles(smiles)

    file_chemutils = write_molecule_str(orig_mol, format=format_name)

    mol_chemutils = read_molecule_str(file_chemutils, format=format_name)

    oms = oechem.oemolostream()
    oms.SetFormat(format)
    oms.openstring()
    oechem.OEWriteMolecule(oms, orig_mol)

    file_openeye = oms.GetString().decode("UTF-8")

    ims = oechem.oemolistream()
    ims.SetFormat(format)
    ims.openstring(file_openeye.encode("UTF-8"))

    mol_openeye = next(ims.GetOEGraphMols())

    assert file_chemutils == file_openeye

    assert_mol_equal(mol_chemutils, mol_openeye)


@pytest.mark.parametrize("smiles", ["CN1C=NC2=C1C(=O)N(C(=O)N2C)C"])
@pytest.mark.parametrize(
    ["format_name", "format"],
    [("mol2", oechem.OEFormat_MOL2), ("SDF", oechem.OEFormat_SDF), ("PDB", oechem.OEFormat_PDB)],
)
@pytest.mark.parametrize("gzip", [True, False])
def test_read_write_bytes(smiles, format_name, format, gzip):
    orig_mol = oemol_from_smiles(smiles)

    file_chemutils = write_molecule_bytes(orig_mol, format=format_name, gzip=gzip)

    mol_chemutils = read_molecule_bytes(file_chemutils, format=format_name, gzip=gzip)

    oms = oechem.oemolostream()
    oms.SetFormat(format)
    oms.Setgz(gzip)
    oms.openstring()
    oechem.OEWriteMolecule(oms, orig_mol)

    file_openeye = oms.GetString()

    ims = oechem.oemolistream()
    ims.SetFormat(format)
    ims.Setgz(gzip)
    ims.openstring(file_openeye)

    mol_openeye = next(ims.GetOEGraphMols())

    assert file_chemutils == file_openeye

    assert_mol_equal(mol_chemutils, mol_openeye)


@pytest.mark.parametrize("format_name", ["mol2", "SDF", "PDB"])
def test_read_invalid_bytes(format_name):
    file_chemutils = b"I'm a garbage file"
    with pytest.raises(ValueError):
        _ = read_molecule_bytes(file_chemutils, format=format_name)


@pytest.mark.parametrize("format_name", ["mol2", "SDF", "PDB"])
def test_read_invalid_str(format_name):
    file_chemutils = "I'm a garbage file"
    with pytest.raises(ValueError):
        _ = read_molecule_str(file_chemutils, format=format_name)


@pytest.mark.parametrize(
    ["smiles1", "smiles2"], [("CN1C=NC2=C1C(=O)N(C(=O)N2C)C", "Cn1cnc2c1c(=O)[nH]c(=O)n2C")]
)
@pytest.mark.parametrize(
    ["format_name", "format"],
    [("mol2", oechem.OEFormat_MOL2), ("SDF", oechem.OEFormat_SDF), ("PDB", oechem.OEFormat_PDB)],
)
def test_read_write_multiple_str(smiles1, smiles2, format_name, format):
    orig_mol1 = oemol_from_smiles(smiles1)
    orig_mol2 = oemol_from_smiles(smiles2)

    file_chemutils = write_molecule_str(orig_mol1, orig_mol2, format=format_name)

    mol_chemutils1, mol_chemutils2 = list(read_molecules_str(file_chemutils, format=format_name))

    oms = oechem.oemolostream()
    oms.SetFormat(format)
    oms.openstring()
    oechem.OEWriteMolecule(oms, orig_mol1)
    oechem.OEWriteMolecule(oms, orig_mol2)

    file_openeye = oms.GetString().decode("UTF-8")

    ims = oechem.oemolistream()
    ims.SetFormat(format)
    ims.openstring(file_openeye.encode("UTF-8"))

    mol_generator = iter(ims.GetOEGraphMols())
    mol_openeye1 = next(mol_generator).CreateCopy()
    mol_openeye2 = next(mol_generator).CreateCopy()

    assert file_chemutils == file_openeye

    assert_mol_equal(mol_chemutils1, mol_openeye1)
    assert_mol_equal(mol_chemutils2, mol_openeye2)
    with pytest.raises(AssertionError):
        assert_mol_equal(mol_chemutils1, mol_chemutils2)


@pytest.mark.parametrize(
    ["smiles1", "smiles2"], [("CN1C=NC2=C1C(=O)N(C(=O)N2C)C", "Cn1cnc2c1c(=O)[nH]c(=O)n2C")]
)
@pytest.mark.parametrize(
    ["format_name", "format"],
    [("mol2", oechem.OEFormat_MOL2), ("SDF", oechem.OEFormat_SDF), ("PDB", oechem.OEFormat_PDB)],
)
@pytest.mark.parametrize("gzip", [True, False])
def test_read_write_multiple_bytes(smiles1, smiles2, format_name, format, gzip):
    orig_mol1 = oemol_from_smiles(smiles1)
    orig_mol2 = oemol_from_smiles(smiles2)

    file_chemutils = write_molecule_bytes(orig_mol1, orig_mol2, format=format_name, gzip=gzip)

    mol_chemutils1, mol_chemutils2 = list(
        read_molecules_bytes(file_chemutils, format=format_name, gzip=gzip)
    )

    oms = oechem.oemolostream()
    oms.SetFormat(format)
    oms.Setgz(gzip)
    oms.openstring()
    oechem.OEWriteMolecule(oms, orig_mol1)
    oechem.OEWriteMolecule(oms, orig_mol2)

    file_openeye = oms.GetString()

    ims = oechem.oemolistream()
    ims.SetFormat(format)
    ims.Setgz(gzip)
    ims.openstring(file_openeye)

    mol_generator = iter(ims.GetOEGraphMols())
    mol_openeye1 = next(mol_generator).CreateCopy()
    mol_openeye2 = next(mol_generator).CreateCopy()

    assert file_chemutils == file_openeye

    assert_mol_equal(mol_chemutils1, mol_openeye1)
    assert_mol_equal(mol_chemutils2, mol_openeye2)
    with pytest.raises(AssertionError):
        assert_mol_equal(mol_chemutils1, mol_chemutils2)


@pytest.mark.parametrize("format_name", ["mol2", "SDF", "PDB"])
def test_read_invalid_multiple_bytes(format_name):
    file_chemutils = b"I'm a garbage file"
    with pytest.raises(ReadMoleculeError):
        for _ in read_molecules_bytes(file_chemutils, format=format_name):
            pass


@pytest.mark.parametrize("format_name", ["mol2", "SDF", "PDB"])
def test_read_invalid_multiple_str(format_name):
    file_chemutils = "I'm a garbage file"
    with pytest.raises(ReadMoleculeError):
        for _ in read_molecules_str(file_chemutils, format=format_name):
            pass


def test_read_molecule_invalid_smiles():
    file = b"c1ccccc2\n"
    with pytest.raises(ReadMoleculeError) as e:
        _ = read_molecule_bytes(file, format="SMI")
    assert str(e.value) == "Failed to read molecule."


def test_read_molecules_invalid_smiles():
    file = b"c1ccccc2\n"
    with pytest.raises(ReadMoleculeError) as e:
        for _ in read_molecules_bytes(file, format="SMI"):
            pass
    assert str(e.value) == "Failed to read molecule at index 0."


def test_read_molecule_invalid_smiles_second_item():
    file = b"c1ccccc1\nc1ccccc2\n"
    # no error, because first SMILES is valid
    _ = read_molecule_bytes(file, format="SMI")


def test_read_molecules_invalid_smiles_second_item():
    file = b"c1ccccc1\nc1ccccc2\n"
    with pytest.raises(ReadMoleculeError) as e:
        for _ in read_molecules_bytes(file, format="SMI"):
            pass
    assert str(e.value) == "Failed to read molecule at index 1."


def test_read_molecule_invalid_smiles_first_item():
    file = b"c1ccccc2\nc1ccccc1\n"
    with pytest.raises(ReadMoleculeError) as e:
        _ = read_molecule_bytes(file, format="SMI")
    assert str(e.value) == "Failed to read molecule."


def test_read_molecules_invalid_smiles_first_item():
    file = b"c1ccccc2\nc1ccccc1\n"
    with pytest.raises(ReadMoleculeError) as e:
        for _ in read_molecules_bytes(file, format="SMI"):
            pass
    assert str(e.value) == "Failed to read molecule at index 0."


def test_read_molecule_invalid_smiles_all_items():
    file = b"c1ccccc3\nc1ccccc2\n"
    with pytest.raises(ReadMoleculeError) as e:
        _ = read_molecule_bytes(file, format="SMI")
    assert str(e.value) == "Failed to read molecule."


def test_read_molecules_invalid_smiles_all_items():
    file = b"c1ccccc3\nc1ccccc2\n"
    with pytest.raises(ReadMoleculeError) as e:
        for _ in read_molecules_bytes(file, format="SMI"):
            pass
    assert str(e.value) == "Failed to read molecule at index 0."
