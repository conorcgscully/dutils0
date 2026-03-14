import pytest
from openeye import oechem

from chemutils.molecule import create_atom, create_residue, modify_residue

CASES = [
    ("residue_name", oechem.OEResidue.GetName, "ALA", "GLY"),
    ("occupancy", oechem.OEResidue.GetOccupancy, 0.5, 1.0),
    ("b_factor", oechem.OEResidue.GetBFactor, 0.5, 1.0),
    ("residue_number", oechem.OEResidue.GetResidueNumber, 1, 2),
    ("sequence_number", oechem.OEResidue.GetSequenceID, 1, 2),
    ("serial_number", oechem.OEResidue.GetSerialNumber, 1, 2),
    ("alternate_location", oechem.OEResidue.GetAlternateLocation, "A", "B"),
    ("chain_id", oechem.OEResidue.GetExtChainID, "A", "B"),
    ("asymmetric_unit_id", oechem.OEResidue.GetSubChainID, "A", "B"),
    ("insertion_code", oechem.OEResidue.GetInsertCode, "A", "B"),
    ("is_hetatm", oechem.OEResidue.IsHetAtom, False, True),
]


@pytest.mark.parametrize(["name", "func", "arg1", "arg2"], CASES)
def test_modify_residue_oeresidue(name, func, arg1, arg2):
    residue = create_residue(**{name: arg1})
    assert func(residue) == arg1
    modify_residue(residue, **{name: arg2})
    assert func(residue) == arg2


@pytest.mark.parametrize(["name", "func", "arg1", "arg2"], CASES)
def test_modify_residue_oeatom(name, func, arg1, arg2):
    mol = oechem.OEGraphMol()
    atom = create_atom(mol, atomic_number=6, **{name: arg1})
    assert func(oechem.OEAtomGetResidue(atom)) == arg1
    modify_residue(atom, **{name: arg2})
    assert func(oechem.OEAtomGetResidue(atom)) == arg2


@pytest.mark.parametrize(["name", "func", "arg1", "arg2"], CASES)
def test_modify_residue_oemol(name, func, arg1, arg2):
    mol = oechem.OEGraphMol()
    create_atom(mol, atomic_number=6, **{name: arg1})
    create_atom(mol, atomic_number=6, **{name: arg1})
    for atom in mol.GetAtoms():
        assert func(oechem.OEAtomGetResidue(atom)) == arg1
    modify_residue(mol, **{name: arg2})
    for atom in mol.GetAtoms():
        assert func(oechem.OEAtomGetResidue(atom)) == arg2


@pytest.mark.parametrize(["name", "func", "arg1", "arg2"], CASES)
def test_modify_residue_oeatombondset(name, func, arg1, arg2):
    mol = oechem.OEGraphMol()
    atom1 = create_atom(mol, atomic_number=6, **{name: arg1})
    atom2 = create_atom(mol, atomic_number=6, **{name: arg1})
    atombondset = oechem.OEAtomBondSet()
    atombondset.AddAtom(atom1)
    assert func(oechem.OEAtomGetResidue(atom1)) == arg1
    assert func(oechem.OEAtomGetResidue(atom2)) == arg1
    modify_residue(atombondset, **{name: arg2})
    assert func(oechem.OEAtomGetResidue(atom1)) == arg2
    assert func(oechem.OEAtomGetResidue(atom2)) == arg1
