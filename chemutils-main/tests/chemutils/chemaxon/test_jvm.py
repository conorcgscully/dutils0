from chemutils.chemaxon.jvm import isjavainstance
from chemutils.chemaxon.jvm.chemaxon.struc import MolAtom, Molecule
from chemutils.chemaxon.molecule import cxmol_from_smiles


def test_isjavainstance():
    cxmol = cxmol_from_smiles("CCCC(=O)O")

    assert isjavainstance(cxmol, Molecule)
    assert not isjavainstance(cxmol, MolAtom)

    cxatom = cxmol.getAtomArray()[0]

    assert not isjavainstance(cxatom, Molecule)
    assert isjavainstance(cxatom, MolAtom)

    # Check Python objects
    assert not isjavainstance("abc", Molecule)
    assert not isjavainstance(1, Molecule)
