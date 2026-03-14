from openeye import oechem

from chemutils.molecule import oemol_from_smiles
from chemutils.molecule.merge import merge_molecules


def test_merge_molecules():
    oemol1 = oemol_from_smiles("CCC")
    oemol2 = oemol_from_smiles("O")
    oemol3 = oemol_from_smiles("CN")

    assert oemol1.NumAtoms() == 3
    assert oemol1.NumBonds() == 2

    assert oemol2.NumAtoms() == 1
    assert oemol2.NumBonds() == 0

    assert oemol3.NumAtoms() == 2
    assert oemol3.NumBonds() == 1

    oemol = merge_molecules(oemol1, oemol2, oemol3)

    assert oemol.NumAtoms() == 6
    assert oemol.NumBonds() == 3


def test_merge_molecules_map():
    oemol1 = oemol_from_smiles("CCC")
    oemol2 = oemol_from_smiles("O")
    oemol3 = oemol_from_smiles("CN")

    atom_map = {}
    bond_map = {}

    oemol = merge_molecules(oemol1, oemol2, oemol3, atom_mapping=atom_map, bond_mapping=bond_map)

    assert len(atom_map) == 6
    assert len(bond_map) == 3

    atom1, atom2, atom3 = oemol1.GetAtoms()
    bond1, bond2 = oemol1.GetBonds()
    (atom4,) = oemol2.GetAtoms()
    atom5, atom6 = oemol3.GetAtoms()
    (bond3,) = oemol3.GetBonds()

    assert atom_map[atom1] == oemol.GetAtom(oechem.OEHasAtomIdx(0))
    assert atom_map[atom2] == oemol.GetAtom(oechem.OEHasAtomIdx(1))
    assert atom_map[atom3] == oemol.GetAtom(oechem.OEHasAtomIdx(2))
    assert atom_map[atom4] == oemol.GetAtom(oechem.OEHasAtomIdx(3))
    assert atom_map[atom5] == oemol.GetAtom(oechem.OEHasAtomIdx(4))
    assert atom_map[atom6] == oemol.GetAtom(oechem.OEHasAtomIdx(5))

    assert bond_map[bond1] == oemol.GetBond(oechem.OEHasBondIdx(0))
    assert bond_map[bond2] == oemol.GetBond(oechem.OEHasBondIdx(1))
    assert bond_map[bond3] == oemol.GetBond(oechem.OEHasBondIdx(2))
