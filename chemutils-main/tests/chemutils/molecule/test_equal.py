import pytest
from openeye import oechem

from chemutils.molecule import create_atom, oemol_from_smiles
from chemutils.molecule.equal import (
    AllAtomProperties,
    AtomProperty,
    assert_atom_equal,
    assert_bond_equal,
    assert_mol_equal,
    atom_equal,
    bond_equal,
    mol_equal,
)


@pytest.mark.parametrize(
    ["smiles1", "smiles2"], [("CCC", "C(C)C"), ("c1cocc1", "c1ccoc1"), ("CCC(C)CCC", "CCCC(CC)C")]
)
def test_equal_smiles_canonical(smiles1, smiles2):
    mol1 = oechem.OEGraphMol()
    oechem.OESmilesToMol(mol1, smiles1)
    mol2 = oechem.OEGraphMol()
    oechem.OESmilesToMol(mol2, smiles2)

    with pytest.raises(AssertionError):
        assert_mol_equal(mol1, mol2)
    assert not mol_equal(mol1, mol2)

    # Canonically order molecules to allow comparison
    oechem.OECanonicalOrderAtoms(mol1)
    oechem.OECanonicalOrderBonds(mol1)
    mol1.Sweep()
    for bond in mol1.GetBonds():
        if bond.GetBgnIdx() > bond.GetEndIdx():
            bond.SwapEnds()

    oechem.OECanonicalOrderAtoms(mol2)
    oechem.OECanonicalOrderBonds(mol2)
    mol2.Sweep()
    for bond in mol2.GetBonds():
        if bond.GetBgnIdx() > bond.GetEndIdx():
            bond.SwapEnds()

    assert_mol_equal(mol1, mol2)
    assert mol_equal(mol1, mol2)


def test_new_molecule_equal():
    mol1 = oechem.OEGraphMol()
    mol1_c = mol1.NewAtom(6)
    mol1_o = mol1.NewAtom(8)
    mol1_bond = mol1.NewBond(mol1_c, mol1_o)

    mol2 = oechem.OEGraphMol()
    mol2_c = mol2.NewAtom(6)
    mol2_o = mol2.NewAtom(8)
    mol2_bond = mol2.NewBond(mol2_c, mol2_o)

    assert_atom_equal(mol1_c, mol2_c)
    assert atom_equal(mol1_c, mol2_c)
    assert_atom_equal(mol1_o, mol2_o)
    assert atom_equal(mol1_o, mol2_o)
    assert_bond_equal(mol1_bond, mol2_bond)
    assert bond_equal(mol1_bond, mol2_bond)

    assert_mol_equal(mol1, mol2)
    assert mol_equal(mol1, mol2)


def test_atom_unequal_index():
    mol1 = oechem.OEGraphMol()
    mol1_c = mol1.NewAtom(6)
    mol1_o = mol1.NewAtom(8)
    mol1.NewBond(mol1_c, mol1_o)

    mol2 = oechem.OEGraphMol()
    mol2_o = mol2.NewAtom(8)
    mol2_c = mol2.NewAtom(6)
    mol2.NewBond(mol2_c, mol2_o)

    with pytest.raises(AssertionError, match="Atom index unequal: 0 != 1"):
        assert_atom_equal(mol1_c, mol2_c)
    assert not atom_equal(mol1_c, mol2_c)

    with pytest.raises(AssertionError, match="Atoms with index 0 are not equal") as exc:
        assert_mol_equal(mol1, mol2)
    assert str(exc.value.__cause__) == "Atomic number unequal: 6 != 8"
    assert not mol_equal(mol1, mol2)


def test_atom_unequal_atomic_number():
    mol1 = oechem.OEGraphMol()
    mol1_c = mol1.NewAtom(6)

    mol2 = oechem.OEGraphMol()
    mol2_o = mol2.NewAtom(8)

    with pytest.raises(AssertionError, match="Atomic number unequal: 6 != 8"):
        assert_atom_equal(mol1_c, mol2_o)
    assert not atom_equal(mol1_c, mol2_o)

    with pytest.raises(AssertionError, match="Atoms with index 0 are not equal") as exc:
        assert_mol_equal(mol1, mol2)
    assert str(exc.value.__cause__) == "Atomic number unequal: 6 != 8"
    assert not mol_equal(mol1, mol2)


def test_atom_unequal_implicit_h():
    mol1 = oechem.OEGraphMol()
    mol1_c = mol1.NewAtom(6)
    mol1_h = mol1.NewAtom(1)
    mol1.NewBond(mol1_c, mol1_h)
    mol1_c.SetImplicitHCount(1)

    mol2 = oechem.OEGraphMol()
    mol2_c = mol2.NewAtom(6)
    mol2_c.SetImplicitHCount(2)

    with pytest.raises(AssertionError, match="Atom implicit H count unequal: 1 != 2"):
        assert_atom_equal(mol1_c, mol2_c)
    assert not atom_equal(mol1_c, mol2_c)

    with pytest.raises(AssertionError, match="Molecule atom count unequal: 2 != 1"):
        assert_mol_equal(mol1, mol2)
    assert not mol_equal(mol1, mol2)


def test_atom_unequal_formal_charge():
    mol1 = oechem.OEGraphMol()
    mol1_c = mol1.NewAtom(6)
    mol1_c.SetFormalCharge(1)

    mol2 = oechem.OEGraphMol()
    mol2_c = mol2.NewAtom(6)
    mol2_c.SetFormalCharge(-1)

    with pytest.raises(AssertionError, match="Atom formal charge unequal: 1 != -1"):
        assert_atom_equal(mol1_c, mol2_c)
    assert not atom_equal(mol1_c, mol2_c)

    with pytest.raises(AssertionError, match="Atoms with index 0 are not equal") as exc:
        assert_mol_equal(mol1, mol2)
    assert str(exc.value.__cause__) == "Atom formal charge unequal: 1 != -1"
    assert not mol_equal(mol1, mol2)


def test_atom_unequal_partial_charge():
    mol1 = oechem.OEGraphMol()
    mol1_c = mol1.NewAtom(6)
    mol1_c.SetPartialCharge(0.5)

    mol2 = oechem.OEGraphMol()
    mol2_c = mol2.NewAtom(6)
    mol2_c.SetPartialCharge(-0.5)

    with pytest.raises(AssertionError, match="Atom partial charge unequal: 0.5 != -0.5"):
        assert_atom_equal(mol1_c, mol2_c)
    assert not atom_equal(mol1_c, mol2_c)

    with pytest.raises(AssertionError, match="Atoms with index 0 are not equal") as exc:
        assert_mol_equal(mol1, mol2)
    assert str(exc.value.__cause__) == "Atom partial charge unequal: 0.5 != -0.5"
    assert not mol_equal(mol1, mol2)


def test_atom_unequal_name():
    mol1 = oechem.OEGraphMol()
    mol1_c = mol1.NewAtom(6)
    mol1_c.SetName("C1")

    mol2 = oechem.OEGraphMol()
    mol2_c = mol2.NewAtom(6)
    mol2_c.SetName("C2")

    with pytest.raises(AssertionError, match="Atom name unequal: C1 != C2"):
        assert_atom_equal(mol1_c, mol2_c)
    assert not atom_equal(mol1_c, mol2_c)
    assert_atom_equal(mol1_c, mol2_c, properties=AllAtomProperties & ~AtomProperty.Name)

    with pytest.raises(AssertionError, match="Atoms with index 0 are not equal") as exc:
        assert_mol_equal(mol1, mol2)
    assert str(exc.value.__cause__) == "Atom name unequal: C1 != C2"
    assert not mol_equal(mol1, mol2)


def test_atom_unequal_aromatic():
    mol1 = oechem.OEGraphMol()
    mol1_c = mol1.NewAtom(6)
    mol1_c.SetAromatic(True)

    mol2 = oechem.OEGraphMol()
    mol2_c = mol2.NewAtom(6)
    mol2_c.SetAromatic(False)

    with pytest.raises(AssertionError, match="Atom aromaticity unequal: True != False"):
        assert_atom_equal(mol1_c, mol2_c)
    assert not atom_equal(mol1_c, mol2_c)

    with pytest.raises(AssertionError, match="Atoms with index 0 are not equal") as exc:
        assert_mol_equal(mol1, mol2)
    assert str(exc.value.__cause__) == "Atom aromaticity unequal: True != False"
    assert not mol_equal(mol1, mol2)


def test_atom_unequal_chiral():
    mol1 = oechem.OEGraphMol()
    mol1_c = mol1.NewAtom(6)
    mol1_c.SetChiral(True)

    mol2 = oechem.OEGraphMol()
    mol2_c = mol2.NewAtom(6)
    mol2_c.SetChiral(False)

    with pytest.raises(AssertionError, match="Atom chirality unequal: True != False"):
        assert_atom_equal(mol1_c, mol2_c)
    assert not atom_equal(mol1_c, mol2_c)
    assert_atom_equal(mol1_c, mol2_c, properties=AllAtomProperties & ~AtomProperty.Chiral)

    with pytest.raises(AssertionError, match="Atoms with index 0 are not equal") as exc:
        assert_mol_equal(mol1, mol2)
    assert str(exc.value.__cause__) == "Atom chirality unequal: True != False"
    assert not mol_equal(mol1, mol2)


def test_atom_unequal_coords():
    mol1 = oechem.OEGraphMol()
    mol1_c = mol1.NewAtom(6)
    mol1.SetCoords(oechem.OEFloatArray([0.0, 0.0, 0.0]))

    mol2 = oechem.OEGraphMol()
    mol2_c = mol2.NewAtom(6)
    mol2.SetCoords(oechem.OEFloatArray([1.0, 1.0, 1.0]))

    with pytest.raises(AssertionError) as exc:
        assert str(exc.value.__cause__) == "Atom coordinates unequal: [0. 0. 0.] != [1. 1. 1.]"
        assert_atom_equal(mol1_c, mol2_c)
    assert not atom_equal(mol1_c, mol2_c)
    assert_atom_equal(mol1_c, mol2_c, properties=AllAtomProperties & ~AtomProperty.Coordinates)

    with pytest.raises(AssertionError, match="Atoms with index 0 are not equal") as exc:
        assert_mol_equal(mol1, mol2)
    assert str(exc.value.__cause__) == "Atom coordinates unequal: [0. 0. 0.] != [1. 1. 1.]"
    assert not mol_equal(mol1, mol2)


def test_atom_unequal_inring():
    mol1 = oechem.OEGraphMol()
    mol1_c = mol1.NewAtom(6)
    mol1_c.SetInRing(True)

    mol2 = oechem.OEGraphMol()
    mol2_c = mol2.NewAtom(6)
    mol2_c.SetInRing(False)

    with pytest.raises(AssertionError, match="Atom ring membership unequal: True != False"):
        assert_atom_equal(mol1_c, mol2_c)
    assert not atom_equal(mol1_c, mol2_c)

    with pytest.raises(AssertionError, match="Atoms with index 0 are not equal") as exc:
        assert_mol_equal(mol1, mol2)
    assert str(exc.value.__cause__) == "Atom ring membership unequal: True != False"
    assert not mol_equal(mol1, mol2)


def test_bond_unequal_order():
    mol1 = oechem.OEGraphMol()
    mol1_c = mol1.NewAtom(6)
    mol1_o = mol1.NewAtom(8)
    mol1_bond = mol1.NewBond(mol1_c, mol1_o, order=1)

    mol2 = oechem.OEGraphMol()
    mol2_c = mol2.NewAtom(6)
    mol2_o = mol2.NewAtom(8)
    mol2_bond = mol2.NewBond(mol2_c, mol2_o, order=2)

    with pytest.raises(AssertionError, match="Bond order unequal: 1 != 2"):
        assert_bond_equal(mol1_bond, mol2_bond)
    assert not bond_equal(mol1_bond, mol2_bond)

    with pytest.raises(AssertionError):
        assert_mol_equal(mol1, mol2)
    assert not mol_equal(mol1, mol2)


def test_bond_unequal_aromatic():
    mol1 = oechem.OEGraphMol()
    mol1_c = mol1.NewAtom(6)
    mol1_o = mol1.NewAtom(8)
    mol1_bond = mol1.NewBond(mol1_c, mol1_o)
    mol1_bond.SetAromatic(True)

    mol2 = oechem.OEGraphMol()
    mol2_c = mol2.NewAtom(6)
    mol2_o = mol2.NewAtom(8)
    mol2_bond = mol2.NewBond(mol2_c, mol2_o)
    mol2_bond.SetAromatic(False)

    with pytest.raises(AssertionError, match="Bond aromaticity unequal: True != False"):
        assert_bond_equal(mol1_bond, mol2_bond)
    assert not bond_equal(mol1_bond, mol2_bond)

    with pytest.raises(AssertionError, match="Bonds with index 0 are not equal") as exc:
        assert_mol_equal(mol1, mol2)
    assert str(exc.value.__cause__) == "Bond aromaticity unequal: True != False"
    assert not mol_equal(mol1, mol2)


def test_bond_unequal_chiral():
    mol1 = oechem.OEGraphMol()
    mol1_c = mol1.NewAtom(6)
    mol1_o = mol1.NewAtom(8)
    mol1_bond = mol1.NewBond(mol1_c, mol1_o)
    mol1_bond.SetChiral(True)

    mol2 = oechem.OEGraphMol()
    mol2_c = mol2.NewAtom(6)
    mol2_o = mol2.NewAtom(8)
    mol2_bond = mol2.NewBond(mol2_c, mol2_o)
    mol2_bond.SetChiral(False)

    with pytest.raises(AssertionError, match="Bond chirality unequal: True != False"):
        assert_bond_equal(mol1_bond, mol2_bond)
    assert not bond_equal(mol1_bond, mol2_bond)

    with pytest.raises(AssertionError, match="Bonds with index 0 are not equal") as exc:
        assert_mol_equal(mol1, mol2)
    assert str(exc.value.__cause__) == "Bond chirality unequal: True != False"
    assert not mol_equal(mol1, mol2)


def test_bond_unequal_ring_membership():
    mol1 = oechem.OEGraphMol()
    mol1_c = mol1.NewAtom(6)
    mol1_o = mol1.NewAtom(8)
    mol1_bond = mol1.NewBond(mol1_c, mol1_o)
    mol1_bond.SetInRing(True)

    mol2 = oechem.OEGraphMol()
    mol2_c = mol2.NewAtom(6)
    mol2_o = mol2.NewAtom(8)
    mol2_bond = mol2.NewBond(mol2_c, mol2_o)
    mol2_bond.SetInRing(False)

    with pytest.raises(AssertionError, match="Bond ring membership unequal: True != False"):
        assert_bond_equal(mol1_bond, mol2_bond)
    assert not bond_equal(mol1_bond, mol2_bond)

    with pytest.raises(AssertionError, match="Bonds with index 0 are not equal") as exc:
        assert_mol_equal(mol1, mol2)
    assert str(exc.value.__cause__) == "Bond ring membership unequal: True != False"
    assert not mol_equal(mol1, mol2)


@pytest.mark.parametrize(
    ["smiles1", "smiles2", "equal"],
    [
        ("[C@](F)(Cl)(Br)I", "[C@](F)(Cl)(Br)I", True),
        ("[C@](F)(Cl)(Br)I", "[C@@](F)(Cl)(Br)I", False),
    ],
)
def test_equal_atom_stereo(smiles1, smiles2, equal):
    assert mol_equal(oemol_from_smiles(smiles1), oemol_from_smiles(smiles2)) == equal


@pytest.mark.parametrize(
    ["smiles1", "smiles2", "equal"],
    [
        (r"C/C=C\C", r"C\C=C/C", True),
        (r"C/C=C/C", r"C\C=C\C", True),
        (r"C/C=C\C", r"C/C=C/C", False),
    ],
)
def test_equal_bond_stereo(smiles1, smiles2, equal):
    assert mol_equal(oemol_from_smiles(smiles1), oemol_from_smiles(smiles2)) == equal


@pytest.mark.parametrize(
    ["setup_residue1", "setup_residue2", "error_msg", "property"],
    [
        (
            lambda res: res.SetName("ALA"),
            lambda res: res.SetName("LEU"),
            "Residue name unequal: ALA != LEU",
            AtomProperty.ResidueName,
        ),
        (
            lambda res: res.SetAlternateLocation(" "),
            lambda res: res.SetAlternateLocation("A"),
            "Alternate location unequal:   != A",
            AtomProperty.AlternateLocation,
        ),
        (
            lambda res: res.SetBFactor(1.0),
            lambda res: res.SetBFactor(0.5),
            "B factor unequal: 1.0 != 0.5",
            AtomProperty.BFactor,
        ),
        (
            lambda res: res.SetExtChainID("A"),
            lambda res: res.SetExtChainID("B"),
            "Chain ID unequal: A != B",
            AtomProperty.ChainID,
        ),
        (
            lambda res: res.SetSubChainID("A"),
            lambda res: res.SetSubChainID("B"),
            "Asymmetric Unit ID unequal: A != B",
            AtomProperty.AsymmetricUnitID,
        ),
        (
            lambda res: res.SetInsertCode(" "),
            lambda res: res.SetInsertCode("B"),
            "Insert code unequal:   != B",
            AtomProperty.InsertionCode,
        ),
        (
            lambda res: res.SetOccupancy(1.0),
            lambda res: res.SetOccupancy(0.5),
            "Occupancy unequal: 1.0 != 0.5",
            AtomProperty.Occupancy,
        ),
        (
            lambda res: res.SetResidueNumber(1),
            lambda res: res.SetResidueNumber(2),
            "Residue number unequal: 1 != 2",
            AtomProperty.ResidueNumber,
        ),
        (
            lambda res: res.SetSerialNumber(1),
            lambda res: res.SetSerialNumber(2),
            "Serial number unequal: 1 != 2",
            AtomProperty.SerialNumber,
        ),
        (
            lambda res: res.SetSequenceID(1),
            lambda res: res.SetSequenceID(2),
            "Sequence ID unequal: 1 != 2",
            AtomProperty.SequenceNumber,
        ),
        (
            lambda res: res.SetHetAtom(True),
            lambda res: res.SetHetAtom(False),
            "HETATM unequal: True != False",
            AtomProperty.HETATM,
        ),
    ],
)
def test_atom_unequal_residue_property(setup_residue1, setup_residue2, error_msg, property):
    mol1 = oechem.OEGraphMol()
    mol1_atom = mol1.NewAtom(6)
    residue1 = oechem.OEResidue()
    setup_residue1(residue1)
    oechem.OEAtomSetResidue(mol1_atom, residue1)

    mol2 = oechem.OEGraphMol()
    mol2_atom = mol2.NewAtom(6)
    residue2 = oechem.OEResidue()
    setup_residue2(residue2)
    oechem.OEAtomSetResidue(mol2_atom, residue2)

    with pytest.raises(AssertionError, match=error_msg):
        assert_atom_equal(mol1_atom, mol2_atom)
    assert not atom_equal(mol1_atom, mol2_atom)
    assert not atom_equal(mol2_atom, mol1_atom)
    assert_atom_equal(mol1_atom, mol2_atom, properties=AllAtomProperties & ~property)
    assert atom_equal(mol1_atom, mol2_atom, properties=AllAtomProperties & ~property)
    assert atom_equal(mol2_atom, mol1_atom, properties=AllAtomProperties & ~property)

    with pytest.raises(AssertionError, match="Atoms with index 0 are not equal") as exc:
        assert_mol_equal(mol1, mol2)
    assert str(exc.value.__cause__) == error_msg
    assert not mol_equal(mol1, mol2)
    assert not mol_equal(mol2, mol1)
    assert_mol_equal(mol1, mol2, atom_properties=AllAtomProperties & ~property)
    assert mol_equal(mol1, mol2, atom_properties=AllAtomProperties & ~property)
    assert mol_equal(mol2, mol1, atom_properties=AllAtomProperties & ~property)


def test_occupancy_atol():
    mol1 = oechem.OEGraphMol()
    atom1 = create_atom(mol1, atomic_number=6, occupancy=0.5)

    mol2 = oechem.OEGraphMol()
    atom2 = create_atom(mol2, atomic_number=6, occupancy=0.505)

    with pytest.raises(AssertionError) as e:
        assert_atom_equal(atom1, atom2)
    assert str(e.value) == "Occupancy unequal: 0.5 != 0.5049999952316284"

    assert_atom_equal(atom1, atom2, occupancy_atol=0.01)
