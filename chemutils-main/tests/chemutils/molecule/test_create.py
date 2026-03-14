import pytest
from openeye import oechem

from chemutils.molecule import (
    assert_atom_equal,
    create_atom,
    create_bond,
    create_residue,
    smiles_from_oemol,
)
from chemutils.molecule.equal import AllAtomProperties, AtomProperty

# Don't compare atom indices
PROPERTIES = AllAtomProperties & ~AtomProperty.Index


def test_create_atom_default_equal():
    mol = oechem.OEGraphMol()
    atom1 = create_atom(mol, atomic_number=6)
    atom2 = create_atom(mol, atomic_number=6)
    assert_atom_equal(atom1, atom2, properties=PROPERTIES)
    assert mol.NumAtoms() == 2


def test_create_atom_copy_equal():
    mol = oechem.OEGraphMol()
    atom1 = create_atom(mol, atomic_number=6, coordinates=(1, -1, 0), partial_charge=0.5)
    atom2 = create_atom(mol, atom=atom1)
    assert_atom_equal(atom1, atom2, properties=PROPERTIES)


def test_create_atom_atomic_number():
    mol = oechem.OEGraphMol()
    atom = create_atom(mol, atomic_number=6)
    assert atom.GetAtomicNum() == 6


def test_create_atom_coordinates():
    mol = oechem.OEGraphMol()
    atom = create_atom(mol, atomic_number=6, coordinates=(1, -2, 0))
    assert atom.GetParent().GetCoords(atom) == (1, -2, 0)


def test_create_atom_coordinates_default():
    mol = oechem.OEGraphMol()
    atom = create_atom(mol, atomic_number=6)
    assert atom.GetParent().GetCoords(atom) == (0, 0, 0)


def test_create_atom_name():
    mol = oechem.OEGraphMol()
    atom = create_atom(mol, atomic_number=6, name="ABC")
    assert atom.GetName() == "ABC"


def test_create_atom_name_default():
    mol = oechem.OEGraphMol()
    atom = create_atom(mol, atomic_number=6)
    assert atom.GetName() == ""


def test_create_atom_isotope():
    mol = oechem.OEGraphMol()
    atom = create_atom(mol, atomic_number=6, isotope=13)
    assert atom.GetIsotope() == 13


def test_create_atom_isotope_default():
    mol = oechem.OEGraphMol()
    atom = create_atom(mol, atomic_number=6)
    assert atom.GetIsotope() == 0


def test_create_atom_formal_charge():
    mol = oechem.OEGraphMol()
    atom = create_atom(mol, atomic_number=6, formal_charge=1)
    assert atom.GetFormalCharge() == 1


def test_create_atom_formal_charge_default():
    mol = oechem.OEGraphMol()
    atom = create_atom(mol, atomic_number=6)
    assert atom.GetFormalCharge() == 0


def test_create_atom_partial_charge():
    mol = oechem.OEGraphMol()
    atom = create_atom(mol, atomic_number=6, partial_charge=0.5)
    assert atom.GetPartialCharge() == 0.5


def test_create_atom_partial_charge_default():
    mol = oechem.OEGraphMol()
    atom = create_atom(mol, atomic_number=6)
    assert atom.GetPartialCharge() == 0.0


def test_create_atom_implicit_h_count():
    mol = oechem.OEGraphMol()
    atom = create_atom(mol, atomic_number=6, implicit_h_count=1)
    assert atom.GetImplicitHCount() == 1


def test_create_atom_implicit_h_count_default():
    mol = oechem.OEGraphMol()
    atom = create_atom(mol, atomic_number=6)
    assert atom.GetImplicitHCount() == 4


def test_create_atom_implicit_h_count_0():
    mol = oechem.OEGraphMol()
    atom = create_atom(mol, atomic_number=6, implicit_h_count=0)
    assert atom.GetImplicitHCount() == 0


def test_create_atom_residue_name():
    mol = oechem.OEGraphMol()
    atom = create_atom(mol, atomic_number=6, residue_name="ABC")
    assert oechem.OEAtomGetResidue(atom).GetName() == "ABC"


def test_create_atom_residue_name_default():
    mol = oechem.OEGraphMol()
    atom = create_atom(mol, atomic_number=6)
    assert oechem.OEAtomGetResidue(atom).GetName() == "UNL"


def test_create_atom_occupancy():
    mol = oechem.OEGraphMol()
    atom = create_atom(mol, atomic_number=6, occupancy=0.5)
    assert oechem.OEAtomGetResidue(atom).GetOccupancy() == 0.5


def test_create_atom_occupancy_default():
    mol = oechem.OEGraphMol()
    atom = create_atom(mol, atomic_number=6)
    assert oechem.OEAtomGetResidue(atom).GetOccupancy() == 1.0


def test_create_atom_b_factor():
    mol = oechem.OEGraphMol()
    atom = create_atom(mol, atomic_number=6, b_factor=0.5)
    assert oechem.OEAtomGetResidue(atom).GetBFactor() == 0.5


def test_create_atom_b_factor_default():
    mol = oechem.OEGraphMol()
    atom = create_atom(mol, atomic_number=6)
    assert oechem.OEAtomGetResidue(atom).GetBFactor() == 20.0


def test_create_atom_residue_number():
    mol = oechem.OEGraphMol()
    atom = create_atom(mol, atomic_number=6, residue_number=5)
    assert oechem.OEAtomGetResidue(atom).GetResidueNumber() == 5


def test_create_atom_residue_number_default():
    mol = oechem.OEGraphMol()
    atom = create_atom(mol, atomic_number=6)
    assert oechem.OEAtomGetResidue(atom).GetResidueNumber() == 1


def test_create_atom_sequence_number():
    mol = oechem.OEGraphMol()
    atom = create_atom(mol, atomic_number=6, sequence_number=1)
    assert oechem.OEAtomGetResidue(atom).GetSequenceID() == 1


def test_create_atom_sequence_number_default():
    mol = oechem.OEGraphMol()
    atom = create_atom(mol, atomic_number=6)
    assert oechem.OEAtomGetResidue(atom).GetSequenceID() == -1


def test_create_atom_serial_number():
    mol = oechem.OEGraphMol()
    atom = create_atom(mol, atomic_number=6, serial_number=4)
    assert oechem.OEAtomGetResidue(atom).GetSerialNumber() == 4


def test_create_atom_serial_number_default():
    mol = oechem.OEGraphMol()
    atom = create_atom(mol, atomic_number=6)
    assert oechem.OEAtomGetResidue(atom).GetSerialNumber() == 0


def test_create_atom_alternate_location():
    mol = oechem.OEGraphMol()
    atom = create_atom(mol, atomic_number=6, alternate_location="A")
    assert oechem.OEAtomGetResidue(atom).GetAlternateLocation() == "A"


def test_create_atom_alternate_location_default():
    mol = oechem.OEGraphMol()
    atom = create_atom(mol, atomic_number=6)
    assert oechem.OEAtomGetResidue(atom).GetAlternateLocation() == " "


def test_create_atom_chain_id():
    mol = oechem.OEGraphMol()
    atom = create_atom(mol, atomic_number=6, chain_id="A")
    assert oechem.OEAtomGetResidue(atom).GetExtChainID() == "A"


def test_create_atom_chain_id_default():
    mol = oechem.OEGraphMol()
    atom = create_atom(mol, atomic_number=6)
    assert oechem.OEAtomGetResidue(atom).GetExtChainID() == " "


def test_create_atom_asymmetric_unit_id():
    mol = oechem.OEGraphMol()
    atom = create_atom(mol, atomic_number=6, asymmetric_unit_id="A")
    assert oechem.OEAtomGetResidue(atom).GetSubChainID() == "A"


def test_create_atom_asymmetric_unit_id_default():
    mol = oechem.OEGraphMol()
    atom = create_atom(mol, atomic_number=6)
    assert oechem.OEAtomGetResidue(atom).GetSubChainID() == " "


def test_create_atom_insertion_code():
    mol = oechem.OEGraphMol()
    atom = create_atom(mol, atomic_number=6, insertion_code="A")
    assert oechem.OEAtomGetResidue(atom).GetInsertCode() == "A"


def test_create_atom_insertion_code_default():
    mol = oechem.OEGraphMol()
    atom = create_atom(mol, atomic_number=6)
    assert oechem.OEAtomGetResidue(atom).GetInsertCode() == " "


def test_create_atom_is_hetatm():
    mol = oechem.OEGraphMol()
    atom = create_atom(mol, atomic_number=6, is_hetatm=False)
    assert not oechem.OEAtomGetResidue(atom).IsHetAtom()


def test_create_atom_is_hetatm_default():
    mol = oechem.OEGraphMol()
    atom = create_atom(mol, atomic_number=6, is_hetatm=True)
    assert oechem.OEAtomGetResidue(atom).IsHetAtom()


def test_create_atom_residue():
    mol = oechem.OEGraphMol()
    residue = create_residue(residue_name="ABC", occupancy=0.5)
    atom = create_atom(mol, atomic_number=6, residue=residue, occupancy=0.8)

    # Property inhertied from the residue
    assert oechem.OEAtomGetResidue(atom).GetName() == "ABC"
    # Property overriden from the residue
    assert oechem.OEAtomGetResidue(atom).GetOccupancy() == pytest.approx(0.8)


def test_create_residue_residue_name():
    residue = create_residue(residue_name="ABC")
    assert residue.GetName() == "ABC"


def test_create_residue_residue_name_default():
    residue = create_residue()
    assert residue.GetName() == "UNL"


def test_create_residue_occupancy():
    residue = create_residue(occupancy=0.5)
    assert residue.GetOccupancy() == 0.5


def test_create_residue_occupancy_default():
    residue = create_residue()
    assert residue.GetOccupancy() == 1.0


def test_create_residue_b_factor():
    residue = create_residue(b_factor=0.5)
    assert residue.GetBFactor() == 0.5


def test_create_residue_b_factor_default():
    residue = create_residue()
    assert residue.GetBFactor() == 20.0


def test_create_residue_residue_number():
    residue = create_residue(residue_number=5)
    assert residue.GetResidueNumber() == 5


def test_create_residue_residue_number_default():
    residue = create_residue()
    assert residue.GetResidueNumber() == 1


def test_create_residue_sequence_number():
    residue = create_residue(sequence_number=1)
    assert residue.GetSequenceID() == 1


def test_create_residue_sequence_number_default():
    residue = create_residue()
    assert residue.GetSequenceID() == -1


def test_create_residue_serial_number():
    residue = create_residue(serial_number=4)
    assert residue.GetSerialNumber() == 4


def test_create_residue_serial_number_default():
    residue = create_residue()
    assert residue.GetSerialNumber() == 0


def test_create_residue_alternate_location():
    residue = create_residue(alternate_location="A")
    assert residue.GetAlternateLocation() == "A"


def test_create_residue_alternate_location_default():
    residue = create_residue()
    assert residue.GetAlternateLocation() == " "


def test_create_residue_chain_id():
    residue = create_residue(chain_id="A")
    assert residue.GetExtChainID() == "A"


def test_create_residue_chain_id_default():
    residue = create_residue()
    assert residue.GetExtChainID() == " "


def test_create_residue_asymmetric_unit_id():
    residue = create_residue(asymmetric_unit_id="A")
    assert residue.GetSubChainID() == "A"


def test_create_residue_asymmetric_unit_id_default():
    residue = create_residue()
    assert residue.GetSubChainID() == " "


def test_create_residue_insertion_code():
    residue = create_residue(insertion_code="A")
    assert residue.GetInsertCode() == "A"


def test_create_residue_insertion_code_default():
    residue = create_residue()
    assert residue.GetInsertCode() == " "


def test_create_residue_is_hetatm():
    residue = create_residue(is_hetatm=False)
    assert not residue.IsHetAtom()


def test_create_residue_is_hetatm_default():
    residue = create_residue(is_hetatm=True)
    assert residue.IsHetAtom()


def test_create_bond():
    mol = oechem.OEGraphMol()
    atom1 = create_atom(mol, atomic_number=6)
    atom2 = create_atom(mol, atomic_number=7)
    bond = create_bond(mol, atom1=atom1, atom2=atom2, order=2)
    assert list(mol.GetBonds()) == [bond]
    assert bond.GetOrder() == 2
    assert_atom_equal(bond.GetBgn(), atom1)
    assert_atom_equal(bond.GetEnd(), atom2)


@pytest.mark.parametrize(
    ["order", "adjust_h_count", "adjust_charge", "smiles"],
    [
        (None, True, True, "CN"),
        (1, True, True, "CN"),
        (2, True, True, "C=N"),
        (3, True, True, "C#N"),
        (4, True, True, "C$[N+]"),
        (None, True, False, "CN"),
        (3, True, False, "C#N"),
        (4, True, False, "error"),
        (None, False, True, "[CH4+][NH3+]"),
        (2, False, True, "[CH4+2]=[NH3+2]"),
        (None, False, False, "[CH4][NH3]"),
        (2, False, False, "[CH4]=[NH3]"),
        (4, False, False, "[CH4]$[NH3]"),
    ],
)
def test_create_bond_adjust(order, adjust_h_count, adjust_charge, smiles):
    mol = oechem.OEGraphMol()
    atom1 = create_atom(mol, atomic_number=6)
    atom2 = create_atom(mol, atomic_number=7)
    if smiles == "error":
        with pytest.raises(ValueError):
            create_bond(
                mol,
                atom1=atom1,
                atom2=atom2,
                order=order,
                adjust_h_count=adjust_h_count,
                adjust_charge=adjust_charge,
            )
    else:
        create_bond(
            mol,
            atom1=atom1,
            atom2=atom2,
            order=order,
            adjust_h_count=adjust_h_count,
            adjust_charge=adjust_charge,
        )
        assert smiles_from_oemol(mol) == smiles
