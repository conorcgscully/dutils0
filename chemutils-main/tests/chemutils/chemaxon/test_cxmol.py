import pytest
from openeye import oechem

from chemutils.chemaxon import cxmol_from_oemol, cxmol_from_smiles
from chemutils.chemaxon.jvm import isjavainstance
from chemutils.chemaxon.jvm.chemaxon.struc import BondType, MolAtom, MolBond, Molecule
from chemutils.chemaxon.molecule import as_cxmol, oemol_from_cxmol, smiles_from_cxmol
from chemutils.molecule import make_hydrogens_implicit, oemol_from_smiles

from ..props.drug_smiles import DRUG_SMILES

SMILES = [*list(DRUG_SMILES.values()), "C[C@@H](C(=O)O)N"]

CXBONDTYPE_TO_ORDER = {
    BondType.SINGLE: 1,
    BondType.DOUBLE: 2,
    BondType.TRIPLE: 3,
}


@pytest.mark.parametrize("smiles", SMILES)
def test_cxmol_from_oemol(smiles):
    oemol = oemol_from_smiles(smiles)
    make_hydrogens_implicit(oemol, remove_isotopic_hydrogens=True)
    atom_mapping: dict[oechem.OEAtomBase, MolAtom] = {}
    bond_mapping: dict[oechem.OEBondBase, MolBond] = {}

    cxmol = cxmol_from_oemol(oemol, atom_mapping=atom_mapping, bond_mapping=bond_mapping)

    for oeatom, cxatom in atom_mapping.items():
        assert oeatom.GetAtomicNum() == cxatom.getAtno()
        assert oeatom.GetValence() == cxatom.getValence()

    for oebond, cxbond in bond_mapping.items():
        oe_aromatic = oebond.IsAromatic()
        cx_aromatic = cxbond.getBondType() == BondType.AROMATIC
        assert oe_aromatic == cx_aromatic
        if not oe_aromatic and not cx_aromatic:
            assert oebond.GetOrder() == CXBONDTYPE_TO_ORDER[cxbond.getBondType()]

    cxmol_smiles = smiles_from_cxmol(cxmol)

    assert oechem.OEMolToSmiles(oemol) == oechem.OEMolToSmiles(oemol_from_smiles(cxmol_smiles))


@pytest.mark.parametrize("smiles", SMILES)
def test_cxmol_from_smiles(smiles):
    oemol = oemol_from_smiles(smiles)
    cxmol = cxmol_from_smiles(smiles)

    cxmol_smiles = smiles_from_cxmol(cxmol)

    assert oechem.OEMolToSmiles(oemol) == oechem.OEMolToSmiles(oemol_from_smiles(cxmol_smiles))


@pytest.mark.parametrize("smiles", SMILES)
def test_oemol_from_cxmol(smiles):
    oemol = oemol_from_smiles(smiles)

    cxmol = cxmol_from_oemol(oemol)
    oemol_round_trip = oemol_from_cxmol(cxmol)

    assert oechem.OEMolToSmiles(oemol) == oechem.OEMolToSmiles(oemol_round_trip)


def test_as_cxmol_cxmol():
    oemol = oemol_from_smiles("CCC")
    cxmol = cxmol_from_oemol(oemol)
    cxmol2 = as_cxmol(cxmol)
    assert cxmol is cxmol2


def test_as_cxmol_smiles():
    cxmol = as_cxmol("CCC")
    assert isjavainstance(cxmol, Molecule)


def test_as_cxmol_oemol():
    oemol = oemol_from_smiles("CCC")
    cxmol = as_cxmol(oemol)
    assert isjavainstance(cxmol, Molecule)
