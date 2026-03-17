from __future__ import annotations

import pytest
from rdkit import Chem

from dutils0.rough.select import select_residue_as_mol


def _tiny_pdb_block() -> str:
    # Two residues in chain A (42, 43) and one residue in chain B (42).
    # Keep minimal but valid columns so RDKit sets AtomPDBResidueInfo.
    return "\n".join(
        [
            "ATOM      1  N   ALA A  42      11.104  13.207   2.100  1.00 20.00           N",
            "ATOM      2  CA  ALA A  42      12.560  13.207   2.100  1.00 20.00           C",
            "ATOM      3  C   ALA A  42      12.960  14.607   2.600  1.00 20.00           C",
            "ATOM      4  O   ALA A  42      12.260  15.607   2.300  1.00 20.00           O",
            "ATOM      5  N   GLY A  43      14.160  14.707   3.200  1.00 20.00           N",
            "ATOM      6  CA  GLY A  43      14.660  16.007   3.600  1.00 20.00           C",
            "ATOM      7  N   SER B  42       9.104  10.207   1.100  1.00 20.00           N",
            "ATOM      8  CA  SER B  42       9.560  11.207   1.500  1.00 20.00           C",
            "TER",
            "END",
            "",
        ]
    )


@pytest.fixture
def pdb_mol() -> Chem.Mol:
    mol = Chem.MolFromPDBBlock(_tiny_pdb_block(), removeHs=False)
    assert mol is not None
    assert mol.GetNumAtoms() > 0
    return mol


def _assert_all_atoms_chain_resi(mol: Chem.Mol, chain: str, resi: int) -> None:
    assert mol.GetNumAtoms() > 0
    for a in mol.GetAtoms():
        info = a.GetPDBResidueInfo()
        assert info is not None
        assert info.GetChainId().strip() == chain
        assert info.GetResidueNumber() == resi


def test_select_residue_as_mol_happy_path(pdb_mol: Chem.Mol) -> None:
    sub = select_residue_as_mol(pdb_mol, "A:42")
    assert isinstance(sub, Chem.Mol)
    assert sub is not pdb_mol
    assert sub.GetNumAtoms() > 0
    assert sub.GetNumAtoms() < pdb_mol.GetNumAtoms()
    _assert_all_atoms_chain_resi(sub, "A", 42)


@pytest.mark.parametrize(
    "bad",
    [
        "A42",
        "A:",
        ":42",
        "A:4.2",
        "A:42:extra",
        ":",
        "",
        "   ",
    ],
)
def test_select_residue_as_mol_invalid_selector_format(pdb_mol: Chem.Mol, bad: str) -> None:
    with pytest.raises((TypeError, ValueError)):
        select_residue_as_mol(pdb_mol, bad)


def test_select_residue_as_mol_non_integer_residue(pdb_mol: Chem.Mol) -> None:
    with pytest.raises(ValueError, match="integer"):
        select_residue_as_mol(pdb_mol, "A:xyz")


def test_select_residue_as_mol_residue_not_found(pdb_mol: Chem.Mol) -> None:
    with pytest.raises(ValueError, match="no atoms match"):
        select_residue_as_mol(pdb_mol, "A:999")


def test_select_residue_as_mol_chain_not_found(pdb_mol: Chem.Mol) -> None:
    with pytest.raises(ValueError, match="no atoms match"):
        select_residue_as_mol(pdb_mol, "Z:42")

