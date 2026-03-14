import pytest
from openeye import oechem

from chemutils.molecule import get_subset, smiles_from_oemol
from chemutils.protein import create_polymer, create_protein
from chemutils.protein.create import SEQUENCE_TYPE_TO_INTERRESIDUE_ATOM_NAMES


@pytest.mark.parametrize(
    ["kwargs", "smiles", "residue_names"],
    [
        pytest.param(
            {"sequence": "G", "sequence_type": "protein"},
            "C(C(=O)O)N",
            ["GLY"],
            id="single-protein",
        ),
        pytest.param(
            {"sequence": "A", "sequence_type": "dna"},
            "c1nc(c2c(n1)n(cn2)[C@H]3C[C@@H]([C@H](O3)COP(=O)(O)O)O)N",
            ["DA"],
            id="single-dna",
        ),
        pytest.param(
            {"sequence": "U", "sequence_type": "rna"},
            "c1cn(c(=O)[nH]c1=O)[C@H]2[C@@H]([C@@H]([C@H](O2)COP(=O)(O)O)O)O",
            ["U"],
            id="single-rna",
        ),
        pytest.param(
            {"sequence": "MP(MSE)", "sequence_type": "protein"},
            "CSCC[C@@H](C(=O)N1CCC[C@H]1C(=O)N[C@@H](CC[Se]C)C(=O)O)N",
            ["MET", "PRO", "MSE"],
            id="protein-nonstandard",
        ),
        pytest.param(
            {"sequence": "G(PSU)C", "sequence_type": "rna"},
            "c1cn(c(=O)nc1N)[C@H]2[C@@H]([C@@H]([C@H](O2)COP(=O)(O)O)OP(=O)(O)OC[C@@H]3[C@H]([C@H]([C@@H](O3)c4c[nH]c(=O)[nH]c4=O)O)OP(=O)(O)OC[C@@H]5[C@H]([C@H]([C@@H](O5)n6cnc7c6nc([nH]c7=O)N)O)O)O",
            ["G", "PSU", "C"],
            id="rna-nonstandard",
        ),
    ],
)
def test_create_polymer(kwargs, smiles, residue_names):
    polymer = create_polymer(**kwargs)

    assert smiles_from_oemol(polymer) == smiles

    assert [residue.GetName() for residue in oechem.OEGetResidues(polymer)] == residue_names

    # Check that the residue numbers and serial numbers are correct
    assert [residue.GetResidueNumber() for residue in oechem.OEGetResidues(polymer)] == list(
        range(1, len(residue_names) + 1)
    )
    assert [oechem.OEAtomGetResidue(atom).GetSerialNumber() for atom in polymer.GetAtoms()] == list(
        range(1, polymer.NumAtoms() + 1)
    )

    # Check that atom names make sense
    assert all(bool(oechem.OEAtomGetResidue(atom).GetName()) for atom in polymer.GetAtoms())

    start_atom_name, end_atom_name, condensed_atom_name = SEQUENCE_TYPE_TO_INTERRESIDUE_ATOM_NAMES[
        kwargs["sequence_type"]
    ]

    # Check that a terminal OXT/OP3 is present
    assert get_subset(polymer, selection=f'name "{condensed_atom_name}"').NumAtoms() == 1
    # Check that there are N-1 peptide/phosphodiester bonds
    assert (
        get_subset(
            polymer, selection=f'name "{start_atom_name}" and bonded to name "{end_atom_name}"'
        ).NumAtoms()
        == len(residue_names) - 1
    )


@pytest.mark.parametrize(
    ["kwargs", "smiles", "residue_names"],
    [
        # Single amino acids
        ({"sequence": "A"}, "C[C@@H](C(=O)O)N", ["ALA"]),
        ({"residue_names": ["LEU"]}, "CC(C)C[C@@H](C(=O)O)N", ["LEU"]),
        # Longer chains
        (
            {"sequence": "HGE"},
            "c1c([nH+]c[nH]1)C[C@@H](C(=O)NCC(=O)N[C@@H](CCC(=O)O)C(=O)O)N",
            ["HIS", "GLY", "GLU"],
        ),
        (
            {"residue_names": ["PHE", "ALA", "VAL"]},
            "C[C@@H](C(=O)N[C@@H](C(C)C)C(=O)O)NC(=O)[C@H](Cc1ccccc1)N",
            ["PHE", "ALA", "VAL"],
        ),
    ],
)
def test_create_protein(kwargs, smiles, residue_names):
    protein = create_protein(**kwargs)

    assert smiles_from_oemol(protein) == smiles

    assert [residue.GetName() for residue in oechem.OEGetResidues(protein)] == residue_names

    # Check that the residue numbers and serial numbers are correct
    assert [residue.GetResidueNumber() for residue in oechem.OEGetResidues(protein)] == list(
        range(1, len(residue_names) + 1)
    )
    assert [oechem.OEAtomGetResidue(atom).GetSerialNumber() for atom in protein.GetAtoms()] == list(
        range(1, protein.NumAtoms() + 1)
    )

    # Check that atom names make sense
    assert all(bool(oechem.OEAtomGetResidue(atom).GetName()) for atom in protein.GetAtoms())
    # Check that a terminal OXT is present
    assert get_subset(protein, selection="name OXT").NumAtoms() == 1
    # Check that there are N-1 peptide bonds
    assert (
        get_subset(protein, selection="name N and bonded to name C").NumAtoms()
        == len(residue_names) - 1
    )


def test_create_protein_residue_numbering():
    protein = create_protein(sequence="ALE")

    assert [residue.GetResidueNumber() for residue in oechem.OEGetResidues(protein)] == [1, 2, 3]

    protein = create_protein(sequence="ALE", start_residue_number=5)

    assert [residue.GetResidueNumber() for residue in oechem.OEGetResidues(protein)] == [5, 6, 7]

    protein = create_protein(sequence="ALE", residue_numbers=[3, 5, 6])

    assert [residue.GetResidueNumber() for residue in oechem.OEGetResidues(protein)] == [3, 5, 6]


def test_create_protein_chain_id():
    protein = create_protein(sequence="ALE")

    assert [residue.GetChainID() for residue in oechem.OEGetResidues(protein)] == [" ", " ", " "]

    protein = create_protein(sequence="ALE", chain_id="A")

    assert [residue.GetChainID() for residue in oechem.OEGetResidues(protein)] == ["A", "A", "A"]
