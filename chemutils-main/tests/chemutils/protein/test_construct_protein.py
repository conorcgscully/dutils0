import io
import os

import numpy as np
import pytest
from Bio.PDB import PDBParser, Structure
from Bio.PDB.PDBIO import PDBIO
from Bio.SeqIO import PdbIO

from chemutils.molecule import read_molecule_str
from chemutils.protein.construct_protein import OPENFOLD_ATOM_NAMES, construct_protein
from chemutils.protein.residues import iterate_residues


def atoms_af2_ordered(structure: Structure):
    # Extract AF2 ordered atoms from a BioPython PDB, which may be in a different order or have additional atoms.
    atoms = []
    for residue in structure.get_residues():
        for atom_name in OPENFOLD_ATOM_NAMES[residue.resname]:
            atom = None
            for res_atom in residue.get_atoms():
                if res_atom.name == atom_name:
                    atom = res_atom
                    break
            if atom is None:
                raise ValueError(f"Residue {residue.resname} does not have atom {atom_name}")
            atoms.append(atom)
    return atoms


@pytest.mark.parametrize(
    "filename",
    ["AF-P82042-F1-model_v3.pdb", "AF-E0CX11-F1-model_v3.pdb", "AF-A2RU14-F1-model_v3.pdb"],
)
@pytest.mark.parametrize("additional_data", [None, "bfactor", "occupancy", "residue_numbers"])
def test_construct_protein(filename, additional_data):
    parser = PDBParser()
    orig_structure = parser.get_structure("", f"{os.path.dirname(__file__)}/data/{filename}")
    orig_atoms = atoms_af2_ordered(orig_structure)
    orig_residues = list(orig_structure.get_residues())
    orig_coords = np.array([atom.coord for atom in orig_atoms])
    orig_sequence = str(next(iter(PdbIO.AtomIterator("", orig_structure))).seq)
    # BioPython inserts residues with 'X' when the residue number jumps
    orig_sequence = orig_sequence.replace("X", "")

    orig_b_factors = np.array(
        [next(iter(residue.get_atoms())).bfactor for residue in orig_structure.get_residues()]
    )
    orig_residue_numbers = np.array([residue.id[1] for residue in orig_structure.get_residues()])

    if additional_data == "bfactor":
        constructed = construct_protein(
            full_atom_coordinates=orig_coords,
            sequence=orig_sequence,
            residue_b_factors=orig_b_factors,
        )
    elif additional_data == "occupancy":
        constructed = construct_protein(
            full_atom_coordinates=orig_coords,
            sequence=orig_sequence,
            residue_occupancies=orig_b_factors,
        )
    elif additional_data == "residue_numbers":
        constructed = construct_protein(
            full_atom_coordinates=orig_coords,
            sequence=orig_sequence,
            residue_numbers=orig_residue_numbers,
        )
    else:
        constructed = construct_protein(full_atom_coordinates=orig_coords, sequence=orig_sequence)

    assert len(orig_atoms) == len(list(constructed.get_atoms()))
    assert len(orig_residues) == len(list(constructed.get_residues()))

    for residue_orig, residue_construct in zip(
        orig_residues, constructed.get_residues(), strict=True
    ):
        assert residue_orig.resname == residue_construct.resname
        if additional_data == "residue_numbers":
            assert residue_orig.id[1] == residue_construct.id[1]
        else:
            assert residue_construct.id[1] == orig_residues.index(residue_orig) + 1

    for atom_orig, atom_construct in zip(orig_atoms, constructed.get_atoms(), strict=True):
        assert atom_orig.name == atom_construct.name
        assert atom_orig.coord == pytest.approx(atom_construct.coord)
        assert atom_orig.fullname.strip() == atom_construct.fullname
        assert atom_orig.element == atom_construct.element
        if additional_data == "bfactor":
            assert atom_orig.bfactor == atom_construct.bfactor
        else:
            assert atom_construct.bfactor == pytest.approx(0.0)
        if additional_data == "occupancy":
            assert atom_orig.bfactor == atom_construct.occupancy
        else:
            assert atom_construct.occupancy == pytest.approx(1.0)

    constructed_sequence = str(next(iter(PdbIO.AtomIterator("", constructed))).seq)
    # BioPython inserts residues with 'X' when the residue number jumps
    constructed_sequence = constructed_sequence.replace("X", "")

    assert orig_sequence == constructed_sequence


@pytest.mark.parametrize(
    "coord_shape",
    [
        (4, 3),  # Too few coords
        (20, 3),  # Too many coords
        (9, 2),  # Too few dimensions
        (9, 4),  # Too many dimensions
        (9,),  # Wrong shape
        (9, 3, 3),  # Wrong shape
    ],
)
def test_construct_protein_mismatch_coords(coord_shape):
    sequence = "GA"
    coords = np.random.random(coord_shape)
    with pytest.raises(ValueError):
        construct_protein(full_atom_coordinates=coords, sequence=sequence)


def test_construct_protein_insertion_codes():
    pdb_str = """ATOM      1  N   ASN B  52      -3.487 -37.610 -22.718  1.00  7.64           N
ATOM      2  CA  ASN B  52      -2.923 -38.940 -22.955  1.00  7.64           C
ATOM      3  C   ASN B  52      -1.518 -38.767 -23.550  1.00  7.64           C
ATOM      4  O   ASN B  52      -1.072 -37.623 -23.713  1.00  7.64           O
ATOM      5  CB  ASN B  52      -3.823 -39.738 -23.883  1.00  7.64           C
ATOM      6  CG  ASN B  52      -3.683 -39.305 -25.298  1.00  7.64           C
ATOM      7  OD1 ASN B  52      -3.846 -38.127 -25.605  1.00  7.64           O
ATOM      8  ND2 ASN B  52      -3.313 -40.237 -26.169  1.00  7.64           N
ATOM      9  N   PRO B  52A     -0.811 -39.882 -23.882  1.00  2.00           N
ATOM     10  CA  PRO B  52A      0.530 -39.717 -24.436  1.00  2.00           C
ATOM     11  C   PRO B  52A      0.551 -38.668 -25.518  1.00  2.00           C
ATOM     12  O   PRO B  52A      0.751 -37.490 -25.226  1.00  2.00           O
ATOM     13  CB  PRO B  52A      0.852 -41.118 -24.971  1.00  2.00           C
ATOM     14  CG  PRO B  52A      0.183 -42.000 -24.023  1.00  2.00           C
ATOM     15  CD  PRO B  52A     -1.161 -41.317 -23.829  1.00  2.00           C
TER      16      PRO B  52A
END
"""
    oemol = read_molecule_str(pdb_str, format="PDB")
    atom_coordinates = np.array([oemol.GetCoords(atom) for atom in oemol.GetAtoms()])
    residue_b_factors = [residue.GetBFactor() for residue in iterate_residues(oemol)]
    residue_occupancies = [residue.GetOccupancy() for residue in iterate_residues(oemol)]
    residue_numbers = [residue.GetResidueNumber() for residue in iterate_residues(oemol)]
    residue_insertion_codes = [residue.GetInsertCode() for residue in iterate_residues(oemol)]
    structure = construct_protein(
        full_atom_coordinates=atom_coordinates,
        sequence="NP",
        residue_b_factors=residue_b_factors,
        residue_occupancies=residue_occupancies,
        residue_numbers=residue_numbers,
        residue_insertion_codes=residue_insertion_codes,
        chain_id="B",
    )
    pdbio = PDBIO()
    pdbio.set_structure(structure)
    out = io.StringIO()
    pdbio.save(out)
    assert [line.strip() for line in out.getvalue().splitlines()] == [
        line.strip() for line in pdb_str.splitlines()
    ]
