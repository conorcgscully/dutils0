import os
from io import StringIO
from pathlib import Path

import numpy as np
import pytest
from Bio.PDB import PDBParser
from openeye import oechem

from chemutils.molecule import read_molecule_str
from chemutils.protein import protein_to_pdbstr


@pytest.mark.parametrize(
    "filename",
    ["AF-P82042-F1-model_v3.pdb", "AF-E0CX11-F1-model_v3.pdb", "AF-A2RU14-F1-model_v3.pdb"],
)
def test_protein_to_pdbstr(filename):
    parser = PDBParser()
    orig_structure = parser.get_structure("", f"{os.path.dirname(__file__)}/../data/{filename}")
    pdbstr = protein_to_pdbstr(orig_structure)
    rec_structure = parser.get_structure("", StringIO(pdbstr))
    assert len(list(orig_structure.get_residues())) == len(list(rec_structure.get_residues()))

    for orig_residue, rec_residue in zip(
        orig_structure.get_residues(), rec_structure.get_residues(), strict=True
    ):
        assert orig_residue.resname == rec_residue.resname
        for orig_atom, rec_atom in zip(
            orig_residue.get_atoms(), rec_residue.get_atoms(), strict=True
        ):
            assert orig_atom.name == rec_atom.name
            assert (orig_atom.coord == rec_atom.coord).all()
            assert orig_atom.bfactor == rec_atom.bfactor
            assert orig_atom.occupancy == rec_atom.occupancy


@pytest.mark.parametrize(
    "filename",
    ["df_output.pdb", "hetatm.pdb", "occupancy.pdb"],
)
def test_read_molecule_str_pdb(filename):
    parser = PDBParser()

    pdb_file = Path(__file__).parent.parent.parent / "data" / "pdb" / filename

    # Parse PDB using Biopython
    structure = parser.get_structure("test", str(pdb_file))
    biopython_b_factors = []
    biopython_occupancies = []
    biopython_atom_names = []
    biopython_coords = []
    biopython_res_names = []
    for model in structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    biopython_b_factors.append(atom.get_bfactor())
                    biopython_occupancies.append(atom.occupancy)
                    biopython_atom_names.append(atom.name)
                    biopython_coords.append(atom.get_coord())
                    biopython_res_names.append(residue.get_resname())

    biopython_b_factors = np.array(biopython_b_factors)
    biopython_occupancies = np.array(biopython_occupancies)
    biopython_atom_names = np.array(biopython_atom_names)
    biopython_coords = np.stack(biopython_coords)
    biopython_res_names = np.array(biopython_res_names)

    # Parse PDB using OpenEye
    with open(str(pdb_file)) as f:
        pdbstr = f.read()
    oemol = read_molecule_str(pdbstr, format="pdb")

    oemol_b_factors = []
    oemol_occupancies = []
    oemol_atom_names = []
    oemol_res_names = []
    oemol_coords = np.stack(list(oemol.GetCoords().values()))
    for atom in oemol.GetAtoms():
        oe_data = oechem.OEAtomGetResidue(atom)
        oemol_b_factors.append(oe_data.GetBFactor())
        oemol_occupancies.append(oe_data.GetOccupancy())
        oemol_atom_names.append(atom.GetName().strip())
        oemol_res_names.append(oe_data.GetName())

    oemol_b_factors = np.array(oemol_b_factors)
    oemol_occupancies = np.array(oemol_occupancies)
    oemol_atom_names = np.array(oemol_atom_names)
    oemol_res_names = np.array(oemol_res_names)

    # Compare outputs
    assert np.allclose(oemol_b_factors, biopython_b_factors), (
        "B factors are inconsistent between Biopython and OpenEye"
    )
    assert np.allclose(oemol_occupancies, biopython_occupancies), (
        "Occupancies are inconsistent between Biopython and OpenEye"
    )
    assert (oemol_atom_names == biopython_atom_names).all(), (
        "Atom names are inconsistent between Biopython and OpenEye"
    )
    assert (oemol_res_names == biopython_res_names).all(), (
        "Residue names are inconsistent between Biopython and OpenEye"
    )
    assert np.allclose(oemol_coords, biopython_coords), (
        "Atom coordinates are inconsistent between Biopython and OpenEye"
    )
