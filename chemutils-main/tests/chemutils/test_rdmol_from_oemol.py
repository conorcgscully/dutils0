import os
from collections import Counter
from io import StringIO
from pathlib import Path

import numpy as np
import pytest
from Bio.PDB import PDBParser
from openeye import oechem
from rdkit import Chem
from rdkit.Chem import AllChem, rdchem, rdmolfiles
from rdkit.Geometry import Point3D

from chemutils.molecule import construct_molecule, oemol_from_smiles, read_molecule_str
from chemutils.rdmol_from_oemol import get_oemol_cistrans, rdmol_from_oemol

from .cip_validation import CIP_VALIDATION_OEMOLS
from .props.drug_smiles import DRUG_SMILES

PROPERTIES = {
    "num_atoms": (lambda oemol: oemol.NumAtoms(), lambda rdmol: rdmol.GetNumAtoms()),
    "num_bonds": (lambda oemol: oemol.NumBonds(), lambda rdmol: rdmol.GetNumBonds()),
    "num_double_bonds": (
        lambda oemol: sum(
            bond.GetOrder() == 2 and not bond.IsAromatic() for bond in oemol.GetBonds()
        ),
        lambda rdmol: sum(bond.GetBondType() == Chem.BondType.DOUBLE for bond in rdmol.GetBonds()),
    ),
    "num_triple_bonds": (
        lambda oemol: sum(bond.GetOrder() == 3 for bond in oemol.GetBonds()),
        lambda rdmol: sum(bond.GetBondType() == Chem.BondType.TRIPLE for bond in rdmol.GetBonds()),
    ),
    "num_aromatic_atoms": (
        lambda oemol: sum(atom.IsAromatic() for atom in oemol.GetAtoms()),
        lambda rdmol: sum(atom.GetIsAromatic() for atom in rdmol.GetAtoms()),
    ),
    "num_aromatic_bonds": (
        lambda oemol: sum(bond.IsAromatic() for bond in oemol.GetBonds()),
        lambda rdmol: sum(bond.GetIsAromatic() for bond in rdmol.GetBonds()),
    ),
    "num_stereo_e_bonds": (
        lambda oemol: sum(
            oechem.OEPerceiveCIPStereo(oemol, bond) == oechem.OECIPBondStereo_E
            for bond in oemol.GetBonds()
        ),
        lambda rdmol: sum(
            bond.GetProp("_CIPCode") == "E" for bond in rdmol.GetBonds() if bond.HasProp("_CIPCode")
        ),
    ),
    "num_stereo_z_bonds": (
        lambda oemol: sum(
            oechem.OEPerceiveCIPStereo(oemol, bond) == oechem.OECIPBondStereo_Z
            for bond in oemol.GetBonds()
        ),
        lambda rdmol: sum(
            bond.GetProp("_CIPCode") == "Z" for bond in rdmol.GetBonds() if bond.HasProp("_CIPCode")
        ),
    ),
    "formal_charge_count": (
        lambda oemol: Counter([atom.GetFormalCharge() for atom in oemol.GetAtoms()]),
        lambda rdmol: Counter([atom.GetFormalCharge() for atom in rdmol.GetAtoms()]),
    ),
    "atomic_numbers": (
        lambda oemol: Counter([atom.GetAtomicNum() for atom in oemol.GetAtoms()]),
        lambda rdmol: Counter(
            [Chem.GetPeriodicTable().GetAtomicNumber(atom.GetSymbol()) for atom in rdmol.GetAtoms()]
        ),
    ),
    "num_sp3": (
        lambda oemol: Counter(
            [atom.GetHyb() == oechem.OEHybridization_sp3 for atom in oemol.GetAtoms()]
        ),
        lambda rdmol: Counter(
            [atom.GetHybridization() == Chem.HybridizationType.SP3 for atom in rdmol.GetAtoms()]
        ),
    ),
    "num_sp2": (
        lambda oemol: Counter(
            [atom.GetHyb() == oechem.OEHybridization_sp2 for atom in oemol.GetAtoms()]
        ),
        lambda rdmol: Counter(
            [atom.GetHybridization() == Chem.HybridizationType.SP2 for atom in rdmol.GetAtoms()]
        ),
    ),
    "num_sp1": (
        lambda oemol: Counter(
            [atom.GetHyb() == oechem.OEHybridization_sp for atom in oemol.GetAtoms()]
        ),
        lambda rdmol: Counter(
            [atom.GetHybridization() == Chem.HybridizationType.SP for atom in rdmol.GetAtoms()]
        ),
    ),
}


@pytest.mark.parametrize(
    "smiles",
    [
        "C",
        "CC",
        "C=C",
        "C#C",
        "C=O",
        "[Na+]",  # Charged
        "c1ccccc1",  # Aromatic ring
        r"F/C=C/F",  # E-stereo bond
        r"F/C=C\F",  # Z-stereo bond
        r"CC2(C)CCCC(\C)=C2\C=C\C(\C)=C\C=C\C(\C)=C\C=C\C=C(/C)\C=C\C=C(/C)\C=C\C1=C(/C)CCCC1(C)C",  # beta-carotene
        *DRUG_SMILES.values(),
        # Crazy ligands from PDB:
        r"CC12C3([Ir]1456(C2(C4(C53C)C)C)[N]7=CC=C(C=C7C(=O)[N]6(CCCc8ccc(cc8)S(=O)(=O)N)Cl)O)C",  # JR8
        # Hmm
        r"c1cc(ccc1[C@@H]2CCNC[C@H]2COc3ccc4c(c3)OCO4)F",  # Has explicit H
    ],
)
def test_rdmol_from_oemol(smiles):
    oemol = oemol_from_smiles(smiles)
    oechem.OEAssignAromaticFlags(oemol)
    oechem.OEAssignHybridization(oemol)
    rdmol = rdmol_from_oemol(oemol)
    for property_name, funcs in PROPERTIES.items():
        oemol_func, rdmol_func = funcs
        assert oemol_func(oemol) == rdmol_func(rdmol), f"{property_name} does not match"


def test_rdmol_from_oemol_coords():
    smiles = "CCC"
    oemol = oemol_from_smiles(smiles)
    coords = np.array([[1.0, 0.0, -1.0], [2.0, -0.5, 2.0], [1.0, 1.0, 0.0]])
    oemol.SetCoords(coords.flatten())
    oemol.SetDimension(3)
    rdmol = rdmol_from_oemol(oemol)
    rdmol.GetConformer()
    rd_coords = np.array(
        [rdmol.GetConformer().GetAtomPosition(i) for i in range(rdmol.GetNumAtoms())]
    )
    np.testing.assert_allclose(coords, rd_coords)


@pytest.mark.parametrize(
    ["smiles", "expected_stereo"],
    [
        ("C\\C=C/C", [(oechem.OEBondStereo_Cis, 0, 3)]),
        ("C\\C=C\\C", [(oechem.OEBondStereo_Trans, 0, 3)]),
        ("C\\C=C(/C)\\C", [(oechem.OEBondStereo_Cis, 0, 3)]),
        ("C\\C=C(\\C)/C", [(oechem.OEBondStereo_Trans, 0, 3)]),
        ("C\\C=C\\C=C/C=CC", [(oechem.OEBondStereo_Trans, 0, 3), (oechem.OEBondStereo_Cis, 2, 5)]),
        ("F\\C=C(/Cl)\\I", [(oechem.OEBondStereo_Cis, 0, 3)]),
        ("F\\C=C(\\I)Cl", [(oechem.OEBondStereo_Trans, 0, 3)]),
    ],
)
def test_get_oemol_cistrans(smiles, expected_stereo):
    mol = oemol_from_smiles(smiles)
    stereo = []
    for bond in mol.GetBonds():
        stereo_class, start, end = get_oemol_cistrans(bond)
        if stereo_class != oechem.OEBondStereo_Undefined:
            stereo.append((stereo_class, start, end))
    assert stereo == expected_stereo


# Molecule where CIP fails, and hence cis/trans stereo is used.
def test_rdmol_from_oemol_cistrans():
    sdf_str = (Path(__file__).parent / "data/cistrans_macrocycle.sdf").read_text()
    mol = read_molecule_str(sdf_str, format="sdf")

    num_bonds_no_cip = sum(
        bond.HasStereoSpecified(oechem.OEBondStereo_CisTrans)
        and not oechem.OEPerceiveCIPStereo(mol, bond)
        for bond in mol.GetBonds()
    )
    assert num_bonds_no_cip == 4

    rdmol = rdmol_from_oemol(mol)

    cistrans_bonds = [
        (bond.GetStereo(), bond.GetStereoAtoms()[0], bond.GetStereoAtoms()[1])
        for bond in rdmol.GetBonds()
        if bond.GetStereo() != rdchem.BondStereo.STEREONONE
    ]

    assert cistrans_bonds == [
        (rdchem.BondStereo.STEREOCIS, 5, 8),
        (rdchem.BondStereo.STEREOTRANS, 9, 13),
        (rdchem.BondStereo.STEREOCIS, 14, 17),
        (rdchem.BondStereo.STEREOCIS, 20, 24),
    ]


def assert_pdb_str_equal(pdb_str_A, pdb_str_B):
    parser = PDBParser()
    structure_A = parser.get_structure("", StringIO(pdb_str_A))
    structure_B = parser.get_structure("", StringIO(pdb_str_B))

    assert len(list(structure_A.get_chains())) == len(list(structure_B.get_chains())), (
        "Chain count mismatch"
    )
    for chain_A, chain_B in zip(structure_A.get_chains(), structure_B.get_chains(), strict=True):
        assert chain_A.get_id() == chain_B.get_id(), "Chain ID mismatch"
        assert len(list(chain_A.get_residues())) == len(list(chain_B.get_residues())), (
            "Residue count mismatch"
        )
        for residue_A, residue_B in zip(
            chain_A.get_residues(), chain_B.get_residues(), strict=True
        ):
            assert residue_A.get_id() == residue_B.get_id(), "Residue ID mismatch"
            assert residue_A.get_resname() == residue_B.get_resname(), "Residue name mismatch"
            assert residue_A.get_segid() == residue_B.get_segid(), "Residue segment ID mismatch"
            assert len(list(residue_A.get_atoms())) == len(list(residue_B.get_atoms())), (
                "Atom count mismatch"
            )
            for atom_A, atom_B in zip(residue_A.get_atoms(), residue_B.get_atoms(), strict=True):
                assert atom_A.name == atom_B.name, "Atom name mismatch"
                assert all(atom_A.coord == atom_B.coord), "Coordinate mismatch"
                assert atom_A.bfactor == atom_B.bfactor, "B Factor mismatch"
                assert atom_A.occupancy == atom_B.occupancy, "Occupancy mismatch"
                assert atom_A.altloc == atom_B.altloc, "AltLoc mismatch"
                assert atom_A.element == atom_B.element, "Element mismatch"


@pytest.mark.parametrize("pdb_filename", os.listdir(Path(__file__).parent / "data" / "pdb"))
def test_rdmol_from_oemol_pdb(pdb_filename):
    ifs = oechem.oemolistream(pdb_filename)
    oe_mol = oechem.OEGraphMol()
    oechem.OEReadMolecule(ifs, oe_mol)

    rd_mol = rdmol_from_oemol(oe_mol)

    ofs = oechem.oemolostream()
    ofs.SetFormat(oechem.OEFormat_PDB)
    ofs.openstring()
    oechem.OEWriteMolecule(ofs, oe_mol)
    oe_pdb_str = ofs.GetString().decode()
    rd_pdb_str = rdmolfiles.MolToPDBBlock(rd_mol)

    assert_pdb_str_equal(oe_pdb_str, rd_pdb_str)


@pytest.mark.parametrize(
    ["value", "serialized_value"],
    [
        ("my_value", "my_value"),
        (1, "1"),
        (1.0, "1.0"),
        (True, "True"),
    ],
)
def test_rdmol_from_oemol_sdtags(value, serialized_value):
    smiles = "CCC"
    oemol = oemol_from_smiles(smiles)
    oechem.OESetSDData(oemol, "my_tag", str(value))
    rdmol = rdmol_from_oemol(oemol)
    assert rdmol.GetProp("my_tag") == serialized_value


def _get_rdmol_coords(rdmol):
    coords = []
    for i in range(rdmol.GetNumAtoms()):
        positions = rdmol.GetConformer().GetAtomPosition(i)
        coords.append([positions.x, positions.y, positions.z])
    return np.array(coords)


@pytest.mark.parametrize(
    ["smiles", "coords"],
    [
        ("C", [[0, 0, 0]]),
        ("C", [[1, 0, 0]]),
        ("C", [[1, 1, 0]]),
        ("C", [[2.2, 1.5, 3.9]]),
        ("CC", [[0.2, -0.5, 0.1], [1.1, -0.4, 0.2]]),
        ("C=C", [[0.2, -0.5, 0.1], [1.1, -0.4, 0.2]]),
        ("C#C", [[0.2, -0.5, 0.1], [1.1, -0.4, 0.2]]),
    ],
)
def test_rdmol_from_oemol_addhs_agree(smiles, coords):
    rdmol = Chem.MolFromSmiles(smiles)
    AllChem.EmbedMolecule(rdmol)
    conf = rdmol.GetConformer()
    for i in range(rdmol.GetNumAtoms()):
        x, y, z = coords[i]
        conf.SetAtomPosition(i, Point3D(x, y, z))
    rdmol = Chem.AddHs(rdmol, addCoords=True)

    oemol = construct_molecule(smiles=smiles, coordinates=coords)
    rdmol2 = rdmol_from_oemol(oemol)
    rdmol2 = Chem.AddHs(rdmol2, addCoords=True)

    np.testing.assert_allclose(_get_rdmol_coords(rdmol), _get_rdmol_coords(rdmol2), atol=1e-3)


@pytest.mark.parametrize("smiles", DRUG_SMILES.values())
def test_rdmol_from_oemol_round_trip_smiles(smiles):
    oemol = oemol_from_smiles(smiles)
    rdmol = rdmol_from_oemol(oemol)
    rd_smiles = Chem.MolToSmiles(rdmol)
    oemol2 = oemol_from_smiles(rd_smiles)
    assert oechem.OEMolToSmiles(oemol) == oechem.OEMolToSmiles(oemol2)


@pytest.mark.parametrize("oemol", CIP_VALIDATION_OEMOLS)
def test_rdmol_from_oemol_round_trip_oemol(oemol):
    if oemol.GetTitle() == "VS132":
        pytest.skip("Unusual stereochemistry that is parsed incorrectly by OpenEye")
    orig_smiles = oechem.OEMolToSmiles(oemol)
    rdmol = rdmol_from_oemol(oemol)
    rd_smiles = Chem.MolToSmiles(rdmol)
    oemol2 = oemol_from_smiles(rd_smiles)
    round_trip_smiles = oechem.OEMolToSmiles(oemol2)

    # Fix for some cases where rotational symmetry means there are two possible canonical smiles
    oemol2_smiles = oechem.OEMolToSmiles(oemol2)
    oemol2 = oemol_from_smiles(oemol2_smiles)
    round_trip_smiles_2 = oechem.OEMolToSmiles(oemol2)

    if round_trip_smiles_2 != round_trip_smiles:
        assert orig_smiles in {round_trip_smiles, round_trip_smiles_2}
    else:
        assert orig_smiles == round_trip_smiles
