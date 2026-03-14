import pytest
from openeye import oechem

from chemutils.molecule import (
    create_atom,
    get_subset,
    oemol_from_smiles,
    smiles_from_oemol,
    subset_molecule,
)


@pytest.mark.parametrize(
    ["smiles", "selection", "bond_handling", "include_internal_bonds", "subset_smiles"],
    [
        ("CCO", oechem.OEHasAtomicNum(6), "dangling", True, "C[CH2]"),
        ("CCO", oechem.OEHasAtomicNum(6), "hydrogens", True, "CC"),
        ("CCO", oechem.OEHasAtomicNum(6), "rgroups", True, "[R1]CC"),
        ("CCO", oechem.OEHasAtomicNum(6), "dangling", False, "[CH3].[CH2]"),
        ("CCO", oechem.OEHasAtomicNum(6), "hydrogens", False, "C.C"),
        ("CCO", oechem.OEHasAtomicNum(6), "rgroups", False, "[R1]C.[R2]C[R3]"),
        ("COC", oechem.OEHasAtomicNum(6), "dangling", True, "[CH3].[CH3]"),
        ("COC", oechem.OEHasAtomicNum(6), "hydrogens", True, "C.C"),
        ("COC", oechem.OEHasAtomicNum(6), "rgroups", True, "[R1]C.[R2]C"),
        ("CCO", lambda atom: atom.GetAtomicNum() == 6, "hydrogens", True, "CC"),
    ],
)
def test_subset_molecule(smiles, subset_smiles, bond_handling, include_internal_bonds, selection):
    mol = oemol_from_smiles(smiles)
    subset_mol = subset_molecule(
        mol,
        selection=selection,
        bond_handling=bond_handling,
        include_internal_bonds=include_internal_bonds,
    )
    assert smiles_from_oemol(subset_mol) == subset_smiles


def test_subset_molecule_atombondset():
    mol = oemol_from_smiles("CCO")
    atombondset = oechem.OEAtomBondSet()
    atombondset.AddAtoms(list(mol.GetAtoms())[:2])

    subset_mol = subset_molecule(mol, selection=atombondset, bond_handling="hydrogens")
    assert smiles_from_oemol(subset_mol) == "CC"

    subset_mol = subset_molecule(
        mol, selection=atombondset, bond_handling="hydrogens", include_internal_bonds=False
    )
    assert smiles_from_oemol(subset_mol) == "C.C"


def test_subset_molecule_atom_list():
    mol = oemol_from_smiles("CCO")
    atoms = list(mol.GetAtoms())[:2]

    subset_mol = subset_molecule(mol, selection=atoms, bond_handling="hydrogens")
    assert smiles_from_oemol(subset_mol) == "CC"

    subset_mol = subset_molecule(
        mol, selection=atoms, bond_handling="hydrogens", include_internal_bonds=False
    )
    assert smiles_from_oemol(subset_mol) == "C.C"


def test_subset_molecule_oeatomiter():
    mol = oemol_from_smiles("CCO")

    subset_mol = subset_molecule(mol, selection=mol.GetAtoms(), bond_handling="hydrogens")
    assert smiles_from_oemol(subset_mol) == "CCO"


def test_subset_molecule_index():
    mol = oemol_from_smiles("CCO")

    subset = subset_molecule(mol, selection=0, bond_handling="hydrogens")
    assert smiles_from_oemol(subset) == "C"

    subset = subset_molecule(mol, selection=[0, 1], bond_handling="hydrogens")
    assert smiles_from_oemol(subset) == "CC"

    subset = subset_molecule(mol, selection={0, 2}, bond_handling="hydrogens")
    assert smiles_from_oemol(subset) == "C.O"


@pytest.mark.parametrize(
    ["smiles", "selection", "indices"],
    [
        # `same molecule as`
        ("CC.CCC.C", "same molecule as index 0", {0, 1}),
        ("CC.CCC.C", "same molecule as index 2", {2, 3, 4}),
        ("CC.CCC.C", "not same molecule as index 0", {2, 3, 4, 5}),
        # `same ring_system as`
        ("C1CCCCC1-c1c2c(nc[nH]2)ncn1", "same ring_system as index 0", {0, 1, 2, 3, 4, 5}),
        (
            "C1CCCCC1-c1c2c(nc[nH]2)ncn1",
            "same ring_system as elem N",
            {6, 7, 8, 9, 10, 11, 12, 13, 14},
        ),
        (
            "CCC",
            "same ring_system as index 0",
            {},
        ),  # Check that non-ring atoms aren't treated as ring systems
        ("c1cccnc1-CN", "same ring_system as elem N", {0, 1, 2, 3, 4, 5}),
        # `same element as`
        ("CCNOCN", "same element as index 0", {0, 1, 4}),
        ("CCNOCN", "not same element as index 2", {0, 1, 3, 4}),
        # `same formal_charge as`
        ("[NH3+]C(=O)[O-]", "same formal_charge as index 0", {0}),
        ("[NH3+]C(=O)[O-]", "not same formal_charge as index 0", {1, 2, 3}),
        # `same degree as`
        ("CNOS", "same degree as index 0", {0}),
        ("CNOS", "same degree as index 2", {2, 3}),
        # `same hvy_degree as`
        ("CC(C)C(C)C", "same hvy_degree as index 0", {0, 2, 4, 5}),
        ("CC(C)C(C)C", "same hvy_degree as index 1", {1, 3}),
        # `same valence as`
        ("CC[SiH3]", "same valence as index 0", {0, 1, 2}),
        # `same hcount as`
        ("CC=CC#C", "same hcount as index 0", {0}),
        ("CC=CC#C", "same hcount as index 1", {1, 2, 4}),
        ("CC=CC#C", "same hcount as index 3", {3}),
        # `all`
        ("CNOS", "all", {0, 1, 2, 3}),
        # `none`
        ("CNOS", "none", set()),
        # Boolean Operators
        # Simple `or`
        ("CNOS", "elem N or elem O", {1, 2}),
        # Multiple `or`
        ("CNOS", "elem N or elem O or elem C", {0, 1, 2}),
        # Simple `and`
        ("CNO[S-]", "degree 2 and formal_charge 0", {2}),
        ("CNOS", "degree 2 and elem N", set()),
        # Precedence
        # `and` is higher priority than `or`
        ("CNOS", "elem C or degree 2 and elem O", {0, 2}),
        # Use parenthesis to make `or` beat `and`
        ("CNOS", "(elem C or degree 2) and elem O", {2}),
        # `not` is higher priority than `or`
        ("CNOS", "not elem C or elem O", {1, 2, 3}),
        # Use parenthesis to make `or` beat `not`
        ("CNOS", "not (elem C or elem O)", {1, 3}),
        # `bonded_to`
        ("CNOS", "bonded to N", {0, 2}),
        ("CNOS", "not bonded to N", {1, 3}),
        ("CNOS", "bonded to index 2", {1, 3}),
        ("CNOS", "bonded to not N", {1}),
        # `index`
        ("CNOS", "index 0", {0}),
        ("CNOS", "index 2", {2}),
        ("CNOS", "index 0 2 3", {0, 2, 3}),
        ("CNOS", "index 1 or index 2", {1, 2}),
        ("CNOS", "index 1:3", {1, 2, 3}),
        ("CNOS", "index 1-3", {1, 2, 3}),
        # `elem`
        ("CCO", "elem C", {0, 1}),
        ("CCO", "elem O", {2}),
        ("CCO", "elem C or elem O", {0, 1, 2}),
        ("CNOS", "elem C O", {0, 2}),
        # `formal_charge`
        ("[NH3+]C(=O)[O-]", "formal_charge 0", {1, 2}),
        ("[NH3+]C(=O)[O-]", "formal_charge -1", {3}),
        ("[NH3+]C(=O)[O-]", "formal_charge 1", {0}),
        ("[NH3+]C(=O)[O-]", "formal_charge +1", {0}),
        # `degree`
        ("CNOS", "degree 4", {0}),
        ("CNOS", "degree 3 2", {1, 2, 3}),
        # `hvy_degree`
        ("CC(C)CC(C)(C)C", "hvy_degree 1", {0, 2, 5, 6, 7}),
        ("CC(C)CC(C)(C)C", "hvy_degree > 2", {1, 4}),
        # `valence`
        ("CNOS", "valence 4", {0}),
        ("CNOS", "valence 2 4", {0, 2, 3}),
        # `hcount`
        ("CNOS", "hcount 1", {1, 3}),
        ("CNOS", "hcount > 0", {0, 1, 3}),
        # `molecule {smarts}` selects full molecules with a given SMARTS
        ("CO.CCO", "CO", {0, 1, 3, 4}),
        ("CO.CCO", "molecule CO", {0, 1}),
        ("CO.CCO", "CO and not molecule CO", {3, 4}),  # CO that isn't carbon monoxide
        # `ring_size`
        ("C1CC1-C1CCC1-C1CCCC1", "ring_size 3", {0, 1, 2}),
        ("C1CC1-C1CCC1-C1CCCC1", "ring_size 3 5", {0, 1, 2, 7, 8, 9, 10, 11}),
        ("C1CC1-C1CCC1-C1CCCC1", "ring_size 3-4", {0, 1, 2, 3, 4, 5, 6}),
        # `aromatic`
        ("c1ccccc1-C1CCCCC1", "aromatic", {0, 1, 2, 3, 4, 5}),
        ("c1ccccc1-C1CCCCC1", "not aromatic", {6, 7, 8, 9, 10, 11}),
        # `aliphatic`
        ("c1ccccc1-C1CCCCC1", "aliphatic", {6, 7, 8, 9, 10, 11}),
        ("c1ccccc1-C1CCCCC1", "not aliphatic", {0, 1, 2, 3, 4, 5}),
        # `halogen`
        ("CNF", "halogen", {2}),
        ("CNF", "not halogen", {0, 1}),
        ("C(Cl)(I)(F)Br", "halogen", {1, 2, 3, 4}),
        # `metal`
        ("[Fe+2]", "metal", {0}),
        ("C", "metal", set()),
        # `ring`
        ("C1CCC1CC", "ring", {0, 1, 2, 3}),
        # SMARTS
        ("CNOS", "C", {0}),
        ("c1ccccc1-C1CCCCC1", "c", {0, 1, 2, 3, 4, 5}),
        ("c1ccccc1-C1CCCCC1", "a", {0, 1, 2, 3, 4, 5}),
        ("c1ccccc1-C1CCCCC1", "C", {6, 7, 8, 9, 10, 11}),
        ("c1ccccc1-C1CCCCC1", "A", {6, 7, 8, 9, 10, 11}),
        ("c1ccccc1-C1CCCCC1", "[#6]", {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11}),
        ("CC(=O)O", "C=O", {1, 2}),
        ("CC(=O)O", "degree 3 and C=O", {1}),
        # SMARTS with branches
        ("CCC(CC)CC", "CC(C)C", {1, 2, 3, 5}),
        # SMARTS with charges
        ("[Fe+2]", "[Fe+2]", {0}),
        ("[Fe+3]", "[Fe+2]", set()),
        # SMARTS with stereo (note that hydrogens are not included)
        ("[C@H](F)(Cl)Br", "[C@H](F)(Cl)Br", {0, 2, 3, 4}),
        # SMARTS with conditions for bonds
        ("C=CCC#C", "C=,#C", {0, 1, 3, 4}),
        ("C=CCC#C", "C!=C", {1, 2, 3, 4}),
    ],
)
def test_selection_language(smiles, selection, indices):
    mol = oemol_from_smiles(smiles)
    subset = get_subset(mol, selection=selection)
    assert {atom.GetIdx() for atom in subset.GetAtoms()} == set(indices)


@pytest.mark.parametrize(
    ["selection", "indices"],
    [
        ("name CA", {0, 1}),
        ('name "O2\'"', {2}),
        ("name", {0, 1, 2}),
        ("not name", {3}),
        ("same name as index 0", {0, 1}),
    ],
)
def test_selection_language_atom_name(selection, indices):
    mol = oechem.OEGraphMol()
    create_atom(mol, atomic_number=6, name="CA")
    create_atom(mol, atomic_number=6, name=" CA")
    create_atom(mol, atomic_number=8, name="O2'")
    create_atom(mol, atomic_number=6)

    assert {atom.GetIdx() for atom in get_subset(mol, selection=selection).GetAtoms()} == indices


@pytest.mark.parametrize(
    ["selection", "indices"],
    [
        ("chain A", {0, 1}),
        ("chain B", {2}),
        ("chain A B", {0, 1, 2}),
        ("chain", {0, 1, 2}),
        ("not chain", {3}),
    ],
)
def test_selection_language_chain(selection, indices):
    mol = oechem.OEGraphMol()
    create_atom(mol, atomic_number=6, chain_id="A")
    create_atom(mol, atomic_number=6, chain_id="A ")
    create_atom(mol, atomic_number=6, chain_id="B")
    create_atom(mol, atomic_number=6, chain_id="")

    assert {atom.GetIdx() for atom in get_subset(mol, selection=selection).GetAtoms()} == indices


@pytest.mark.parametrize(
    ["selection", "indices"],
    [
        ("resname ALA", {0}),
        ("resname LYS", {1, 2}),
        ("resname 'LYS'", {1}),
        ('resname " LYS"', {2}),
        ("resname ' LYS'", {2}),
        ("resname ALA LYS", {0, 1, 2}),
        ("resname", {0, 1, 2}),
        ("not resname", {3}),
        ("same resname as index 0", {0}),
        ("same resname as index 1", {1, 2}),
    ],
)
def test_selection_language_residue_name(selection, indices):
    mol = oechem.OEGraphMol()
    create_atom(mol, atomic_number=6, residue_name="ALA")
    create_atom(mol, atomic_number=6, residue_name="LYS")
    create_atom(mol, atomic_number=6, residue_name=" LYS")
    create_atom(mol, atomic_number=6, residue_name="")

    assert {atom.GetIdx() for atom in get_subset(mol, selection=selection).GetAtoms()} == indices


@pytest.mark.parametrize(
    ["selection", "indices"],
    [
        ("resname 57V", {0}),
    ],
)
def test_selection_language_residue_name_starts_number(selection, indices):
    mol = oechem.OEGraphMol()
    create_atom(mol, atomic_number=6, residue_name="57V")

    assert {atom.GetIdx() for atom in get_subset(mol, selection=selection).GetAtoms()} == indices


@pytest.mark.parametrize(
    ["selection", "indices"],
    [
        ("resname 584", {0}),
    ],
)
def test_selection_language_residue_name_is_number(selection, indices):
    mol = oechem.OEGraphMol()
    create_atom(mol, atomic_number=6, residue_name="584")

    assert {atom.GetIdx() for atom in get_subset(mol, selection=selection).GetAtoms()} == indices


@pytest.mark.parametrize(
    ["selection", "indices"],
    [
        ("resnum -2", {0}),
        ("resnum 1", {1}),
        ("resnum 1 2", {1, 2, 3}),
        ("resnum = 1", {1}),
        ("resnum > 0", {1, 2, 3}),
        ("resnum < 0", {0}),
        ("resnum >= 1", {1, 2, 3}),
        ("resnum <= 1", {0, 1}),
        ("resnum != 1", {0, 2, 3}),
        ("same resnum as index 2", {2, 3}),
    ],
)
def test_selection_language_residue_number(selection, indices):
    mol = oechem.OEGraphMol()
    create_atom(mol, atomic_number=6, chain_id="A", residue_number=-2)
    create_atom(mol, atomic_number=6, chain_id="A", residue_number=1)
    create_atom(mol, atomic_number=6, chain_id="A", residue_number=2)
    create_atom(mol, atomic_number=6, chain_id="B", residue_number=2)

    assert {atom.GetIdx() for atom in get_subset(mol, selection=selection).GetAtoms()} == indices


@pytest.mark.parametrize(
    ["selection", "indices"],
    [
        ("altloc A", {0}),
        ("altloc B", {1, 2}),
        ("altloc A B", {0, 1, 2}),
        ("altloc", {0, 1, 2}),
        ("not altloc", {3, 4}),
        ("same altloc as index 0", {0}),
        ("same altloc as index 1", {1, 2}),
    ],
)
def test_selection_language_altloc(selection, indices):
    mol = oechem.OEGraphMol()
    create_atom(mol, atomic_number=6, alternate_location="A")
    create_atom(mol, atomic_number=6, alternate_location="B")
    create_atom(mol, atomic_number=6, alternate_location="B")
    create_atom(mol, atomic_number=6, alternate_location=" ")
    create_atom(mol, atomic_number=6, alternate_location=None)

    assert {atom.GetIdx() for atom in get_subset(mol, selection=selection).GetAtoms()} == indices


@pytest.mark.parametrize(
    ["selection", "indices"],
    [
        ("serial 1", {0}),
        ("serial 20", {1, 3}),
        ("serial 214", {2}),
        ("serial 1 20", {0, 1, 3}),
        ("serial > 15", {1, 2, 3}),
        ("serial 0:25", {0, 1, 3}),
    ],
)
def test_selection_language_atom_serial(selection, indices):
    mol = oechem.OEGraphMol()
    create_atom(mol, atomic_number=6, serial_number=1)
    create_atom(mol, atomic_number=6, serial_number=20)
    create_atom(mol, atomic_number=6, serial_number=214)
    create_atom(mol, atomic_number=6, serial_number=20)

    assert {atom.GetIdx() for atom in get_subset(mol, selection=selection).GetAtoms()} == indices


@pytest.mark.parametrize(
    ["selection", "indices"],
    [
        ("b_factor 56.3", {0}),
        ("b_factor 20.1", {1, 3}),
        ("b_factor 96.9", {2}),
        ("b_factor != 96.9", {0, 1, 3}),
        ("b_factor 56.3 20.1", {0, 1, 3}),
        ("b_factor > 30", {0, 2}),
        ("b_factor 0:60", {0, 1, 3}),
        ("same b_factor as index 1", {1, 3}),
        ("not same b_factor as index 1", {0, 2}),
    ],
)
def test_selection_language_b_factor(selection, indices):
    mol = oechem.OEGraphMol()
    create_atom(mol, atomic_number=6, b_factor=56.3)
    create_atom(mol, atomic_number=6, b_factor=20.1)
    create_atom(mol, atomic_number=6, b_factor=96.9)
    create_atom(mol, atomic_number=6, b_factor=20.1)

    assert {atom.GetIdx() for atom in get_subset(mol, selection=selection).GetAtoms()} == indices


@pytest.mark.parametrize(
    ["selection", "indices"],
    [
        ("occ 0.24", {0}),
        ("occ 0.64", {1, 3}),
        ("occ 1.00", {2}),
        ("occ != 1.00", {0, 1, 3}),
        ("occ 0.24 0.64", {0, 1, 3}),
        ("occ > 0.5", {1, 2, 3}),
        ("occ 0:0.7", {0, 1, 3}),
        ("same occ as index 1", {1, 3}),
    ],
)
def test_selection_language_occupancy(selection, indices):
    mol = oechem.OEGraphMol()
    create_atom(mol, atomic_number=6, occupancy=0.24)
    create_atom(mol, atomic_number=6, occupancy=0.64)
    create_atom(mol, atomic_number=6, occupancy=1.00)
    create_atom(mol, atomic_number=6, occupancy=0.64)

    assert {atom.GetIdx() for atom in get_subset(mol, selection=selection).GetAtoms()} == indices


@pytest.mark.parametrize(
    ["selection", "indices"],
    [
        ("same chain as index 0", {0, 1}),
        ("same chain as index 2", {2, 3}),
        ("same chain as index 4", {4}),
        ("not same chain as index 0", {2, 3, 4}),
    ],
)
def test_same_chain_as(selection, indices):
    mol = oechem.OEGraphMol()
    create_atom(mol, atomic_number=6, chain_id="A")
    create_atom(mol, atomic_number=6, chain_id="A")
    create_atom(mol, atomic_number=6, chain_id="B")
    create_atom(mol, atomic_number=6, chain_id="B")
    create_atom(mol, atomic_number=6, chain_id="C")

    assert {atom.GetIdx() for atom in get_subset(mol, selection=selection).GetAtoms()} == indices


@pytest.mark.parametrize(
    ["selection", "indices"],
    [
        ("same residue as index 0", {0, 1}),
        ("same residue as index 2", {2, 3}),
        ("same residue as index 4", {4}),
        ("same residue as index 5", {5}),
        ("not same residue as index 0", {2, 3, 4, 5}),
    ],
)
def test_same_residue_as(selection, indices):
    mol = oechem.OEGraphMol()
    create_atom(mol, atomic_number=6, chain_id="A", residue_name="ALA", residue_number=1)
    create_atom(mol, atomic_number=6, chain_id="A", residue_name="ALA", residue_number=1)
    create_atom(mol, atomic_number=6, chain_id="A", residue_name="LYS", residue_number=2)
    create_atom(mol, atomic_number=6, chain_id="A", residue_name="LYS", residue_number=2)
    create_atom(mol, atomic_number=6, chain_id="A", residue_name="ALA", residue_number=3)
    create_atom(mol, atomic_number=6, chain_id="B", residue_name="ALA", residue_number=1)

    assert {atom.GetIdx() for atom in get_subset(mol, selection=selection).GetAtoms()} == indices


@pytest.mark.parametrize(
    ["selection", "indices"],
    [
        ("hetatm", {0, 2}),
        ("not hetatm", {1}),
    ],
)
def test_hetatm(selection, indices):
    mol = oechem.OEGraphMol()
    create_atom(mol, atomic_number=6, is_hetatm=True)
    create_atom(mol, atomic_number=6)
    create_atom(mol, atomic_number=6, is_hetatm=True)

    assert {atom.GetIdx() for atom in get_subset(mol, selection=selection).GetAtoms()} == indices


@pytest.mark.parametrize(
    ["selection", "indices"],
    [("water", {1, 2}), ("not water", {0, 3})],
)
def test_water(selection, indices):
    mol = oechem.OEGraphMol()
    create_atom(mol, atomic_number=8)  # Oxygen by itself doesn't count
    create_atom(mol, atomic_number=8, residue_name="HOH")
    create_atom(mol, atomic_number=8, residue_name="TIP")
    create_atom(mol, atomic_number=6)

    assert {atom.GetIdx() for atom in get_subset(mol, selection=selection).GetAtoms()} == indices


@pytest.mark.parametrize(
    ["selection", "indices"],
    [("backbone", {0}), ("not backbone", {1})],
)
def test_backbone(selection, indices):
    mol = oechem.OEGraphMol()
    create_atom(mol, atomic_number=6, residue_name="ALA", name=" CA ")
    create_atom(mol, atomic_number=6, residue_name="ALA", name=" CB ")

    assert {atom.GetIdx() for atom in get_subset(mol, selection=selection).GetAtoms()} == indices


@pytest.mark.parametrize(
    ["selection", "indices"],
    [("alpha_carbon", {1}), ("not alpha_carbon", {0, 2})],
)
def test_alpha_carbon(selection, indices):
    mol = oechem.OEGraphMol()
    create_atom(mol, atomic_number=6, residue_name="ALA", name="CA")
    create_atom(mol, atomic_number=6, residue_name="ALA", name=" CA ")
    create_atom(mol, atomic_number=6, residue_name="ALA", name="CB")

    assert {atom.GetIdx() for atom in get_subset(mol, selection=selection).GetAtoms()} == indices


@pytest.mark.parametrize(
    ["selection", "indices"],
    [
        ("within 1.5 of index 0", {0, 1}),
        ("around 1.5 of index 0", {1}),
        ("beyond 1.5 of index 0", {2, 3, 4, 5}),
        ("within 2 of index 2", {1, 2, 3, 4}),
        ("beyond 1.3 of index 2", {0, 4, 5}),
        ("around 2.7 of index 2", {0, 1, 3, 4, 5}),
    ],
)
def test_distances(selection, indices):
    mol = oechem.OEGraphMol()
    mol.SetDimension(3)
    create_atom(mol, atomic_number=6, coordinates=[0, 0, 0])
    create_atom(mol, atomic_number=6, coordinates=[1, 0, 0])
    create_atom(mol, atomic_number=6, coordinates=[2, 0, 0])
    create_atom(mol, atomic_number=6, coordinates=[3, 0, 0])
    create_atom(mol, atomic_number=6, coordinates=[3, 1, 0])
    create_atom(mol, atomic_number=6, coordinates=[3, 2, 0])

    assert {atom.GetIdx() for atom in get_subset(mol, selection=selection).GetAtoms()} == indices


def test_atom_interpolation():
    oemol = oemol_from_smiles("CCO.CN")
    atoms = list(oemol.GetAtoms())

    assert {atom.GetIdx() for atom in get_subset(oemol, selection=f"{atoms[0]}").GetAtoms()} == {0}
    assert {
        atom.GetIdx() for atom in get_subset(oemol, selection=f"bonded to {atoms[1]}").GetAtoms()
    } == {0, 2}
    assert {
        atom.GetIdx()
        for atom in get_subset(oemol, selection=f"same molecule as {atoms[3]}").GetAtoms()
    } == {3, 4}


def test_atombondset_interpolation():
    oemol = oemol_from_smiles("CCO.CN")

    sel1 = get_subset(oemol, selection="elem O")
    assert {atom.GetIdx() for atom in get_subset(oemol, selection=f"{sel1}").GetAtoms()} == {2}
    assert {
        atom.GetIdx() for atom in get_subset(oemol, selection=f"bonded to {sel1}").GetAtoms()
    } == {1}
    assert {
        atom.GetIdx()
        for atom in get_subset(oemol, selection=f"not same molecule as {sel1}").GetAtoms()
    } == {3, 4}


@pytest.mark.parametrize(
    ["smiles", "selection", "bond_handling", "expected_atom_mapping", "expected_bond_mapping"],
    [
        ("CCC", "index 0", "hydrogens", {0: 0}, {}),
        ("CCC", "index <= 1", "hydrogens", {0: 0, 1: 1}, {0: 0}),
        ("CCC", "index >= 1", "hydrogens", {1: 0, 2: 1}, {1: 0}),
    ],
)
def test_mappings(smiles, selection, bond_handling, expected_atom_mapping, expected_bond_mapping):
    oemol = oemol_from_smiles(smiles)

    atom_mapping = {}
    bond_mapping = {}
    # Must keep reference to subset molecule, or else atoms will segfault
    _ = subset_molecule(
        oemol,
        selection=selection,
        bond_handling=bond_handling,
        atom_mapping=atom_mapping,
        bond_mapping=bond_mapping,
    )
    assert {k.GetIdx(): v.GetIdx() for k, v in atom_mapping.items()} == expected_atom_mapping
    assert {k.GetIdx(): v.GetIdx() for k, v in bond_mapping.items()} == expected_bond_mapping
