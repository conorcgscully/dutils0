import numpy as np
import pytest
from openeye import oechem

from chemutils.molecule.embed import embed_molecule_3d
from chemutils.molecule.enhanced_stereo import reflect_stereocenter
from chemutils.molecule.equal import assert_mol_equal
from chemutils.molecule.smiles import oemol_from_smiles, smiles_from_oemol, std_oemol_from_smiles
from chemutils.molecule.standardise import standardise_oemol

STANDARDISE_CASES = [
    ("C", 1, "C"),
    ("[CH4]", 1, "C"),
    ("C(H)(H)(H)H", 1, "C"),
    ("C([2H])([2H])([2H])[2H]", 1, "C"),
    ("CC.C", 2, "CC"),
    ("[NH4+]", 1, "N"),
    ("N[C@@H](C)C(=O)O", 6, "C[C@@H](C(=O)O)N"),
]


def oemol_from_smiles_then_standardise(smiles):
    mol = oemol_from_smiles(smiles)
    standardise_oemol(mol)
    return mol


SMILES_TO_MOL_FUNCS = [std_oemol_from_smiles, oemol_from_smiles_then_standardise]


@pytest.mark.parametrize("smiles_to_oemol", SMILES_TO_MOL_FUNCS)
@pytest.mark.parametrize(["smiles", "num_atoms", "expected_smiles"], STANDARDISE_CASES)
def test_standardise_oemol(smiles_to_oemol, smiles, num_atoms, expected_smiles):
    mol = smiles_to_oemol(smiles)
    assert oechem.OEMolToSmiles(mol) == expected_smiles
    assert mol.NumAtoms() == num_atoms
    assert not any(atom.GetAtomicNum() == 1 for atom in mol.GetAtoms())
    assert np.array(list(mol.GetCoords().values())).shape == (num_atoms, 3)


@pytest.mark.parametrize(
    "group",
    [
        # Atom ordering
        ("CONS", "SNOC", "O(C)NS", "N(OC)S"),
        # Atom ordering and bond orders
        ("O=PC#N", "P(=O)C#N", "N#CP=O", "C(P=O)#N"),
        # Kekulization and aromaticity tests (benzene)
        ("c1ccccc1", "C1=CC=CC=C1", "C1C=CC=CC=1"),
        # Tetrahedral isomerisation
        (
            "F[C@](Cl)(I)Br",
            "F[C@](Br)(Cl)I",
            "F[C@@](Cl)(Br)I",
            "Cl[C@](F)(Br)I",
        ),
        # Charge removal and ordering (gylcine)
        ("C(C(=O)O)N", "C(C(=O)[O-])[NH3+]", "NCC(O)=O", "[NH3+]CC(=O)O", "[O-]C(=O)CN"),
    ],
)
def test_standardise_oemol_isomeric(group):
    mols = []
    for smiles in group:
        mol = oemol_from_smiles(smiles)
        standardise_oemol(mol)
        mols.append(mol)
    for i in range(1, len(mols)):
        a = mols[0]
        b = mols[i]
        try:
            assert_mol_equal(a, b)
        except AssertionError as e:
            raise AssertionError(
                f"Standardised molecule mismatch: {group[0]} does not match {group[i]}"
            ) from e


@pytest.mark.parametrize(
    ["smiles", "hybridisations"],
    [
        ("C", [oechem.OEHybridization_sp3]),
        ("C=C", [oechem.OEHybridization_sp2, oechem.OEHybridization_sp2]),
        ("C#C", [oechem.OEHybridization_sp, oechem.OEHybridization_sp]),
    ],
)
def test_standardise_oemol_hybridisation(smiles, hybridisations):
    oemol = oemol_from_smiles(smiles)
    standardise_oemol(oemol)
    assert [atom.GetHyb() for atom in oemol.GetAtoms()] == hybridisations


@pytest.mark.parametrize("smiles", ["N(O)(Cl)C", "CCCN(C)CC"])
def test_standardise_oemol_chirality(smiles):
    oemol = std_oemol_from_smiles(smiles)
    oemol2 = oemol.CreateCopy()
    oechem.OEClearChiralPerception(oemol2)
    oechem.OEPerceiveChiral(oemol2, carbonOnly=True)
    with pytest.raises(AssertionError):
        assert_mol_equal(oemol, oemol2)
    oemol3 = oemol2.CreateCopy()
    standardise_oemol(oemol3)
    assert_mol_equal(oemol, oemol3)


def test_standardise_fixes_chirality():
    # Regression test for standardise_oemol not correcting chirality when centers are flipped without
    # changing the relevant coordinates

    smiles = "[C@H](F)(Cl)Br"
    mol = oemol_from_smiles(smiles)
    embed_molecule_3d(mol, seed=0)
    assert smiles_from_oemol(mol) == "[C@H](F)(Cl)Br"
    standardise_oemol(mol)
    assert smiles_from_oemol(mol) == "[C@H](F)(Cl)Br"
    reflect_stereocenter(next(iter(mol.GetAtoms())))
    assert smiles_from_oemol(mol) == "[C@@H](F)(Cl)Br"
    standardise_oemol(mol)
    # Standardise should flip back based to the label based on coordinates
    assert smiles_from_oemol(mol) == "[C@H](F)(Cl)Br"
