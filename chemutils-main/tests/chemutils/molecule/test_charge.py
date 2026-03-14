import pytest
from openeye import oechem
from rdkit import Chem

from chemutils.molecule import oemol_from_smiles, smiles_from_oemol
from chemutils.molecule.charge import remove_formal_charges
from chemutils.molecule.neutralize import neutralize_molecule

CASES = [
    ("[H+]", "[H+]"),
    ("[Li+]", "[Li+]"),
    ("[Na+]", "[Na+]"),
    ("[K+]", "[K+]"),
    ("[Be+2]", "[Be+2]"),
    ("[Mg+2]", "[Mg+2]"),
    ("[Ca+2]", "[Ca+2]"),
    ("[Al+3]", "[Al+3]"),
    ("[Zn+2]", "[Zn+2]"),
    ("[Fe+2]", "[Fe+2]"),
    ("[Fe+3]", "[Fe+3]"),
    ("[Cu+3]", "[Cu+3]"),
    ("[BH4+]", "[BH4+]"),
    ("B", "B"),
    ("[BH2-]", "[BH2-]"),
    ("[CH5+]", "C"),
    ("C", "C"),
    ("[CH3-]", "C"),
    ("[SiH5+]", "[SiH4]"),
    ("[SiH4]", "[SiH4]"),
    ("[SiH3-]", "[SiH4]"),
    ("[GeH5+]", "[GeH5+]"),
    ("[GeH4]", "[GeH4]"),
    ("[GeH3-]", "[GeH3-]"),
    ("[NH4+]", "N"),
    ("N", "N"),
    ("[NH2-]", "N"),
    ("[NH-2]", "N"),
    ("[N-3]", "N"),
    ("[PH4+]", "P"),
    ("P", "P"),
    ("[PH2-]", "P"),
    ("[PH-2]", "P"),
    ("[P-3]", "P"),
    ("[AsH4+][AsH3]", "[AsH3][AsH3]"),
    ("[AsH2-]", "[AsH3]"),
    ("[AsH-2]", "[AsH3]"),
    ("[As-3]", "[AsH3]"),
    ("[SbH4+]", "[SbH4+]"),
    ("[SbH3]", "[SbH3]"),
    ("[SbH2-]", "[SbH2-]"),
    ("[SbH-2]", "[SbH-2]"),
    ("[Sb-3]", "[Sb-3]"),
    ("[OH3+]", "O"),
    ("O", "O"),
    ("[OH-]", "O"),
    ("[O-2]", "O"),
    ("[SH3+]", "S"),
    ("S", "S"),
    ("[SH-]", "S"),
    ("[S-2]", "S"),
    ("[SeH3+]", "[SeH2]"),
    ("[SeH2]", "[SeH2]"),
    ("[SeH-]", "[SeH2]"),
    ("[Se-2]", "[SeH2]"),
    ("[TeH3+]", "[TeH3+]"),
    ("[TeH2]", "[TeH2]"),
    ("[TeH-]", "[TeH-]"),
    ("[Te-2]", "[Te-2]"),
    ("[FH2+]", "F"),
    ("F", "F"),
    ("[F-]", "F"),
    ("[ClH2+]", "Cl"),
    ("Cl", "Cl"),
    ("[Cl-]", "Cl"),
    ("[IH2+]", "I"),
    ("I", "I"),
    ("[I-]", "I"),
    ("[BrH2+]", "Br"),
    ("Br", "Br"),
    ("[Br-]", "Br"),
    ("[AtH2+]", "[AtH2+]"),
    ("[At]", "[At]"),
    ("[At-]", "[At-]"),
    ("[HeH+]", "[HeH+]"),
    ("[He]", "[He]"),
    # Multi-atom ions
    ("C[NH3+]", "CN"),
    ("C(=O)[O-]", "C(=O)O"),
    ("c1cc[nH+]cc1", "c1ccncc1"),
    ("[N+](C)(C)(C)C", "C[N+](C)(C)C"),
    ("O=N[O-]", "N(=O)O"),
    ("O=N([O-])[O-]", "N(=O)(O)O"),
    ("O=S([O-])[O-]", "OS(=O)O"),
    ("O=S(=O)([O-])[O-]", "OS(=O)(=O)O"),
    ("O=P([O-])([O-])[O-]", "OP(=O)(O)O"),
    ("P([O-])([O-])[O-]", "OP(O)O"),
    # Clearly wrong valence, but simple rules still apply
    ("[NH7+]", "[NH6]"),
    ("[NH6+]", "[NH5]"),
    ("[NH-]", "[NH2]"),
    ("[CH3+2]", "[CH]"),
    ("[O-]", "[OH]"),
    # If we can't remove enough hydrogens to get to 0 charge, don't do anything
    ("[OH+2]", "[OH+2]"),
    # Stereochemistry cases
    ("C[N@@H+](N)O", "C[N@@](N)O"),
    ("[C@@-](F)(Cl)Br", "[C@H](F)(Cl)Br"),
    ("[C@@H](F)(Cl)Br", "[C@@H](F)(Cl)Br"),
    ("[C@@](F)(Cl)(Br)I", "[C@@](F)(Cl)(Br)I"),
    ("[C@@H+](F)(Cl)(Br)I", "C(F)(Cl)(Br)I"),
    # Isotope behaviour
    ("[N+](H)(H)(H)H", "N"),
    ("[N+](H)([2H])(H)H", "[2H]N"),
    # Deletes the first hydrogen, even though it's an isotope
    ("[N+]([2H])(H)(H)H", "N"),
    # Behaviour when it can't fix a single atom
    ("[OH-].[Ca+]", "[OH-].[Ca+]"),
    ("[O-]C[N+](C)(C)C", "C[N+](C)(C)CO"),
]


@pytest.mark.parametrize(["smiles", "expected"], CASES)
def test_remove_formal_charges(smiles, expected):
    oemol = oemol_from_smiles(smiles)
    remove_formal_charges(oemol)
    assert oechem.OEMolToSmiles(oemol) == expected


@pytest.mark.parametrize(["smiles", "expected"], CASES)
def test_remove_formal_charges_explicit_H(smiles, expected):
    oemol = oemol_from_smiles(smiles)
    oechem.OEAddExplicitHydrogens(oemol)
    remove_formal_charges(oemol)
    assert oechem.OEMolToSmiles(oemol) == expected


@pytest.mark.parametrize(
    ["smiles", "expected"],
    [
        # Nitro
        ("C[N+](=O)[O-]", "C[N+](=O)[O-]"),
        ("C[N+](=O)O", "C[N+](=O)[O-]"),
        ("C[N](=O)[O-]", "C[N+](=O)[O-]"),
        ("C[N](=O)O", "C[N+](=O)[O-]"),
        ("CN(=O)[O-]", "C[N+](=O)[O-]"),
        ("CN(=O)O", "C[N+](=O)[O-]"),
        # Nitrone
        ("C=[N+](O)C", "C[N+](=C)[O-]"),
        ("C=[N](O)C", "C[N+](=C)[O-]"),
        ("C=N(O)C", "C[N+](=C)[O-]"),
        # Aromatic N-Oxide
        ("c1cc[n+](cc1)O", "c1cc[n+](cc1)[O-]"),
        ("c1cc[n+](cc1)[O-]", "c1cc[n+](cc1)[O-]"),
        ("O[N]1=CC=CC=C1", "c1cc[n+](cc1)[O-]"),
        # Aromatic nitro
        ("O[N](C1=CC=CC=C1)=O", "c1ccc(cc1)[N+](=O)[O-]"),
        # Azido
        ("CN=[N+]=[N-]", "CN=[N+]=[N-]"),
        ("CN=[N+]=N", "CN=[N+]=[N-]"),
        # Isonitrile
        ("C[N+]#[C-]", "C[N+]#[C-]"),
        ("C[N+]#C", "C[N+]#[C-]"),
        # Ammonium
        ("[NH4+]", "N"),
        # Tertiary ammonia
        ("C[NH+](C)C", "CN(C)C"),
        # Quaternary ammonia
        ("C[N+](C)(C)C", "C[N+](C)(C)C"),
        # Not limited to copying OpenEye! Can actually neutralize properly.
        ("[OH-].[Ca+]", "O.[Ca+]"),
        # Multi-atom ions
        ("C[NH3+]", "CN"),
        ("C(=O)[O-]", "C(=O)O"),
        ("c1cc[nH+]cc1", "c1ccncc1"),
        ("[N+](C)(C)(C)C", "C[N+](C)(C)C"),
        ("O=N[O-]", "N(=O)O"),
        ("O=N([O-])[O-]", "[N+](=O)(O)[O-]"),
        ("O=S([O-])[O-]", "OS(=O)O"),
        ("O=S(=O)([O-])[O-]", "OS(=O)(=O)O"),
        ("O=P([O-])([O-])[O-]", "OP(=O)(O)O"),
        ("P([O-])([O-])[O-]", "OP(O)O"),
        # Fixes wacky valencies
        ("[NH7+]", "N"),
        ("[NH6+]", "N"),
        ("[NH-]", "N"),
        ("[CH3+2]", "C"),
        ("[O-]", "O"),
        ("[OH+2]", "O"),
        # Other cases
        ("B(C)(C)C", "B(C)(C)C"),
        ("[B](C)(C)(C)C", "[B-](C)(C)(C)C"),
        ("FCl(F)F", "FCl(F)F"),
        ("OCl(=O)(=O)=O", "OCl(=O)(=O)=O"),
        ("O=Br(=O)F", "O=Br(=O)F"),
    ],
)
def test_neutralize_molecule(smiles, expected):
    oemol = oemol_from_smiles(smiles)
    neutralize_molecule(oemol)
    assert smiles_from_oemol(oemol) == expected


def test_neutralize_molecule_helps_rdkit():
    # Common smiles for nitro, with tetravalent nitrogen
    # RDKit can't handle this
    smiles = "CN(=O)O"

    # Failed to read
    mol = Chem.MolFromSmiles(smiles)
    assert mol is None

    # Fix using neutralize_molecule
    oemol = oemol_from_smiles(smiles)
    neutralize_molecule(oemol)
    fixed_smiles = smiles_from_oemol(oemol)

    assert fixed_smiles == "C[N+](=O)[O-]"

    # Now fixed
    mol = Chem.MolFromSmiles(fixed_smiles)
    assert mol is not None
