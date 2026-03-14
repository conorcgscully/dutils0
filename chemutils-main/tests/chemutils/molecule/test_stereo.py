from dataclasses import dataclass

import numpy as np
import pytest
from openeye import oechem
from rdkit import Chem
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers

from chemutils.molecule import canonicalize_smiles, cxsmiles_from_oemol, oemol_from_smiles
from chemutils.molecule.coordinates import get_coordinates
from chemutils.molecule.graph import get_adjacency_matrix
from chemutils.molecule.mask import get_atom_mask
from chemutils.molecule.standardise import standardise_oemol
from chemutils.molecule.stereo import (
    UnderspecifiedStereoError,
    are_diastereomers,
    are_enantiomers,
    are_mirror_images,
    are_stereo_related,
    are_stereoisomers,
    check_tetrahedral_coordinates,
    enumerate_stereoisomers,
    get_all_tetrahedral_stereo_from_coordinates,
    get_atom_tetrahedral_stereo,
    get_bond_cistrans_stereo,
    is_chiral_molecule,
    is_meso_isomer,
    is_specific_stereoisomer,
    reflect_molecule_stereo,
    reflect_tetrahedral_stereo,
    set_atom_tetrahedral_stereo,
    set_bond_cistrans_stereo,
)

from ..cip_validation import CIP_VALIDATION_OEMOLS_3D


@dataclass
class ReflectedStereoCase:
    smiles: str
    smiles_reflected: str
    is_stereoisomer: bool
    is_mirror: bool
    is_chiral: bool
    is_enantiomer: bool
    is_meso: bool


# smiles, reflected_smiles, is_chiral_molecule, is_enantiomer, is_meso
REFLECTED_CASES = [
    ReflectedStereoCase("[C@](F)(Cl)(Br)I", "[C@@](F)(Cl)(Br)I", True, True, True, True, False),
    ReflectedStereoCase(
        "[C](F)(Cl)(Br)I",
        "[C](F)(Cl)(Br)I",
        UnderspecifiedStereoError,
        UnderspecifiedStereoError,
        UnderspecifiedStereoError,
        UnderspecifiedStereoError,
        False,
    ),
    # Tartaric acid
    ReflectedStereoCase(
        "[C@@H]([C@H](C(=O)O)O)(C(=O)O)O",
        "[C@H]([C@@H](C(=O)O)O)(C(=O)O)O",
        True,
        True,
        True,
        True,
        False,
    ),
    ReflectedStereoCase(
        "[C@H]([C@@H](C(=O)O)O)(C(=O)O)O",
        "[C@@H]([C@H](C(=O)O)O)(C(=O)O)O",
        True,
        True,
        True,
        True,
        False,
    ),
    ReflectedStereoCase(
        "[C@@H]([C@@H](C(=O)O)O)(C(=O)O)O",
        "[C@@H]([C@@H](C(=O)O)O)(C(=O)O)O",
        False,
        True,
        False,
        False,
        True,
    ),
    # 1,2-dimethylcyclopropane
    ReflectedStereoCase("C[C@H]1C[C@@H]1C", "C[C@@H]1C[C@H]1C", True, True, True, True, False),
    ReflectedStereoCase("C[C@@H]1C[C@H]1C", "C[C@H]1C[C@@H]1C", True, True, True, True, False),
    ReflectedStereoCase("C[C@@H]1C[C@@H]1C", "C[C@@H]1C[C@@H]1C", False, True, False, False, True),
]


@dataclass
class StereoPairCase:
    smiles1: str
    smiles2: str
    is_stereo_related: bool
    is_stereoisomer: bool
    is_mirror_image: bool
    is_enantiomer: bool
    is_diastereomer: bool


PAIR_CASES = [
    StereoPairCase("[C@](F)(Cl)(Br)I", "[C@@](F)(Cl)(Br)I", True, True, True, True, False),
    StereoPairCase("[C@](F)(Cl)(Br)I", "[C@](F)(Cl)(Br)I", True, False, False, False, False),
    StereoPairCase(
        "[C](F)(Cl)(Br)I",
        "[C@@](F)(Cl)(Br)I",
        True,
        UnderspecifiedStereoError,
        UnderspecifiedStereoError,
        UnderspecifiedStereoError,
        UnderspecifiedStereoError,
    ),
    StereoPairCase("[C@](F)(Cl)(Br)I", "[C@](N)(Cl)(Br)I", False, False, False, False, False),
    # Tartaric acid
    StereoPairCase(
        "[C@@H]([C@H](C(=O)O)O)(C(=O)O)O",
        "[C@H]([C@@H](C(=O)O)O)(C(=O)O)O",
        True,
        True,
        True,
        True,
        False,
    ),
    StereoPairCase(
        "[C@H]([CH](C(=O)O)O)(C(=O)O)O",
        "[CH]([C@H](C(=O)O)O)(C(=O)O)O",
        True,
        UnderspecifiedStereoError,
        UnderspecifiedStereoError,
        UnderspecifiedStereoError,
        UnderspecifiedStereoError,
    ),
    StereoPairCase(
        "[CH]([C@@H](C(=O)O)O)(C(=O)O)O",
        "[C@@H]([C@@H](C(=O)O)O)(C(=O)O)O",
        True,
        UnderspecifiedStereoError,
        UnderspecifiedStereoError,
        UnderspecifiedStereoError,
        UnderspecifiedStereoError,
    ),
]


@pytest.mark.parametrize(
    ["smiles1", "smiles2", "is_stereo"],
    [
        *[(case.smiles1, case.smiles2, case.is_stereo_related) for case in PAIR_CASES],
        *[(case.smiles, case.smiles_reflected, True) for case in REFLECTED_CASES],
    ],
)
def test_are_stereo_related(smiles1, smiles2, is_stereo):
    oemol1 = oemol_from_smiles(smiles1)
    oemol2 = oemol_from_smiles(smiles2)
    assert are_stereo_related(oemol1, oemol2) == is_stereo
    assert are_stereo_related(oemol2, oemol1) == is_stereo


@pytest.mark.parametrize(
    ["smiles1", "smiles2", "is_stereo"],
    [
        *[(case.smiles1, case.smiles2, case.is_stereoisomer) for case in PAIR_CASES],
        *[(case.smiles, case.smiles_reflected, case.is_stereoisomer) for case in REFLECTED_CASES],
    ],
)
def test_are_stereoisomers(smiles1, smiles2, is_stereo):
    oemol1 = oemol_from_smiles(smiles1)
    oemol2 = oemol_from_smiles(smiles2)
    if isinstance(is_stereo, bool):
        assert are_stereoisomers(oemol1, oemol2) == is_stereo
        assert are_stereoisomers(oemol2, oemol1) == is_stereo
    else:
        with pytest.raises(is_stereo):
            _ = are_stereoisomers(oemol1, oemol2)
        with pytest.raises(is_stereo):
            _ = are_stereoisomers(oemol2, oemol1)


@pytest.mark.parametrize(
    ["smiles1", "smiles2", "is_diastereomer"],
    [(case.smiles1, case.smiles2, case.is_diastereomer) for case in PAIR_CASES],
)
def test_are_diastereomers(smiles1, smiles2, is_diastereomer):
    oemol1 = oemol_from_smiles(smiles1)
    oemol2 = oemol_from_smiles(smiles2)
    if isinstance(is_diastereomer, bool):
        assert are_diastereomers(oemol1, oemol2) == is_diastereomer
        assert are_diastereomers(oemol2, oemol1) == is_diastereomer
    else:
        with pytest.raises(UnderspecifiedStereoError):
            _ = are_diastereomers(oemol1, oemol2)
        with pytest.raises(UnderspecifiedStereoError):
            _ = are_diastereomers(oemol2, oemol1)


@pytest.mark.parametrize(
    ["smiles", "smiles_reflected"],
    [(case.smiles, case.smiles_reflected) for case in REFLECTED_CASES],
)
def test_reflect_stereoisomer(smiles, smiles_reflected):
    oemol = oemol_from_smiles(smiles)
    oemol_reflected = oemol_from_smiles(smiles_reflected)
    assert oechem.OEMolToSmiles(reflect_molecule_stereo(oemol)) == oechem.OEMolToSmiles(
        oemol_reflected
    )


@pytest.mark.parametrize(
    ["smiles1", "smiles2", "is_enantiomer"],
    [
        *[(case.smiles1, case.smiles2, case.is_enantiomer) for case in PAIR_CASES],
        *[(case.smiles, case.smiles_reflected, case.is_enantiomer) for case in REFLECTED_CASES],
    ],
)
def test_are_enantiomers(smiles1, smiles2, is_enantiomer):
    oemol1 = oemol_from_smiles(smiles1)
    oemol2 = oemol_from_smiles(smiles2)
    if isinstance(is_enantiomer, bool):
        assert are_enantiomers(oemol1, oemol2) == is_enantiomer
        assert are_enantiomers(oemol2, oemol1) == is_enantiomer
    else:
        with pytest.raises(is_enantiomer):
            _ = are_enantiomers(oemol1, oemol2)
        with pytest.raises(is_enantiomer):
            _ = are_enantiomers(oemol2, oemol1)


@pytest.mark.parametrize(
    ["smiles1", "smiles2", "is_mirror"],
    [
        *[(case.smiles1, case.smiles2, case.is_mirror_image) for case in PAIR_CASES],
        *[(case.smiles, case.smiles_reflected, case.is_mirror) for case in REFLECTED_CASES],
    ],
)
def test_are_mirror_images(smiles1, smiles2, is_mirror):
    oemol1 = oemol_from_smiles(smiles1)
    oemol2 = oemol_from_smiles(smiles2)
    if isinstance(is_mirror, bool):
        assert are_mirror_images(oemol1, oemol2) == is_mirror
        assert are_mirror_images(oemol2, oemol1) == is_mirror
    else:
        with pytest.raises(is_mirror):
            assert are_mirror_images(oemol1, oemol2) == is_mirror
        with pytest.raises(is_mirror):
            assert are_mirror_images(oemol2, oemol1) == is_mirror


@pytest.mark.parametrize(
    ["smiles", "is_meso"], [(case.smiles, case.is_meso) for case in REFLECTED_CASES]
)
def test_is_meso_isomer(smiles, is_meso):
    oemol = oemol_from_smiles(smiles)
    assert is_meso_isomer(oemol) == is_meso


@pytest.mark.parametrize(
    ["smiles", "is_chiral"], [(case.smiles, case.is_chiral) for case in REFLECTED_CASES]
)
def test_is_chiral_molecule(smiles, is_chiral):
    oemol = oemol_from_smiles(smiles)
    if isinstance(is_chiral, bool):
        assert is_chiral_molecule(oemol) == is_chiral
    else:
        with pytest.raises(is_chiral):
            _ = is_chiral_molecule(oemol)


SPECIFIC_STEREOISOMER_CASES = [
    ("[C@](F)(Cl)(Br)I", True, True),
    ("[C@@](F)(Cl)(Br)I", True, True),
    ("[C](F)(Cl)(Br)I", False, False),
    ("[C@@H]([C@H](C(=O)O)O)(C(=O)O)O", True, True),
    ("[CH]([C@H](C(=O)O)O)(C(=O)O)O", False, False),
    ("[CH]([CH](C(=O)O)O)(C(=O)O)O", False, False),
    ("CC[N@]1CCC[C@H](C)C1", True, True),
    ("CC[N@]1CCCC(C)C1", False, False),
    ("CCN1CCC[C@H](C)C1", False, True),
    ("CCN1CCCC(C)C1", False, False),
    ("C1(C)C(C)C1", False, False),
    ("CC1CCC(C)CC1", False, False),
    ("C[C@H]1CCC(C)CC1", False, False),
    ("C[C@H]1CC[C@@H](C)CC1", True, True),
    ("C(F)(Cl)=C(F)(Cl)", False, False),
    (r"C(\F)(/Cl)=C(/F)(\Cl)", True, True),
    ("C1=CCCCC1", True, True),
]


@pytest.mark.parametrize(
    ["smiles", "is_specific_with_nitrogens"],
    [
        (smiles, is_specific_with_nitrogens)
        for smiles, is_specific_with_nitrogens, _ in SPECIFIC_STEREOISOMER_CASES
    ],
)
def test_is_specific_stereoisomer(smiles, is_specific_with_nitrogens):
    oemol = oemol_from_smiles(smiles)
    assert is_specific_stereoisomer(oemol) == is_specific_with_nitrogens


@pytest.mark.parametrize(
    ["smiles", "is_specific_no_nitrogens"],
    [
        (smiles, is_specific_no_nitrogens)
        for smiles, _, is_specific_no_nitrogens in SPECIFIC_STEREOISOMER_CASES
    ],
)
def test_is_specific_stereoisomer_no_nitogens(smiles, is_specific_no_nitrogens):
    oemol = oemol_from_smiles(smiles)
    assert is_specific_stereoisomer(oemol, include_nitrogens=False) == is_specific_no_nitrogens


@pytest.mark.parametrize(
    ["stereo", "reflected_stereo"],
    [
        (oechem.OEAtomStereo_LeftHanded, oechem.OEAtomStereo_RightHanded),
        (oechem.OEAtomStereo_RightHanded, oechem.OEAtomStereo_LeftHanded),
        (oechem.OEAtomStereo_Undefined, oechem.OEAtomStereo_Undefined),
    ],
)
def test_reflect_tetrahedral_stereo(stereo, reflected_stereo):
    assert reflect_tetrahedral_stereo(stereo) == reflected_stereo


def test_reflect_tetrahedral_stereo_error():
    with pytest.raises(ValueError):
        _ = reflect_tetrahedral_stereo(-1)


@pytest.mark.parametrize(
    ["smiles", "expected"],
    [
        ("[C@](F)(Cl)(Br)I", oechem.OEAtomStereo_Left),
        ("[C@@](F)(Cl)(Br)I", oechem.OEAtomStereo_Right),
        ("C(F)(Cl)(Br)I", oechem.OEAtomStereo_Undefined),
    ],
)
def test_get_atom_tetrahedral_stereo(smiles, expected):
    oemol = oemol_from_smiles(smiles)
    assert get_atom_tetrahedral_stereo(next(iter(oemol.GetAtoms()))) == expected


@pytest.mark.parametrize(
    ["smiles", "expected"],
    [
        (r"C/C=C/C", oechem.OEBondStereo_Trans),
        (r"C\C=C\C", oechem.OEBondStereo_Trans),
        (r"C/C=C\C", oechem.OEBondStereo_Cis),
        (r"C\C=C/C", oechem.OEBondStereo_Cis),
        ("CC=CC", oechem.OEBondStereo_Undefined),
    ],
)
def test_get_bond_cistrans_stereo(smiles, expected):
    oemol = oemol_from_smiles(smiles)
    assert get_bond_cistrans_stereo(list(oemol.GetBonds())[1]) == expected


def test_set_atom_tetrahedral_stereo():
    oemol = oemol_from_smiles("C(F)(Cl)(Br)I")
    atom = next(iter(oemol.GetAtoms()))

    assert get_atom_tetrahedral_stereo(atom) == oechem.OEAtomStereo_Undefined

    set_atom_tetrahedral_stereo(atom, stereo_label=oechem.OEAtomStereo_Right)
    assert oechem.OEMolToSmiles(oemol) == "[C@@](F)(Cl)(Br)I"
    assert get_atom_tetrahedral_stereo(atom) == oechem.OEAtomStereo_Right

    set_atom_tetrahedral_stereo(atom, stereo_label=oechem.OEAtomStereo_Left)
    assert oechem.OEMolToSmiles(oemol) == "[C@](F)(Cl)(Br)I"
    assert get_atom_tetrahedral_stereo(atom) == oechem.OEAtomStereo_Left

    set_atom_tetrahedral_stereo(atom, stereo_label=oechem.OEAtomStereo_Undefined)
    assert oechem.OEMolToSmiles(oemol) == "C(F)(Cl)(Br)I"
    assert get_atom_tetrahedral_stereo(atom) == oechem.OEAtomStereo_Undefined


def test_set_bond_cistrans_stereo():
    oemol = oemol_from_smiles("CC=CC")
    bond = list(oemol.GetBonds())[1]

    assert get_bond_cistrans_stereo(bond) == oechem.OEBondStereo_Undefined

    set_bond_cistrans_stereo(bond, stereo_label=oechem.OEBondStereo_Cis)
    assert oechem.OEMolToSmiles(oemol) == r"C/C=C\C"
    assert get_bond_cistrans_stereo(bond) == oechem.OEBondStereo_Cis

    set_bond_cistrans_stereo(bond, stereo_label=oechem.OEBondStereo_Trans)
    assert oechem.OEMolToSmiles(oemol) == "C/C=C/C"
    assert get_bond_cistrans_stereo(bond) == oechem.OEBondStereo_Trans

    set_bond_cistrans_stereo(bond, stereo_label=oechem.OEBondStereo_Undefined)
    assert oechem.OEMolToSmiles(oemol) == "CC=CC"
    assert get_bond_cistrans_stereo(bond) == oechem.OEBondStereo_Undefined


@pytest.mark.parametrize(
    "oemol",
    [oemol.CreateCopy() for oemol in CIP_VALIDATION_OEMOLS_3D],
)
def test_get_all_tetrahedral_stereo_from_coordinates(oemol):
    standardise_oemol(oemol)

    coordinates = get_coordinates(oemol)
    adjacency_matrix = get_adjacency_matrix(oemol)
    chiral_mask = get_atom_mask(oemol, lambda atom: atom.IsChiral())

    assert get_all_tetrahedral_stereo_from_coordinates(
        coordinates=coordinates, adjacency_matrix=adjacency_matrix, chiral_mask=chiral_mask
    ) == [
        get_atom_tetrahedral_stereo(atom) if atom.IsChiral() else oechem.OEAtomStereo_Undefined
        for atom in oemol.GetAtoms()
    ]


@pytest.mark.parametrize("neighbour_size", [(3, 3), (4, 3)])
def test_check_tetrahedral_coordinates(neighbour_size):
    check_tetrahedral_coordinates(
        atom_coordinates=np.random.randn(3), neighbour_coordinates=np.random.randn(*neighbour_size)
    )


@pytest.mark.parametrize(
    ["size", "neighbour_size"],
    [
        (3, (2, 3)),
        (3, (5, 3)),
        (2, (3, 3)),
        (4, (3, 3)),
        (3, (3, 4)),
        (3, (3, 2)),
    ],
)
def test_check_tetrahedral_coordinates_fail(size, neighbour_size):
    with pytest.raises(ValueError):
        check_tetrahedral_coordinates(
            atom_coordinates=np.random.randn(size),
            neighbour_coordinates=np.random.randn(*neighbour_size),
        )


ENUMERATE_CASES = [
    ("[C@H](Cl)(Br)I", {"[C@H](Cl)(Br)I"}),
    ("[C@@H](Cl)(Br)I", {"[C@@H](Cl)(Br)I"}),
    ("[CH](Cl)(Br)I", {"[C@H](Cl)(Br)I", "[C@@H](Cl)(Br)I"}),
    ("[C@H](Cl)(Br)I |a:0|", {"[C@H](Cl)(Br)I"}),
    ("[C@H](Cl)(Br)I |&1:0|", {"[C@H](Cl)(Br)I", "[C@@H](Cl)(Br)I"}),
    ("[C@H](Cl)(Br)I |o1:0|", {"[C@H](Cl)(Br)I", "[C@@H](Cl)(Br)I"}),
    ("[C@H]([C@H](Br)I)(F)Cl", {"[C@H]([C@H](Br)I)(F)Cl"}),
    ("[CH]([C@H](Br)I)(F)Cl", {"[C@H]([C@H](Br)I)(F)Cl", "[C@@H]([C@H](Br)I)(F)Cl"}),
    ("[C@H]([CH](Br)I)(F)Cl", {"[C@H]([C@H](Br)I)(F)Cl", "[C@H]([C@@H](Br)I)(F)Cl"}),
    (
        "[CH]([CH](Br)I)(F)Cl",
        {
            "[C@H]([C@H](Br)I)(F)Cl",
            "[C@H]([C@@H](Br)I)(F)Cl",
            "[C@@H]([C@H](Br)I)(F)Cl",
            "[C@@H]([C@@H](Br)I)(F)Cl",
        },
    ),
    ("[C@H]([C@H](Br)I)(F)Cl |&1:0,1|", {"[C@H]([C@H](Br)I)(F)Cl", "[C@@H]([C@@H](Br)I)(F)Cl"}),
    # Meso compound
    (
        "[CH]([CH](F)Cl)(F)Cl",
        {"[C@H]([C@@H](F)Cl)(F)Cl", "[C@@H]([C@H](F)Cl)(F)Cl", "[C@@H]([C@@H](F)Cl)(F)Cl"},
    ),
    # Chiral nitrogen, unspecified stereo
    ("CCCN(C)CC", {"CCCN(C)CC"}),
    # Chiral ammonium, unspecified stereo
    ("CCCC[N+](C)(CC)CCC", {"CCCC[N@+](C)(CC)CCC", "CCCC[N@@+](C)(CC)CCC"}),
    # Chiral ammonium, specified stereo
    ("CCCC[N@+](C)(CC)CCC", {"CCCC[N@+](C)(CC)CCC"}),
    # Meso compound
    # Double substitued cyclohexane
    ("CC1CCC(O)CC1", {"C[C@H]1CC[C@@H](CC1)O", "C[C@H]1CC[C@H](CC1)O"}),
    ("C1CC(O)CC[C@@H]1C", {"C[C@H]1CC[C@@H](CC1)O", "C[C@H]1CC[C@H](CC1)O"}),
    # Single substituted cyclohexane
    ("CC1CCCCC1", {"CC1CCCCC1"}),
    # Combination of double and singly substituted cyclohexane
    ("C1CC(C2CCCCC2)CCC1C", {"C[C@H]1CC[C@@H](CC1)C2CCCCC2", "C[C@H]1CC[C@H](CC1)C2CCCCC2"}),
    (
        "C1CC(C)CCC1C1CCC(C)CC1",
        {
            "C[C@H]1CC[C@@H](CC1)[C@H]2CC[C@@H](CC2)C",
            "C[C@H]1CC[C@@H](CC1)[C@H]2CC[C@H](CC2)C",
            "C[C@H]1CC[C@H](CC1)[C@H]2CC[C@H](CC2)C",
        },
    ),
    # Triple substituted cyclohexane
    ("C1C(C)CC(C)CC1C", {"C[C@H]1C[C@@H](C[C@@H](C1)C)C", "C[C@H]1C[C@H](C[C@H](C1)C)C"}),
    # Double substitued cyclopropane
    ("C1(C)CC1C", {"C[C@@H]1C[C@@H]1C", "C[C@@H]1C[C@H]1C", "C[C@H]1C[C@@H]1C"}),
    # Double bond
    ("C(F)(Cl)=C(F)(Cl)", {r"C(=C(/F)\Cl)(\F)/Cl", r"C(=C(\F)/Cl)(\F)/Cl"}),
    # Oxime
    ("CC(CC)=NO", {r"CC/C(=N/O)/C", r"CC/C(=N\O)/C"}),
]

# Cases that only apply to chemutils and not rdkit
ENUMERATE_CASES_CHEMUTILS = [
    # Chiral nitrogen, specified stereo
    ("CCC[N@](C)CC", {"CCC[N@](C)CC"}),
    # Double bond across from a stereocenter in a ring
    ("CC=C1CCC(C)CC1", {r"C/C=C/1\CC[C@@H](CC1)C", r"C/C=C/1\CC[C@H](CC1)C"}),
    # Xylene-like
    ("CC=C1C=CC(=CC)C=C1", {r"C/C=c/1\cc/c(=C/C)/cc1", r"C/C=c/1\cc/c(=C\C)/cc1"}),
]


@pytest.mark.parametrize(["smiles", "expected"], ENUMERATE_CASES + ENUMERATE_CASES_CHEMUTILS)
def test_enumerate_stereoisomers(smiles, expected):
    assert sorted(
        cxsmiles_from_oemol(oemol) for oemol in enumerate_stereoisomers(oemol_from_smiles(smiles))
    ) == sorted(expected)


@pytest.mark.parametrize(["smiles", "expected"], ENUMERATE_CASES)
def test_enumerate_stereoisomers_rdkit(smiles, expected):
    m = Chem.MolFromSmiles(smiles)
    isomers = tuple(EnumerateStereoisomers(m))

    assert sorted(canonicalize_smiles(Chem.MolToSmiles(isomer)) for isomer in isomers) == sorted(
        expected
    )
