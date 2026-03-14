from pathlib import Path

import pytest
from openeye import oechem

from chemutils import fs
from chemutils.molecule import (
    cxsmiles_from_oemol,
    oemol_from_smiles,
    read_molecules_str,
    write_molecule_str,
)
from chemutils.molecule.enhanced_stereo import (
    convert_enhanced_stereo_to_unspecified,
    get_canonical_stereoisomer,
    get_stereoisomer_mixtures_from_oemol,
    remove_redundant_enhanced_stereo,
)

CASES = [
    # Examples from
    # https://support.collaborativedrug.com/hc/en-us/articles/360020872171
    # Mixture of 4 diastereoisomers
    (
        "C[C@H]([C@@H](N)Cl)O |&1:1,&2:2|",
        [
            [
                "C[C@@H]([C@@H](N)Cl)O",
                "C[C@@H]([C@H](N)Cl)O",
                "C[C@H]([C@@H](N)Cl)O",
                "C[C@H]([C@H](N)Cl)O",
            ]
        ],
        "C[C@@H]([C@@H](N)Cl)O |&1:1,&2:2|",
    ),
    # Mixture of 2 enantiomers
    (
        "C[C@H]([C@@H](N)Cl)O |&1:1,2|",
        [["C[C@@H]([C@H](N)Cl)O", "C[C@H]([C@@H](N)Cl)O"]],
        "C[C@@H]([C@H](N)Cl)O |&1:1,2|",
    ),
    # Mixture of 2 enantiomers
    (
        "C[C@H]([C@H](N)Cl)O |&3:1,2|",
        [["C[C@@H]([C@@H](N)Cl)O", "C[C@H]([C@H](N)Cl)O"]],
        "C[C@@H]([C@@H](N)Cl)O |&1:1,2|",
    ),
    # Pure sample of unknown configuration
    (
        "C[C@H]([C@@H](N)Cl)O |o1:1,o2:2|",
        [
            ["C[C@@H]([C@@H](N)Cl)O"],
            ["C[C@@H]([C@H](N)Cl)O"],
            ["C[C@H]([C@@H](N)Cl)O"],
            ["C[C@H]([C@H](N)Cl)O"],
        ],
        "C[C@@H]([C@@H](N)Cl)O |o1:1,o2:2|",
    ),
    # Pure sample of unknown absolute configuration
    (
        "C[C@H]([C@@H](N)Cl)O |o1:1,2|",
        [["C[C@@H]([C@H](N)Cl)O"], ["C[C@H]([C@@H](N)Cl)O"]],
        "C[C@@H]([C@H](N)Cl)O |o1:1,2|",
    ),
    # Pure sample of unknown absolute configuration
    (
        "C[C@H]([C@H](N)Cl)O |o1:1,2|",
        [["C[C@@H]([C@@H](N)Cl)O"], ["C[C@H]([C@H](N)Cl)O"]],
        "C[C@@H]([C@@H](N)Cl)O |o1:1,2|",
    ),
]


@pytest.mark.parametrize(["smiles", "mixture"], [(case[0], case[1]) for case in CASES])
def test_get_stereoisomer_mixtures_from_oemol(smiles, mixture):
    oemol = oemol_from_smiles(smiles)
    assert [
        [oechem.OEMolToSmiles(oemol) for oemol in group]
        for group in get_stereoisomer_mixtures_from_oemol(oemol)
    ] == mixture


@pytest.mark.parametrize(["smiles", "canonical_smiles"], [(case[0], case[2]) for case in CASES])
def test_get_canonical_stereoisomer(smiles, canonical_smiles):
    oemol = oemol_from_smiles(smiles)
    assert cxsmiles_from_oemol(get_canonical_stereoisomer(oemol)) == canonical_smiles


@pytest.mark.parametrize(
    ["canonical_smiles", "equivalent_smiles_list"],
    [
        # Ensure order/ids of the options are not important
        (
            "C[C@@H]([C@@H](N)Cl)O |&1:1,&2:2|",
            [
                "C[C@@H]([C@@H](N)Cl)O |&1:2,&2:1|",
                "C[C@@H]([C@@H](N)Cl)O |&2:2,&1:1|",
                "C[C@@H]([C@@H](N)Cl)O |&2:1,&1:2|",
                "C[C@@H]([C@@H](N)Cl)O |&4:1,&7:2|",
            ],
        )
    ],
)
def test_equivalent_canonical_stereoisomers(canonical_smiles, equivalent_smiles_list):
    for equivalent_smiles in equivalent_smiles_list:
        oemol = oemol_from_smiles(equivalent_smiles)
        equivalent_canonical_smiles = cxsmiles_from_oemol(get_canonical_stereoisomer(oemol))
        assert equivalent_canonical_smiles == canonical_smiles


def test_enhanced_stereo_real_example():
    molecules = list(fs.read_molecules(Path(__file__).parent / "data" / "enhanced_stereo.sdf"))

    assert [cxsmiles_from_oemol(mol) for mol in molecules] == [
        "CC1([C@@H]([C@@H]1C(=O)O[C@@H](C#N)c2cccc(c2)Oc3ccccc3)C=C(Cl)Cl)C |&1:2,3,&2:7| C22H19Cl2NO3",
        "CC1([C@H]([C@H]1C(=O)O[C@H](C#N)c2cccc(c2)Oc3ccccc3)C=C(Cl)Cl)C |&1:2,3,&2:7| C22H19Cl2NO3",
        "CC1([C@@H]([C@@H]1C(=O)O[C@@H](C#N)c2cccc(c2)Oc3ccccc3)C=C(Cl)Cl)C |&1:2,3,&2:7| C22H19Cl2NO3",
        "CC1([C@H]([C@H]1C(=O)O[C@@H](C#N)c2cccc(c2)Oc3ccccc3)C=C(Cl)Cl)C |&1:2,3,&2:7| C22H19Cl2NO3",
        "CC1([C@@H]([C@@H]1C(=O)O[C@H](C#N)c2cccc(c2)Oc3ccccc3)C=C(Cl)Cl)C |&1:2,3,&2:7| C22H19Cl2NO3",
    ]

    molecules = [get_canonical_stereoisomer(mol) for mol in molecules]

    sdf_str = write_molecule_str(*molecules, format="SDF")
    molecules = list(read_molecules_str(sdf_str, format="SDF"))

    assert [cxsmiles_from_oemol(mol) for mol in molecules] == [
        "CC1([C@@H]([C@@H]1C(=O)O[C@@H](C#N)c2cccc(c2)Oc3ccccc3)C=C(Cl)Cl)C |&1:2,3,&2:7| C22H19Cl2NO3",
        "CC1([C@@H]([C@@H]1C(=O)O[C@@H](C#N)c2cccc(c2)Oc3ccccc3)C=C(Cl)Cl)C |&1:2,3,&2:7| C22H19Cl2NO3",
        "CC1([C@@H]([C@@H]1C(=O)O[C@@H](C#N)c2cccc(c2)Oc3ccccc3)C=C(Cl)Cl)C |&1:2,3,&2:7| C22H19Cl2NO3",
        "CC1([C@@H]([C@@H]1C(=O)O[C@@H](C#N)c2cccc(c2)Oc3ccccc3)C=C(Cl)Cl)C |&1:2,3,&2:7| C22H19Cl2NO3",
        "CC1([C@@H]([C@@H]1C(=O)O[C@@H](C#N)c2cccc(c2)Oc3ccccc3)C=C(Cl)Cl)C |&1:2,3,&2:7| C22H19Cl2NO3",
    ]

    # Check that SDFs are *identical*
    assert len({sdf.strip() for sdf in sdf_str.split("$$$$") if sdf.strip()}) == 1


@pytest.mark.parametrize(
    ["smiles", "expected_smiles"],
    [
        ("[C@](F)(Cl)(Br)I", "[C@](F)(Cl)(Br)I"),
        ("[C@](F)(Cl)(Br)I |a:0|", "[C@](F)(Cl)(Br)I"),
        ("[C@](F)(Cl)(Br)I |&1:0|", "C(F)(Cl)(Br)I"),
        ("[C@](F)(Cl)(Br)I |o1:0|", "C(F)(Cl)(Br)I"),
        ("[C@]([C@@](F)(Cl)Br)(N)(O)I", "[C@]([C@@](F)(Cl)Br)(N)(O)I"),
        ("[C@]([C@@](F)(Cl)Br)(N)(O)I |a:0,&1:1|", "[C@](C(F)(Cl)Br)(N)(O)I"),
        ("[C@]([C@@](F)(Cl)Br)(N)(O)I |a:0,o1:1|", "[C@](C(F)(Cl)Br)(N)(O)I"),
        ("[C@]([C@@](F)(Cl)Br)(N)(O)I |&1:0|", "[C@@](C(N)(O)I)(F)(Cl)Br"),
        ("[C@]([C@@](F)(Cl)Br)(N)(O)I |o1:0,1|", "C(C(F)(Cl)Br)(N)(O)I"),
    ],
)
def test_convert_enhanced_stereo_to_unspecified(smiles, expected_smiles):
    mol = oemol_from_smiles(smiles)
    assert cxsmiles_from_oemol(mol) == smiles
    convert_enhanced_stereo_to_unspecified(mol)
    assert cxsmiles_from_oemol(mol) == expected_smiles


@pytest.mark.parametrize(
    ["smiles", "expected"],
    [
        # Remove redundant absolute stereo labels
        ("[C@](F)(Cl)(Br)I |a:0|", "[C@](F)(Cl)(Br)I"),
        ("[C@@](F)(Cl)(Br)I |a:0|", "[C@@](F)(Cl)(Br)I"),
        # Remove stereo labels for things that aren't actually chiral
        ("[C@](F)(F)(Br)I |&1:0|", "[C@@](I)(Br)(F)F"),
        ("[C@](F)(F)(Br)I |o1:0|", "[C@@](I)(Br)(F)F"),
        # Don't remove and's and or's that mean something
        ("[C@](F)(Cl)(Br)I |&1:0|", "[C@](F)(Cl)(Br)I |&1:0|"),
        ("[C@@](F)(Cl)(Br)I |&1:0|", "[C@@](F)(Cl)(Br)I |&1:0|"),
        ("[C@](F)(Cl)(Br)I |o1:0|", "[C@](F)(Cl)(Br)I |o1:0|"),
        ("[C@@](F)(Cl)(Br)I |o1:0|", "[C@@](F)(Cl)(Br)I |o1:0|"),
        # Meso compounds, remove redunant labels
        ("[C@H](F)(Cl)[C@H](F)(Cl) |&1:0,3|", "[C@@H]([C@@H](F)Cl)(F)Cl"),
        ("[C@@H](F)(Cl)[C@@H](F)(Cl) |&1:0,3|", "[C@@H]([C@@H](F)Cl)(F)Cl"),
        ("[C@H](F)(Cl)[C@H](F)(Cl) |o1:0,3|", "[C@@H]([C@@H](F)Cl)(F)Cl"),
        ("[C@@H](F)(Cl)[C@@H](F)(Cl) |o1:0,3|", "[C@@H]([C@@H](F)Cl)(F)Cl"),
        # Non meso compounds, keep stereo
        ("[C@H](F)(Cl)[C@@H](F)(Cl) |&1:0,3|", "[C@H]([C@@H](F)Cl)(F)Cl |&1:0,1|"),
        ("[C@@H](F)(Cl)[C@H](F)(Cl) |&1:0,3|", "[C@@H]([C@H](F)Cl)(F)Cl |&1:0,1|"),
        ("[C@H](F)(Cl)[C@@H](F)(Cl) |o1:0,3|", "[C@H]([C@@H](F)Cl)(F)Cl |o1:0,1|"),
        ("[C@@H](F)(Cl)[C@H](F)(Cl) |o1:0,3|", "[C@@H]([C@H](F)Cl)(F)Cl |o1:0,1|"),
        ("[C@H](F)(Cl)[C@H](F)(Cl) |&1:0,&2:3|", "[C@@H]([C@@H](F)Cl)(F)Cl |&1:0,&2:1|"),
        ("[C@H](F)(Cl)[C@H](F)(Cl) |o1:0,o2:3|", "[C@@H]([C@@H](F)Cl)(F)Cl |o1:0,o2:1|"),
        # Cyclic meso compounds
        ("C[C@H]1C[C@H](C=C1)C", "C[C@@H]1C[C@@H](C=C1)C"),
        ("C[C@@H]1C[C@@H](C=C1)C", "C[C@@H]1C[C@@H](C=C1)C"),
        ("C[C@H]1C[C@H](C=C1)C |&1:1,3|", "C[C@@H]1C[C@@H](C=C1)C"),
        ("C[C@H]1C[C@H](C=C1)C |o1:1,3|", "C[C@@H]1C[C@@H](C=C1)C"),
        # Non-meso compounds
        ("[C@@H]([C@H](F)Cl)([C@@H](F)Cl)F |&1:1,4|", "[C@@H]([C@H](F)Cl)([C@@H](F)Cl)F |&1:1,4|"),
        ("[C@@H]([C@H](F)Cl)([C@H](F)Cl)F |&1:1,4|", "[C@H](F)([C@@H](Cl)F)[C@@H](Cl)F |&1:2,5|"),
    ],
)
def test_remove_redundant_enhanced_stereo(smiles, expected):
    oemol = oemol_from_smiles(smiles)
    remove_redundant_enhanced_stereo(oemol)
    assert cxsmiles_from_oemol(oemol) == expected
