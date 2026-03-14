from pathlib import Path

import numpy as np
import pytest
from openeye import oechem

from chemutils import fs
from chemutils.molecule import (
    InvalidSMARTSError,
    get_substructure_counts,
    get_substructure_indices,
    get_substructure_mask,
    get_target_mask_from_match,
    has_substructure,
    iterate_substructure_matches,
    oemol_from_smiles,
)

SMARTS_CASES = fs.read_yaml(Path(__file__).parent / "smarts.yaml")


@pytest.mark.parametrize(
    ["smiles", "smarts", "has_match"],
    (
        ("CCC", "[#6]", True),
        ("c1c2ccccc2ccc1", "c1ccccc1", True),
        ("c1ccccc1", "C", False),
        *[(case["smiles"], case["smarts"], True) for case in SMARTS_CASES],
    ),
)
def test_has_substructure_match(smiles, smarts, has_match):
    mol = oemol_from_smiles(smiles)
    assert has_substructure(pattern=smarts, target=mol) == has_match


@pytest.mark.parametrize(
    ["pattern", "target", "expected_unique_matches", "expected_nonunique_matches"],
    [
        ("c1ccccc1", "c1ccccc1", 1, 12),
        ("c1ccccc1", "c1ccccc1C", 1, 12),
        ("c1ccccc1", "c1ccccc1Cc1ccccc1", 2, 24),
        ("CCC", "CC(C)C", 3, 6),
        ("CCC", "CC(C)(C)C", 6, 12),
    ],
)
def test_num_substructure_matches(
    pattern, target, expected_unique_matches, expected_nonunique_matches
):
    pattern = oemol_from_smiles(pattern)
    target = oemol_from_smiles(target)
    unique_matches = list(iterate_substructure_matches(pattern=pattern, target=target, unique=True))
    nonunique_matches = list(
        iterate_substructure_matches(pattern=pattern, target=target, unique=False)
    )
    assert len(unique_matches) == expected_unique_matches
    assert len(nonunique_matches) == expected_nonunique_matches


@pytest.mark.parametrize(
    ["smiles", "smarts", "expected_matches"],
    [(case["smiles"], case["smarts"], case["matches"]) for case in SMARTS_CASES],
)
def test_get_substructure_matches(smiles, smarts, expected_matches):
    mol = oemol_from_smiles(smiles)
    matches = iterate_substructure_matches(pattern=smarts, target=mol, unique=True)
    matches = [
        {atom.pattern.GetIdx(): atom.target.GetIdx() for atom in match.GetAtoms()}
        for match in matches
    ]
    assert matches == expected_matches


@pytest.mark.parametrize(
    ["smiles", "smarts", "expected_matches"],
    [(case["smiles"], case["smarts"], case["matches"]) for case in SMARTS_CASES],
)
def test_get_target_mask(smiles, smarts, expected_matches):
    mol = oechem.OEGraphMol()
    oechem.OESmilesToMol(mol, smiles)
    for match, expected_match in zip(
        iterate_substructure_matches(pattern=smarts, target=mol, unique=True),
        expected_matches,
        strict=True,
    ):
        mask = get_target_mask_from_match(match)
        expected_mask = np.zeros(mol.NumAtoms())
        for atom_idx in expected_match.values():
            expected_mask[atom_idx] = True
        assert np.all(mask == expected_mask)


ACCEPTOR_SMARTS = (
    "[$([O;H1;v2]),"
    "$([O;H0;v2;!$(O=N-*),"
    "$([O;-;!$(*-N=O)]),"
    "$([o;+0])]),"
    "$([n;+0;!X3;!$([n;H1](cc)cc),"
    "$([$([N;H0]#[C&v4])]),"
    "$([N&v3;H0;$(Nc)])]),"
    "$([F;$(F-[#6]);!$(FC[F,Cl,Br,I])])]"
)

DONOR_SMARTS = (
    "[$([N&!H0&v3,N&!H0&+1&v4,n&H1&+0,$([$([Nv3](-C)(-C)-C)]),"
    "$([$(n[n;H1]),"
    "$(nc[n;H1])])]),"
    "$([NX3,NX2]([!O,!S])!@C(!@[NX3,NX2]([!O,!S]))!@[NX3,NX2]([!O,!S])),"
    "$([O,S;H1;+0])]"
)


ERROR_SMARTS = ["Q", "not a smarts", "[E2]"]


@pytest.mark.parametrize("smarts", ERROR_SMARTS)
def test_get_substructure_matches_error(smarts):
    mol = oemol_from_smiles("CCC")
    with pytest.raises(InvalidSMARTSError):
        for _ in iterate_substructure_matches(pattern=smarts, target=mol, unique=True):
            pass


@pytest.mark.parametrize(
    ["smiles", "smarts", "matches"],
    (*[(case["smiles"], case["smarts"], case["matches"]) for case in SMARTS_CASES],),
)
def test_get_substructure_counts(smiles, smarts, matches):
    oemol = oemol_from_smiles(smiles)

    counts = [0] * oemol.NumAtoms()
    for match in matches:
        for target_idx in match.values():
            counts[target_idx] += 1

    assert np.all(counts == get_substructure_counts(pattern=smarts, target=smiles, unique=True))


@pytest.mark.parametrize(
    ["smiles", "smarts", "matches"],
    (*[(case["smiles"], case["smarts"], case["matches"]) for case in SMARTS_CASES],),
)
def test_get_substructure_indices(smiles, smarts, matches):
    indices = []
    for match in matches:
        for target_idx in match.values():
            indices.append(target_idx)

    assert np.all(indices == list(get_substructure_indices(pattern=smarts, target=smiles)))


@pytest.mark.parametrize(
    ["smiles", "smarts", "matches"],
    (*[(case["smiles"], case["smarts"], case["matches"]) for case in SMARTS_CASES],),
)
def test_get_substructure_mask(smiles, smarts, matches):
    oemol = oemol_from_smiles(smiles)

    mask = np.zeros(oemol.NumAtoms(), dtype=np.bool_)
    for match in matches:
        for target_idx in match.values():
            mask[target_idx] = True

    assert np.all(mask == get_substructure_mask(pattern=smarts, target=smiles, unique=True))
