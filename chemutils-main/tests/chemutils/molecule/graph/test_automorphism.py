import pytest
from openeye import oechem

from chemutils.molecule import (
    has_nontrivial_automorphism,
    iterate_automorphism_matches,
    make_hydrogens_implicit,
    oemol_from_smiles,
)

from ...props.drug_smiles import DRUG_SMILES


@pytest.mark.parametrize(
    ["smiles", "expected"],
    [
        ["C", 1],
        ["CO", 1],
        ["CC", 2],
        ["CCC", 2],
        ["c1ccccc1", 12],
        ["S(=O)(=O)(O)O", 24],
        # Equal to 24 * 24 * 2!
        ["S(=O)(=O)(O)O.S(=O)(=O)(O)O", 1152],
        # Equal to 24 * 24 * 24 * 3!
        ["S(=O)(=O)(O)O.S(=O)(=O)(O)O.S(=O)(=O)(O)O", 82944],
    ],
)
def test_num_automorphisms(smiles, expected):
    oemol = oemol_from_smiles(smiles)
    assert len(list(iterate_automorphism_matches(oemol))) == expected
    assert has_nontrivial_automorphism(oemol) == (expected > 1)


# Ensure that automorphism behaviour matches OpenEye
@pytest.mark.parametrize(
    "smiles",
    ["S(=O)(=O)(O)O", *DRUG_SMILES.values()],
)
def test_automorphisms_match_openeye(smiles):
    oemol = oemol_from_smiles(smiles)
    make_hydrogens_implicit(oemol, remove_isotopic_hydrogens=True)

    automorphs_chemutils = {
        tuple(atom.GetIdx() for atom in match.GetTargetAtoms())
        for match in iterate_automorphism_matches(oemol)
    }
    automorphs_openeye = {
        tuple(atom.GetIdx() for atom in match.GetTargetAtoms())
        for match in oechem.OEGetAutomorphs(oemol)
    }
    assert automorphs_chemutils == automorphs_openeye
