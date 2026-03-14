import pytest

from chemutils.molecule import (
    oemol_from_smiles,
    representative_cxsmiles_from_oemol,
    representative_cxsmiles_from_smiles,
    representative_smiles_from_oemol,
    representative_smiles_from_smiles,
)
from chemutils.props import (
    calculate_properties_for_representative_smiles,
    calculate_properties_for_smiles,
)
from chemutils.props.identifier import CXSMILES, SMILES

CASES = [
    # Basic example
    (
        "CC(C)(C)C",
        "CC(C)(C)C",
        "CC(C)(C)C",
        "CC(C)(C)C",
        "CC(C)(C)C",
    ),
    # Uncharged carboxylic acid
    (
        "CC(=O)O",
        "CC(=O)O",
        "CC(=O)O",
        "CC(=O)O",
        "CC(=O)O",
    ),
    # Charged carboxylic acid
    (
        "CC(=O)[O-]",
        "CC(=O)[O-]",
        "CC(=O)[O-]",
        "CC(=O)O",
        "CC(=O)O",
    ),
    # Nitro group
    (
        "C[N+](=O)O",
        "C[N+](=O)O",
        "C[N+](=O)O",
        "C[N+](=O)[O-]",
        "C[N+](=O)[O-]",
    ),
    # CXSMILES, absolute
    (
        "[C@](F)(Cl)(Br)I |a:0|",
        "[C@](F)(Cl)(Br)I",
        "[C@](F)(Cl)(Br)I |a:0|",
        "[C@](F)(Cl)(Br)I",
        "[C@](F)(Cl)(Br)I",
    ),
    # CXSMILES, 'and' enhanced stereo
    (
        "[C@](F)(Cl)(Br)I |&1:0|",
        "C(F)(Cl)(Br)I",
        "[C@](F)(Cl)(Br)I |&1:0|",
        "C(F)(Cl)(Br)I",
        "[C@@](F)(Cl)(Br)I |&1:0|",
    ),
    # CXSMILES, 'or' enhanced stereo
    (
        "[C@](F)(Cl)(Br)I |o1:0|",
        "C(F)(Cl)(Br)I",
        "[C@](F)(Cl)(Br)I |o1:0|",
        "C(F)(Cl)(Br)I",
        "[C@@](F)(Cl)(Br)I |o1:0|",
    ),
    # CXSMILES, meso compound
    (
        "[C@H](F)(Cl)[C@H](F)(Cl) |&1:0,3|",
        "C(C(F)Cl)(F)Cl",
        "[C@@H]([C@@H](F)Cl)(F)Cl |&1:0,1|",
        "[C@@H]([C@@H](F)Cl)(F)Cl",
        "[C@@H]([C@@H](F)Cl)(F)Cl",
    ),
]


@pytest.mark.parametrize(["input", "smiles", "cxsmiles", "rep_smiles", "rep_cxsmiles"], CASES)
def test_rep_smiles(input, smiles, cxsmiles, rep_smiles, rep_cxsmiles):
    assert calculate_properties_for_smiles(input, [SMILES, CXSMILES]) == {
        "smiles": smiles,
        "cxsmiles": cxsmiles,
    }
    assert calculate_properties_for_representative_smiles(input, [SMILES, CXSMILES]) == {
        "smiles": rep_smiles,
        "cxsmiles": rep_cxsmiles,
    }


@pytest.mark.parametrize(["input", "smiles", "cxsmiles", "rep_smiles", "rep_cxsmiles"], CASES)
def test_rep_smiles_from(input, smiles, cxsmiles, rep_smiles, rep_cxsmiles):
    assert representative_smiles_from_smiles(input) == rep_smiles
    assert representative_cxsmiles_from_smiles(input) == rep_cxsmiles
    oemol = oemol_from_smiles(input)
    assert representative_smiles_from_oemol(oemol) == rep_smiles
    assert representative_cxsmiles_from_oemol(oemol) == rep_cxsmiles
