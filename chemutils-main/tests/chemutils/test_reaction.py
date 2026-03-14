import pytest

from chemutils.molecule import (
    InvalidReactionError,
    MultipleReactionMatchesError,
    ReactionNotApplicableError,
    apply_reaction,
    oemol_from_smiles,
    smiles_from_oemol,
)


@pytest.mark.parametrize(
    ["smiles", "smirks", "correct_valency", "expected"],
    [
        ("C=C", "[C:1]=[C:2]>>[C:1][C:2]", False, "[CH2][CH2]"),
        ("C=C", "[C:1]=[C:2]>>[C:1][C:2]", True, "CC"),
        # Alchemical reaction
        ("CCN", "[C:1][N:2]>>[C:1][C:2]", False, "CC[CH2]"),
        ("CCN", "[C:1][N:2]>>[C:1][C:2]", True, "CC[CH2]"),
        # Protonation
        ("CCN", "[NH2:1]>>[NH3+:1]", True, "CC[NH3+]"),
        ("CCN", "[NH2:1]>>[NH3+:1]", True, "CC[NH3+]"),
        # Acrylamide Addition
        ("CC=CC(=O)C", "[C:1]=[C:2][C:3]=[O:4]>>[C:1][C:2][C:3]=[O:4]", False, "C[CH][CH]C(=O)C"),
        ("CC=CC(=O)C", "[C:1]=[C:2][C:3]=[O:4]>>[C:1][C:2][C:3]=[O:4]", True, "CCCC(=O)C"),
    ],
)
def test_apply_reaction(smiles, smirks, correct_valency, expected):
    oemol = oemol_from_smiles(smiles)
    oemol2 = apply_reaction(reaction=smirks, mol=oemol, correct_valency=correct_valency)
    assert smiles_from_oemol(oemol2) == expected


@pytest.mark.parametrize(
    ["smiles", "smirks"],
    [
        ("C=C", "[C:1]=[C:2>>[C:1][C:2]"),
        ("C=C", "[C:1]=[C:2]>[C:1][C:2]"),
        ("C=C", "[C:1]=[C:3]>>[C:1][C:2]"),
    ],
)
def test_apply_reaction_invalid_smirks(smiles, smirks):
    with pytest.raises(InvalidReactionError) as exc:
        oemol = oemol_from_smiles(smiles)
        _ = apply_reaction(reaction=smirks, mol=oemol)

    assert str(exc.value) == f"Invalid reaction `{smirks}`."


@pytest.mark.parametrize(["smiles", "smirks"], [("C#C", "[C:1]=[C:2]>>[C:1][C:2]")])
def test_apply_reaction_no_match(smiles, smirks):
    with pytest.raises(ReactionNotApplicableError) as exc:
        oemol = oemol_from_smiles(smiles)
        _ = apply_reaction(reaction=smirks, mol=oemol)

    assert str(exc.value) == f"Reaction `{smirks}` not applicable to molecule `{smiles}`."


@pytest.mark.parametrize(["smiles", "smirks"], [("C=CC=C", "[C:1]=[C:2]>>[C:1][C:2]")])
def test_apply_reaction_multiple_matches(smiles, smirks):
    with pytest.raises(MultipleReactionMatchesError) as exc:
        oemol = oemol_from_smiles(smiles)
        _ = apply_reaction(reaction=smirks, mol=oemol)

    assert str(exc.value) == f"Reaction `{smirks}` matches more than one part of `{smiles}`."


def test_apply_reaction_copy_map_idx():
    smiles = "CC=CC"
    smirks = "[C:1]=[C:2]>>[C:1][C:2]"

    oemol1 = oemol_from_smiles(smiles)
    oemol1_product = apply_reaction(reaction=smirks, mol=oemol1, correct_valency=True)

    assert [atom.GetMapIdx() for atom in oemol1_product.GetAtoms()] == [0, 0, 0, 0]

    oemol2 = oemol_from_smiles(smiles)
    oemol2_product = apply_reaction(
        reaction=smirks, mol=oemol2, correct_valency=True, copy_map_idx=True
    )

    assert [atom.GetMapIdx() for atom in oemol2_product.GetAtoms()] == [0, 1, 2, 0]
