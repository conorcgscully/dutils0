import pytest

from chemutils.molecule import (
    has_maximum_common_substructure,
    iterate_maximum_common_substructure_matches,
    oemol_from_smiles,
)


@pytest.mark.parametrize(
    ["smiles1", "smiles2", "min_atoms", "num_unique_matches", "num_nonunique_matches"],
    [
        ["CC", "CCC", None, 2, 4],
        ["c1ccccc1", "c1ccccc1C", None, 1, 12],
        ["CCCNN", "CCNNN", None, 1, 1],
        ["CCCNN", "CCNNN", 5, 0, 0],
        ["CO", "NF", None, 0, 0],
    ],
)
def test_maximum_substructure(
    smiles1, smiles2, min_atoms, num_unique_matches, num_nonunique_matches
):
    oemol1 = oemol_from_smiles(smiles1)
    oemol2 = oemol_from_smiles(smiles2)

    matches_unique = list(
        iterate_maximum_common_substructure_matches(
            pattern=oemol1, target=oemol2, unique=True, min_atoms=min_atoms
        )
    )
    assert len(matches_unique) == num_unique_matches

    matches_nonunique = list(
        iterate_maximum_common_substructure_matches(
            pattern=oemol1, target=oemol2, unique=False, min_atoms=min_atoms
        )
    )
    assert len(matches_nonunique) == num_nonunique_matches

    assert has_maximum_common_substructure(pattern=oemol1, target=oemol2, min_atoms=min_atoms) == (
        len(matches_unique) > 0
    )
