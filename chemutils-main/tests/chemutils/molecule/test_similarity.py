from unittest.mock import patch

import numpy as np
import pytest

from chemutils.molecule.similarity import (
    InvalidNumResultsError,
    InvalidSmilesError,
    SimilarMolecule,
    _fingerprint_from_mol,
    calculate_similarity_in_parallel,
    get_most_similar_molecules,
)
from chemutils.molecule.smiles import oemol_from_smiles


def test_calculate_similarity_single_target() -> None:
    query_smiles = "CCO"
    target_smiles_list = ["CCO"]
    expected_similarity = {
        "CCO": 1.0,
    }
    result: dict[str, np.float16 | None] = calculate_similarity_in_parallel(
        query=query_smiles, targets=target_smiles_list
    )
    assert result == expected_similarity


def test_calculate_similarity_multiple_targets() -> None:
    query_smiles = "CCO"
    target_smiles_list = ["CCO", "CCC", "CCN"]
    result = calculate_similarity_in_parallel(query=query_smiles, targets=target_smiles_list)
    # Convert to float16 for comparison
    assert result == {
        "CCO": np.float16(1.0),
        "CCC": np.float16(0.4285),
        "CCN": np.float16(0.3333),
    }


def test_calculate_similarity_no_targets() -> None:
    query_smiles = "CCO"
    target_smiles_list: list[str] = []
    result: dict[str, np.float16 | None] = calculate_similarity_in_parallel(
        query=query_smiles, targets=target_smiles_list
    )
    assert result == {}


def test_calculate_similarity_invalid_smiles() -> None:
    query_smiles = "invalid_smiles"
    target_smiles_list = ["CCCC"]
    with pytest.raises(InvalidSmilesError):
        calculate_similarity_in_parallel(query=query_smiles, targets=target_smiles_list)


def test_calculate_similarity_invalid_target_smiles() -> None:
    query_smiles = "CCCC"
    target_smiles_list = ["invalid_smiles"]
    result = calculate_similarity_in_parallel(query=query_smiles, targets=target_smiles_list)
    assert result == {"invalid_smiles": None}


def test_calculate_similarity_duplicate_targets() -> None:
    query_smiles = "CCO"
    target_smiles_list = ["CCO", "CCO", "CCO"]
    result = calculate_similarity_in_parallel(query=query_smiles, targets=target_smiles_list)
    assert result == {"CCO": np.float16(1.0)}


def test_calculate_similarity_with_valid_and_invalid_targets() -> None:
    query_smiles = "CCO"
    target_smiles_list = ["CCO", "invalid_smiles", "CCC"]
    result = calculate_similarity_in_parallel(query=query_smiles, targets=target_smiles_list)
    assert result == {
        "CCO": np.float16(1.0),
        "invalid_smiles": None,
        "CCC": np.float16(0.4285),
    }


def test_calculate_similarity_with_oe_molecules() -> None:
    query = oemol_from_smiles("CCO")
    target = oemol_from_smiles("CCO")
    result = calculate_similarity_in_parallel(query=query, targets=[target])
    assert result == {target: np.float16(1.0)}


def test_calculate_similarity_with_oe_molecules_and_smiles() -> None:
    query = oemol_from_smiles("CCO")
    target_smiles_list = ["CCO", "CCC", "CCN"]
    result = calculate_similarity_in_parallel(query=query, targets=target_smiles_list)
    assert result == {
        "CCO": np.float16(1.0),
        "CCC": np.float16(0.4285),
        "CCN": np.float16(0.3333),
    }


def test_calculate_similarity_with_oe_molecules_and_invalid_smiles() -> None:
    query = oemol_from_smiles("CCO")
    target_smiles_list = ["invalid_smiles"]
    result = calculate_similarity_in_parallel(query=query, targets=target_smiles_list)
    assert result == {"invalid_smiles": None}


def test_get_most_similar_molecules() -> None:
    query_smiles = "CCO"
    target_smiles_list = ["CCO", "CCC", "CCN"]
    result = get_most_similar_molecules(
        query=query_smiles, targets=target_smiles_list, num_results=1
    )
    assert result == [SimilarMolecule(molecule="CCO", similarity=np.float16(1.0))]


def test_get_most_similar_molecules_multiple() -> None:
    query_smiles = "CCO"
    target_smiles_list = ["CCO", "CCC", "CCN"]
    result = get_most_similar_molecules(
        query=query_smiles, targets=target_smiles_list, num_results=2
    )
    assert result == [
        SimilarMolecule(molecule="CCO", similarity=np.float16(1.0)),
        SimilarMolecule(molecule="CCC", similarity=np.float16(0.4285)),
    ]


def test_get_most_similar_molecules_no_targets() -> None:
    query_smiles = "CCO"
    target_smiles_list: list[str] = []
    result: list[SimilarMolecule[str]] = get_most_similar_molecules(
        query=query_smiles, targets=target_smiles_list, num_results=1
    )
    assert result == []


def test_get_most_similar_molecules_invalid_smiles() -> None:
    query_smiles = "invalid_smiles"
    target_smiles_list = ["CCCC"]
    with pytest.raises(InvalidSmilesError):
        get_most_similar_molecules(query=query_smiles, targets=target_smiles_list, num_results=1)


def test_get_most_similar_molecules_invalid_target_smiles() -> None:
    query_smiles = "CCCC"
    target_smiles_list = ["invalid_smiles"]
    result = get_most_similar_molecules(
        query=query_smiles, targets=target_smiles_list, num_results=1
    )
    assert result == []


def test_get_most_similar_molecules_duplicate_targets() -> None:
    query_smiles = "CCO"
    target_smiles_list = ["CCO", "CCO", "CCO"]
    result = get_most_similar_molecules(
        query=query_smiles, targets=target_smiles_list, num_results=1
    )
    assert result == [SimilarMolecule(molecule="CCO", similarity=np.float16(1.0))]


def test_get_most_similar_molecules_top_n_results_returns_min_of_n_and_length() -> None:
    query_smiles = "CCO"
    target_smiles_list = ["CCO", "CCC", "CCN"]
    result = get_most_similar_molecules(
        query=query_smiles, targets=target_smiles_list, num_results=100
    )
    assert result == [
        SimilarMolecule(molecule="CCO", similarity=np.float16(1.0)),
        SimilarMolecule(molecule="CCC", similarity=np.float16(0.4285)),
        SimilarMolecule(molecule="CCN", similarity=np.float16(0.3333)),
    ]


def test_get_most_similar_molecules_with_negative_top_n_results() -> None:
    query_smiles = "CCO"
    target_smiles_list = ["CCO", "CCC", "CCN"]
    with pytest.raises(InvalidNumResultsError):
        get_most_similar_molecules(query=query_smiles, targets=target_smiles_list, num_results=-1)


def test_get_most_similar_molecules_with_oe_molecules() -> None:
    query = oemol_from_smiles("CCO")
    target = oemol_from_smiles("CCO")
    result = get_most_similar_molecules(query=query, targets=[target], num_results=1)
    assert result == [SimilarMolecule(molecule=target, similarity=np.float16(1.0))]


def test_get_most_similar_molecules_with_oe_molecules_and_smiles() -> None:
    query = oemol_from_smiles("CCO")
    target_smiles_list = ["CCO", "CCC", "CCN"]
    result = get_most_similar_molecules(query=query, targets=target_smiles_list, num_results=1)
    assert result == [SimilarMolecule(molecule="CCO", similarity=np.float16(1.0))]


def test_get_most_similar_molecules_with_oe_molecules_and_invalid_smiles() -> None:
    query = oemol_from_smiles("CCO")
    target_smiles_list = ["invalid_smiles"]
    result = get_most_similar_molecules(query=query, targets=target_smiles_list, num_results=1)
    assert result == []


def test_get_most_similar_molecules_with_all_precalculated_fps() -> None:
    query_smiles = "CCO"
    target_smiles_list = ["CCO", "CCC", "CCN"]
    target_fps = {
        "CCO": _fingerprint_from_mol("CCO"),
        "CCC": _fingerprint_from_mol("CCC"),
        "CCN": _fingerprint_from_mol("CCN"),
    }

    with patch(
        "chemutils.molecule.similarity._fingerprint_from_mol", wraps=_fingerprint_from_mol
    ) as mock_fingerprint_fn:
        result = get_most_similar_molecules(
            query=query_smiles, targets=target_smiles_list, num_results=1, fingerprints=target_fps
        )

        assert result == [SimilarMolecule(molecule="CCO", similarity=np.float16(1.0))]

        mock_fingerprint_fn.assert_called_once()  # for the query molecule


def test_get_most_similar_molecules_with_some_precalculated_fps() -> None:
    query_smiles = "CCO"
    target_smiles_list = ["CCO", "CCC", "CCN"]
    target_fps = {
        "CCO": _fingerprint_from_mol("CCO"),
        "CCN": _fingerprint_from_mol("CCN"),
    }

    result = get_most_similar_molecules(
        query=query_smiles, targets=target_smiles_list, num_results=1, fingerprints=target_fps
    )

    assert result == [SimilarMolecule(molecule="CCO", similarity=np.float16(1.0))]


def test_get_most_similar_molecules_with_flipped_fingerprints() -> None:
    # Prove the provided fingeprints are being used
    query_smiles = "CCO"
    target_smiles_list = ["CCO", "CCC", "CCN"]
    target_fps = {
        "CCO": _fingerprint_from_mol("CCC"),
        "CCC": _fingerprint_from_mol("CCN"),
        "CCN": _fingerprint_from_mol("CCO"),
    }

    result = get_most_similar_molecules(
        query=query_smiles, targets=target_smiles_list, num_results=1, fingerprints=target_fps
    )

    assert result == [SimilarMolecule(molecule="CCN", similarity=np.float16(1.0))]
