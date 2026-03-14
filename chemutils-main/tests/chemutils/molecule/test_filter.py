import pytest

from chemutils.molecule.filter import get_medchem_filter_alerts


@pytest.mark.parametrize(
    ["smiles", "expected"],
    [
        (
            "CC=CC=C",
            {
                "num_red_alerts": 1,
                "num_amber_alerts": 1,
                "num_pains_alerts": 0,
                "alerts": [
                    {"rule": "polyene", "num_matches": 1, "category": "red"},
                    {"rule": "terminal_vinyl", "num_matches": 1, "category": "amber"},
                ],
            },
        ),
        (
            "C(O)C(O)C(O)",
            {
                "num_red_alerts": 1,
                "num_amber_alerts": 0,
                "num_pains_alerts": 0,
                "alerts": [{"rule": "3_aliphatic_OH", "num_matches": 1, "category": "red"}],
            },
        ),
        (
            "C(O)C(O)C(O)C(O)",
            {
                "num_red_alerts": 1,
                "num_amber_alerts": 0,
                "num_pains_alerts": 0,
                "alerts": [{"rule": "3_aliphatic_OH", "num_matches": 1, "category": "red"}],
            },
        ),
        (
            "Nc(cc1)ccc1Oc1cnccc1",
            {
                "num_red_alerts": 0,
                "num_amber_alerts": 1,
                "num_pains_alerts": 1,
                "alerts": [
                    {"rule": "aniline", "num_matches": 1, "category": "amber"},
                    {"rule": "anil_no_alk", "num_matches": 1, "category": "pains"},
                ],
            },
        ),
    ],
)
def test_get_medchem_filter_alerts(smiles, expected):
    result = get_medchem_filter_alerts(smiles)
    assert {
        "num_red_alerts": result["num_red_alerts"],
        "num_amber_alerts": result["num_amber_alerts"],
        "num_pains_alerts": result["num_pains_alerts"],
        "alerts": [
            {
                "rule": alert["rule"],
                "num_matches": len(alert["matches"]),
                "category": alert["category"],
            }
            for alert in result["alerts"]
        ],
    } == expected
