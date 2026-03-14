import pytest

from chemutils.molecule import oemol_from_smiles
from chemutils.molecule.graph.automorphism_group.nauty import (
    generate_dreadnaut_input_from_oemol,
    parse_cycles,
)


def test_generate_dreadnaut_input_from_oemol():
    assert generate_dreadnaut_input_from_oemol(
        oemol_from_smiles("c1ccccc1CCNC(=O)O")
    ).splitlines() == [
        "As",
        "l=0",
        "n=12 g",
        " 0 : 5 1 ;",
        " 1 : 2 ;",
        " 2 : 3 ;",
        " 3 : 4 ;",
        " 4 : 5 ;",
        " 5 : 6 ;",
        " 6 : 7 ;",
        " 7 : 8 ;",
        " 8 : 9 ;",
        " 9 : 10 11 ;",
        " 10 :  ;",
        " 11 :  ;",
        "f=[0,1,2,3,4,5|8|10,11|6,7,9]",
        "x",
    ]


@pytest.mark.parametrize(
    ["cycle_str", "cycles"],
    [
        ("(1 2 3)", ((1, 2, 3),)),
        ("(1 2 3)(4 5)", ((1, 2, 3), (4, 5))),
        ("(1 2 3)(4 5)(6 7 8)", ((1, 2, 3), (4, 5), (6, 7, 8))),
    ],
)
def test_parse_cycles(cycle_str, cycles):
    assert parse_cycles(cycle_str) == cycles
