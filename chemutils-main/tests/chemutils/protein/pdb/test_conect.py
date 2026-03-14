import pytest

from chemutils.protein.pdb.conect import (
    CONECTLine,
    get_bonds_from_conect_records,
    parse_conect_line,
)


@pytest.mark.parametrize(
    ["line", "expected"],
    [
        (
            "CONECT 1179  746 1184 1195 1203                    ",
            CONECTLine(atom_serial=1179, bonded_serials=[746, 1184, 1195, 1203]),
        ),
        (
            "CONECT 1179 1211 1222                              ",
            CONECTLine(atom_serial=1179, bonded_serials=[1211, 1222]),
        ),
        (
            "CONECT 1021  544 1017 1020 1022",
            CONECTLine(atom_serial=1021, bonded_serials=[544, 1017, 1020, 1022]),
        ),
        (
            "CONECT2008220078200882010320106",
            CONECTLine(atom_serial=20082, bonded_serials=[20078, 20088, 20103, 20106]),
        ),
    ],
)
def test_parse_conect_line(line, expected):
    assert parse_conect_line(line) == expected


@pytest.mark.parametrize(
    ["lines", "expected"],
    [
        (
            [
                "CONECT 1179  746 1184 1195 1203",
                "CONECT 1179 1211 1222",
                "CONECT 1021  544 1017 1020 1022",
            ],
            {
                (746, 1179): 1,
                (1179, 1184): 1,
                (1179, 1195): 1,
                (1179, 1203): 1,
                (1179, 1211): 1,
                (1179, 1222): 1,
                (544, 1021): 1,
                (1017, 1021): 1,
                (1020, 1021): 1,
                (1021, 1022): 1,
            },
        ),
        (
            [
                "CONECT 3582 3596",
                "CONECT 3579 3582",
                "CONECT 3579 3581 3581",
                "CONECT 3579 3580 3580",
                "CONECT 3583 3597",
            ],
            {
                (3582, 3596): 1,
                (3579, 3582): 1,
                (3579, 3581): 2,
                (3579, 3580): 2,
                (3583, 3597): 1,
            },
        ),
        (
            [
                "CONECT2008220078200882010320106",
                "CONECT2008320080200842010020101",
                "CONECT200842007520083",
            ],
            {
                (20078, 20082): 1,
                (20082, 20088): 1,
                (20082, 20103): 1,
                (20082, 20106): 1,
                (20080, 20083): 1,
                (20083, 20084): 1,
                (20083, 20100): 1,
                (20083, 20101): 1,
                (20075, 20084): 1,
            },
        ),
    ],
)
def test_get_bonds_from_conect_records(lines, expected):
    assert get_bonds_from_conect_records(lines) == expected
