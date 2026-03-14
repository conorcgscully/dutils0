import pytest

from chemutils.protein.pdb.seqadv import get_seqadv_line, parse_seqadv_line


@pytest.mark.parametrize(
    ["kwargs", "expected"],
    [
        # Examples from https://www.wwpdb.org/documentation/file-format-content/format33/sect3.html#SEQADV
        (
            {
                "id_code": "3ABC",
                "res_name": "MET",
                "chain_id": "A",
                "seq_num": -1,
                "database": "UNP",
                "db_accession": "P10725",
                "conflict": "EXPRESSION TAG",
            },
            "SEQADV 3ABC MET A   -1  UNP  P10725              EXPRESSION TAG                 ",
        ),
        (
            {
                "id_code": "2QLE",
                "res_name": "CRO",
                "chain_id": "A",
                "seq_num": 66,
                "database": "UNP",
                "db_accession": "P42212",
                "db_res_name": "SER",
                "db_seq_num": 65,
                "conflict": "CHROMOPHORE",
            },
            "SEQADV 2QLE CRO A   66  UNP  P42212    SER    65 CHROMOPHORE                    ",
        ),
        (
            {
                "id_code": "2OKW",
                "res_name": "LEU",
                "chain_id": "A",
                "seq_num": 64,
                "database": "NOR",
                "db_accession": "NOR00669",
                "db_res_name": "PHE",
                "db_seq_num": 14,
                "conflict": "SEE REMARK 999",
            },
            "SEQADV 2OKW LEU A   64  NOR  NOR00669  PHE    14 SEE REMARK 999                 ",
        ),
        # Overflow
        (
            {
                "id_code": "ABCDEFG",
                "res_name": "LEUCINE",
                "chain_id": "ABC",
                "seq_num": 123456,
                "database": "DATABASE",
                "db_accession": "ACCESSIONDATABASE",
                "db_res_name": "ALANINE",
                "db_seq_num": 4567890,
                "conflict": "A VERY VERY VERY LONG COMMENT THAT DOES NOT FIT",
            },
            "SEQADV ABCD LEU A 1234  DATA ACCESSION ALA 45678 A VERY VERY VERY LONG          ",
        ),
    ],
)
def test_get_seqadv_line(kwargs, expected):
    line = get_seqadv_line(**kwargs)
    assert line == expected
    assert len(line) == 80


@pytest.mark.parametrize(
    ["parts", "line"],
    [
        # Examples from https://www.wwpdb.org/documentation/file-format-content/format33/sect3.html#SEQADV
        (
            {
                "id_code": "3ABC",
                "res_name": "MET",
                "chain_id": "A",
                "seq_num": -1,
                "insert_code": None,
                "database": "UNP",
                "db_res_name": None,
                "db_seq_num": None,
                "db_accession": "P10725",
                "conflict": "EXPRESSION TAG",
            },
            "SEQADV 3ABC MET A   -1  UNP  P10725              EXPRESSION TAG                 ",
        ),
        (
            {
                "id_code": "2QLE",
                "res_name": "CRO",
                "chain_id": "A",
                "seq_num": 66,
                "insert_code": None,
                "database": "UNP",
                "db_accession": "P42212",
                "db_res_name": "SER",
                "db_seq_num": 65,
                "conflict": "CHROMOPHORE",
            },
            "SEQADV 2QLE CRO A   66  UNP  P42212    SER    65 CHROMOPHORE                    ",
        ),
        (
            {
                "id_code": "2OKW",
                "res_name": "LEU",
                "chain_id": "A",
                "seq_num": 64,
                "insert_code": None,
                "database": "NOR",
                "db_accession": "NOR00669",
                "db_res_name": "PHE",
                "db_seq_num": 14,
                "conflict": "SEE REMARK 999",
            },
            "SEQADV 2OKW LEU A   64  NOR  NOR00669  PHE    14 SEE REMARK 999                 ",
        ),
    ],
)
def test_parse_seqadv_line(line, parts):
    assert parse_seqadv_line(line) == parts
