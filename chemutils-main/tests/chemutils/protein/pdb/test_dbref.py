from pathlib import Path

import fsutils as fs
import pytest

from chemutils.protein.pdb.dbref import (
    DBREF1Line,
    DBREF2Line,
    DBREFLine,
    get_dbref1_line,
    get_dbref2_line,
    get_dbref_line,
    get_dbref_lines,
    parse_dbref1_line,
    parse_dbref2_line,
    parse_dbref_line,
    parse_dbref_lines,
)

CASES = fs.read_yaml(Path(__file__).parent / "dbref.yaml")


DBREF_ENTRY = DBREFLine(
    id_code="ABCD",
    chain_id="A",
    seq_begin=1234,
    insert_begin="A",
    seq_end=6789,
    insert_end="B",
    database="ABCDEF",
    db_accession="ABCDEFGH",
    db_id_code="ABCDEFGHIJKL",
    db_seq_begin=12345,
    db_insert_begin="A",
    db_seq_end=67890,
    db_insert_end="B",
)
DBREF_LINE = "DBREF  ABCD A 1234A 6789B ABCDEF ABCDEFGH ABCDEFGHIJKL 12345A 67890B            "

DBREF_ENTRY_NO_INS = DBREFLine(
    id_code="ABCD",
    chain_id="A",
    seq_begin=1234,
    insert_begin="A",
    seq_end=6789,
    insert_end="B",
    database="ABCDEF",
    db_accession="ABCDEFGH",
    db_id_code="ABCDEFGHIJKL",
    db_seq_begin=12345,
    db_insert_begin=None,
    db_seq_end=67890,
    db_insert_end=None,
)

DBREF1_ENTRY = DBREF1Line(
    id_code="ABCD",
    chain_id="A",
    seq_begin=1234,
    insert_begin="A",
    seq_end=6789,
    insert_end="B",
    database="ABCDEF",
    db_id_code="ABCDEFGHIJKLMNOPQRST",
)

DBREF2_ENTRY = DBREF2Line(
    id_code="ABCD",
    chain_id="A",
    db_accession="ABCDEFGHIJKLMNOPQRSTUV",
    db_seq_begin=1234567890,
    db_seq_end=6789012345,
)


@pytest.mark.parametrize(
    ["lines", "entries"],
    [(case["lines"], case["entries"]) for case in CASES],
)
def test_parse_dbref_line(lines, entries):
    assert [line for entry in entries for line in get_dbref_lines(**entry)] == lines


@pytest.mark.parametrize(
    ["lines", "entries"],
    [(case["lines"], case["entries"]) for case in CASES],
)
def test_parse_dbref_lines(lines, entries):
    assert parse_dbref_lines(lines) == entries


def test_get_dbref_line():
    assert get_dbref_line(**DBREF_ENTRY) == DBREF_LINE


def test_parse_dbref_line_full():
    assert parse_dbref_line(DBREF_LINE) == DBREF_ENTRY


def test_get_dbref_lines_long_but_has_db_inserts():
    with pytest.raises(ValueError) as exc:
        _ = get_dbref_lines(**(DBREF_ENTRY | {"db_accession": "ABCDEFGHIJ"}))
    assert (
        str(exc.value) == "Cannot have `db_insert_begin` or `db_insert_end` if line must be DBREF1."
    )


def test_get_dbref_lines_long_db_accession():
    # Check that `db_accession` longer than 8 gives DBREF1/DBREF2
    assert get_dbref_lines(**(DBREF_ENTRY_NO_INS | {"db_accession": "ABCDEFGHIJ"})) == [
        "DBREF1 ABCD A 1234A 6789B ABCDEF               ABCDEFGHIJKL                     ",
        "DBREF2 ABCD A     ABCDEFGHIJ                      12345       67890             ",
    ]


def test_get_dbref_lines_long_db_id_code():
    # Check that `db_id_code` longer than 12 gives DBREF1/DBREF2
    assert get_dbref_lines(**(DBREF_ENTRY_NO_INS | {"db_id_code": "ABCDEFGHIJKLMNOPQRST"})) == [
        "DBREF1 ABCD A 1234A 6789B ABCDEF               ABCDEFGHIJKLMNOPQRST             ",
        "DBREF2 ABCD A     ABCDEFGH                        12345       67890             ",
    ]


def test_get_dbref_lines_long_db_seq_begin():
    # Check that `db_seq_begin` longer than 5 gives DBREF1/DBREF2
    assert get_dbref_lines(**(DBREF_ENTRY_NO_INS | {"db_seq_begin": 1234567890})) == [
        "DBREF1 ABCD A 1234A 6789B ABCDEF               ABCDEFGHIJKL                     ",
        "DBREF2 ABCD A     ABCDEFGH                   1234567890       67890             ",
    ]


def test_get_dbref_lines_long_db_seq_end():
    # Check that `db_seq_end` longer than 5 gives DBREF1/DBREF2
    assert get_dbref_lines(**(DBREF_ENTRY_NO_INS | {"db_seq_end": 5679012345})) == [
        "DBREF1 ABCD A 1234A 6789B ABCDEF               ABCDEFGHIJKL                     ",
        "DBREF2 ABCD A     ABCDEFGH                        12345  5679012345             ",
    ]


@pytest.mark.parametrize(
    ["long_args", "exception_message"],
    [
        ({"id_code": "ABCDEF"}, "`id_code` must be less than 5 characters."),
        ({"chain_id": "ABC"}, "`chain_id` must be a single character."),
        ({"seq_begin": 1234567890}, "`seq_begin` must be less than 5 characters."),
        ({"insert_begin": "ABC"}, "`insert_begin` must be a single character."),
        ({"seq_end": 1234567890}, "`seq_end` must be less than 5 characters."),
        ({"insert_end": "ABC"}, "`insert_end` must be a single character."),
        ({"database": "ABCDEFGHIJK"}, "`database` must be less than 7 characters."),
        ({"db_accession": "ABCDEFGHIJK"}, "`db_accession` must be less than 9 characters."),
        ({"db_id_code": "ABCDEFGHIJKLM"}, "`db_id_code` must be less than 13 characters."),
        ({"db_seq_begin": 12345678901}, "`db_seq_begin` must be less than 6 characters."),
        ({"db_insert_begin": "ABC"}, "`db_insert_begin` must be a single character."),
        ({"db_seq_end": 12345678901}, "`db_seq_end` must be less than 6 characters."),
        ({"db_insert_end": "ABC"}, "`db_insert_end` must be a single character."),
    ],
)
def test_get_dbref_line_too_long(long_args, exception_message):
    with pytest.raises(ValueError) as exc:
        _ = get_dbref_line(**(DBREF_ENTRY_NO_INS | long_args))
    assert str(exc.value) == exception_message


@pytest.mark.parametrize(
    ["long_args", "exception_message"],
    [
        ({"id_code": "ABCDEF"}, "`id_code` must be less than 5 characters."),
        ({"chain_id": "ABC"}, "`chain_id` must be a single character."),
        ({"seq_begin": 1234567890}, "`seq_begin` must be less than 5 characters."),
        ({"insert_begin": "ABC"}, "`insert_begin` must be a single character."),
        ({"seq_end": 1234567890}, "`seq_end` must be less than 5 characters."),
        ({"insert_end": "ABC"}, "`insert_end` must be a single character."),
        ({"database": "ABCDEFGHIJK"}, "`database` must be less than 7 characters."),
        (
            {"db_id_code": "ABCDEFGHIJKLMNOPQRSTUVWX"},
            "`db_id_code` must be less than 21 characters.",
        ),
    ],
)
def test_get_dbref1_line_too_long(long_args, exception_message):
    with pytest.raises(ValueError) as exc:
        _ = get_dbref1_line(**(DBREF1_ENTRY | long_args))
    assert str(exc.value) == exception_message


@pytest.mark.parametrize(
    ["long_args", "exception_message"],
    [
        ({"id_code": "ABCDEF"}, "`id_code` must be less than 5 characters."),
        ({"chain_id": "ABC"}, "`chain_id` must be a single character."),
        (
            {"db_accession": "ABCDEFGHIJKLMOPQRSTUVWXYZ"},
            "`db_accession` must be less than 23 characters.",
        ),
        ({"db_seq_begin": 12345678901234}, "`db_seq_begin` must be less than 11 characters."),
        ({"db_seq_end": 12345678901234}, "`db_seq_end` must be less than 11 characters."),
    ],
)
def test_get_dbref2_line_too_long(long_args, exception_message):
    with pytest.raises(ValueError) as exc:
        _ = get_dbref2_line(**(DBREF2_ENTRY | long_args))
    assert str(exc.value) == exception_message


def test_parse_dbref_line_dbref1_is_not_dbref():
    with pytest.raises(ValueError) as exc:
        _ = parse_dbref_line("DBREF1      ")
    assert str(exc.value) == "Line is not `DBREF` record."


def test_parse_dbref_line_not_dbref():
    with pytest.raises(ValueError) as exc:
        _ = parse_dbref_line("NOTDBR      ")
    assert str(exc.value) == "Line is not `DBREF` record."


def test_parse_dbref_line_not_dbref1():
    with pytest.raises(ValueError) as exc:
        _ = parse_dbref1_line("NOTDBR      ")
    assert str(exc.value) == "Line is not `DBREF1` record."


def test_parse_dbref2_line_not_dbref1():
    with pytest.raises(ValueError) as exc:
        _ = parse_dbref2_line("NOTDBR      ")
    assert str(exc.value) == "Line is not `DBREF2` record."


def test_blank_dbref_line():
    assert parse_dbref_line("DBREF ") == DBREFLine(
        id_code=None,
        chain_id=None,
        seq_begin=None,
        insert_begin=None,
        seq_end=None,
        insert_end=None,
        database=None,
        db_accession=None,
        db_id_code=None,
        db_seq_begin=None,
        db_insert_begin=None,
        db_seq_end=None,
        db_insert_end=None,
    )


def test_blank_dbref1_line():
    assert parse_dbref1_line("DBREF1") == DBREF1Line(
        id_code=None,
        chain_id=None,
        seq_begin=None,
        insert_begin=None,
        seq_end=None,
        insert_end=None,
        database=None,
        db_id_code=None,
    )


def test_blank_dbref2_line():
    assert parse_dbref2_line("DBREF2") == DBREF2Line(
        id_code=None,
        chain_id=None,
        db_accession=None,
        db_seq_begin=None,
        db_seq_end=None,
    )
