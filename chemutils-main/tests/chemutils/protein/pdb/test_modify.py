import pytest

from chemutils.protein.pdb import filter_lines_in_pdb, insert_lines_into_pdb
from chemutils.protein.pdb.modify import get_line_order

CASES = [
    # No new lines
    (
        ["HEADER    MY HEADER"],
        [],
        ["HEADER    MY HEADER"],
    ),
    # Insert HEADER before title
    (
        ["TITLE     SOME_TITLE"],
        ["HEADER    MY HEADER"],
        ["HEADER    MY HEADER", "TITLE     SOME_TITLE"],
    ),
    # Insert SEQRES between SEQADV and MODRES
    (
        ["SEQADV    ABCD", "MODRES    EFGH"],
        ["SEQRES    IJKL"],
        ["SEQADV    ABCD", "SEQRES    IJKL", "MODRES    EFGH"],
    ),
    # Insert SOURCE after existing SOURCE line but before KEYWDS
    (
        ["COMPND    ABCD", "SOURCE    EFGH", "KEYWDS    IJKL"],
        ["SOURCE    MNOP"],
        ["COMPND    ABCD", "SOURCE    EFGH", "SOURCE    MNOP", "KEYWDS    IJKL"],
    ),
    # Maintain order of DBREF1 and DBREF2 records
    (
        ["DBREF1    FIRST", "DBREF2    FIRST", "SEQADV    FIRST"],
        [
            "DBREF1    SECOND",
            "DBREF2    SECOND",
            "SEQADV    SECOND",
            "DBREF1    THIRD",
            "DBREF2    THIRD",
            "SEQADV    THIRD",
        ],
        [
            "DBREF1    FIRST",
            "DBREF2    FIRST",
            "DBREF1    SECOND",
            "DBREF2    SECOND",
            "DBREF1    THIRD",
            "DBREF2    THIRD",
            "SEQADV    FIRST",
            "SEQADV    SECOND",
            "SEQADV    THIRD",
        ],
    ),
    # Insert lines before the body but at the end of the header
    (
        ["MODRES    ABCD", "ATOM    ATOM_1"],
        ["LINK      EFGH"],
        ["MODRES    ABCD", "LINK      EFGH", "ATOM    ATOM_1"],
    ),
    # REMARK are sorted by number
    (
        ["REMARK   2 SECOND", "REMARK   4 FOURTH", "REMARK 100 SIXTH"],
        ["REMARK   1 FIRST", "REMARK   3 THIRD", "REMARK   4 FIFTH"],
        [
            "REMARK   1 FIRST",
            "REMARK   2 SECOND",
            "REMARK   3 THIRD",
            "REMARK   4 FOURTH",
            "REMARK   4 FIFTH",
            "REMARK 100 SIXTH",
        ],
    ),
]


@pytest.mark.parametrize("existing_lines, new_lines, expected_lines", CASES)
def test_insert_lines_into_pdb_list_str(existing_lines, new_lines, expected_lines):
    assert insert_lines_into_pdb(pdb=existing_lines, lines=new_lines) == expected_lines


@pytest.mark.parametrize("existing_lines, new_lines, expected_lines", CASES)
def test_insert_lines_into_pdb_str(existing_lines, new_lines, expected_lines):
    assert insert_lines_into_pdb(pdb="\n".join(existing_lines), lines=new_lines) == "\n".join(
        expected_lines
    )


@pytest.mark.parametrize("existing_lines, new_lines, expected_lines", CASES)
def test_insert_lines_into_pdb_bytes(existing_lines, new_lines, expected_lines):
    assert (
        insert_lines_into_pdb(pdb="\n".join(existing_lines).encode(), lines=new_lines)
        == "\n".join(expected_lines).encode()
    )


def test_insert_lines_into_pdb_not_real_record():
    with pytest.raises(ValueError) as exc:
        insert_lines_into_pdb(pdb=["HEADER    MY HEADER"], lines=["NOTREC    NOT_A_LINE"])

    assert str(exc.value) == "Cannot insert invalid PDB record `NOTREC`."


def test_insert_lines_into_pdb_header_after_body():
    with pytest.raises(ValueError) as exc:
        insert_lines_into_pdb(
            pdb=["HEADER    MY HEADER", "ATOM  SOME_ATOM", "MODRES   TO_FAR_AFTER"],
            lines=["TITLE    MY_TITLE"],
        )

    assert str(exc.value) == "PDB has header record `MODRES` following body records. Cannot insert."


def test_insert_lines_into_pdb_pdb_has_invalid_record():
    with pytest.raises(ValueError) as exc:
        insert_lines_into_pdb(pdb=["NOTREC    NOT_A_LINE"], lines=["TITLE    MY_TITLE"])

    assert str(exc.value) == "PDB has unknown record `NOTREC`."


def test_get_line_order_dbref():
    assert get_line_order("DBREF1") == get_line_order("DBREF2")


def test_get_line_order_lt():
    assert get_line_order("AUTHOR") < get_line_order("JRNL")


def test_filter_lines_in_pdb():
    pdb = ["MODRES  ABCD", "ATOM    ABCD", "END     ABCD"]
    pdb = filter_lines_in_pdb(pdb=pdb, filter=lambda line: not line.startswith("ATOM  "))
    assert pdb == ["MODRES  ABCD", "END     ABCD"]
