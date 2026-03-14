"""
Handles DBREF, DBREF1 and DBREF2 records for cross-referencing sequences with reference databases.

https://www.wwpdb.org/documentation/file-format-content/format33/sect3.html
"""

from typing import TypedDict

DBREF_FORMAT = "DBREF  {id_code:<4.4s} {chain_id:>1.1s} {seq_begin!s:>4.4s}{insert_begin:>1.1s} {seq_end!s:>4.4s}{insert_end:>1.1s} {database:<6.6s} {db_accession:<8.8s} {db_id_code:<12.12s} {db_seq_begin!s:>5.5s}{db_insert_begin:>1.1s} {db_seq_end!s:>5.5s}{db_insert_end:>1.1s}            "
DBREF1_FORMAT = "DBREF1 {id_code:<4.4s} {chain_id:>1.1s} {seq_begin!s:>4.4s}{insert_begin:>1.1s} {seq_end!s:>4.4s}{insert_end:>1.1s} {database:<6.6s}               {db_id_code:<20.20}             "
DBREF2_FORMAT = "DBREF2 {id_code:<4.4s} {chain_id:>1.1s}     {db_accession:<22.22s}     {db_seq_begin!s:>10.10s}  {db_seq_end!s:>10.10s}             "


class DBREFLine(TypedDict):
    """Representation of a single DBREF line, or a pair of DBREF1/DBREF2 lines."""

    id_code: str | None
    chain_id: str | None
    seq_begin: int | None
    insert_begin: str | None
    seq_end: int | None
    insert_end: str | None
    database: str | None
    db_accession: str | None
    db_id_code: str | None
    db_seq_begin: int | None
    db_insert_begin: str | None
    db_seq_end: int | None
    db_insert_end: str | None


class DBREF1Line(TypedDict):
    """Representation of a single DBREF1 line."""

    id_code: str | None
    chain_id: str | None
    seq_begin: int | None
    insert_begin: str | None
    seq_end: int | None
    insert_end: str | None
    database: str | None
    db_id_code: str | None


class DBREF2Line(TypedDict):
    """Representation of a single DBREF2 line."""

    id_code: str | None
    chain_id: str | None
    db_accession: str | None
    db_seq_begin: int | None
    db_seq_end: int | None


def get_dbref_line(
    *,
    id_code: str | None = None,
    chain_id: str | None = None,
    seq_begin: int | None = None,
    insert_begin: str | None = None,
    seq_end: int | None = None,
    insert_end: str | None = None,
    database: str | None = None,
    db_accession: str | None = None,
    db_id_code: str | None = None,
    db_seq_begin: int | None = None,
    db_insert_begin: str | None = None,
    db_seq_end: int | None = None,
    db_insert_end: str | None = None,
) -> str:
    """Generate a DBREF line from its various components."""
    if id_code and len(id_code) > 4:
        raise ValueError("`id_code` must be less than 5 characters.")
    if chain_id and len(chain_id) > 1:
        raise ValueError("`chain_id` must be a single character.")
    if seq_begin and len(str(seq_begin)) > 4:
        raise ValueError("`seq_begin` must be less than 5 characters.")
    if insert_begin and len(insert_begin) > 1:
        raise ValueError("`insert_begin` must be a single character.")
    if seq_end and len(str(seq_end)) > 4:
        raise ValueError("`seq_end` must be less than 5 characters.")
    if insert_end and len(insert_end) > 1:
        raise ValueError("`insert_end` must be a single character.")
    if database and len(database) > 6:
        raise ValueError("`database` must be less than 7 characters.")
    if db_accession and len(db_accession) > 8:
        raise ValueError("`db_accession` must be less than 9 characters.")
    if db_id_code and len(db_id_code) > 12:
        raise ValueError("`db_id_code` must be less than 13 characters.")
    if db_seq_begin and len(str(db_seq_begin)) > 5:
        raise ValueError("`db_seq_begin` must be less than 6 characters.")
    if db_insert_begin and len(db_insert_begin) > 1:
        raise ValueError("`db_insert_begin` must be a single character.")
    if db_seq_end and len(str(db_seq_end)) > 5:
        raise ValueError("`db_seq_end` must be less than 6 characters.")
    if db_insert_end and len(db_insert_end) > 1:
        raise ValueError("`db_insert_end` must be a single character.")
    return DBREF_FORMAT.format(
        id_code=id_code or "",
        chain_id=chain_id or "",
        seq_begin=seq_begin if seq_begin is not None else "",
        insert_begin=insert_begin or "",
        seq_end=seq_end if seq_end is not None else "",
        insert_end=insert_end or "",
        database=database or "",
        db_accession=db_accession or "",
        db_id_code=db_id_code or "",
        db_seq_begin=db_seq_begin if db_seq_begin is not None else "",
        db_insert_begin=db_insert_begin or "",
        db_seq_end=db_seq_end if db_seq_end is not None else "",
        db_insert_end=db_insert_end or "",
    )


def get_dbref1_line(
    *,
    id_code: str | None = None,
    chain_id: str | None = None,
    seq_begin: int | None = None,
    insert_begin: str | None = None,
    seq_end: int | None = None,
    insert_end: str | None = None,
    database: str | None = None,
    db_id_code: str | None = None,
) -> str:
    """Generate a DBREF1 line from its various components."""
    if id_code and len(id_code) > 4:
        raise ValueError("`id_code` must be less than 5 characters.")
    if chain_id and len(chain_id) > 1:
        raise ValueError("`chain_id` must be a single character.")
    if seq_begin and len(str(seq_begin)) > 4:
        raise ValueError("`seq_begin` must be less than 5 characters.")
    if insert_begin and len(insert_begin) > 1:
        raise ValueError("`insert_begin` must be a single character.")
    if seq_end and len(str(seq_end)) > 4:
        raise ValueError("`seq_end` must be less than 5 characters.")
    if insert_end and len(insert_end) > 1:
        raise ValueError("`insert_end` must be a single character.")
    if database and len(database) > 6:
        raise ValueError("`database` must be less than 7 characters.")
    if db_id_code and len(db_id_code) > 20:
        raise ValueError("`db_id_code` must be less than 21 characters.")
    return DBREF1_FORMAT.format(
        id_code=id_code or "",
        chain_id=chain_id or "",
        seq_begin=seq_begin if seq_begin is not None else "",
        insert_begin=insert_begin or "",
        seq_end=seq_end if seq_end is not None else "",
        insert_end=insert_end or "",
        database=database or "",
        db_id_code=db_id_code or "",
    )


def get_dbref2_line(
    *,
    id_code: str | None = None,
    chain_id: str | None = None,
    db_accession: str | None = None,
    db_seq_begin: int | None = None,
    db_seq_end: int | None = None,
) -> str:
    """Generate a DBREF2 line from its various components."""
    if id_code and len(id_code) > 4:
        raise ValueError("`id_code` must be less than 5 characters.")
    if chain_id and len(chain_id) > 1:
        raise ValueError("`chain_id` must be a single character.")
    if db_accession and len(db_accession) > 22:
        raise ValueError("`db_accession` must be less than 23 characters.")
    if db_seq_begin and len(str(db_seq_begin)) > 10:
        raise ValueError("`db_seq_begin` must be less than 11 characters.")
    if db_seq_end and len(str(db_seq_end)) > 10:
        raise ValueError("`db_seq_end` must be less than 11 characters.")
    return DBREF2_FORMAT.format(
        id_code=id_code or "",
        chain_id=chain_id or "",
        db_accession=db_accession or "",
        db_seq_begin=db_seq_begin if db_seq_begin is not None else "",
        db_seq_end=db_seq_end if db_seq_end is not None else "",
    )


def get_dbref_lines(
    *,
    id_code: str | None = None,
    chain_id: str | None = None,
    seq_begin: int | None = None,
    insert_begin: str | None = None,
    seq_end: int | None = None,
    insert_end: str | None = None,
    database: str | None = None,
    db_accession: str | None = None,
    db_id_code: str | None = None,
    db_seq_begin: int | None = None,
    db_insert_begin: str | None = None,
    db_seq_end: int | None = None,
    db_insert_end: str | None = None,
) -> list[str]:
    """Generate either a DBREF line or pair of DBREF1/DBREF2 lines from its various components."""
    db_accession_long = db_accession and len(db_accession) > 8
    db_id_code_long = db_id_code and len(db_id_code) > 12
    db_seq_begin_long = db_seq_begin and len(str(db_seq_begin)) > 5
    db_seq_end_long = db_seq_end and len(str(db_seq_end)) > 5
    if db_accession_long or db_id_code_long or db_seq_begin_long or db_seq_end_long:
        if db_insert_begin or db_insert_end:
            raise ValueError(
                "Cannot have `db_insert_begin` or `db_insert_end` if line must be DBREF1."
            )
        return [
            get_dbref1_line(
                id_code=id_code,
                chain_id=chain_id,
                seq_begin=seq_begin,
                insert_begin=insert_begin,
                seq_end=seq_end,
                insert_end=insert_end,
                database=database,
                db_id_code=db_id_code,
            ),
            get_dbref2_line(
                id_code=id_code,
                chain_id=chain_id,
                db_accession=db_accession,
                db_seq_begin=db_seq_begin,
                db_seq_end=db_seq_end,
            ),
        ]
    else:
        return [
            get_dbref_line(
                id_code=id_code,
                chain_id=chain_id,
                seq_begin=seq_begin,
                insert_begin=insert_begin,
                seq_end=seq_end,
                insert_end=insert_end,
                database=database,
                db_accession=db_accession,
                db_id_code=db_id_code,
                db_seq_begin=db_seq_begin,
                db_insert_begin=db_insert_begin,
                db_seq_end=db_seq_end,
                db_insert_end=db_insert_end,
            )
        ]


def parse_dbref_line(line: str, /) -> DBREFLine:
    if len(line) != 80:
        line = line.ljust(80)
    if line[:6] != "DBREF ":
        raise ValueError("Line is not `DBREF` record.")
    return {
        "id_code": line[7:11].strip() or None,
        "chain_id": line[12].strip() or None,
        "seq_begin": int(num) if (num := line[14:18].strip()) else None,
        "insert_begin": line[18].strip() or None,
        "seq_end": int(num) if (num := line[20:24].strip()) else None,
        "insert_end": line[24].strip() or None,
        "database": line[26:32].strip() or None,
        "db_accession": line[33:41].strip() or None,
        "db_id_code": line[42:54].strip() or None,
        "db_seq_begin": int(num) if (num := line[55:60].strip()) else None,
        "db_insert_begin": line[60].strip() or None,
        "db_seq_end": int(num) if (num := line[62:67].strip()) else None,
        "db_insert_end": line[67].strip() or None,
    }


def parse_dbref1_line(line: str, /) -> DBREF1Line:
    if len(line) != 80:
        line = line.ljust(80)
    if line[:6] != "DBREF1":
        raise ValueError("Line is not `DBREF1` record.")
    return {
        "id_code": line[7:11].strip() or None,
        "chain_id": line[12].strip() or None,
        "seq_begin": int(num) if (num := line[14:18].strip()) else None,
        "insert_begin": line[18].strip() or None,
        "seq_end": int(num) if (num := line[20:24].strip()) else None,
        "insert_end": line[24].strip() or None,
        "database": line[26:32].strip() or None,
        "db_id_code": line[47:67].strip() or None,
    }


def parse_dbref2_line(line: str, /) -> DBREF2Line:
    if len(line) != 80:
        line = line.ljust(80)
    if line[:6] != "DBREF2":
        raise ValueError("Line is not `DBREF2` record.")
    return {
        "id_code": line[7:11].strip() or None,
        "chain_id": line[12].strip() or None,
        "db_accession": line[18:40].strip() or None,
        "db_seq_begin": int(num) if (num := line[45:55].strip()) else None,
        "db_seq_end": int(num) if (num := line[57:67].strip()) else None,
    }


def parse_dbref_lines(lines: list[str], /) -> list[DBREFLine]:
    """
    Parse `DBREF` records (and `DBREF1/DBREF2` pairs) from a list of lines.

    This also supports parsing PDBs directly line-by-line, as lines that aren't DBREF are ignored.

    Args:
        lines: List of lines.

    Returns:
        List of `DBREFLine` objects describing the cross references.
    """
    dbref_lines = [line for line in lines if line[:6] in {"DBREF ", "DBREF1", "DBREF2"}]
    result = []
    i = 0
    while i < len(dbref_lines):
        line = dbref_lines[i]
        if line[:6] == "DBREF ":
            result.append(parse_dbref_line(line))
        elif line[:6] == "DBREF1":
            if i == len(dbref_lines) - 1:
                raise ValueError("DBREF1 not followed by DBREF2.")
            line2 = dbref_lines[i + 1]
            i += 1
            if not line2[:6] == "DBREF2":
                raise ValueError("DBREF1 not followed by DBREF2.")
            result.append(
                {
                    **parse_dbref1_line(line),
                    **parse_dbref2_line(line2),
                    "db_insert_begin": None,
                    "db_insert_end": None,
                }
            )
        i += 1
    return result
