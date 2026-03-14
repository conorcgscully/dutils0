from typing import TypedDict


class SEQADVLine(TypedDict):
    id_code: str | None
    res_name: str | None
    chain_id: str | None
    seq_num: int | None
    insert_code: str | None
    database: str | None
    db_accession: str | None
    db_res_name: str | None
    db_seq_num: int | None
    conflict: str | None


def get_seqadv_line(
    *,
    id_code: str | None = None,
    res_name: str | None = None,
    chain_id: str | None = None,
    seq_num: int | None = None,
    insert_code: str | None = None,
    database: str | None = None,
    db_accession: str | None = None,
    db_res_name: str | None = None,
    db_seq_num: int | None = None,
    conflict: str | None = None,
) -> str:
    return f"SEQADV {id_code or '':<4.4s} {res_name or '':<3.3s} {chain_id or '':<1.1s} {seq_num if seq_num is not None else ''!s:>4.4s}{insert_code or '':<1.1s} {database or '':<4.4s} {db_accession or '':<9.9s} {db_res_name or ''!s:<3.3s} {db_seq_num if db_seq_num is not None else ''!s:>5.5s} {conflict or '':<21.21s}          "


def parse_seqadv_line(line: str, /) -> SEQADVLine:
    if len(line) < 80:
        line = line.ljust(80)
    return {
        "id_code": line[7:11].strip() or None,
        "res_name": line[12:15].strip() or None,
        "chain_id": line[16].strip() or None,
        "seq_num": int(num) if (num := line[18:22].strip()) else None,
        "insert_code": line[22].strip() or None,
        "database": line[24:28].strip() or None,
        "db_accession": line[29:38].strip() or None,
        "db_res_name": line[39:42].strip() or None,
        "db_seq_num": int(num) if (num := line[43:48].strip()) else None,
        "conflict": line[49:70].strip() or None,
    }
