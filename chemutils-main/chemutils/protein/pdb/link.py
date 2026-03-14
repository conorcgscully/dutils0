from typing import TypedDict


class LINKLine(TypedDict):
    atom_name1: str | None
    alt_loc1: str | None
    res_name1: str | None
    chain_id1: str | None
    seq_num1: int | None
    insertion_code1: str | None
    atom_name2: str | None
    alt_loc2: str | None
    res_name2: str | None
    chain_id2: str | None
    seq_num2: int | None
    insertion_code2: str | None
    sym1: str | None
    sym2: str | None
    length: float | None


def parse_link_line(line: str, /) -> LINKLine:
    """
    Parse a PDB line starting with LINK.

    The format of a LINK line is described at
    https://www.wwpdb.org/documentation/file-format-content/format33/sect6.html#LINK.

    Args:
        line: PDB line starting LINK.

    Returns:
        Dictionary with the fields parsed from the line.
    """
    if len(line) < 80:
        line = line.ljust(80)
    if line[:6] != "LINK  ":
        raise ValueError("Line is not `LINK` record.")
    return {
        "atom_name1": line[12:16].strip() or None,
        "alt_loc1": line[16].strip() or None,
        "res_name1": line[17:20].strip() or None,
        "chain_id1": line[21].strip() or None,
        "seq_num1": int(num) if (num := line[22:26].strip()) else None,
        "insertion_code1": line[26].strip() or None,
        "atom_name2": line[42:46].strip() or None,
        "alt_loc2": line[46].strip() or None,
        "res_name2": line[47:50].strip() or None,
        "chain_id2": line[51].strip() or None,
        "seq_num2": int(num) if (num := line[52:56].strip()) else None,
        "insertion_code2": line[56].strip() or None,
        "sym1": line[59:65].strip() or None,
        "sym2": line[66:72].strip() or None,
        "length": float(line[73:78].strip()) if line[73:78].strip() else None,
    }
