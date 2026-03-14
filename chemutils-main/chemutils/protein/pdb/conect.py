import itertools
from collections import Counter
from typing import TypedDict


class CONECTLine(TypedDict):
    atom_serial: int
    bonded_serials: list[int]


def parse_conect_line(line: str, /) -> CONECTLine:
    """
    Parse a PDB line starting with CONECT.

    The format of a CONECT line is described at
    https://www.wwpdb.org/documentation/file-format-content/format33/sect10.html#CONECT.

    Args:
        line: PDB line starting CONECT.

    Returns:
        Dictionary with the atom serial number and a list of bonded serial numbers.
    """
    if len(line) < 31:
        line = line.ljust(31)
    if line[:6] != "CONECT":
        raise ValueError("Line is not `CONECT` record.")
    return {
        "atom_serial": int(line[6:11]),
        "bonded_serials": [int(num) for i in range(11, 27, 5) if (num := line[i : i + 5].strip())],
    }


def get_bonds_from_conect_records(lines: str | list[str], /) -> dict[tuple[int, int], int]:
    """
    Get the bond pairs and bond orders from a list of CONECT records.

    Repeated specification of a bonded serial number is used to indicate a bond order.

    Args:
        lines: Either a PDB string or a list of CONECT records. Non-CONECT records are ignored.

    Returns:
        Dictionary mapping bond pairs (as sorted tuples of the two serial numbers) to the bond order.
    """
    if isinstance(lines, str):
        lines = lines.splitlines()
    conects = [parse_conect_line(line) for line in lines if line.startswith("CONECT")]
    bonds = {}
    for atom_serial, conect_grouped in itertools.groupby(conects, key=lambda x: x["atom_serial"]):
        bonded_serials = Counter(
            bonded_serial for conect in conect_grouped for bonded_serial in conect["bonded_serials"]
        )

        for bonded_serial, bond_order in bonded_serials.items():
            bond_name = tuple(sorted([atom_serial, bonded_serial]))
            bonds[bond_name] = bond_order
    return bonds
