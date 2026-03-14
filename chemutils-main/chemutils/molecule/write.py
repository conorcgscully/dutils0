from openeye import oechem

from .format import setup_stream_format


def write_molecule_bytes(
    *oemols: oechem.OEGraphMol,
    format: str | int,
    flavor: int | None = None,
    gzip: bool = False,
) -> bytes:
    """
    Write one or more OpenEye molecules to bytes in memory.

    Args:
        oemols: One or more OpenEye molecules to write.
        format: Either an `OEFlavor_` from `oechem`, or a string such as `PDB`.
        flavor: Optional flavor from OpenEye, modifying the format.
        gzip: Should the file be written in a gzip format.

    Returns:
        Bytes representing the molecule in the specified file format.
    """
    oms = oechem.oemolostream()
    setup_stream_format(oms, format=format, flavor=flavor)
    oms.Setgz(gzip)
    oms.openstring()
    for oemol in oemols:
        oechem.OEWriteConstMolecule(oms, oemol)
    return oms.GetString()  # type: ignore


def write_molecule_str(
    *oemols: oechem.OEGraphMol, format: str | int, flavor: int | None = None
) -> str:
    """
    Write one or more OpenEye molecules to string in memory using UTF-8.

    Args:
        oemols: One or more OpenEye molecules to write.
        format: Either an `OEFlavor_` from `oechem`, or a string such as `PDB`.
        flavor: Optional flavor from OpenEye, modifying the format.

    Returns:
        String representing the molecule in the specified file format.
    """
    return write_molecule_bytes(*oemols, format=format, flavor=flavor).decode("UTF-8")
