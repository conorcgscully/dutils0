from __future__ import annotations

from collections.abc import Generator

from openeye import oechem

from .format import setup_stream_format


class ReadMoleculeError(ValueError):
    index: int | None = None

    def __init__(self, *args: object, index: int | None = None) -> None:
        super().__init__(*args)
        self.index = index


def read_molecules_bytes(
    *bytes_list: bytes, format: str | int, flavor: int | None = None, gzip: bool = False
) -> Generator[oechem.OEGraphMol, None, None]:
    """
    Read OpenEye molecules from one or more byte strings in memory.

    This can both read multiple molecules from a byte string, and read from multiple
    byte strings sequentially.

    If no molecules can be read of a file, a `ValueError` is raised. This is because
    OpenEye provides no way of knowing if a file is invalid.

    Args:
        bytes_list: One more more raw bytes representing molecular files.
        format: Either an `OEFlavor_` from `oechem`, or a string such as `PDB`.
        flavor: Optional flavor from OpenEye, modifying the format.
        gzip: Do the bytes represent a gzipped file.

    Yields:
        OpenEye molecules.

    Raises:
        ValueError: Failed to read molecule.
    """
    for bytes in bytes_list:
        ifs = oechem.oemolistream()
        setup_stream_format(ifs, format=format, flavor=flavor)
        ifs.Setgz(gzip)
        ifs.openstring(bytes)

        last_successful_index = -1

        for index, mol in enumerate(ifs.GetOEGraphMols()):
            if ifs.GetReadFailureCt() > 0:
                raise ReadMoleculeError(
                    f"Failed to read molecule at index {last_successful_index + 1}.",
                    index=last_successful_index + 1,
                )
            last_successful_index = index
            # Copy because OpenEye reuses the same `mol` object in the generator.
            yield mol.CreateCopy()

        if ifs.GetReadFailureCt() > 0:
            raise ReadMoleculeError(
                f"Failed to read molecule at index {last_successful_index + 1}.",
                index=last_successful_index + 1,
            )

        # Sometimes, OpenEye just yields no molecules when the string is not blank, with no errors
        if len(bytes) > 0 and last_successful_index == -1:
            raise ReadMoleculeError("Failed to read molecule.")


def read_molecule_bytes(
    bytes: bytes, /, *, format: str | int, flavor: int | None = None, gzip: bool = False
) -> oechem.OEGraphMol:
    """
    Read an OpenEye molecule from bytes in memory.

    Args:
        bytes: Raw bytes representing a molecular file.
        format: Either an `OEFlavor_` from `oechem`, or a string such as `PDB`.
        flavor: Optional flavor from OpenEye, modifying the format.
        gzip: Do the bytes represent a gzipped file.

    Returns:
        Parsed OpenEye molecule.

    Raises:
        ValueError: Failed to read molecule.
    """
    try:
        return next(read_molecules_bytes(bytes, format=format, flavor=flavor, gzip=gzip))
    except ReadMoleculeError:
        raise ReadMoleculeError("Failed to read molecule.") from None
    except StopIteration:
        raise ReadMoleculeError("Failed to read molecule.") from None


def read_molecules_str(
    *strings: str, format: str | int, flavor: int | None = None
) -> Generator[oechem.OEGraphMol, None, None]:
    """
    Read OpenEye molecules from one or more strings in memory.

    This can both read multiple molecules from a string, and read from multiple strings
    sequentially.

    Args:
        strings: One or more strings representing molecular files.
        format: Either an `OEFlavor_` from `oechem`, or a string such as `PDB`.
        flavor: Optional flavor from OpenEye, modifying the format.

    Yields:
        OpenEye molecules.
    """
    for string in strings:
        yield from read_molecules_bytes(string.encode("UTF-8"), format=format, flavor=flavor)


def read_molecule_str(
    str: str, /, *, format: str | int, flavor: int | None = None
) -> oechem.OEGraphMol:
    """
    Read an OpenEye molecule from a string in memory.

    Args:
        str: String representing a molecular file.
        format: Either an `OEFlavor_` from `oechem`, or a string such as `PDB`.
        flavor: Optional flavor from OpenEye, modifying the format.

    Returns:
        Parsed OpenEye molecule.
    """
    return read_molecule_bytes(str.encode("UTF-8"), format=format, flavor=flavor)
