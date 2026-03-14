from __future__ import annotations

from collections.abc import Generator
from os import PathLike
from urllib.parse import urlparse

from fsutils import read_bytes
from openeye import oechem

from chemutils.molecule import read_molecule_bytes, read_molecules_bytes
from chemutils.molecule.format import get_format_and_gzip_from_filename


def read_molecule(
    path: str | PathLike[str],
    *,
    format: int | str | None = None,
    flavor: int | None = None,
    gzip: bool | None = None,
) -> oechem.OEGraphMol:
    """
    Read an OpenEye molecule from a path.

    This path could be a local filepath, a path on S3 or a web URL.

    Args:
        path: Path with relevant schema.
        format: Optional format, either an `OEFlavor_` from `oechem`, or a string such as `PDB`.
            This is infered from the path if not provided.
        flavor: Optional flavor from OpenEye, modifying the format.
        gzip: Optional indication of if the file is gzipped. This is infered from the path if
            not provided.

    Returns:
        Parsed OpenEye molecule.

    Raises:
        ValueError: Failed to read molecule.
    """
    filename = urlparse(str(path)).path.split("/")[-1]

    if format is None or gzip is None:
        inferred_format, inferred_gzip = get_format_and_gzip_from_filename(filename)

    format = format if format is not None else inferred_format
    gzip = gzip if gzip is not None else inferred_gzip

    bytes = read_bytes(path)
    return read_molecule_bytes(bytes, format=format, flavor=flavor, gzip=gzip)


def read_molecules(
    path: str | PathLike[str],
    *,
    format: int | str | None = None,
    flavor: int | None = None,
    gzip: bool | None = None,
) -> Generator[oechem.OEGraphMol, None, None]:
    """
    Iterate over multiple OpenEye molecules from a path.

    This path could be a local filepath, a path on S3 or a web URL.

    Args:
        path: Path with relevant schema.
        format: Optional format, either an `OEFlavor_` from `oechem`, or a string such as `PDB`.
            This is infered from the path if not provided.
        flavor: Optional flavor from OpenEye, modifying the format.
        gzip: Optional indication of if the file is gzipped. This is infered from the path if
            not provided.

    Yields:
        Parsed OpenEye molecules.

    Raises:
        ValueError: Failed to read molecule.
    """
    filename = urlparse(str(path)).path.split("/")[-1]

    if format is None or gzip is None:
        inferred_format, inferred_gzip = get_format_and_gzip_from_filename(filename)

    format = format if format is not None else inferred_format
    gzip = gzip if gzip is not None else inferred_gzip

    bytes = read_bytes(path)
    yield from read_molecules_bytes(bytes, format=format, flavor=flavor, gzip=gzip)
