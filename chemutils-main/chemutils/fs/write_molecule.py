from __future__ import annotations

from collections.abc import Collection
from os import PathLike
from urllib.parse import urlparse

from fsutils import write_bytes
from openeye import oechem

from chemutils.molecule import write_molecule_bytes
from chemutils.molecule.format import get_format_and_gzip_from_filename


def write_molecule(
    path: str | PathLike[str],
    *,
    mol: oechem.OEGraphMol,
    format: int | str | None = None,
    flavor: int | None = None,
    gzip: bool | None = None,
) -> oechem.OEGraphMol:
    """
    Write an OpenEye molecule to a path.

    This path could be a local filepath, a path on S3 or a web URL.

    Args:
        path: Path with relevant schema.
        mol: OpenEye molecule.
        format: Optional format, either an `OEFlavor_` from `oechem`, or a string such as `PDB`.
            This is infered from the path if not provided.
        flavor: Optional flavor from OpenEye, modifying the format.
        gzip: Optional indication of if the file is gzipped. This is infered from the path if
            not provided.

    Raises:
        ValueError: Cannot infer file format.
    """
    filename = urlparse(str(path)).path.split("/")[-1]

    if format is None or gzip is None:
        inferred_format, inferred_gzip = get_format_and_gzip_from_filename(filename)

    format = format if format is not None else inferred_format
    gzip = gzip if gzip is not None else inferred_gzip

    bytes = write_molecule_bytes(mol, format=format, flavor=flavor, gzip=gzip)
    write_bytes(path, contents=bytes)


def write_molecules(
    path: str | PathLike[str],
    *,
    mols: Collection[oechem.OEGraphMol],
    format: int | str | None = None,
    flavor: int | None = None,
    gzip: bool | None = None,
) -> oechem.OEGraphMol:
    """
    Write one or more OpenEye molecules to a single path.

    This path could be a local filepath, a path on S3 or a web URL.

    Args:
        path: Path with relevant schema.
        mols: Collection of one or more OpenEye molecules.
        format: Optional format, either an `OEFlavor_` from `oechem`, or a string such as `PDB`.
            This is infered from the path if not provided.
        flavor: Optional flavor from OpenEye, modifying the format.
        gzip: Optional indication of if the file is gzipped. This is infered from the path if
            not provided.

    Raises:
        ValueError: Cannot infer file format.
    """
    filename = urlparse(str(path)).path.split("/")[-1]

    if format is None or gzip is None:
        inferred_format, inferred_gzip = get_format_and_gzip_from_filename(filename)

    format = format if format is not None else inferred_format
    gzip = gzip if gzip is not None else inferred_gzip

    bytes = write_molecule_bytes(*mols, format=format, flavor=flavor, gzip=gzip)
    write_bytes(path, contents=bytes)
