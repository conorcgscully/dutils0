from __future__ import annotations

from functools import lru_cache

import fsutils as fs
from openeye import oechem


@lru_cache
def _get_openeye_formats() -> dict[str, int]:
    formats: dict[str, int] = {}

    for attr in oechem.__dict__:
        if attr.startswith("OEFormat_"):
            name = attr.split("OEFormat_")[1]
            formats[name.upper()] = getattr(oechem, attr)
    return formats


def get_format(name: str) -> int:
    """Get the OpenEye format associated with a given name."""
    return _get_openeye_formats()[name.upper()]


EXTENSION_TO_OEFORMAT = {
    ext: format
    for format in range(1, oechem.OEFormat_MAXFORMAT)
    for ext in oechem.OEGetFormatExtension(format).split(",")
}


def get_oeformat_from_extension(ext: str) -> int:
    """Get the OpenEye format associated with a given extension."""
    return EXTENSION_TO_OEFORMAT[ext.lower()]


def setup_stream_format(
    stream: oechem.oestream, format: str | int, flavor: int | None = None
) -> None:
    """
    Configure an OpenEye stream in-place with a given format and optional flavor.

    Args:
        stream: OpenEye molecule stream.
        format: Either an `OEFlavor_` from `oechem`, or a string such as `PDB`.
        flavor: Optional flavor from OpenEye, modifying the format.

    Raises:
        ValueError: Unrecognized format.
    """
    if isinstance(format, str):
        try:
            format = get_format(format)
        except KeyError:
            raise ValueError(f"Unrecognized format {format}") from None
    stream.SetFormat(format)
    if flavor:
        stream.SetFlavor(format, flavor)


def get_format_and_gzip_from_filename(filename: str, /) -> tuple[int, bool]:
    """Infer the OpenEye file format and whether GZip compression is used from a filename."""
    # Openeye does not expose how they infer formats
    # Instead, when opening a stream to some temporary file with the same suffix
    # (without reading anything), the stream autodetects the format and if it is
    # gzipped.
    extension = fs.utils.get_extension(filename).removeprefix(".")

    if not extension:
        raise ValueError(f"Filename {filename} does not have an extension")

    parts = extension.split(".")

    if parts[-1] == "gz" or parts[-1] == "GZ":
        gzip = True
        parts.pop()
    else:
        gzip = False

    if len(parts) == 0:
        raise ValueError("Filename only has a .gz extension")

    try:
        ext = get_oeformat_from_extension(parts[-1])
    except KeyError:
        raise ValueError(f"Unrecognized extension .{parts[-1]}") from None

    return ext, gzip
