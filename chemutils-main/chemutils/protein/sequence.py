"""
Module for code reading and writing sequences using BioPython.

BioPython expects a file handle that has already been opened in either
binary or text mode, depending on the format. The `read_sequences_str`
and `read_sequences_bytes` removes the need to know what kind of handle
required.
"""

import io
from collections.abc import Generator
from typing import Final

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

BINARY_FORMATS: Final = ("abi", "abi-trim", "gck", "seqxml", "sff", "sff-trim", "snapgene", "xdna")
"""
BioPython formats that must be read in binary mode.
"""


def read_sequences_str(content: str, /, *, format: str) -> Generator[SeqRecord, None, None]:
    """
    Iterate over BioPython sequence records in the contents of a file.

    Args:
        content: File contents.
        format: File format as understood by BioPython.

    Yields:
        BioPython `SeqRecord`s read from the file.
    """
    if format in BINARY_FORMATS:
        yield from SeqIO.parse(io.BytesIO(content.encode()), format=format)
    else:
        yield from SeqIO.parse(io.StringIO(content), format=format)


def read_sequences_bytes(content: bytes, /, *, format: str) -> Generator[SeqRecord, None, None]:
    """
    Iterate over BioPython sequence records in the contents of a file.

    Args:
        content: File contents.
        format: File format as understood by BioPython.

    Yields:
        BioPython `SeqRecord`s read from the file.
    """
    if format in BINARY_FORMATS:
        yield from SeqIO.parse(io.BytesIO(content), format=format)
    else:
        yield from SeqIO.parse(io.StringIO(content.decode()), format=format)


def read_sequence_str(content: str, /, *, format: str) -> SeqRecord:
    """
    Read a single BioPython sequence record from the contents of a file.

    Args:
        content: File contents.
        format: File format as understood by BioPython.

    Returns:
        BioPython `SeqRecord` read from the file.
    """
    if format in BINARY_FORMATS:
        return SeqIO.read(io.BytesIO(content.encode()), format=format)
    else:
        return SeqIO.read(io.StringIO(content), format=format)


def read_sequence_bytes(content: bytes, /, *, format: str) -> SeqRecord:
    """
    Read a single BioPython sequence records from the contents of a file.

    Args:
        content: File contents.
        format: File format as understood by BioPython.

    Returns:
        BioPython `SeqRecord` read from the file.
    """
    if format in BINARY_FORMATS:
        return SeqIO.read(io.BytesIO(content), format=format)
    else:
        return SeqIO.read(io.StringIO(content.decode()), format=format)


def write_sequences_str(*sequences: SeqRecord, format: str) -> str:
    """
    Write one or more BioPython sequence to a specific file format as a string.

    Args:
        sequences: BioPython sequences to write to file.
        format: File format as understood by BioPython.

    Returns:
        File contents as string.
    """
    if format in BINARY_FORMATS:
        bytes_stream = io.BytesIO()
        SeqIO.write(sequences, bytes_stream, format=format)
        return bytes_stream.getvalue().decode()
    else:
        str_stream = io.StringIO()
        SeqIO.write(sequences, str_stream, format=format)
        return str_stream.getvalue()


def write_sequences_bytes(*sequences: SeqRecord, format: str) -> bytes:
    """
    Write one or more BioPython sequence to a specific file format as bytes.

    Args:
        sequences: BioPython sequences to write to file.
        format: File format as understood by BioPython.

    Returns:
        File contents as string.
    """
    if format in BINARY_FORMATS:
        bytes_stream = io.BytesIO()
        SeqIO.write(sequences, bytes_stream, format=format)
        return bytes_stream.getvalue()
    else:
        str_stream = io.StringIO()
        SeqIO.write(sequences, str_stream, format=format)
        return str_stream.getvalue().encode()
