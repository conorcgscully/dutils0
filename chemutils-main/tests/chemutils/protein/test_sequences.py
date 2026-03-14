import io
from typing import Final

import pytest
from Bio import SeqIO

from chemutils.protein.sequence import (
    BINARY_FORMATS,
    read_sequence_bytes,
    read_sequence_str,
    read_sequences_bytes,
    read_sequences_str,
    write_sequences_bytes,
    write_sequences_str,
)

CASES: Final = [
    ("3100.ab1", "abi"),
    ("3100.ab1", "abi-trim"),
    ("contig1.ace", "ace"),
    ("1A7G.cif", "cif-atom"),
    ("1A7G.cif", "cif-seqres"),
    ("muscle.aln", "clustal"),
    ("SC10H5.embl", "embl"),
    ("protein_lib.fa", "fasta"),
    ("example.fastq", "fastq"),
    ("artificial.gck", "gck"),
    ("KF527485.gbk", "genbank"),
    ("VIF_mase-pro.txt", "ig"),
    ("A04195.imgt", "imgt"),
    ("1A8O.pdb", "pdb-seqres"),
    ("1A8O.pdb", "pdb-atom"),
    ("phd1.txt", "phd"),
    ("four.dat", "phylip"),
    ("dna_example.xml", "seqxml"),
    ("greek.sff", "sff"),
    ("greek.sff", "sff-trim"),
    ("sample-e.dna", "snapgene"),
    ("cath3.sth", "stockholm"),
    ("sample-a.xdna", "xdna"),
]

READ_CASES: Final = [
    ("3100.ab1", "abi"),
    ("3100.ab1", "abi-trim"),
    ("1A7G.cif", "cif-atom"),
    ("1A7G.cif", "cif-seqres"),
    ("SC10H5.embl", "embl"),
    ("artificial.gck", "gck"),
    ("KF527485.gbk", "genbank"),
    ("A04195.imgt", "imgt"),
    ("1A8O.pdb", "pdb-seqres"),
    ("1A8O.pdb", "pdb-atom"),
    ("sample-e.dna", "snapgene"),
    ("sample-a.xdna", "xdna"),
]

OUT_FORMATS = [
    "clustal",
    "embl",
    "fasta",
    "genbank",
    "imgt",
    "phd",
    "phylip",
    "xml-seqxml",
    "sff",
    "stockholm",
    "xdna",
]


@pytest.mark.parametrize(["filename", "format"], CASES)
def test_read_sequences_bytes(filename, format):
    with open(f"tests/chemutils/protein/data/{filename}", "rb") as f:
        content = f.read()

    chemutils_seqs = list(read_sequences_bytes(content, format=format))

    if format in BINARY_FORMATS:
        biopython_seqs = list(SeqIO.parse(io.BytesIO(content), format=format))
    else:
        biopython_seqs = list(SeqIO.parse(io.StringIO(content.decode()), format=format))

    # Biopython does not support comparison of sequences
    assert repr(chemutils_seqs) == repr(biopython_seqs)


@pytest.mark.parametrize(["filename", "format"], CASES)
def test_read_sequences_str(filename, format):
    if format in BINARY_FORMATS and format != "seqxml":
        return

    with open(f"tests/chemutils/protein/data/{filename}") as f:
        content = f.read()

    chemutils_seqs = list(read_sequences_str(content, format=format))

    if format in BINARY_FORMATS:
        biopython_seqs = list(SeqIO.parse(io.BytesIO(content.encode()), format=format))
    else:
        biopython_seqs = list(SeqIO.parse(io.StringIO(content), format=format))

    # Biopython does not support comparison of sequences
    assert repr(chemutils_seqs) == repr(biopython_seqs)


@pytest.mark.parametrize(["filename", "format"], READ_CASES)
def test_read_sequence_bytes(filename, format):
    with open(f"tests/chemutils/protein/data/{filename}", "rb") as f:
        content = f.read()

    chemutils_seqs = read_sequence_bytes(content, format=format)

    if format in BINARY_FORMATS:
        biopython_seqs = SeqIO.read(io.BytesIO(content), format=format)
    else:
        biopython_seqs = SeqIO.read(io.StringIO(content.decode()), format=format)

    # Biopython does not support comparison of sequences
    assert repr(chemutils_seqs) == repr(biopython_seqs)


@pytest.mark.parametrize(["filename", "format"], READ_CASES)
def test_read_sequence_str(filename, format):
    if format in BINARY_FORMATS and format != "seqxml":
        return

    with open(f"tests/chemutils/protein/data/{filename}") as f:
        content = f.read()

    chemutils_seqs = read_sequence_str(content, format=format)

    if format in BINARY_FORMATS:
        biopython_seqs = SeqIO.read(io.BytesIO(content.encode()), format=format)
    else:
        biopython_seqs = SeqIO.read(io.StringIO(content), format=format)

    # Biopython does not support comparison of sequences
    assert repr(chemutils_seqs) == repr(biopython_seqs)


@pytest.mark.parametrize(["filename", "format"], [case for case in CASES if case[1] in OUT_FORMATS])
def test_write_sequences_bytes(filename, format):
    filename = f"tests/chemutils/protein/data/{filename}"
    with open(filename, "rb" if format in BINARY_FORMATS else "r") as f:
        sequences = list(SeqIO.parse(f, format=format))

    chemutils_out = write_sequences_bytes(*sequences, format=format)

    if format in BINARY_FORMATS:
        bytes_io = io.BytesIO()
        SeqIO.write(sequences, bytes_io, format=format)
        biopython_out = bytes_io.getvalue()
    else:
        string_io = io.StringIO()
        SeqIO.write(sequences, string_io, format=format)
        biopython_out = string_io.getvalue().encode()

    assert chemutils_out == biopython_out


@pytest.mark.parametrize(["filename", "format"], [case for case in CASES if case[1] in OUT_FORMATS])
def test_write_sequences_str(filename, format):
    # `seqxml` is a binary format that can be written as a string
    if format in BINARY_FORMATS and format != "seqxml":
        return

    filename = f"tests/chemutils/protein/data/{filename}"
    with open(filename, "rb" if format in BINARY_FORMATS else "r") as f:
        sequences = list(SeqIO.parse(f, format=format))

    chemutils_out = write_sequences_str(*sequences, format=format)

    if format in BINARY_FORMATS:
        bytes_io = io.BytesIO()
        SeqIO.write(sequences, bytes_io, format=format)
        biopython_out = bytes_io.getvalue().decode()
    else:
        string_io = io.StringIO()
        SeqIO.write(sequences, string_io, format=format)
        biopython_out = string_io.getvalue()

    assert chemutils_out == biopython_out
