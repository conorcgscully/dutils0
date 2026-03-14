import inspect

import pytest

from chemutils.molecule import read_molecule_str
from chemutils.protein.residues import (
    NonstandardResidueHandling,
    SequenceType,
    UnknownCodeHandling,
    enumerate_sequence,
    iterate_residues,
    residue_names_from_sequence,
    sequence_from_residue_names,
)


@pytest.mark.parametrize(
    ["sequence", "sequence_type", "unknown_code_handling", "allow_ambiguous", "expected"],
    [
        ("", SequenceType.Protein, UnknownCodeHandling.Exception, False, []),
        ("A", "protein", UnknownCodeHandling.Exception, False, ["ALA"]),
        ("HPE", SequenceType.Protein, UnknownCodeHandling.Exception, False, ["HIS", "PRO", "GLU"]),
        (
            "TYCGK",
            SequenceType.Protein,
            UnknownCodeHandling.Exception,
            False,
            ["THR", "TYR", "CYS", "GLY", "LYS"],
        ),
        (
            "TY--K",
            SequenceType.Protein,
            UnknownCodeHandling.Exception,
            False,
            ["THR", "TYR", "-", "-", "LYS"],
        ),
        (
            "T(PCA)K",
            SequenceType.Protein,
            UnknownCodeHandling.Exception,
            False,
            ["THR", "PCA", "LYS"],
        ),
        ("X", SequenceType.Protein, UnknownCodeHandling.Exception, False, ValueError),
        ("X", SequenceType.Protein, UnknownCodeHandling.ReturnNone, False, [None]),
        ("BE", "protein", UnknownCodeHandling.Exception, False, ValueError),
        ("BE", SequenceType.Protein, UnknownCodeHandling.ReturnNone, False, [None, "GLU"]),
        (
            "(PCA)BE",
            SequenceType.Protein,
            UnknownCodeHandling.ReturnNone,
            False,
            ["PCA", None, "GLU"],
        ),
        ("BE", SequenceType.Protein, UnknownCodeHandling.ReturnNone, True, ["ASX", "GLU"]),
        ("BE", SequenceType.Protein, UnknownCodeHandling.Exception, True, ["ASX", "GLU"]),
        ("", "rna", UnknownCodeHandling.Exception, False, []),
        ("A", SequenceType.RNA, UnknownCodeHandling.Exception, False, ["A"]),
        ("UC", SequenceType.RNA, UnknownCodeHandling.Exception, False, ["U", "C"]),
        ("UAGC", SequenceType.RNA, UnknownCodeHandling.Exception, False, ["U", "A", "G", "C"]),
        ("U-G-", "rna", UnknownCodeHandling.Exception, False, ["U", "-", "G", "-"]),
        ("XA", SequenceType.RNA, UnknownCodeHandling.Exception, False, ValueError),
        ("XA", SequenceType.RNA, UnknownCodeHandling.ReturnNone, False, [None, "A"]),
        ("", SequenceType.DNA, UnknownCodeHandling.Exception, False, []),
        ("A", "dna", UnknownCodeHandling.Exception, False, ["DA"]),
        ("TC", SequenceType.DNA, UnknownCodeHandling.Exception, False, ["DT", "DC"]),
        ("TAGC", SequenceType.DNA, UnknownCodeHandling.Exception, False, ["DT", "DA", "DG", "DC"]),
        ("XA", "dna", UnknownCodeHandling.Exception, False, ValueError),
        ("XA", SequenceType.DNA, UnknownCodeHandling.ReturnNone, False, [None, "DA"]),
    ],
)
def test_residue_names_from_sequence(
    sequence, sequence_type, unknown_code_handling, allow_ambiguous, expected
):
    if inspect.isclass(expected) and issubclass(expected, Exception):
        with pytest.raises(expected):
            _ = residue_names_from_sequence(
                sequence,
                sequence_type=sequence_type,
                unknown_code_handling=unknown_code_handling,
                allow_ambiguous=allow_ambiguous,
            )
    else:
        assert (
            residue_names_from_sequence(
                sequence,
                sequence_type=sequence_type,
                unknown_code_handling=unknown_code_handling,
                allow_ambiguous=allow_ambiguous,
            )
            == expected
        )


@pytest.mark.parametrize(
    ["residue_names", "nonstandard_handling", "allow_ambiguous", "expected"],
    [
        ([], NonstandardResidueHandling.Exception, False, ""),
        (["ALA"], NonstandardResidueHandling.Exception, False, "A"),
        (["HIS", "PRO", "GLU"], NonstandardResidueHandling.Exception, False, "HPE"),
        (["-", "PRO", "-"], NonstandardResidueHandling.Exception, False, "-P-"),
        (["THR", "TYR", "CYS", "GLY", "LYS"], NonstandardResidueHandling.Exception, False, "TYCGK"),
        (["PCA"], NonstandardResidueHandling.Exception, False, ValueError),
        (["PCA"], NonstandardResidueHandling.AnyCode, False, "X"),
        (["PCA"], NonstandardResidueHandling.Parentheses, False, "(PCA)"),
        ([None], NonstandardResidueHandling.Exception, False, ValueError),
        ([None], NonstandardResidueHandling.AnyCode, False, "X"),
        ([None], NonstandardResidueHandling.Parentheses, False, ValueError),
        (["ASX", "GLU"], NonstandardResidueHandling.Exception, False, ValueError),
        (["ASX", "GLU"], NonstandardResidueHandling.AnyCode, False, "XE"),
        (["ASX", "GLU"], NonstandardResidueHandling.AnyCode, True, "BE"),
        (["ASX", "GLU"], NonstandardResidueHandling.Exception, True, "BE"),
    ],
)
def test_sequence_from_residue_names(
    residue_names, nonstandard_handling, allow_ambiguous, expected
):
    if inspect.isclass(expected) and issubclass(expected, Exception):
        with pytest.raises(expected):
            _ = sequence_from_residue_names(
                residue_names,
                nonstandard_handling=nonstandard_handling,
                allow_ambiguous=allow_ambiguous,
            )
    else:
        assert (
            sequence_from_residue_names(
                residue_names,
                nonstandard_handling=nonstandard_handling,
                allow_ambiguous=allow_ambiguous,
            )
            == expected
        )


@pytest.mark.parametrize(
    ["residue_names", "sequence_type", "expected"],
    [
        (["DG", "G"], None, "GG"),
        (["DG", "G"], "dna", "G(G)"),
        (["DG", "G"], "rna", "(DG)G"),
        (["DG", "G"], "protein", "(DG)(G)"),
    ],
)
def test_sequence_from_residue_names_sequence_type(residue_names, sequence_type, expected):
    assert (
        sequence_from_residue_names(
            residue_names,
            sequence_type=sequence_type,
            nonstandard_handling=NonstandardResidueHandling.Parentheses,
        )
        == expected
    )


@pytest.mark.parametrize(
    ["residue_names", "expected"],
    [(["PCA", "ALA"], "QA"), (["G", "8AH"], "GA")],
)
def test_sequence_from_residue_names_replace_nonstd_with_parent(residue_names, expected):
    assert (
        sequence_from_residue_names(
            residue_names,
            replace_nonstd_with_parent=True,
            nonstandard_handling=NonstandardResidueHandling.Parentheses,
        )
        == expected
    )


@pytest.mark.parametrize("sequence", ["AA(B", "AA(B)C)"])
def test_enumerate_sequence_error(sequence):
    with pytest.raises(ValueError):
        for _ in enumerate_sequence(sequence):
            pass


def test_iterate_residues():
    PDB = """
ATOM      2  CA  THR A   2     -16.753 -13.175 -12.762  1.00 30.00           C
ATOM      9  CA  GLU A   3     -14.487 -11.863 -15.565  1.00 37.14           C
ATOM     18  CA  TYR A   4     -12.690  -8.531 -14.835  1.00 34.56           C
ATOM     30  CA  LYS A   5     -10.120  -7.211 -17.369  1.00 24.73           C
ATOM     37  NZ  LYS A   5      -9.468  -8.399 -21.633  1.00 52.17           N
ATOM     39  CA  LEU A   6      -7.479  -5.407 -15.247  1.00 22.66           C  """
    oemol = read_molecule_str(PDB, format="PDB")
    residues = list(iterate_residues(oemol))
    assert len(residues) == 5
    assert [residue.GetName().strip() for residue in residues] == [
        "THR",
        "GLU",
        "TYR",
        "LYS",
        "LEU",
    ]
