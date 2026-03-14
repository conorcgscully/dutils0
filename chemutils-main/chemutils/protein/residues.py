"""
Various utilities for dealing with conversions between sequences and residue names.

The full description of a sequence of residues is the list of residue names. These are defined
by the RCSB, and include amino acids (ALA for alanine), nucleotides (DA for deoxyadenosine) and
other nonstandard residues.

A common shorthand are the single letter codes assigned to each standard amino acid and nucleotide
by IUPAC. For example, 'AGPT' could encode the amino acids ['ALA', 'GLY', 'PRO', 'THR']. For
nonstandard residues, there are two approaches we could use:

* Use the 'X' code to represent any residue.
* Encode the nonstandard names in full, using parentheses. For example, 'AG(PCA)' encodes
  ['ALA', 'GLY', 'PCA'].

For protein sequences, we can also consider the three special codes 'B', 'Z' and 'J' (and their three
letter equivalents 'ASX', 'GLX' and 'XLE') for indicating ambiguity between related pairs of amino acids.
"""

from collections.abc import Generator, Sequence
from enum import Enum, auto
from typing import Literal, overload

from openeye import oechem

from ._chemical_component_dictionary import _get_chemical_component_dictionary

PROTEIN_CODES = (
    ("ALA", "A"),
    ("CYS", "C"),
    ("ASP", "D"),
    ("GLU", "E"),
    ("PHE", "F"),
    ("GLY", "G"),
    ("HIS", "H"),
    ("ILE", "I"),
    ("LYS", "K"),
    ("LEU", "L"),
    ("MET", "M"),
    ("ASN", "N"),
    ("PRO", "P"),
    ("GLN", "Q"),
    ("ARG", "R"),
    ("SER", "S"),
    ("THR", "T"),
    ("VAL", "V"),
    ("TRP", "W"),
    ("TYR", "Y"),
)
"""Standard IUPAC codes for amino acids."""

AMBIGUOUS_PROTEIN_CODES = (
    ("ASX", "B"),
    ("GLX", "Z"),
    ("XLE", "J"),
)
"""Standard IUPAC codes for ambiguous amino acids."""

DNA_CODES = (
    ("DA", "A"),
    ("DT", "T"),
    ("DU", "U"),
    ("DG", "G"),
    ("DC", "C"),
)
"""Standard IUPAC codes for DNA nucleotides."""

RNA_CODES = (
    ("A", "A"),
    ("U", "U"),
    ("G", "G"),
    ("C", "C"),
)
"""Standard IUPAC codes for RNA nucleotides."""

SEQUENCE_TYPE_TO_UNKNOWN_RESIDUE_NAME = {"protein": "UNK", "dna": "DN", "rna": "N"}

PROTEIN_CODE_ONE_TO_THREE = {oneletter: threeletter for threeletter, oneletter in PROTEIN_CODES}
AMBIGIOUS_PROTEIN_CODE_ONE_TO_THREE = {
    oneletter: threeletter for threeletter, oneletter in AMBIGUOUS_PROTEIN_CODES
}
DNA_CODE_ONE_TO_THREE = {oneletter: threeletter for threeletter, oneletter in DNA_CODES}
RNA_CODE_ONE_TO_THREE = {oneletter: threeletter for threeletter, oneletter in RNA_CODES}

PROTEIN_CODE_THREE_TO_ONE = dict(PROTEIN_CODES)
AMBIGIOUS_PROTEIN_CODE_THREE_TO_ONE = dict(AMBIGUOUS_PROTEIN_CODES)
DNA_CODE_THREE_TO_ONE = dict(DNA_CODES)
RNA_CODE_THREE_TO_ONE = dict(RNA_CODES)


class SequenceType(str, Enum):
    """Identify the type of sequence, to distinguish 'AGUC' for example."""

    Protein = "protein"
    DNA = "dna"
    RNA = "rna"


SequenceTypeArg = SequenceType | Literal["protein", "dna", "rna"]


SEQUENCE_TYPE_ONE_TO_THREE: dict[SequenceTypeArg, dict[str, str]] = {
    SequenceType.Protein: PROTEIN_CODE_ONE_TO_THREE,
    SequenceType.DNA: DNA_CODE_ONE_TO_THREE,
    SequenceType.RNA: RNA_CODE_ONE_TO_THREE,
}

SEQUENCE_TYPE_THREE_TO_ONE: dict[SequenceTypeArg, dict[str, str]] = {
    SequenceType.Protein: PROTEIN_CODE_THREE_TO_ONE,
    SequenceType.DNA: DNA_CODE_THREE_TO_ONE,
    SequenceType.RNA: RNA_CODE_THREE_TO_ONE,
}


class NonstandardResidueHandling(Enum):
    """How to encode nonstandard residues into a sequence."""

    Exception = auto()
    """Raise an exception if a nonstandard residue name is encountered."""
    AnyCode = auto()
    """Use an X to indicate a nonstandard residue."""
    Parentheses = auto()
    """Use parentheses (such as `(PCA)`) to indicate nonstandard residues."""


class UnknownCodeHandling(Enum):
    """How to handle an 'X' when converting a sequence to residue names."""

    Exception = auto()
    """Raise an exception if a nonstandard character is encountered."""
    ReturnNone = auto()
    """Return None for the residue name of an 'X' character."""
    ResidueName = auto()
    """Use the standard residue name for an unknown residue."""


def enumerate_sequence(sequence: str, /) -> Generator[tuple[str, bool], None, None]:
    """
    Enumertate a sequence, supporting bracket notation for custom residue names.

    For example, 'AK(PCA)L' yields ('A', False), ('K', False), ('PCA', True) and ('L', False).

    The reason why you cannot use the length to distinguish single-letter sequence codes and residue
    names is that there are single letter residue names (such as 'E').

    Args:
        sequence: Sequence, consisting of single letters.

    Yields:
        Tuples of the form (code, is_nonstd). If is_nonstd is True, then the code may be more than one
        letter and represents a residue name. If it is False, then it represents a single character.

    Raises:
        ValueError: Sequence is malformed.
    """
    current_nonstd = None
    for char in sequence:
        if char == "(":
            if current_nonstd is not None:
                raise ValueError(f"Invalid sequence {sequence}")
            current_nonstd = ""
        elif char == ")":
            if current_nonstd is None:
                raise ValueError(f"Invalid sequence {sequence}")
            yield current_nonstd, True
            current_nonstd = None
        elif current_nonstd is not None:
            current_nonstd += char
        else:
            yield char, False
    if current_nonstd is not None:
        raise ValueError(f"Invalid sequence {sequence}")


@overload
def residue_names_from_sequence(
    sequence: str,
    /,
    *,
    sequence_type: SequenceTypeArg,
    unknown_code_handling: Literal[UnknownCodeHandling.Exception]
    | Literal[UnknownCodeHandling.ResidueName],
    nonstandard_code_handling: Literal[UnknownCodeHandling.Exception]
    | Literal[UnknownCodeHandling.ResidueName]
    | None,
    allow_ambiguous: bool = ...,
) -> list[str]:
    pass


@overload
def residue_names_from_sequence(
    sequence: str,
    /,
    *,
    sequence_type: SequenceTypeArg,
    unknown_code_handling: UnknownCodeHandling = ...,
    nonstandard_code_handling: UnknownCodeHandling | None = ...,
    allow_ambiguous: bool = ...,
) -> list[str] | list[str | None]:
    pass


def residue_names_from_sequence(
    sequence: str,
    /,
    *,
    sequence_type: SequenceTypeArg,
    unknown_code_handling: UnknownCodeHandling = UnknownCodeHandling.Exception,
    nonstandard_code_handling: UnknownCodeHandling | None = None,
    allow_ambiguous: bool = False,
) -> list[str] | list[str | None]:
    """
    Convert a sequence of characters to their corresponding standard residue names.

    Args:
        sequence: Sequence consisting of upper case characters.
        sequence_type: Sequence type defining which residue name-character mapping to use.
        unknown_code_handling: How to handle 'X' in the sequence.
        nonstandard_code_handling: How to handle nonstandard characters (other than 'X') in the sequence.
        allow_ambiguous: If True, allow the use of ambiguous amino acids ASX, GLX and XLE
            (and their corresponding single characters). If False, raise an exception if
            one of these characters is encountered.

    Returns:
        List of residue names, or a list of optional residue names if `nonstandard_handling`
        is True.

    Raises:
        ValueError: Sequence contains nonstandard characters and `allow_nonstandard` is not True.
    """
    if nonstandard_code_handling is None:
        nonstandard_code_handling = unknown_code_handling
    one_to_three = SEQUENCE_TYPE_ONE_TO_THREE[sequence_type] | {"-": "-"}
    if allow_ambiguous and sequence_type == SequenceType.Protein:
        one_to_three |= AMBIGIOUS_PROTEIN_CODE_ONE_TO_THREE

    residue_names: list[str | None] = []
    for char, nonstd in enumerate_sequence(sequence):
        if nonstd:
            residue_names.append(char)
        elif char in one_to_three:
            residue_names.append(one_to_three[char])
        else:
            handling = unknown_code_handling if char == "X" else nonstandard_code_handling
            match handling:
                case UnknownCodeHandling.Exception:
                    raise ValueError(f"Sequence contains non-standard character {char}.")
                case UnknownCodeHandling.ReturnNone:
                    residue_names.append(None)
                case UnknownCodeHandling.ResidueName:
                    residue_names.append(SEQUENCE_TYPE_TO_UNKNOWN_RESIDUE_NAME[sequence_type])

    return residue_names


@overload
def sequence_from_residue_names(
    residue_names: list[str],
    /,
    *,
    sequence_type: SequenceTypeArg | None = None,
    nonstandard_handling: Literal[NonstandardResidueHandling.Parentheses] = ...,
    allow_ambiguous: bool = ...,
    replace_nonstd_with_parent: bool = False,
) -> str:
    pass


@overload
def sequence_from_residue_names(
    residue_names: Sequence[str] | Sequence[str | None],
    /,
    *,
    sequence_type: SequenceTypeArg | None = None,
    nonstandard_handling: NonstandardResidueHandling = ...,
    allow_ambiguous: bool = ...,
    replace_nonstd_with_parent: bool = False,
) -> str:
    pass


def sequence_from_residue_names(
    residue_names: Sequence[str] | Sequence[str | None],
    /,
    *,
    sequence_type: SequenceTypeArg | None = None,
    nonstandard_handling: NonstandardResidueHandling = NonstandardResidueHandling.Exception,
    allow_ambiguous: bool = False,
    replace_nonstd_with_parent: bool = False,
) -> str:
    """
    Convert a list of residue names to the corresponding one-letter sequence.

    Args:
        residue_names: List of residue names to convert to a sequence.
        nonstandard_handling: How to encode nonstandard residue names.
        allow_ambiguous: If True, allow the parsing of ambiguous amino acids ASX, GLX and XLE
            to their corresponding single letters. If False, raise an exception if one of these
            characters is encountered.
        replace_nonstd_with_parent: If True, replace nonstandard residues that have a standard
            parent with the parent's name. For example, 'PCA' would be replaced with 'GLN'.
        sequence_type: Type of sequence to use for the conversion.

    Returns:
        Sequence of one-letter codes as a string.

    Raises:
        ValueError: Residue names contains nonstandard residue names and `allow_nonstandard` is not True.
    """
    if sequence_type:
        three_to_one = SEQUENCE_TYPE_THREE_TO_ONE[sequence_type] | {"-": "-"}
    else:
        three_to_one = (
            PROTEIN_CODE_THREE_TO_ONE | RNA_CODE_THREE_TO_ONE | DNA_CODE_THREE_TO_ONE | {"-": "-"}
        )
    if allow_ambiguous and sequence_type in ["protein", None]:
        three_to_one |= AMBIGIOUS_PROTEIN_CODE_THREE_TO_ONE

    if replace_nonstd_with_parent:
        residue_names = [
            replace_residue_name_with_monomer_parent(residue_name) if residue_name else None
            for residue_name in residue_names
        ]
    match nonstandard_handling:
        case NonstandardResidueHandling.AnyCode:
            return "".join(three_to_one.get(residue_name, "X") for residue_name in residue_names)  # type: ignore
        case NonstandardResidueHandling.Parentheses:
            if None in residue_names:
                raise ValueError("Cannot encode residue name `None`.")
            return "".join(
                three_to_one.get(residue_name, f"({residue_name})")  # type: ignore
                for residue_name in residue_names
            )
        case NonstandardResidueHandling.Exception:
            try:
                return "".join(three_to_one[residue_name] for residue_name in residue_names)  # type: ignore
            except KeyError as e:
                raise ValueError("Nonstandard residue name encountered.") from e
        case _:
            raise ValueError(f"Unknown nonstandard handling {nonstandard_handling}.")


def iterate_residues(oemol: oechem.OEMolBase, /) -> Generator[oechem.OEResidue, None, None]:
    """
    Iterate over the residues in the order they appear in an OpenEye molecule.

    Args:
        oemol: OpenEye molecule to iterate over.

    Yields:
        `OEResidue` objects in order that they appear in the molecule.
    """
    hv = oechem.OEHierView(oemol)
    for residue in hv.GetResidues():
        yield residue.GetOEResidue()


def replace_residue_name_with_monomer_parent(residue_name: str) -> str:
    ccd = _get_chemical_component_dictionary()
    if residue_name not in ccd:
        return residue_name
    return ccd[residue_name].get("monomer_parent_id", residue_name)
