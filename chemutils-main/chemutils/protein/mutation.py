import re
from dataclasses import dataclass

from .residues import PROTEIN_CODE_ONE_TO_THREE, PROTEIN_CODE_THREE_TO_ONE

INSERTION_RE = re.compile(
    r"(?:(?P<residue_preceding>\w)(?P<resnum_preceding>\d+))?_(?:(?P<residue_subsequent>\w)(?P<resnum_subsequent>\d+))?ins(?P<inserted_resnames>\w+)"
)
"""Match for insertion patterns, such as `_M1insGP`, `A24_Y25insGP` and `A24_insK`."""

DELETION_RE = re.compile(
    r"(?P<residue_start>\w)(?P<resnum_start>\d+)(?:_(?P<residue_end>\w)(?P<resnum_end>\d+))?del"
)
"""Match for deletion patterns, such as `A24del` or `B72_T102del`."""

MISSENSE_RE = re.compile(r"(?P<residue>\w)(?P<resnum>\d+)(?P<replacement_resname>\w)")
"""Match for mutation such as `G12D`."""


@dataclass
class InsertionMutation:
    preceding_resname: str | None
    preceding_resnum: int | None
    subsequent_resname: str | None
    subsequent_resnum: int | None
    inserted_resnames: list[str]

    @property
    def nomenclature(self) -> str:
        preceding = None
        if self.preceding_resname is not None and self.preceding_resnum is not None:
            preceding = (
                f"{PROTEIN_CODE_THREE_TO_ONE[self.preceding_resname]}{self.preceding_resnum}"
            )
        subsequent = None
        if self.subsequent_resname is not None and self.subsequent_resnum is not None:
            subsequent = (
                f"{PROTEIN_CODE_THREE_TO_ONE[self.subsequent_resname]}{self.subsequent_resnum}"
            )
        inserted = "".join(PROTEIN_CODE_THREE_TO_ONE[resname] for resname in self.inserted_resnames)
        return f"{preceding if preceding else ''}_{subsequent if subsequent else ''}ins{inserted}"


@dataclass
class DeletionMutation:
    from_resname: str
    from_resnum: int
    to_resname: str
    to_resnum: int

    @property
    def nomenclature(self) -> str:
        if self.from_resnum != self.to_resnum:
            from_ = f"{PROTEIN_CODE_THREE_TO_ONE[self.from_resname]}{self.from_resnum}"
            to_ = f"{PROTEIN_CODE_THREE_TO_ONE[self.to_resname]}{self.to_resnum}"
            return f"{from_}_{to_}del"
        else:
            from_ = f"{PROTEIN_CODE_THREE_TO_ONE[self.from_resname]}{self.from_resnum}"
            return f"{from_}del"


@dataclass
class MissenseMutation:
    resname: str
    resnum: int
    replacement_resname: str

    @property
    def nomenclature(self) -> str:
        return f"{PROTEIN_CODE_THREE_TO_ONE[self.resname]}{self.resnum}{PROTEIN_CODE_THREE_TO_ONE[self.replacement_resname]}"


Mutation = InsertionMutation | DeletionMutation | MissenseMutation


class InvalidMutationError(ValueError):
    pass


def parse_mutation_nomenclature(nomenclature: str, /) -> list[Mutation]:
    return [_parse_mutation_nomenclature_part(part) for part in nomenclature.split()]


def _parse_mutation_nomenclature_part(nomenclature: str, /) -> Mutation:
    if match := INSERTION_RE.match(nomenclature):
        # Seqnum of residue in reference that precedes the insertion, if applicable
        resnum_preceding = (
            int(match.group("resnum_preceding")) if match.group("resnum_preceding") else None
        )
        # Residue in reference that precedes the insertion, if applicable
        residue_preceding = match.group("residue_preceding")
        # Seqnum of residue in reference that follows the insertion, if applicable
        resnum_subsequent = (
            int(match.group("resnum_subsequent")) if match.group("resnum_subsequent") else None
        )
        # Residue in reference that follows the insertion, if applicable
        residue_subsequent = match.group("residue_subsequent")
        # Sequence of residues inserted into the reference.
        inserted_resnames = match.group("inserted_resnames")

        return InsertionMutation(
            preceding_resname=PROTEIN_CODE_ONE_TO_THREE[residue_preceding]
            if residue_preceding is not None
            else None,
            preceding_resnum=resnum_preceding,
            subsequent_resname=PROTEIN_CODE_ONE_TO_THREE[residue_subsequent]
            if residue_subsequent is not None
            else None,
            subsequent_resnum=resnum_subsequent,
            inserted_resnames=[PROTEIN_CODE_ONE_TO_THREE[code] for code in inserted_resnames],
        )
    elif match := DELETION_RE.match(nomenclature):
        resnum_start = int(match.group("resnum_start"))
        residue_start = PROTEIN_CODE_ONE_TO_THREE[match.group("residue_start")]

        resnum_end = (
            int(match.group("resnum_end")) if match.group("resnum_end") is not None else None
        )
        residue_end = (
            PROTEIN_CODE_ONE_TO_THREE[match.group("residue_end")]
            if match.group("residue_end") is not None
            else None
        )

        return DeletionMutation(
            from_resname=residue_start,
            from_resnum=resnum_start,
            to_resname=residue_end if residue_end is not None else residue_start,
            to_resnum=resnum_end if resnum_end is not None else resnum_start,
        )

    elif match := MISSENSE_RE.match(nomenclature):
        resnum = int(match.group("resnum"))
        resname = PROTEIN_CODE_ONE_TO_THREE[match.group("residue")]
        replacement = PROTEIN_CODE_ONE_TO_THREE[match.group("replacement_resname")]

        return MissenseMutation(resnum=resnum, resname=resname, replacement_resname=replacement)

    raise InvalidMutationError(f"Invalid mutation `{nomenclature}`.")
