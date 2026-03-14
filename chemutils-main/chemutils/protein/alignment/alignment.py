from __future__ import annotations

from dataclasses import dataclass
from typing import Protocol

from .from_nomenclature import parse_alignment_from_mutation_nomenclature
from .to_nomenclature import get_mutation_nomenclature_from_alignment
from .types import SequencePairs


# BioPython's pairwise2.Alignment does not actually exist as a real type
class BioPythonAlignment(Protocol):
    seqA: str
    seqB: str


@dataclass(frozen=True)
class SequenceAlignment:
    """
    Object describing the alignment of a sequence to a reference sequence.

    The internal representation of an alignment is as a list of 'pairs', where
    each pair represents the alignment of a part of the sequence to a part of the
    reference. For example, the sequence and reference:

    ```
    ACE-  # sequence
    -CDG  # reference
    ```

    would be represented as four pairs:

    * `('A', None)`, an insertion of something into the reference
    * `('C', 'C')`, a simple alignment
    * `('E', 'D')`, a mutation of D of the reference to E of the sequence
    * `(None, 'G')`, a deletion of something from the reference

    Each pair can therefore represent the four possible states of an alignment.
    """

    pairs: SequencePairs

    @property
    def sequence(self) -> str:
        """Sequence involved in the alignment."""
        return "".join([char for char, _ in self.pairs if char is not None])

    @property
    def reference(self) -> str:
        """Reference sequence involved in the alignment."""
        return "".join([char for _, char in self.pairs if char is not None])

    @property
    def aligned_sequence(self) -> str:
        """Sequence padded with `-` to align with `aligned_reference."""
        return "".join([char if char is not None else "-" for char, _ in self.pairs])

    @property
    def aligned_reference(self) -> str:
        """Reference sequence padded with `-` to align with `aligned_sequence."""
        return "".join([char if char is not None else "-" for _, char in self.pairs])

    @property
    def nomenclature(self) -> str:
        return get_mutation_nomenclature_from_alignment(self.pairs)

    @classmethod
    def from_aligned_sequences(
        cls, *, aligned_sequence: str, aligned_reference: str
    ) -> SequenceAlignment:
        """Create an alignment from a pair of aligned sequences."""
        if len(aligned_sequence) != len(aligned_reference):
            raise ValueError("Aligned sequence and reference must be the same length")
        pairs: SequencePairs = tuple(
            (char_seq if char_seq != "-" else None, char_ref if char_ref != "-" else None)
            for char_seq, char_ref in zip(aligned_sequence, aligned_reference, strict=True)
        )
        return cls(pairs=pairs)

    @classmethod
    def from_nomenclature(
        cls, *, sequence: str, reference: str, nomenclature: str
    ) -> SequenceAlignment:
        return cls(
            pairs=parse_alignment_from_mutation_nomenclature(
                sequence=sequence, reference=reference, mutation=nomenclature
            )
        )

    @classmethod
    def from_biopython(cls, obj: BioPythonAlignment, /) -> SequenceAlignment:
        """Create an alignment from a pair of aligned sequences."""
        return cls.from_aligned_sequences(aligned_sequence=obj.seqA, aligned_reference=obj.seqB)
