"""
Get alignment from standardised nomenclature for mutations.

See the docstring of `chemutils.protein.alignment` for more details.
"""

from ..mutation import (
    DeletionMutation,
    InsertionMutation,
    MissenseMutation,
    parse_mutation_nomenclature,
)
from ..residues import PROTEIN_CODE_THREE_TO_ONE
from .types import SequencePairs


def parse_alignment_from_mutation_nomenclature(
    *, sequence: str, reference: str, mutation: str
) -> SequencePairs:
    """
    Parse an alignment given a sequence, a reference sequence and a string describing their differences in standardised notation.

    Args:
        sequence: Sequence given by mutating the reference sequence.
        reference: Original reference sequence.
        mutation: Mutation between the two, in standardised notation.

    Returns:
        Sequence pairs, as a list of alignments. See the docstring of the `SequenceAlignment` class for more details.
    """
    pairs: list[tuple[str | None, str | None]] = []
    index = 0  # Current index in the sequence
    index_ref = 0  # Current index in the reference sequence
    length = len(sequence)
    length_ref = len(reference)

    def add_sequence_pair_to_output(in_sequence: bool, in_reference: bool, change: bool) -> None:
        """
        Extend the sequence alignment, given information about the next difference.

        Args:
            in_sequence: Does the next pairing involve the next element of `sequence`.
            in_reference: Does the next pairing involve the next element of `reference`.
            change: Is a change involved.
        """
        nonlocal index, index_ref
        if in_sequence and in_reference and index < len(sequence) and index_ref < len(reference):
            if not change and sequence[index] != reference[index_ref]:
                raise ValueError(
                    f"Sequence with {sequence[index]} at position {index + 1} does not match reference with {reference[index_ref]} at position {index_ref + 1}"
                )
            if change and sequence[index] == reference[index_ref]:
                raise ValueError(
                    f"Sequence with {sequence[index]} at position {index + 1} should not match reference with {reference[index_ref]} at position {index_ref + 1}"
                )
        if in_sequence:
            index += 1
        if in_reference:
            index_ref += 1
        pairs.append(
            (
                sequence[index - 1] if in_sequence else None,
                reference[index_ref - 1] if in_reference else None,
            )
        )

    mutations = parse_mutation_nomenclature(mutation)

    # Iterate over all parts of the mutation
    for mut in mutations:
        if isinstance(mut, InsertionMutation):
            # Index of residue in reference that precedes the insertion, if applicable
            reference_index_preceding = (
                mut.preceding_resnum - 1 if mut.preceding_resnum is not None else None
            )
            # Residue in reference that precedes the insertion, if applicable
            reference_residue_preceding = (
                PROTEIN_CODE_THREE_TO_ONE[mut.preceding_resname]
                if mut.preceding_resname is not None
                else None
            )
            # Index of residue in reference that follows the insertion, if applicable
            reference_index_subsequent = (
                mut.subsequent_resnum - 1 if mut.subsequent_resnum is not None else None
            )
            # Residue in reference that follows the insertion, if applicable
            reference_residue_subsequent = (
                PROTEIN_CODE_THREE_TO_ONE[mut.subsequent_resname]
                if mut.subsequent_resname is not None
                else None
            )
            # Sequence of residues inserted into the reference.
            inserted_residues = "".join(
                PROTEIN_CODE_THREE_TO_ONE[resname] for resname in mut.inserted_resnames
            )

            # Assume up to this insertion, everything is aligned
            if reference_index_preceding is not None:
                for _ in range(reference_index_preceding - index_ref + 1):
                    add_sequence_pair_to_output(True, True, False)

            # Check that if there is a preceding residue, that it matches the reference
            if (
                reference_index_preceding is not None
                and reference_residue_preceding != reference[index_ref - 1]
            ):
                raise ValueError(
                    f"Mutation {mut} with {reference_residue_preceding} at {reference_index_preceding + 1} does not match reference with {reference[index_ref - 1]} at position {index_ref}"
                )

            # Check that if there is a subsequent residue, that it matches the reference
            if (
                reference_index_subsequent is not None
                and reference_residue_subsequent != reference[index_ref]
            ):
                raise ValueError(
                    f"Mutation {mut} with {reference_residue_subsequent} at {reference_index_subsequent + 1} does not match reference with {reference[index_ref]} at position {index_ref + 1}"
                )

            # Check that the inserted residues match the sequence
            if inserted_residues != sequence[index : index + len(inserted_residues)]:
                raise ValueError(
                    f"Mutation {mut} inserts {inserted_residues} but sequence has {sequence[index : index + len(inserted_residues)]}"
                )

            for _ in inserted_residues:
                add_sequence_pair_to_output(True, False, True)

        elif isinstance(mut, DeletionMutation):
            # Index of the start of the deletion, relative to the reference sequence
            reference_index_start = int(mut.from_resnum) - 1
            # Index of the end of the deletion (inclusive), relative to the reference sequence, if applicable
            reference_index_end = (
                int(mut.to_resnum) - 1 if mut.to_resnum != mut.from_resnum else None
            )
            # Assume up to this deletion, everything is aligned
            for _ in range(reference_index_start - index_ref):
                add_sequence_pair_to_output(True, True, False)

            reference_residue_start = PROTEIN_CODE_THREE_TO_ONE[mut.from_resname]

            # Check that the residue at the start of the deletion matches the reference
            if reference_residue_start != reference[index_ref]:
                raise ValueError(
                    f"Mutation {mut} with {reference_residue_start} at {index_ref + 1} does not match reference with {reference[index_ref]} at position {index_ref + 1}"
                )

            if reference_index_end is None:
                add_sequence_pair_to_output(False, True, True)
            else:
                # Total number of deleted residues
                deletion_length = reference_index_end - reference_index_start + 1

                # Check that the residue at the end of the deletion matches the reference
                reference_residue_end = PROTEIN_CODE_THREE_TO_ONE[mut.to_resname]
                if (
                    reference_residue_end is not None
                    and reference_residue_end != reference[index_ref + deletion_length - 1]
                ):
                    raise ValueError(
                        f"Mutation {mut} with {reference_residue_end} at {index_ref + deletion_length} does not match reference with {reference[index_ref + deletion_length - 1]} at position {index_ref + deletion_length}"
                    )

                for _ in range(deletion_length):
                    add_sequence_pair_to_output(False, True, True)

        elif isinstance(mut, MissenseMutation):
            # 0-based index in reference where the mutation occured
            reference_mutation_index = mut.resnum - 1

            # Align until we get to this mutation
            for _ in range(reference_mutation_index - index_ref):
                add_sequence_pair_to_output(True, True, False)

            # Check if mutation's residue labels actually make sense
            reference_residue = PROTEIN_CODE_THREE_TO_ONE[mut.resname]
            if reference_residue != reference[index_ref]:
                raise ValueError(
                    f"Mutation {mut} with {reference_residue} at {index_ref + 1} does not match reference with {reference[index_ref]} at position {index_ref + 1}"
                )

            sequence_residue = PROTEIN_CODE_THREE_TO_ONE[mut.replacement_resname]
            if sequence_residue != sequence[index]:
                raise ValueError(
                    f"Mutation {mut} with {sequence_residue} at {index + 1} does not match sequence with {sequence[index]} at position {index + 1}"
                )

            # Add mutation
            add_sequence_pair_to_output(True, True, True)

    # Align end of sequence
    if length - index != length_ref - index_ref:
        raise ValueError(
            f"Sequence and reference have mismatched lengths: sequence has {length - index} residues remaining after all mutations resolved, and reference has {length_ref - index_ref}."
        )

    for _ in range(length - index):
        add_sequence_pair_to_output(True, True, False)

    return tuple(pairs)
