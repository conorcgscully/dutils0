"""
Generate standardised nomenclature for protein alignments & mutations.

See the docstring of `chemutils.protein.alignment` for more details.
"""

from itertools import groupby

from .types import SequencePairs


def format_sequence_region(*, sequence: str, start: int, end: int) -> str:
    """
    Format a region of a sequence as a string, e.g. "A1_B2", "_A1" or "A1_".

    Args:
        sequence: Sequence of single letter characters.
        start: Start index of region (inclusive). `-1` indicates the region starts before the sequence.
        end: End index of region (inclusive)
    """
    if start == -1 and end == 0:
        return f"_{sequence[0]}{1}"
    elif start == end:
        return f"{sequence[start]}{start + 1}"
    elif start < len(sequence) and end < len(sequence):
        return f"{sequence[start]}{start + 1}_{sequence[end]}{end + 1}"
    elif start == len(sequence) - 1 and end == len(sequence):
        return f"{sequence[start]}{start + 1}_"
    raise ValueError(f"Cannot format sequence region of {sequence} from {start} to {end}")


def get_nomenclature_for_insertion(*, sequence: str, start: int, inserted: str) -> str:
    """
    Get the standard nomenclature for an insertion.

    Args:
        sequence: Sequence being inserted into.
        start: 0-based index of the residue preceding the insertion.
        inserted: Sequence being inserted.
    """
    return f"{format_sequence_region(sequence=sequence, start=start, end=start + 1)}ins{inserted}"


def get_nomenclature_for_deletion(*, sequence: str, start: int, end: int) -> str:
    """
    Get the standard nomenclature for a deletion.

    Args:
        sequence: Sequence being deleted from.
        start: 0-based index of the first residue being deleted.
        end: 0-based index of the last residue being deleted.
    """
    return f"{format_sequence_region(sequence=sequence, start=start, end=end)}del"


def get_mutation_nomenclature_from_alignment(pairs: SequencePairs, /) -> str:
    """Get the protein mutation nomenclature that describes this difference between a sequence and a reference sequence."""
    nomenclature = []

    last_seq_index = -1
    last_ref_index = -1

    sequence = "".join([char for char, _ in pairs if char is not None])
    reference = "".join([char for _, char in pairs if char is not None])

    pair_infos = [
        (sequence is not None, reference is not None, sequence != reference)
        for sequence, reference in pairs
    ]

    # Using `groupby` works on loci info as it groups together adjacent mutation/insertions/deletions
    for pair_info, group in groupby(pair_infos):
        in_sequence, in_reference, is_change = pair_info
        group_size = len(list(group))

        next_seq_index = last_seq_index + 1 + group_size if in_sequence else last_seq_index + 1
        next_ref_index = last_ref_index + 1 + group_size if in_reference else last_ref_index + 1

        if in_sequence and not in_reference:
            nomenclature.append(
                get_nomenclature_for_insertion(
                    sequence=reference,
                    start=last_ref_index,
                    inserted=sequence[last_seq_index + 1 : last_seq_index + 1 + group_size],
                )
            )
        elif in_reference and not in_sequence:
            nomenclature.append(
                get_nomenclature_for_deletion(
                    sequence=reference,
                    start=last_ref_index + 1,
                    end=last_ref_index + group_size,
                )
            )
        elif is_change:
            for i in range(group_size):
                nomenclature.append(
                    f"{reference[last_ref_index + 1 + i]}{last_ref_index + 2 + i}{sequence[last_seq_index + 1 + i]}"
                )

        last_seq_index = next_seq_index - 1
        last_ref_index = next_ref_index - 1

    return " ".join(nomenclature)
