from collections import namedtuple

import pytest

from chemutils.protein.alignment import SequenceAlignment
from chemutils.protein.alignment.to_nomenclature import get_mutation_nomenclature_from_alignment

Case = namedtuple(
    "Point", "sequence reference aligned_sequence aligned_reference nomenclature pairs"
)

CASES = [
    # Simple 1-to-1 match
    Case(
        sequence="ACE",
        reference="ACE",
        aligned_sequence="ACE",
        aligned_reference="ACE",
        nomenclature="",
        pairs=(("A", "A"), ("C", "C"), ("E", "E")),
    ),
    # Deletion at start
    Case(
        sequence="CE",
        reference="ACE",
        aligned_sequence="-CE",
        aligned_reference="ACE",
        nomenclature="A1del",
        pairs=((None, "A"), ("C", "C"), ("E", "E")),
    ),
    # Deletion in middle
    Case(
        sequence="AE",
        reference="ACE",
        aligned_sequence="A-E",
        aligned_reference="ACE",
        nomenclature="C2del",
        pairs=(("A", "A"), (None, "C"), ("E", "E")),
    ),
    # Deletion at end
    Case(
        sequence="AC",
        reference="ACE",
        aligned_sequence="AC-",
        aligned_reference="ACE",
        nomenclature="E3del",
        pairs=(("A", "A"), ("C", "C"), (None, "E")),
    ),
    # Deletion of two residues at start
    Case(
        sequence="EG",
        reference="ACEG",
        aligned_sequence="--EG",
        aligned_reference="ACEG",
        nomenclature="A1_C2del",
        pairs=((None, "A"), (None, "C"), ("E", "E"), ("G", "G")),
    ),
    # Deletion of two residues in middle
    Case(
        sequence="AG",
        reference="ACEG",
        aligned_sequence="A--G",
        aligned_reference="ACEG",
        nomenclature="C2_E3del",
        pairs=(("A", "A"), (None, "C"), (None, "E"), ("G", "G")),
    ),
    # Deletion of two residues at end
    Case(
        sequence="AC",
        reference="ACEG",
        aligned_sequence="AC--",
        aligned_reference="ACEG",
        nomenclature="E3_G4del",
        pairs=(("A", "A"), ("C", "C"), (None, "E"), (None, "G")),
    ),
    # Deletion of three residues at start
    Case(
        sequence="GI",
        reference="ACEGI",
        aligned_sequence="---GI",
        aligned_reference="ACEGI",
        nomenclature="A1_E3del",
        pairs=((None, "A"), (None, "C"), (None, "E"), ("G", "G"), ("I", "I")),
    ),
    # Deletion of three residues in middle
    Case(
        sequence="AI",
        reference="ACEGI",
        aligned_sequence="A---I",
        aligned_reference="ACEGI",
        nomenclature="C2_G4del",
        pairs=(("A", "A"), (None, "C"), (None, "E"), (None, "G"), ("I", "I")),
    ),
    # Deletion of three residues at end
    Case(
        sequence="AC",
        reference="ACEGI",
        aligned_sequence="AC---",
        aligned_reference="ACEGI",
        nomenclature="E3_I5del",
        pairs=(("A", "A"), ("C", "C"), (None, "E"), (None, "G"), (None, "I")),
    ),
    # Insertion at start
    Case(
        sequence="ACE",
        reference="CE",
        aligned_sequence="ACE",
        aligned_reference="-CE",
        nomenclature="_C1insA",
        pairs=(("A", None), ("C", "C"), ("E", "E")),
    ),
    # Insertion in middle
    Case(
        sequence="ACE",
        reference="AE",
        aligned_sequence="ACE",
        aligned_reference="A-E",
        nomenclature="A1_E2insC",
        pairs=(("A", "A"), ("C", None), ("E", "E")),
    ),
    # Insertion at end
    Case(
        sequence="ACE",
        reference="AC",
        aligned_sequence="ACE",
        aligned_reference="AC-",
        nomenclature="C2_insE",
        pairs=(("A", "A"), ("C", "C"), ("E", None)),
    ),
    # Insertion of two residues at start
    Case(
        sequence="ACEG",
        reference="EG",
        aligned_sequence="ACEG",
        aligned_reference="--EG",
        nomenclature="_E1insAC",
        pairs=(("A", None), ("C", None), ("E", "E"), ("G", "G")),
    ),
    # Insertion of two residues in middle
    Case(
        sequence="ACEG",
        reference="AG",
        aligned_sequence="ACEG",
        aligned_reference="A--G",
        nomenclature="A1_G2insCE",
        pairs=(("A", "A"), ("C", None), ("E", None), ("G", "G")),
    ),
    # Insertion of two residues at end
    Case(
        sequence="ACEG",
        reference="AC",
        aligned_sequence="ACEG",
        aligned_reference="AC--",
        nomenclature="C2_insEG",
        pairs=(("A", "A"), ("C", "C"), ("E", None), ("G", None)),
    ),
    # Insertion of three residues at start
    Case(
        sequence="ACEGI",
        reference="GI",
        aligned_sequence="ACEGI",
        aligned_reference="---GI",
        nomenclature="_G1insACE",
        pairs=(("A", None), ("C", None), ("E", None), ("G", "G"), ("I", "I")),
    ),
    # Insertion of three residues in middle
    Case(
        sequence="ACEGI",
        reference="AI",
        aligned_sequence="ACEGI",
        aligned_reference="A---I",
        nomenclature="A1_I2insCEG",
        pairs=(("A", "A"), ("C", None), ("E", None), ("G", None), ("I", "I")),
    ),
    # Insertion of three residues at end
    Case(
        sequence="ACEGI",
        reference="AC",
        aligned_sequence="ACEGI",
        aligned_reference="AC---",
        nomenclature="C2_insEGI",
        pairs=(("A", "A"), ("C", "C"), ("E", None), ("G", None), ("I", None)),
    ),
    # Mutation at start
    Case(
        sequence="KCE",
        reference="ACE",
        aligned_sequence="KCE",
        aligned_reference="ACE",
        nomenclature="A1K",
        pairs=(("K", "A"), ("C", "C"), ("E", "E")),
    ),
    # Mutation in middle
    Case(
        sequence="AKE",
        reference="ACE",
        aligned_sequence="AKE",
        aligned_reference="ACE",
        nomenclature="C2K",
        pairs=(("A", "A"), ("K", "C"), ("E", "E")),
    ),
    # Mutation at end
    Case(
        sequence="ACK",
        reference="ACE",
        aligned_sequence="ACK",
        aligned_reference="ACE",
        nomenclature="E3K",
        pairs=(("A", "A"), ("C", "C"), ("K", "E")),
    ),
    # Two mutations at start
    Case(
        sequence="KPE",
        reference="ACE",
        aligned_sequence="KPE",
        aligned_reference="ACE",
        nomenclature="A1K C2P",
        pairs=(("K", "A"), ("P", "C"), ("E", "E")),
    ),
    # Two mutations at end
    Case(
        sequence="AKP",
        reference="ACE",
        aligned_sequence="AKP",
        aligned_reference="ACE",
        nomenclature="C2K E3P",
        pairs=(("A", "A"), ("K", "C"), ("P", "E")),
    ),
    # Two mutations at either end
    Case(
        sequence="KCP",
        reference="ACE",
        aligned_sequence="KCP",
        aligned_reference="ACE",
        nomenclature="A1K E3P",
        pairs=(("K", "A"), ("C", "C"), ("P", "E")),
    ),
    # Deletion in middle, all same letter (checks regression of bug)
    Case(
        sequence="TT",
        reference="TTT",
        aligned_sequence="T-T",
        aligned_reference="TTT",
        nomenclature="T2del",
        pairs=(("T", "T"), (None, "T"), ("T", "T")),
    ),
    # Insetion in middle, all same letter (checks regression of bug)
    Case(
        sequence="TTT",
        reference="TT",
        aligned_sequence="TTT",
        aligned_reference="T-T",
        nomenclature="T1_T2insT",
        pairs=(("T", "T"), ("T", None), ("T", "T")),
    ),
]


@pytest.mark.parametrize(
    ["sequence", "reference", "aligned_sequence", "aligned_reference", "nomenclature", "pairs"],
    CASES,
)
def test_sequence_alignment_from_pairs(
    sequence, reference, aligned_sequence, aligned_reference, nomenclature, pairs
):
    alignment = SequenceAlignment(pairs)
    assert alignment.pairs == pairs
    assert alignment.nomenclature == nomenclature
    assert alignment.sequence == sequence
    assert alignment.reference == reference
    assert alignment.aligned_sequence == aligned_sequence
    assert alignment.aligned_reference == aligned_reference


@pytest.mark.parametrize(
    ["sequence", "reference", "aligned_sequence", "aligned_reference", "nomenclature", "pairs"],
    CASES,
)
def test_sequence_alignment_from_aligned_sequences(
    sequence, reference, aligned_sequence, aligned_reference, nomenclature, pairs
):
    alignment = SequenceAlignment.from_aligned_sequences(
        aligned_sequence=aligned_sequence, aligned_reference=aligned_reference
    )
    assert alignment.pairs == pairs
    assert alignment.nomenclature == nomenclature
    assert alignment.sequence == sequence
    assert alignment.reference == reference
    assert alignment.aligned_sequence == aligned_sequence
    assert alignment.aligned_reference == aligned_reference


@pytest.mark.parametrize(
    ["sequence", "reference", "aligned_sequence", "aligned_reference", "nomenclature", "pairs"],
    CASES,
)
def test_sequence_alignment_from_nomenclature(
    sequence, reference, aligned_sequence, aligned_reference, nomenclature, pairs
):
    alignment = SequenceAlignment.from_nomenclature(
        sequence=sequence, reference=reference, nomenclature=nomenclature
    )
    assert alignment.pairs == pairs
    assert alignment.nomenclature == nomenclature
    assert alignment.sequence == sequence
    assert alignment.reference == reference
    assert alignment.aligned_sequence == aligned_sequence
    assert alignment.aligned_reference == aligned_reference


@pytest.mark.parametrize(
    ["nomenclature", "pairs"], [(case.nomenclature, case.pairs) for case in CASES]
)
def test_get_mutation_nomenclature_from_alignment(nomenclature, pairs):
    assert get_mutation_nomenclature_from_alignment(pairs) == nomenclature


def test_sequence_alignment_hash():
    alignment = SequenceAlignment.from_aligned_sequences(
        aligned_sequence="AB-C", aligned_reference="-BCE"
    )
    assert isinstance(hash(alignment), int)
