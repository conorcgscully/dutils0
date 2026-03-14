import numpy as np
import pytest

from chemutils.fs import read_molecule
from chemutils.protein.alignment import SequenceAlignment
from chemutils.tmalign import TMAlignError, TMAlignResult, tmalign


@pytest.fixture
def proteins():
    return {
        "4OQG": read_molecule("tests/data/4OQG.cif.gz"),
        "8IWV": read_molecule("tests/data/8IWV.cif.gz"),
    }


@pytest.mark.parametrize(
    ["protein", "reference_protein", "protein_chain_id", "reference_protein_chain_id", "expected"],
    [
        (
            "4OQG",
            "8IWV",
            "A",
            None,
            TMAlignResult(
                aligned_length=258,
                rmsd=1.61,
                sequence_identity=0.372,
                tm_score=0.94053,
                alignment=SequenceAlignment.from_aligned_sequences(
                    aligned_sequence="HPETLVKVKDAEDQLGARVGYIELDLNSGKILESFRPEERFPMMSTFKVLLCGAVLSRIDAGQEQLGRRIHYSQNDLVEYSPVTEKHLTDGMTVRELCSAAITMSDNTAANLLLTTIGGPKELTAFLHNMGDHVTRLDRWEPELNEAIPNDERDTTMPVAMATTLRKLLTGELLTLASRQQLIDWMEADKVAGPLLRSALPAGWFIADKSGAGERGSRGIIAALGPDGKPSRIVVIYTTGSQ--ATMDERNRQIAEIGASLIKHW",
                    aligned_reference="-SV-DKQLAELERNANGRLGVAMINTGN-GTKILYRAAQRFPFCSTFKFMLAAAVLDQSQSQPNLLNKHINYHESDLLSYAPITRKNLAHGMTVSELCAATIQYSDNTAANLLIKELGGLAAVNQFARSIGDQMFRLDRWEPDLNTARPNDPRDTTTPAAMAASMNKLVLGDALRPAQRSQLAVWLKGNTTGDATIRAGAPTDWIVGDKTGSGDYGTTNDIAVLWPTKGAPIVLVVYFTQREKDAKP--RRDVLASVTKIILSQI",
                ),
                transformation=np.array(
                    [
                        [-0.4939, -0.8417, -0.2183, 34.8621],
                        [0.8514, -0.4171, -0.3180, 17.1313],
                        [0.1767, -0.3429, 0.9226, 9.2706],
                        [0, 0, 0, 1],
                    ]
                ),
            ),
        ),
        # Align two chains in the same protein
        (
            "4OQG",
            "4OQG",
            "A",
            "B",
            TMAlignResult(
                aligned_length=263,
                rmsd=0.4,
                sequence_identity=1.0,
                tm_score=0.99571,
                alignment=SequenceAlignment.from_aligned_sequences(
                    aligned_sequence="HPETLVKVKDAEDQLGARVGYIELDLNSGKILESFRPEERFPMMSTFKVLLCGAVLSRIDAGQEQLGRRIHYSQNDLVEYSPVTEKHLTDGMTVRELCSAAITMSDNTAANLLLTTIGGPKELTAFLHNMGDHVTRLDRWEPELNEAIPNDERDTTMPVAMATTLRKLLTGELLTLASRQQLIDWMEADKVAGPLLRSALPAGWFIADKSGAGERGSRGIIAALGPDGKPSRIVVIYTTGSQATMDERNRQIAEIGASLIKHW",
                    aligned_reference="HPETLVKVKDAEDQLGARVGYIELDLNSGKILESFRPEERFPMMSTFKVLLCGAVLSRIDAGQEQLGRRIHYSQNDLVEYSPVTEKHLTDGMTVRELCSAAITMSDNTAANLLLTTIGGPKELTAFLHNMGDHVTRLDRWEPELNEAIPNDERDTTMPVAMATTLRKLLTGELLTLASRQQLIDWMEADKVAGPLLRSALPAGWFIADKSGAGERGSRGIIAALGPDGKPSRIVVIYTTGSQATMDERNRQIAEIGASLIKHW",
                ),
                transformation=np.array(
                    [
                        [0.9996, 0.0204, 0.0213, -8.7645],
                        [-0.0208, 0.9996, 0.0188, 37.7466],
                        [-0.0209, -0.0193, 0.9996, 1.0437],
                        [0, 0, 0, 1],
                    ]
                ),
            ),
        ),
        # Test self-alignment
        (
            "8IWV",
            "8IWV",
            None,
            None,
            TMAlignResult(
                aligned_length=260,
                rmsd=0.0,
                sequence_identity=1.0,
                tm_score=1.0,
                alignment=SequenceAlignment.from_aligned_sequences(
                    aligned_sequence="SVDKQLAELERNANGRLGVAMINTGNGTKILYRAAQRFPFCSTFKFMLAAAVLDQSQSQPNLLNKHINYHESDLLSYAPITRKNLAHGMTVSELCAATIQYSDNTAANLLIKELGGLAAVNQFARSIGDQMFRLDRWEPDLNTARPNDPRDTTTPAAMAASMNKLVLGDALRPAQRSQLAVWLKGNTTGDATIRAGAPTDWIVGDKTGSGDYGTTNDIAVLWPTKGAPIVLVVYFTQREKDAKPRRDVLASVTKIILSQI",
                    aligned_reference="SVDKQLAELERNANGRLGVAMINTGNGTKILYRAAQRFPFCSTFKFMLAAAVLDQSQSQPNLLNKHINYHESDLLSYAPITRKNLAHGMTVSELCAATIQYSDNTAANLLIKELGGLAAVNQFARSIGDQMFRLDRWEPDLNTARPNDPRDTTTPAAMAASMNKLVLGDALRPAQRSQLAVWLKGNTTGDATIRAGAPTDWIVGDKTGSGDYGTTNDIAVLWPTKGAPIVLVVYFTQREKDAKPRRDVLASVTKIILSQI",
                ),
                transformation=np.eye(4),
            ),
        ),
    ],
)
def test_tmalign(
    proteins,
    protein,
    reference_protein,
    protein_chain_id,
    reference_protein_chain_id,
    expected: TMAlignResult,
):
    protein = proteins[protein]
    reference_protein = proteins[reference_protein]

    result = tmalign(
        protein=protein,
        reference_protein=reference_protein,
        chain_id=protein_chain_id,
        reference_chain_id=reference_protein_chain_id,
    )

    assert result.aligned_length == expected.aligned_length
    assert result.rmsd == pytest.approx(expected.rmsd, abs=1e-2)
    assert result.sequence_identity == pytest.approx(expected.sequence_identity, abs=1e-3)
    assert result.tm_score == pytest.approx(expected.tm_score, abs=1e-4)
    assert result.alignment == expected.alignment
    np.testing.assert_allclose(result.transformation, expected.transformation, atol=1e-1)


def test_tmalign_no_chain_id(proteins):
    protein = proteins["4OQG"]
    reference_protein = proteins["8IWV"]

    with pytest.raises(TMAlignError) as e:
        tmalign(protein=protein, reference_protein=reference_protein)

    assert (
        str(e.value)
        == "TMAlign error: Multiple protein chains found in protein. Specify a chain ID with `chain_id`."
    )


def test_tmalign_no_reference_chain_id(proteins):
    protein = proteins["8IWV"]
    reference_protein = proteins["4OQG"]

    with pytest.raises(TMAlignError) as e:
        tmalign(protein=protein, reference_protein=reference_protein)

    assert (
        str(e.value)
        == "TMAlign error: Multiple protein chains found in reference protein. Specify a chain ID with `reference_chain_id`."
    )


def test_tmalign_missing_chain_id(proteins):
    protein = proteins["4OQG"]
    reference_protein = proteins["8IWV"]

    with pytest.raises(TMAlignError) as e:
        tmalign(protein=protein, reference_protein=reference_protein, chain_id="Z")

    assert str(e.value) == "TMAlign error: Chain ID `Z` is not a protein or peptide chain."


def test_tmalign_missing_reference_chain_id(proteins):
    protein = proteins["8IWV"]
    reference_protein = proteins["4OQG"]

    with pytest.raises(TMAlignError) as e:
        tmalign(protein=protein, reference_protein=reference_protein, reference_chain_id="Z")

    assert (
        str(e.value) == "TMAlign error: Reference chain ID `Z` is not a protein or peptide chain."
    )
