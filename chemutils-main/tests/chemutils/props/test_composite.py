import pytest

from chemutils.props import (
    BloodBrainBarrierScore,
    MolecularWeight,
    NumberAromaticRings,
    NumberHeavyAtoms,
    NumberHydrogenAcceptors,
    NumberHydrogenDonors,
    PKaConjugateAcid,
    TwoDimensionalTPSA,
    get_properties_to_calculate,
)


@pytest.mark.parametrize(
    ["properties", "required_properties"],
    [
        (
            [BloodBrainBarrierScore],
            [
                PKaConjugateAcid,
                TwoDimensionalTPSA,
                NumberHydrogenDonors,
                NumberHydrogenAcceptors,
                MolecularWeight,
                NumberHeavyAtoms,
                NumberAromaticRings,
                BloodBrainBarrierScore,
            ],
        ),
        # Moves MolecularWeight to before BloodBrainBarrierScore
        (
            [BloodBrainBarrierScore, MolecularWeight],
            [
                PKaConjugateAcid,
                TwoDimensionalTPSA,
                NumberHydrogenDonors,
                NumberHydrogenAcceptors,
                MolecularWeight,
                NumberHeavyAtoms,
                NumberAromaticRings,
                BloodBrainBarrierScore,
            ],
        ),
    ],
)
def test_get_properties_to_calculate(properties, required_properties):
    props = get_properties_to_calculate(properties=properties)
    assert props == required_properties
