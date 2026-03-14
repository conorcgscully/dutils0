import pytest

from chemutils.props import (
    MEDCHEM_PROPERTIES,
    calculate_properties_for_representative_smiles,
)

EXPECTED = {
    "% Uncharged at pH 6": 100.0,
    "% Uncharged at pH 7.4": 100.0,
    "2D TPSA": 61.82,
    "BBB Score": 4.52,
    "CNS MPO": 6.0,
    "CXSMILES": "Cn1cnc2c1c(=O)n(c(=O)n2C)C",
    "Chemaxon LogD": -0.55,
    "Chemaxon LogP": -0.55,
    "Conjugate Acid pKa": -1.15,
    "Fraction C sp3": 0.38,
    "Is Macrocycle": False,
    "Molecular Weight": 194.2,
    "Num Aromatic Rings": 2,
    "Num Chiral Atoms": 0,
    "Num HBA": 3,
    "Num HBA Lone Pairs": 5,
    "Num HBD": 0,
    "Num HBD Hydrogens": 0,
    "Num Heavy Atoms": 14,
    "Num Lipinski HBA": 6,
    "Num Lipinski HBD": 0,
    "Num Rotatable Bonds": 0,
    "Num Spiro Atoms": 0,
    "OpenEye XLogP": -0.61,
    "QED": 0.54,
    "SA Score": 2.30,
    "SMILES": "Cn1cnc2c1c(=O)n(c(=O)n2C)C",
    "Spacial Score": 11.07,
    "Wildman-Crippen ALogP": -1.03,
    "pKa": None,
    "Num Red Alerts": 0,
    "Num Amber Alerts": 0,
    "Num PAINS Alerts": 0,
}


def test_medchem_properties():
    smiles = "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"
    props = calculate_properties_for_representative_smiles(smiles, MEDCHEM_PROPERTIES)

    for key, calculated in props.items():
        assert calculated == pytest.approx(EXPECTED[key], abs=0.01), f"Property {key} did not match"
