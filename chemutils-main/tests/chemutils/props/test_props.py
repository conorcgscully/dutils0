import fsutils as fs
import pytest
from rdkit import Chem

from chemutils.molecule.smiles import oemol_from_smiles
from chemutils.props import (
    ALL_PROPERTIES,
    calculate_properties_for_representative_smiles,
    calculate_properties_for_smiles,
)
from chemutils.props.atomic import NumberHeavyAtoms
from chemutils.props.logp import OEXLogP, WildmanCrippenLogP
from chemutils.rdmol_from_oemol import rdmol_from_oemol

DRUG_CASES = fs.read_yaml("tests/data/drugs.yaml")


def test_custom_props_name():
    props = calculate_properties_for_smiles("CC=O", {"my_prop": NumberHeavyAtoms})
    assert props == {"my_prop": 3}


@pytest.mark.parametrize(
    "smiles",
    [
        r"CC12C3([Ir]1456(C2(C4(C53C)C)C)[N]7=CC=C(C=C7C(=O)[N]6(CCCc8ccc(cc8)S(=O)(=O)N)Cl)O)C",
    ],
)
def test_all_props_no_exceptions(smiles):
    _ = calculate_properties_for_smiles(smiles, ALL_PROPERTIES)


@pytest.mark.parametrize(
    ("smiles", "name", "expected"),
    [(case["smiles"], case["name"], case) for case in DRUG_CASES],
)
def test_all_props(smiles, name, expected):
    result = calculate_properties_for_representative_smiles(smiles, ALL_PROPERTIES)
    result["morgan_fingerprint"] = "".join(f"{x:02x}" for x in result["morgan_fingerprint"])
    result["name"] = name
    assert result.keys() == expected.keys()
    for key in result:
        assert expected[key] == pytest.approx(result[key], abs=0.01), f"Property {key} not equal"


def test_prop_no_rdkit():
    smiles = "CC=O"
    oemol = oemol_from_smiles(smiles)
    rdmol = rdmol_from_oemol(oemol)
    with pytest.raises(NotImplementedError):
        _ = OEXLogP.get_rdkit(rdmol)


def test_prop_no_openeye_fallback(monkeypatch):
    smiles = "CC=O"
    oemol = oemol_from_smiles(smiles)
    count = 0

    def counter(*args, **kwargs):
        nonlocal count
        count += 1
        return 1

    monkeypatch.setattr("chemutils.props.logp.WildmanCrippenLogP.get_rdkit", counter)
    assert WildmanCrippenLogP.get_openeye(oemol) == 1
    assert count == 1


def test_prop_get_rdkit_prop():
    smiles = "CC=O"
    oemol = oemol_from_smiles(smiles)
    rdmol = Chem.MolFromSmiles(smiles)

    assert (
        len(
            {
                WildmanCrippenLogP.get(smiles),
                WildmanCrippenLogP.get(oemol),
                WildmanCrippenLogP.get(rdmol),
                WildmanCrippenLogP.get_openeye(oemol),
                WildmanCrippenLogP.get_rdkit(rdmol),
            }
        )
        == 1
    )
