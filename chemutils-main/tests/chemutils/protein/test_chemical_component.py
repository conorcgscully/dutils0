import pytest

from chemutils.molecule import smiles_from_oemol
from chemutils.protein import get_chemical_component


@pytest.mark.parametrize(
    ["chemical_component", "smiles", "atom_names"],
    [
        ("ALA", "C[C@@H](C(=O)O)N", ["CB", "CA", "C", "O", "OXT", "N"]),
        (
            "DA",
            "c1nc(c2c(n1)n(cn2)[C@H]3C[C@@H]([C@H](O3)COP(=O)(O)O)O)N",
            [
                "C2",
                "N1",
                "C6",
                "C5",
                "C4",
                "N3",
                "N9",
                "C8",
                "N7",
                "C1'",
                "C2'",
                "C3'",
                "C4'",
                "O4'",
                "C5'",
                "O5'",
                "P",
                "OP1",
                "OP3",
                "OP2",
                "O3'",
                "N6",
            ],
        ),
        ("PEG", "C(COCCO)O", ["C1", "C2", "O2", "C3", "C4", "O4", "O1"]),
    ],
)
def test_get_chemical_component(chemical_component, smiles, atom_names):
    oemol = get_chemical_component(chemical_component)
    assert smiles_from_oemol(oemol) == smiles
    assert [atom.GetName().strip() for atom in oemol.GetAtoms()] == atom_names


def test_get_chemical_component_missing():
    with pytest.raises(KeyError):
        get_chemical_component("MISSING")
